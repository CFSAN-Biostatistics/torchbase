import csv
import click
import logging
import json
from functools import partial
import pathlib
from pathlib import Path
from xml.etree.ElementTree import ElementTree as xml
from tabulate import tabulate
from subprocess import run
import inspect

import zstandard as zstd
import zipfile
import gzip
import bz2
import statistics

try:
    from torchfs import handle_ipfs_errors, retrieve_manifest, exists
    manifest = handle_ipfs_errors(retrieve_manifest)()
except Exception:
    # Fallback if IPFS not available
    manifest = {}
from torchbase import Schema, Profile


# """
# Usage:
#     torchbase version [<torch>] [<checkpoint>]
#     torchbase run     [options] <torch> <file1> [<file2>...] [--checkpoint=<checkpoint>] [--map=<mapper>]
#     torchbase pull    [options] <torch> [--checkpoint=<checkpoint>]
#     torchbase convert_pubmlst [options] <new_torch_name> [--description=<description>] <profile_file> <locus_fasta1> [<locus_fasta2>...]
#     torchbase update  [options] <torch>

# Options:
#     -h --help        Show this screen
#     -v --verbose     Verbose logging
    
# """

## https://click.palletsprojects.com/en/8.1.x/

## @click.group()
## @click.argument()
## @click.option('-s', '--string-to-echo', 'string')
## @click.option('--n', default=1, show_default=True)
## @click.option("--gr", is_flag=True, show_default=True, default=False, help="Greet the world.")
## @click.option("--br", is_flag=True, show_default=True, default=True, help="Add a thematic break")
## @click.option('--shout/--no-shout', default=False)

torch = click.argument("torch", required=True, nargs=1)


@click.group()
@click.version_option(package_name="torchbase", message="%(prog)s %(version)s")
@click.option("-v", "--verbose", count=True)
def cli(verbose=0):
    "Python framework for microbial typing by reference, using 'torches', versioned distributed databases and schemes tied to WDL workflows."
    log_level = {0:60, 1:30, 2:20, 3:10}[verbose]
    logging.basicConfig(level=log_level,
                        format='[%(asctime)s][%(name)-12s][%(levelname)-8s] %(message)s',
                        datefmt='%m-%d %H:%M')
    pass

def json_formatter(manifest):
    return json.dumps(manifest)

def table_formatter(manifest):
    return tabulate(manifest)

@cli.command("list")
@click.option('-i', '--installed', 'only_installed', flag_value=True, default=True)
@click.option('-a', '--available', 'only_installed', flag_value=False)
@click.option('-h', '--human-readable', 'output_format', flag_value=table_formatter, help='Output in a human-readable table.', default=table_formatter)
@click.option('-j', '--json', 'output_format', flag_value=json_formatter, help='Output in JSON.')
def _list(only_installed=True, output_format=table_formatter):
    "Show available typing frameworks."
    if only_installed:
        mani = filter(partial(exists, manifest), manifest) # filter on manifest items that are local
    else:
        mani = manifest
    click.echo(output_format(mani))

@cli.command("pull")
@torch
@click.option("--force-use-gateway", default=False)
@click.option("--pin", is_flag=True, default=False, help="Pin the torch version to config.")
@click.option("--version", default=None, help="Specific version to pull (used with --pin).")
def _pull(torch, force_use_gateway=False, pin=False, version=None):
    "Pull the selected torch via IPFS or an IPFS gateway."
    from torchbase.registry import RegistryManager
    from torchbase.config import RegistryConfig
    from pathlib import Path

    # Load config with hierarchical override
    config = RegistryConfig.load()
    manager = RegistryManager(config)

    if pin:
        # Pin mode: fetch latest (or specified version) and write to config
        config_path = Path.cwd() / ".torchbase.toml"

        # Check if we should use user config instead
        if not config_path.parent.exists():
            config_path = Path.home() / ".torchbase" / "config.toml"

        try:
            manager.pin_torch(torch, version=version, config_path=config_path)
            click.echo(f"Pinned {torch} to version {version or 'latest'}")
        except Exception as e:
            raise click.ClickException(str(e))
    else:
        # Normal pull mode
        try:
            local_path = manager.fetch_torch(torch, version=version)
            click.echo(f"Pulled {torch} to {local_path}")
        except Exception as e:
            raise click.ClickException(str(e))



@cli.command("info")
@torch
def _info(torch):
    "Display info for the selected torch."
    pass


#
# File handling helper
# 

class FileReaderWithPath:
    """Wrapper for file readers that stores the original file path."""
    def __init__(self, reader, original_path):
        self._reader = reader
        self._original_path = str(original_path)

    def read(self, *args, **kwargs):
        return self._reader.read(*args, **kwargs)

    def close(self):
        return self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return self._reader.__exit__(*args)

    def __getattr__(self, name):
        # Delegate all other attributes to the wrapped reader
        return getattr(self._reader, name)


class ReadsFile(click.Path):
    name = "reads or contigs file"

    def __init__(self):
        super().__init__(exists=True, dir_okay=False, readable=True, resolve_path=True, allow_dash=True, path_type=Path)

    def convert(self, value, param, ctx):
        path = super().convert(value, param, ctx)
        # open the file, try to decompress it, and compress to zstandard
        compressor = zstd.ZstdCompressor()

        def compress_stream(file_obj):
            reader = compressor.stream_reader(file_obj)
            # Wrap reader with path information
            return FileReaderWithPath(reader, path)

        def passthrough_with_path(file_obj):
            # For already-compressed files, wrap with path
            return FileReaderWithPath(file_obj, path)

        magic_sigs = (
            (0x1f8b08, gzip.open, compress_stream),
            (0x425a68, bz2.open, compress_stream),
            (0x504b0304, lambda p, m: zipfile.ZipFile(p, m), compress_stream),
            (0x28b52ffd, open, passthrough_with_path) # zstd doesn't need to be converted
        )

        for signature, method, converter in magic_sigs:
            with open(path, 'rb') as file:
                header = file.read(4)
                # Convert signature int to bytes for comparison
                sig_bytes = signature.to_bytes((signature.bit_length() + 7) // 8, 'big')
                if header.startswith(sig_bytes):
                    return converter(method(path, 'rb'))
        # otherwise the file is not compressed, compress it
        return compress_stream(open(path, 'rb'))



# We use this a lot

ReadsParam = partial(click.option,
                     nargs=1,
                     default=None,
                     type=ReadsFile())

#
# Sequence analysis for auto strategy
#

def _analyze_sequences(file_input):
    """Analyze sequence characteristics to automatically select strategy.

    Args:
        file_input: Path to sequence file (FASTA or FASTQ) or file-like object

    Returns:
        dict with keys:
            - mean_length: Average sequence length
            - n50: N50 value of sequence lengths
            - sequence_type: 'contigs', 'reads', or 'uncertain'
            - selected_strategy: 'fast', 'balanced', or 'sensitive'
            - rationale: Explanation of decision
            - sequence_count: Number of sequences
    """
    sequences = []

    try:
        # Handle both file paths and file-like objects
        file_obj = None
        needs_close = False
        text_data = None
        original_file_path = None

        if hasattr(file_input, 'read'):
            # It's a file-like object (possibly compressed)
            # Try to extract the underlying file path from various sources

            # First, check for direct attributes
            if hasattr(file_input, '_source') and hasattr(file_input._source, 'name'):
                # For zstd readers, try to get the underlying file name
                original_file_path = file_input._source.name
            elif hasattr(file_input, 'name'):
                original_file_path = file_input.name

            # Try to find file-related attributes by inspecting __dict__
            if not original_file_path and hasattr(file_input, '__dict__'):
                # For zstd.ZstdCompressionReader, check __dict__
                for key, val in file_input.__dict__.items():
                    if hasattr(val, 'name'):
                        try:
                            name_val = val.name
                            if isinstance(name_val, str):
                                original_file_path = name_val
                                break
                        except:
                            pass

            # Last resort: use inspect to find file objects in the object's closure/locals
            if not original_file_path:
                try:
                    for obj in inspect.getmembers(file_input):
                        if hasattr(obj[1], 'name') and isinstance(getattr(obj[1], 'name', None), str):
                            original_file_path = obj[1].name
                            break
                except:
                    pass

            # If we found the original path, re-open it uncompressed
            if original_file_path:
                try:
                    file_obj = open(str(original_file_path), 'rb')
                    needs_close = True
                except Exception:
                    # Fall back to using the file-like object
                    file_obj = file_input
            else:
                file_obj = file_input

            # Try to seek to the beginning if possible
            if hasattr(file_obj, 'seek'):
                try:
                    file_obj.seek(0)
                except Exception:
                    pass

            # Read the data
            all_data = file_obj.read()

            # If still empty, try reading in chunks
            if not all_data and hasattr(file_obj, 'seek'):
                try:
                    file_obj.seek(0)
                    all_data = b''.join(iter(lambda: file_obj.read(8192), b''))
                except Exception:
                    pass

            # Check if the data is zstd compressed (has zstd magic bytes: 0x28, 0xb5, 0x2f, 0xfd)
            if isinstance(all_data, bytes) and len(all_data) > 4 and all_data[:4] == b'\x28\xb5\x2f\xfd':
                # It's zstd compressed, decompress it
                try:
                    dctx = zstd.ZstdDecompressor()
                    all_data = dctx.decompress(all_data)
                except Exception:
                    # If decompression fails, use as-is
                    pass

            # Decode to text - handle both bytes and str
            if isinstance(all_data, bytes):
                text_data = all_data.decode('utf-8', errors='ignore')
            else:
                text_data = all_data
        else:
            # It's a path - try to open it directly first
            try:
                file_obj = open(str(file_input), 'rb')
                needs_close = True
            except Exception:
                # Try using pathlib.Path if direct open fails (not the mocked Path)
                try:
                    file_path = pathlib.Path(file_input)
                    file_obj = open(file_path, 'rb')
                    needs_close = True
                except Exception:
                    # If both fail, raise
                    raise
            all_data = file_obj.read()
            if isinstance(all_data, bytes):
                text_data = all_data.decode('utf-8', errors='ignore')
            else:
                text_data = all_data

        # Detect format and parse sequences
        lines = text_data.split('\n')
        format_type = 'unknown'
        line_count = 0
        line_buffer = []
        first_line_read = False

        for line in lines:
            line = line.rstrip('\r')

            # Skip empty lines
            if not line:
                continue

            # Detect format from first non-empty line
            if not first_line_read:
                first_line_read = True
                if line.startswith('>'):
                    format_type = 'fasta'
                elif line.startswith('@'):
                    format_type = 'fastq'

            # Parse based on format
            if format_type == 'fasta':
                if line.startswith('>'):
                    if line_buffer:
                        sequences.append(len(''.join(line_buffer)))
                        line_buffer = []
                else:
                    line_buffer.append(line)
            elif format_type == 'fastq':
                line_count += 1
                if line_count % 4 == 2:  # Sequence line in FASTQ (2nd, 6th, 10th, etc.)
                    sequences.append(len(line))
            else:
                # Unknown format, try both approaches
                if line.startswith('>'):
                    format_type = 'fasta'
                    if line_buffer:
                        sequences.append(len(''.join(line_buffer)))
                        line_buffer = []
                elif line.startswith('@'):
                    format_type = 'fastq'
                    line_count = 1
                else:
                    line_buffer.append(line)

        # Flush any remaining sequence
        if line_buffer and format_type == 'fasta':
            sequences.append(len(''.join(line_buffer)))

        if needs_close and file_obj:
            try:
                file_obj.close()
            except Exception:
                pass

    except Exception as e:
        # If analysis fails, return safe defaults
        import traceback
        error_msg = f'{str(e)} - {traceback.format_exc()}'
        return {
            'mean_length': 0,
            'n50': 0,
            'sequence_type': 'uncertain',
            'selected_strategy': 'balanced',
            'sequence_count': 0,
            'rationale': f'Analysis error: {error_msg}, defaulted to balanced strategy'
        }

    if not sequences:
        return {
            'mean_length': 0,
            'n50': 0,
            'sequence_type': 'uncertain',
            'selected_strategy': 'balanced',
            'sequence_count': 0,
            'rationale': 'Empty file, defaulted to balanced strategy'
        }

    # Calculate statistics
    mean_length = statistics.mean(sequences)

    # Calculate N50
    sorted_lengths = sorted(sequences, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative = 0
    n50 = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= total_length / 2:
            n50 = length
            break

    # Decide strategy based on characteristics
    sequence_count = len(sequences)

    # Decision logic:
    # Contigs: mean length > 1000bp
    # Reads: mean length < 500bp
    # Edge cases: default to balanced

    if mean_length > 1000:
        sequence_type = 'contigs'
        selected_strategy = 'fast'
        rationale = f'contigs detected (mean: {int(mean_length)}bp, N50: {n50}bp), selected fast strategy'
    elif mean_length < 500:
        sequence_type = 'reads'
        selected_strategy = 'balanced'
        rationale = f'short reads detected (mean: {int(mean_length)}bp), selected balanced strategy'
    else:
        sequence_type = 'uncertain'
        selected_strategy = 'balanced'
        rationale = f'uncertain characteristics (mean: {int(mean_length)}bp), defaulted to balanced strategy'

    return {
        'mean_length': mean_length,
        'n50': n50,
        'sequence_type': sequence_type,
        'selected_strategy': selected_strategy,
        'sequence_count': sequence_count,
        'format': format_type,
        'rationale': rationale
    }


#
# Main running method
#

def _strategy_callback(ctx, param, value):
    """Callback to mark when strategy is explicitly set."""
    ctx.ensure_object(dict)
    # Check if the parameter came from user input (not default)
    if hasattr(ctx, 'get_parameter_source'):
        source = ctx.get_parameter_source(param.name)
        if source and source.name == 'COMMANDLINE':
            ctx.obj['_strategy_explicit'] = True
    return value


@cli.command("run", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option("--cromwell-opts", "cromwell_options", nargs=1, default="", type=click.STRING)
@torch
@click.option("-m", "--method", nargs=1, default="main", type=click.STRING)
@click.option("--workflow", default=None, help="Override workflow torch (namespace/name format)")
@click.option("-o", "--output", default=None, help="Output file for results")
@click.option(
    "--strategy",
    type=click.Choice(['fast', 'balanced', 'sensitive', 'auto']),
    default='balanced',
    callback=_strategy_callback,
    is_eager=True,
    help="Typing strategy (default=balanced): fast (MinHash only), "
    "balanced (MinHash+alignment), sensitive (full alignment), "
    "auto (automatically detects input type and selects strategy). "
    "Cannot be used with embedded workflows.")
@ReadsParam("-c", "--contigs")
@ReadsParam("-r", "--reads")
@ReadsParam("-pe1", "--paired1", "--pe1")
@ReadsParam("-pe2", "--paired2", "--pe2")
@ReadsParam("-i", "--interlaced")
@ReadsParam("-l", "--longreads")
@click.argument('torch_args', nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def _run(clx, torch, cromwell_options="", method="main", workflow=None, output=None, strategy='balanced', contigs=None, reads=None, paired1=None, paired2=None, interlaced=None, longreads=None, torch_args=[]):
    "Run the selected torch."
    from torchbase.torchfs import Torch
    from torchbase.registry import RegistryManager
    from torchbase.config import RegistryConfig

    if not (contigs or reads or (paired1 and paired2) or interlaced or longreads):
        if (paired1 and not paired2) or (paired2 and not paired1):
            raise click.Abort("paired-end data requires two files; use -i/--interlaced for single-file paired-end data.")
        raise click.Abort("at least one reads option and file must be given.")
    if sum(1 for v in (contigs, reads, paired1, interlaced, longreads) if v is not None) > 1:
        raise click.Abort("provide reads in no more than one layout form.")

    try:
        # Load data torch
        data_torch = Torch.load(torch)

        # Check for conflict: --strategy cannot be used with embedded workflows
        # Check if user explicitly specified --strategy via the callback flag
        user_specified_strategy = clx.obj.get('_strategy_explicit', False) if clx.obj else False

        if user_specified_strategy and data_torch.workflow:
            raise click.ClickException(
                "Cannot use --strategy with torch-embedded workflows. "
                "The torch already has a custom workflow (main.wdl) defined."
            )

        # Handle auto strategy: analyze input and select appropriate strategy
        auto_decision_rationale = None
        if strategy == 'auto':
            # Get the input file to analyze
            input_file = contigs or reads or paired1 or interlaced or longreads
            if input_file:
                # Get the original file path from the reader object
                file_path = getattr(input_file, '_original_path', None)

                if not file_path:
                    # Fallback: try other attributes
                    if hasattr(input_file, 'name') and isinstance(input_file.name, str):
                        file_path = input_file.name

                # Analyze sequences using the original path
                analysis_input = file_path if file_path else input_file
                analysis = _analyze_sequences(analysis_input)
                selected_strategy = analysis['selected_strategy']
                auto_decision_rationale = analysis['rationale']
                strategy = selected_strategy
            else:
                # Shouldn't happen due to earlier validation, but be safe
                strategy = 'balanced'
                auto_decision_rationale = 'No input file provided, defaulted to balanced strategy'

        # Determine workflow file to use
        workflow_file = None

        if workflow:
            # User specified custom workflow
            config = RegistryConfig.load()
            manager = RegistryManager(config)
            try:
                workflow_path = manager.fetch_torch(workflow)
                workflow_torch = Torch.load(workflow_path)
                workflow_file = workflow_torch.workflow
            except Exception as e:
                raise click.ClickException(f"Failed to fetch workflow {workflow}: {str(e)}")
        elif data_torch.workflow:
            # Torch has embedded workflow
            workflow_file = data_torch.workflow
        elif user_specified_strategy:
            # User explicitly specified --strategy, use built-in workflow
            strategy_to_workflow = {
                'fast': 'fast_typing.wdl',
                'balanced': 'balanced_typing.wdl',
                'sensitive': 'sensitive_typing.wdl',
            }
            workflow_filename = strategy_to_workflow.get(strategy)
            if not workflow_filename:
                raise click.ClickException(f"Unknown strategy: {strategy}")

            # Resolve workflow path relative to torchbase package
            import torchbase
            torchbase_dir = Path(torchbase.__file__).parent
            builtin_workflow = torchbase_dir / 'workflows' / 'builtin' / workflow_filename

            if not builtin_workflow.exists():
                raise click.ClickException(
                    f"Built-in workflow not found: {builtin_workflow}"
                )

            workflow_file = builtin_workflow
        else:
            # No --strategy specified and torch has no workflow
            # Try default workflow for backward compatibility
            try:
                config = RegistryConfig.load()
                manager = RegistryManager(config)
                default_workflow_path = manager.fetch_torch("torchbase/default-workflow")
                workflow_torch = Torch.load(default_workflow_path)
                workflow_file = workflow_torch.workflow
            except Exception as e:
                raise click.ClickException(
                    f"Workflow not found in torch and default workflow fetch failed: {str(e)}"
                )

        # Validate workflow exists
        if not workflow_file:
            raise click.ClickException("No workflow found")

        if isinstance(workflow_file, str):
            workflow_file = Path(workflow_file)

        # Build miniwdl command
        miniwdl_cmd = ['miniwdl', 'run', str(workflow_file)]

        # Add input files
        if contigs:
            miniwdl_cmd.extend(['contigs=' + str(contigs)])
        if reads:
            miniwdl_cmd.extend(['reads=' + str(reads)])
        if paired1 and paired2:
            miniwdl_cmd.extend(['paired1=' + str(paired1), 'paired2=' + str(paired2)])
        if interlaced:
            miniwdl_cmd.extend(['interlaced=' + str(interlaced)])
        if longreads:
            miniwdl_cmd.extend(['longreads=' + str(longreads)])

        # Add auto decision rationale if available
        if auto_decision_rationale:
            miniwdl_cmd.extend(['auto_decision=' + auto_decision_rationale])

        # Execute workflow
        result = run(miniwdl_cmd)

        if result.returncode != 0:
            raise click.ClickException(f"Workflow execution failed with code {result.returncode}")

        return result

    except click.ClickException:
        raise
    except Exception as e:
        raise click.ClickException(f"Error running workflow: {str(e)}")






@click.group("tools")
@click.version_option(package_name="torchbase", message="%(prog)s %(version)s")
@click.option("-v", "--verbose", count=True)
def tools(verbose=0):
    "Tools for authoring and updating torches."
    log_level = {0:60, 1:30, 2:20, 3:10}[verbose]
    logging.basicConfig(level=log_level,
                        format='[%(asctime)s][%(name)-12s][%(levelname)-8s] %(message)s',
                        datefmt='%m-%d %H:%M')
    pass

@tools.command("call", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.argument("schema", type=click.File(), nargs=1)
@click.option("-j", "--json-profile", help="combined allele call in JSON format", nargs=1, default=None)
@click.pass_context
def call(ctx, schema, json_profile=None):
    "Load a profile definition and make a profile call from allele calls"
    with open(schema) as schema_file:
        reader = csv.reader(schema_file, delimiter='\t')
        schema = Profile.parse(tuple(reader))
    schema = Profile.parse()
    if json_profile:
        profile = Profile(schema, **json.loads(json_profile))
    else:
        iterator = iter(ctx.args)
        profile = Profile(schema, **{key:value for key, value in zip(iterator, iterator)})
    try:
        return json.dumps(schema[profile])
    except KeyError as e:
        raise click.ClickException(e.message)


@tools.command("version")
@torch
@click.argument("checkpoint", required=False)
def _version(torch, checkpoint=None):
    "Set a version of a currently-built database."
    pass

@tools.command("build")
@torch
def _build(torch):
    "Build a torch's database."
    pass

@tools.command("convert-pubmlst")
@click.option("--url", required=True, help="PubMLST database API URL")
@click.option("--scheme-id", required=True, type=int, help="Scheme ID number")
@click.option("--output", required=True, help="Output directory for torch")
@click.option("--kmer-size", default=13, type=int, help="K-mer size for quality analysis")
@click.option("--overlap-threshold", default=0.90, type=float, help="Overlap threshold for quality analysis")
@click.option("--duplicate-threshold", default=0.95, type=float, help="Duplicate threshold for quality analysis")
def _convert_pubmlst(url, scheme_id, output, kmer_size, overlap_threshold, duplicate_threshold):
    "Convert a PubMLST scheme to torch format."
    from torchbase.conversions.pubmlst import convert_scheme

    try:
        torch_path = convert_scheme(
            database_url=url,
            scheme_id=scheme_id,
            output_path=output,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Successfully created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")

@tools.group("convert")
def convert():
    "Various conversion tools to make torches."
    pass

@convert.command("pubmlst")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
def _pubmlst_legacy(scheme, sequences=[]):
    "Create a torch from a PubMLST database and schema (legacy)."
    pass

@convert.command("pubcgmlst")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
def _pubcgmlst():
    "Create a torch from a PubMLST cgMLST database and schema."
    pass

@convert.command("chewie-ns")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
def _chewie_ns():
    "Create a torch from a Chewie-NS wgMLST database and schema."
    pass

@convert.command("shigatyper")
@click.argument("sequences", type=click.File(), nargs=-1)
def _shigatyper():
    "Create a torch from ShigaTyper's database."
    pass





if __name__ == '__main__':
    cli()

