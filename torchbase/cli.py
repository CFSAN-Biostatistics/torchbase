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
import shutil
import tempfile

import requests
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
    from torchbase.torchfs import Torch
    from torchbase.registry import RegistryManager
    from torchbase.config import RegistryConfig
    import toml

    config = RegistryConfig.load()
    manager = RegistryManager(config)
    try:
        torch_path = manager.fetch_torch(torch)
    except Exception as e:
        raise click.ClickException(str(e))

    metadata_path = torch_path / "metadata.toml"
    if not metadata_path.exists():
        raise click.ClickException(f"metadata.toml not found at {torch_path}")

    with open(metadata_path, encoding="utf-8") as f:
        metadata = toml.load(f)

    rows = []
    for key in ("namespace", "name", "version"):
        if key in metadata:
            rows.append((key, metadata[key]))
    for section in ("provenance", "typing"):
        for k, v in metadata.get(section, {}).items():
            rows.append((f"{section}.{k}", v))
    click.echo(tabulate(rows, headers=["Field", "Value"], tablefmt="simple"))


@cli.group("workflow")
def workflow():
    "Workflow management and inspection commands."
    pass


@workflow.command("inspect")
@click.argument("workflow_spec", required=True)
@click.option("--verbose", is_flag=True, default=False, help="Show detailed parameter information.")
def inspect(workflow_spec, verbose=False):
    "Inspect a workflow and display its structure as an ASCII diagram."
    from torchbase.workflow_inspect import inspect_workflow

    try:
        diagram = inspect_workflow(workflow_spec, verbose=verbose)
        click.echo(diagram)
    except FileNotFoundError as e:
        raise click.ClickException(f"Workflow not found: {e}")
    except ValueError as e:
        raise click.ClickException(f"WDL parsing error: {e}")
    except Exception as e:
        raise click.ClickException(f"Error inspecting workflow: {e}")


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


def _input_path(reads_param) -> str:
    """Filesystem path a reads/contigs option refers to.

    `ReadsFile.convert` hands back an open (possibly recompressing) reader, so
    the raw option value stringifies to a repr, not a path. miniwdl needs the
    path on disk.
    """
    return str(getattr(reads_param, "_original_path", reads_param))


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

    if mean_length > 10000:
        sequence_type = 'long_reads'
        selected_strategy = 'sensitive'
        rationale = f'long reads detected (mean: {int(mean_length)}bp), selected sensitive strategy'
    elif mean_length > 1000:
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
@click.option("--quality-json", type=click.Path(exists=True), default=None, help="Quality JSON file for suspect data filtering")
@click.option("--include-suspect-alleles", "allele_filter", flag_value="include", default=True, help="Include suspect alleles (default)")
@click.option("--exclude-suspect-alleles", "allele_filter", flag_value="exclude", help="Exclude suspect alleles")
@click.option("--exclude-suspect-loci", is_flag=True, default=False, help="Exclude suspect loci")
@click.option("--exclude-suspect-profiles", is_flag=True, default=False, help="Exclude suspect profiles")
@click.argument('torch_args', nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def _run(clx, torch, cromwell_options="", method="main", workflow=None, output=None, strategy='balanced', contigs=None, reads=None, paired1=None, paired2=None, interlaced=None, longreads=None, quality_json=None, allele_filter="include", exclude_suspect_loci=False, exclude_suspect_profiles=False, torch_args=[]):
    "Run the selected torch."
    from torchbase.torchfs import Torch
    from torchbase.registry import RegistryManager
    from torchbase.config import RegistryConfig
    from torchbase import runner

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

        if user_specified_strategy and data_torch.typing_model == "operon":
            raise click.ClickException(
                "Cannot use --strategy with operon torches. "
                "The operon typing model has exactly one built-in workflow "
                "(operon_typing.wdl); --strategy is a speed/accuracy tier "
                "over the allelic typing model and does not apply."
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

        # A torch may supply its own workflow, or the user may point at one.
        # Those produce whatever their author decided, so they are run and
        # reported as-is. Everything else is a built-in typing model, which
        # this layer dispatches and then interprets.
        external_workflow = None
        if workflow:
            config = RegistryConfig.load()
            manager = RegistryManager(config)
            try:
                workflow_path = manager.fetch_torch(workflow)
                workflow_torch = Torch.load(workflow_path)
                external_workflow = workflow_torch.workflow
            except Exception as e:
                raise click.ClickException(f"Failed to fetch workflow {workflow}: {str(e)}")
            if not external_workflow:
                raise click.ClickException(f"Workflow torch {workflow} has no workflow file")
        elif data_torch.workflow:
            external_workflow = data_torch.workflow

        if external_workflow is None:
            from torchbase import typing_run

            try:
                if data_torch.typing_model == "operon":
                    if not contigs:
                        raise click.ClickException(
                            "Operon typing runs on an assembly; pass -c/--contigs."
                        )
                    result = typing_run.type_operon(data_torch, _input_path(contigs))
                else:
                    query = contigs or reads or paired1 or interlaced or longreads
                    input_type = (
                        "contigs" if contigs else
                        "reads" if reads else
                        "paired" if paired1 else
                        "interlaced" if interlaced else
                        "longreads"
                    )
                    result = typing_run.type_allelic(
                        data_torch,
                        _input_path(query),
                        input_type=input_type,
                        strategy=strategy,
                        quality_json=quality_json,
                        exclude_suspect_alleles=(allele_filter == "exclude"),
                        exclude_suspect_loci=exclude_suspect_loci,
                        exclude_suspect_profiles=exclude_suspect_profiles,
                    )
                    if auto_decision_rationale and isinstance(result.get("method"), dict):
                        result["method"]["auto_decision"] = auto_decision_rationale
            except (typing_run.TypingError, runner.WorkflowError) as e:
                raise click.ClickException(str(e))

            runner.emit(result, output)
            return result

        workflow_file = Path(external_workflow)
        try:
            allele_fasta_path, profiles_table_path = data_torch.get_unified_files()
        except Exception as e:
            raise click.ClickException(f"Failed to generate torch data files: {str(e)}")

        # Parse the workflow for its parameter schema and auto-provision the
        # torch's data files into whatever names it declares.
        if True:
            from torchbase.workflow_inspect import WDLParser
            from torchbase.workflow_params import (
                parse_torch_args,
                auto_provision_torch_parameters,
                validate_required_parameters
            )

            # Parse the embedded WDL to extract parameter schema
            try:
                with open(workflow_file, encoding="utf-8") as wf:
                    parser = WDLParser(wf.read(), wdl_dir=workflow_file.parent)
            except (FileNotFoundError, IOError) as e:
                raise click.ClickException(f"Failed to read workflow file {workflow_file}: {str(e)}")
            except ValueError as e:
                # WDL parsing failed - workflow may have invalid syntax
                # Fall back to non-parametric invocation for backward compatibility
                click.echo(f"Warning: Could not parse workflow parameters: {str(e)}", err=True)
                click.echo(f"Falling back to basic workflow invocation.", err=True)

                # Build basic miniwdl command without parameter validation
                miniwdl_cmd = ['miniwdl', 'run', str(workflow_file)]

                # Add all input files and torch args as-is
                if contigs:
                    miniwdl_cmd.append('contigs=' + _input_path(contigs))
                if reads:
                    miniwdl_cmd.append('reads=' + _input_path(reads))
                if paired1 and paired2:
                    miniwdl_cmd.extend(['paired1=' + _input_path(paired1),
                                        'paired2=' + _input_path(paired2)])
                if interlaced:
                    miniwdl_cmd.append('interlaced=' + _input_path(interlaced))
                if longreads:
                    miniwdl_cmd.append('longreads=' + _input_path(longreads))

                # Execute and return early
                result = run(miniwdl_cmd)
                if result.returncode != 0:
                    raise click.ClickException(f"Workflow execution failed with code {result.returncode}")
                return result

            # Auto-provision torch data files
            auto_params = auto_provision_torch_parameters(
                parser.workflow_inputs,
                allele_fasta_path,
                profiles_table_path
            )

            # Parse user-provided parameters from torch_args
            try:
                user_params = parse_torch_args(torch_args)
            except ValueError as e:
                raise click.ClickException(str(e))

            # Merge: user parameters override auto-provisioned
            all_params = {**auto_params, **user_params}

            # Add sequence input (user provided via -c/-r/etc)
            first_input = contigs or reads or paired1 or interlaced or longreads
            if first_input:
                all_params['query_sequences'] = _input_path(first_input)

            # Add optional parameters if provided via existing flags
            if quality_json:
                all_params['quality_json'] = str(quality_json)
            if allele_filter == "exclude":
                all_params['exclude_suspect_alleles'] = 'true'
            if exclude_suspect_loci:
                all_params['exclude_suspect_loci'] = 'true'
            if exclude_suspect_profiles:
                all_params['exclude_suspect_profiles'] = 'true'

            # Validate required parameters
            errors = validate_required_parameters(parser.workflow_inputs, all_params)
            if errors:
                raise click.ClickException(
                    f"Missing required workflow parameters:\n" + "\n".join(errors) +
                    f"\n\nRun 'torchbase workflow inspect {torch}' to see available parameters."
                )

            # Build miniwdl command with all parameters
            miniwdl_cmd = ['miniwdl', 'run', str(workflow_file)]
            for key, value in all_params.items():
                miniwdl_cmd.append(f'{key}={value}')

        # Execute a torch-supplied workflow. Its output shape is the torch
        # author's business, so it is reported as-is rather than interpreted.
        result = run(miniwdl_cmd)

        if result.returncode != 0:
            raise click.ClickException(f"Workflow execution failed with code {result.returncode}")

        return result

    except click.ClickException:
        raise
    except Exception as e:
        raise click.ClickException(f"Error running workflow: {str(e)}")
    finally:
        # Clean up temporary torch files
        for temp_path in (
            locals().get('allele_fasta_path'),
            locals().get('profiles_table_path'),
        ):
            if temp_path is not None and temp_path.exists():
                temp_path.unlink()






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

@tools.command("compress")
@click.argument("torch_path", type=click.Path(exists=True, dir_okay=True))
@click.option("--level", default=3, type=click.IntRange(1, 22), show_default=True,
              help="Zstandard compression level (higher = better compression, slower)")
@click.option("--keep-original", is_flag=True,
              help="Keep original .fasta files after compressing")
def compress_torch_alleles(torch_path, level, keep_original):
    """Compress FASTA allele files in torch to .fasta.zst format.

    Compresses all .fasta files in _resources/ or schemes/*/alleles/.
    Skips already-compressed files. Idempotent.
    """
    torch_path = Path(torch_path)
    cctx = zstd.ZstdCompressor(level=level)
    compressed_count = 0

    # Find all FASTA files
    fasta_files = []

    # Single-scheme format
    resources = torch_path / "_resources"
    if resources.exists():
        fasta_files.extend(resources.glob("*.fasta"))

    # Multi-scheme format
    schemes_dir = torch_path / "schemes"
    if schemes_dir.exists():
        fasta_files.extend(schemes_dir.glob("*/alleles/*.fasta"))

    for fasta_file in fasta_files:
        output_file = fasta_file.with_suffix(fasta_file.suffix + '.zst')

        if output_file.exists():
            click.echo(f"Skip {fasta_file.name} (already compressed)")
            continue

        with open(fasta_file, 'rb') as in_f:
            with open(output_file, 'wb') as out_f:
                cctx.copy_stream(in_f, out_f)

        if not keep_original:
            fasta_file.unlink()

        compressed_count += 1
        click.echo(f"Compressed {fasta_file.name} → {output_file.name}")

    click.echo(f"\nCompressed {compressed_count} files at level {level}")


@tools.command("call", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.argument("schema", type=click.File(), nargs=1)
@click.option("-j", "--json-profile", help="combined allele call in JSON format", nargs=1, default=None)
@click.pass_context
def call(ctx, schema, json_profile=None):
    "Load a profile definition and make a profile call from allele calls"
    with open(schema, newline="", encoding="utf-8") as schema_file:
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
@click.option("--bump", type=click.Choice(["patch", "minor", "major"]), default="patch", show_default=True,
              help="Version component to increment.")
@click.argument("checkpoint", required=False)
def _version(torch, bump="patch", checkpoint=None):
    "Bump the version of a currently-built torch."
    import toml
    torch_path = Path(torch)
    metadata_path = torch_path / "metadata.toml"
    if not metadata_path.exists():
        raise click.ClickException(f"metadata.toml not found at {torch_path}")

    with open(metadata_path, encoding="utf-8") as f:
        metadata = toml.load(f)

    version_str = str(metadata.get("version", "1.0.0"))
    try:
        parts = [int(x) for x in version_str.split(".")]
        while len(parts) < 3:
            parts.append(0)
    except ValueError:
        raise click.ClickException(f"Cannot parse version: {version_str}")

    if bump == "major":
        parts = [parts[0] + 1, 0, 0]
    elif bump == "minor":
        parts = [parts[0], parts[1] + 1, 0]
    else:
        parts = [parts[0], parts[1], parts[2] + 1]

    new_version = ".".join(str(p) for p in parts)
    metadata["version"] = new_version

    with open(metadata_path, "w", encoding="utf-8") as f:
        toml.dump(metadata, f)

    click.echo(f"Bumped version {version_str} -> {new_version}")


@tools.command("build")
@torch
def _build(torch):
    "Validate a torch's database structure."
    from torchbase.torchfs import Torch

    torch_path = Path(torch)
    try:
        result = Torch.load(torch_path)
    except (FileNotFoundError, ValueError) as e:
        raise click.ClickException(f"Validation failed: {e}")

    rows = [
        ("path", str(result.path)),
        ("schemes", len(result.schemes) if result.schemes else 1),
        ("references", len(result.references) if result.references else sum(
            len(v) for v in result.scheme_references.values())),
        ("workflow", str(result.workflow) if result.workflow else "built-in"),
    ]
    click.echo(tabulate(rows, headers=["Field", "Value"], tablefmt="simple"))
    click.echo("Torch structure is valid.")

@tools.group("convert")
@click.option("-v", "--verbose", count=True,
              help="Increase output verbosity. Use -v for progress, -vv for per-locus detail, -vvv for diagnostic trace.")
@click.pass_context
def convert(ctx, verbose):
    "Various conversion tools to make torches."
    ctx.ensure_object(dict)
    ctx.obj["verbosity"] = verbose

@convert.command("pubmlst")
@click.option("--url", default="https://rest.pubmlst.org", show_default=True,
              help="PubMLST/BIGSdb REST API root URL")
@click.option("--scheme-id", multiple=True, type=int,
              help="Scheme ID to include; repeat for multi-scheme torches. "
                   "Mutually exclusive with --all.")
@click.option("--all", "fetch_all", is_flag=True, default=False,
              help="Fetch every scheme from every database and produce one "
                   "'pubmlst' torch.  Mutually exclusive with --scheme-id.")
@click.option("--database", "database_names", multiple=True,
              help="With --all: restrict to these seqdef database names "
                   "(e.g. pubmlst_salmonella_seqdef). Repeat to include "
                   "several. Unset fetches every database.")
@click.option("--scheme-description", default=None, metavar="REGEX",
              help="With --all: only include schemes whose description "
                   "matches this regex (case-insensitive, anchored at the "
                   "start). e.g. '^MLST\\\\b' for classic 7-gene-style MLST, "
                   "skipping a database's cgMLST/wgMLST schemes.")
@click.option("--output", required=True, help="Output directory for torch")
@click.option("--name", default=None,
              help="Torch name (default: derived from first scheme, or 'pubmlst' with --all)")
@click.option("--namespace", default="pubmlst", show_default=True, help="Torch namespace")
@click.option("--cutoff-date", default="2024-12-31", show_default=True,
              help="Exclude alleles entered after this date (YYYY-MM-DD). "
                   "Defaults to 2024-12-31, the last day of freely-redistributable "
                   "PubMLST data before the 2025 licensing change.")
@click.option("--no-skip-errors", is_flag=True, default=False,
              help="With --all: abort on the first scheme that fails instead of skipping it.")
@click.option("--kmer-size", default=13, type=int, help="K-mer size for quality analysis")
@click.option("--overlap-threshold", default=0.90, type=float, help="Overlap threshold for quality analysis")
@click.option("--duplicate-threshold", default=0.95, type=float, help="Duplicate threshold for quality analysis")
@click.option("--max-alleles-for-quality", default=800, type=int, show_default=True,
              help="Skip pairwise similarity analysis for a locus with more "
                   "alleles than this. That analysis is O(n^2); PubMLST loci "
                   "commonly hold thousands of alleles, where it can run for "
                   "hours. Skipped loci are still fully present in the torch, "
                   "just without a quality.json suspect-pair analysis.")
@click.option("--ca-bundle", default=None, metavar="PATH",
              help="Path to a CA certificate bundle file (useful when a VPN or corporate proxy "
                   "performs SSL inspection). The REQUESTS_CA_BUNDLE environment variable is "
                   "also respected.")
@click.option("--no-ssl-verify", is_flag=True, default=False,
              help="Disable SSL certificate verification. Use only as a last resort in "
                   "environments where certificate inspection cannot be bypassed; prefer "
                   "--ca-bundle instead.")
@click.pass_context
def _convert_pubmlst(ctx, url, scheme_id, fetch_all, database_names, scheme_description,
                     output, name, namespace, cutoff_date,
                     no_skip_errors, kmer_size, overlap_threshold, duplicate_threshold,
                     max_alleles_for_quality, ca_bundle, no_ssl_verify):
    """Convert PubMLST schemes into a multi-scheme torch.

    Pass --scheme-id multiple times to bundle specific schemes (e.g., MLST and
    cgMLST for one organism) into a single torch.

    Pass --all to enumerate schemes across PubMLST databases into one merged
    torch. Unrestricted, that is every scheme of every database (tens of
    gigabytes); scope it with --database (repeatable) and/or
    --scheme-description (e.g. '^MLST\\b' for classic 7-gene MLST, excluding
    a database's larger cgMLST/wgMLST schemes) -- e.g. the foodborne-pathogen
    "pubmlst/mlst" torch this project ships was built with a --database list
    and --scheme-description '^MLST\\b'.

    The default cutoff date (2024-12-31) restricts alleles to those with
    freely-redistributable terms, matching the dataset snapshot bundled in
    tseemann/mlst v2.33.0.
    """
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    from datetime import date as _date
    from torchbase.conversions.pubmlst import convert_schemes, convert_all

    if fetch_all and scheme_id:
        raise click.UsageError("--all and --scheme-id are mutually exclusive")
    if not fetch_all and not scheme_id:
        raise click.UsageError("Specify at least one --scheme-id or use --all")

    try:
        cutoff = _date.fromisoformat(cutoff_date)
    except ValueError:
        raise click.ClickException(f"Invalid cutoff-date '{cutoff_date}': expected YYYY-MM-DD")

    if no_ssl_verify:
        verify = False
        click.echo(
            "WARNING: SSL certificate verification disabled. "
            "This is insecure; use --ca-bundle if possible.",
            err=True,
        )
    elif ca_bundle:
        verify = ca_bundle
    else:
        verify = True

    try:
        if fetch_all:
            torch_path = convert_all(
                base_url=url,
                output_path=output,
                namespace=namespace,
                torch_name=name or "pubmlst",
                kmer_size=kmer_size,
                overlap_threshold=overlap_threshold,
                duplicate_threshold=duplicate_threshold,
                max_alleles_for_quality=max_alleles_for_quality,
                cutoff_date=cutoff,
                skip_errors=not no_skip_errors,
                verify=verify,
                database_names=list(database_names) or None,
                scheme_description_pattern=scheme_description,
            )
        else:
            torch_path = convert_schemes(
                database_url=url,
                scheme_ids=list(scheme_id),
                output_path=output,
                namespace=namespace,
                torch_name=name,
                kmer_size=kmer_size,
                overlap_threshold=overlap_threshold,
                duplicate_threshold=duplicate_threshold,
                max_alleles_for_quality=max_alleles_for_quality,
                cutoff_date=cutoff,
                verify=verify,
            )
        click.echo(f"Successfully created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")

@convert.command("pubcgmlst")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--namespace", default="pubcgmlst", show_default=True)
@click.option("--name", default=None, help="Torch name (default: scheme filename stem)")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _pubcgmlst(ctx, scheme, sequences, output, namespace, name, version, kmer_size, overlap_threshold, duplicate_threshold):
    "Create a torch from a PubMLST cgMLST database and schema."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    from torchbase.conversions.pubcgmlst import convert_local

    try:
        torch_path = convert_local(
            scheme_file=scheme,
            sequence_files=list(sequences),
            output_path=output,
            namespace=namespace,
            name=name,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")

@convert.command("seqsero2")
@click.argument("sequences", type=click.File(), nargs=-1)
@click.option("--profiles", type=click.File(), default=None, help="Serotype definitions TSV (Serotype, O, H1, H2)")
@click.option("--download", is_flag=True, default=False, help="Download FASTA files from denglab/SeqSero2")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--name", default="seqsero2", show_default=True, help="Torch name")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _seqsero2(ctx, sequences, profiles, output, name, version, kmer_size, overlap_threshold, duplicate_threshold, download):
    "Create a torch from SeqSero2 Salmonella serotyping database files."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.seqsero2 import convert_local, DOWNLOAD_SOURCES

    if download:
        if sequences:
            raise click.UsageError("Cannot use --download with explicit SEQUENCES arguments.")
        from torchbase.conversions.seqsero2 import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_seqsero2_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        sequences = sources["sequences"]

    try:
        torch_path = convert_local(
            sequence_files=list(sequences),
            profiles_file=profiles,
            output_path=output,
            name=name,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")


@convert.command("seqsero2s")
@click.argument("sequences", type=click.File(), nargs=-1,
                metavar="[MLST_FASTAS]...",)
@click.option("--antigen-db", type=click.File(), default=None,
              help="Combined antigen FASTA (H_and_O_and_specific_genes.fasta)")
@click.option("--mlst-profiles", type=click.File(), default=None,
              help="MLST allele-to-ST table (salmonella_profile.txt)")
@click.option("--serotype-profiles", type=click.File(), default=None,
              help="Serotype definitions TSV (Serotype, O, H1, H2)")
@click.option("--download", is_flag=True, default=False,
              help="Download antigen DB, MLST locus FASTAs, and MLST profiles from LSTUGA/SeqSero2S")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--name", default="seqsero2s", show_default=True, help="Torch name")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _seqsero2s(ctx, sequences, antigen_db, mlst_profiles, serotype_profiles,
               output, name, version, kmer_size, overlap_threshold, duplicate_threshold, download):
    "Create a torch from SeqSero2S (LSTUGA) Salmonella serotyping + MLST database files."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.seqsero2s import convert_local, DOWNLOAD_SOURCES

    if download:
        if sequences or antigen_db or mlst_profiles:
            raise click.UsageError(
                "Cannot use --download with explicit MLST_FASTAS, --antigen-db, or --mlst-profiles."
            )
        from torchbase.conversions.seqsero2s import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_seqsero2s_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        sequences = sources["sequences"]
        antigen_db = sources["antigen_db"]
        mlst_profiles = sources["mlst_profiles"]

    try:
        torch_path = convert_local(
            sequence_files=list(sequences),
            antigen_db=antigen_db,
            mlst_profiles=mlst_profiles,
            serotype_profiles=serotype_profiles,
            output_path=output,
            name=name,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")


@convert.command("ectyper")
@click.argument("sequences", type=click.File(), nargs=-1)
@click.option("--profiles", type=click.File(), default=None, help="Serotype definitions TSV (Serotype, O, H)")
@click.option("--db", "db_fasta", type=click.File(), default=None,
              help="Combined ECTyper database FASTA (split automatically into O/H files)")
@click.option("--download", is_flag=True, default=False,
              help="Download ECTyperDB.fasta and allele profiles from phac-nml/ectyper")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--name", default="ectyper", show_default=True, help="Torch name")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _ectyper(ctx, sequences, profiles, db_fasta, output, name, version, kmer_size, overlap_threshold, duplicate_threshold, download):
    "Create a torch from ECTyper E. coli / Shigella serotyping database files."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.ectyper import convert_local, DOWNLOAD_SOURCES

    if download:
        if sequences or db_fasta:
            raise click.UsageError("Cannot use --download with explicit SEQUENCES or --db arguments.")
        from torchbase.conversions.ectyper import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_ectyper_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        db_fasta = sources["db_fasta"]
        if profiles is None:
            profiles = sources["profiles"]
    elif not sequences and db_fasta is None:
        raise click.UsageError("Provide SEQUENCES files, --db, or --download.")

    try:
        torch_path = convert_local(
            sequence_files=list(sequences),
            profiles_file=profiles,
            db_fasta=db_fasta,
            output_path=output,
            name=name,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")


@convert.command("lissero")
@click.argument("sequences", type=click.File(), nargs=-1)
@click.option("--profiles", type=click.File(), default=None,
              help="Serogroup definitions TSV (Serogroup, prs, LMOSA, LMOSB, ORF2110, ORF2819, ldh, lin0764, lin1118)")
@click.option("--download", is_flag=True, default=False,
              help="Download locus FASTAs and serogroup profiles from MDU-PHL/LisSero")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--name", default="lissero", show_default=True, help="Torch name")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _lissero(ctx, sequences, profiles, output, name, version, kmer_size, overlap_threshold, duplicate_threshold, download):
    "Create a torch from LisSero Listeria monocytogenes serogroup database files."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.lissero import convert_local, DOWNLOAD_SOURCES

    if download:
        if sequences:
            raise click.UsageError("Cannot use --download with explicit SEQUENCES arguments.")
        from torchbase.conversions.lissero import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_lissero_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        sequences = sources["sequences"]
        if profiles is None:
            profiles = sources["profiles"]

    try:
        torch_path = convert_local(
            sequence_files=list(sequences),
            profiles_file=profiles,
            output_path=output,
            name=name,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")


@convert.command("chewie-ns")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
def _chewie_ns():
    "Create a torch from a Chewie-NS wgMLST database and schema."
    pass

@convert.command("shigatyper")
@click.argument("sequences", type=click.File(), nargs=-1)
@click.option("--profiles", type=click.File(), default=None, help="Serotype profiles TSV")
@click.option("--download", is_flag=True, default=False,
              help="Download FASTA files from CFSAN-Biostatistics/ShigaTyper")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--name", default="shigatyper", show_default=True, help="Torch name")
@click.option("--namespace", default="hfp", show_default=True,
              help="Torch namespace. Names the authority for the data: 'hfp' "
                   "because ShigaTyper is an FDA Human Foods Program product.")
@click.option("--version", default="1.0.0", show_default=True, help="Torch version")
@click.option("--kmer-size", default=13, type=int, show_default=True)
@click.option("--overlap-threshold", default=0.90, type=float, show_default=True)
@click.option("--duplicate-threshold", default=0.95, type=float, show_default=True)
@click.pass_context
def _shigatyper(ctx, sequences, profiles, output, name, namespace, version, kmer_size, overlap_threshold, duplicate_threshold, download):
    "Create a torch from ShigaTyper's database."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.shigatyper import convert_local, DOWNLOAD_SOURCES

    if download:
        if sequences:
            raise click.UsageError("Cannot use --download with explicit SEQUENCES arguments.")
        from torchbase.conversions.shigatyper import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_shigatyper_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        sequences = sources["sequences"]

    try:
        torch_path = convert_local(
            sequence_files=list(sequences),
            profiles_file=profiles,
            output_path=output,
            name=name,
            namespace=namespace,
            version=version,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")


@convert.command("stxtyper")
@click.argument("stx_prot", type=click.File(), required=False)
@click.option("--download", is_flag=True, default=False,
              help="Download stx.prot + version.txt from ncbi/stxtyper")
@click.option("--stxtyper-version", default=None,
              help="Upstream StxTyper version this torch claims parity against")
@click.option("--output", default=".", show_default=True, help="Output directory")
@click.option("--namespace", default="ncbi", show_default=True, help="Torch namespace")
@click.option("--name", default="stxtyper", show_default=True, help="Torch name")
@click.option("--version", default=None, help="Torch version (default: stxtyper version)")
@click.pass_context
def _stxtyper(ctx, stx_prot, download, stxtyper_version, output, namespace, name, version):
    "Create an operon torch from NCBI StxTyper's stx.prot reference set."
    from torchbase.conversions.log import setup_logging
    setup_logging(ctx.obj.get("verbosity", 0))
    import tempfile
    from torchbase.conversions.stxtyper import convert_local, DOWNLOAD_SOURCES

    if download:
        if stx_prot:
            raise click.UsageError("Cannot use --download with an explicit STX_PROT argument.")
        from torchbase.conversions.stxtyper import download_sources
        dest = Path(tempfile.mkdtemp(prefix="torchtools_stxtyper_"))
        click.echo(f"Downloading from {DOWNLOAD_SOURCES['repo']} → {dest}")
        sources = download_sources(dest)
        stx_prot = sources["stx_prot"]
        stxtyper_version = stxtyper_version or sources["stxtyper_version"]
    elif not stx_prot:
        raise click.UsageError("Provide STX_PROT or pass --download.")

    try:
        torch_path = convert_local(
            stx_prot_file=stx_prot,
            output_path=output,
            namespace=namespace,
            name=name,
            version=version,
            stxtyper_version=stxtyper_version,
        )
        click.echo(f"Created torch at: {torch_path}")
    except Exception as e:
        raise click.ClickException(f"Conversion failed: {str(e)}")




def _yubikey_signer_options(func):
    """Decorator that adds --yubikey/--slot/--pin options."""
    func = click.option("--pin", default=None, help="YubiKey PIV PIN (prompted if omitted)")(func)
    func = click.option("--slot", default="9c", type=click.Choice(["9a", "9c", "9d", "9e"]),
                        show_default=True, help="YubiKey PIV slot")(func)
    func = click.option("--yubikey", is_flag=True, default=False,
                        help="Use YubiKey for signing/key storage")(func)
    return func


def _build_signer(namespace, yubikey=False, slot="9c", pin=None):
    """Return a signer instance for the given backend."""
    from torchbase.signing import FileKeySigner, YubiKeySigner
    if yubikey:
        return YubiKeySigner(slot=slot, pin=pin)
    key_path = Path.home() / ".torchbase" / "keys" / f"{namespace}.key"
    if not key_path.exists():
        raise click.ClickException(
            f"No key found at {key_path}. Run: torchtools keygen --namespace {namespace}"
        )
    return FileKeySigner(key_path)


@tools.command("keygen")
@click.option("--namespace", required=True, help="Namespace to generate a key for")
@click.option("--output", default=None, help="Key output directory (default: ~/.torchbase/keys/)")
@_yubikey_signer_options
def _keygen(namespace, output, yubikey, slot, pin):
    """Generate a signing key for a namespace."""
    from torchbase.signing import generate_software_keypair, setup_yubikey_slot

    key_dir = Path(output) if output else Path.home() / ".torchbase" / "keys"

    if yubikey:
        pub_path = setup_yubikey_slot(slot, namespace, key_dir, pin=pin)
        click.echo(f"YubiKey key generated in slot {slot}")
        click.echo(f"Public key written to: {pub_path}")
    else:
        priv_path, pub_path = generate_software_keypair(namespace, key_dir)
        click.echo(f"Private key: {priv_path}")
        click.echo(f"Public key:  {pub_path}")

    pub_text = pub_path.read_text(encoding="utf-8").strip()
    click.echo(f"\nPublic key (add to key registry):\n{namespace} = \"{pub_text}\"")


@tools.command("pubkey")
@click.argument("namespace", required=True)
@click.option("--key-dir", default=None, help="Key directory (default: ~/.torchbase/keys/)")
def _pubkey(namespace, key_dir):
    """Print the public key for a namespace."""
    key_dir = Path(key_dir) if key_dir else Path.home() / ".torchbase" / "keys"
    pub_path = key_dir / f"{namespace}.pub"
    if not pub_path.exists():
        raise click.ClickException(f"No public key found at {pub_path}")
    click.echo(pub_path.read_text(encoding="utf-8").strip())


@tools.command("sign")
@torch
@_yubikey_signer_options
def _sign(torch, yubikey, slot, pin):
    """Sign a torch, writing signature.toml."""
    from torchbase.signing import sign_torch
    from torchbase import torchfs

    t = torchfs.Torch.load(torch)
    meta_path = Path(torch) / "metadata.toml"
    import toml as _toml
    meta = _toml.load(meta_path)
    namespace = meta["namespace"]

    signer = _build_signer(namespace, yubikey=yubikey, slot=slot, pin=pin)
    sig_path = sign_torch(Path(torch), signer)
    click.echo(f"Signed: {sig_path}")


@tools.command("publish")
@torch
@_yubikey_signer_options
@click.pass_context
def _publish(ctx, torch, yubikey, slot, pin):
    """Deprecated: use 'torchtools manifest add' instead."""
    click.echo(
        "Warning: 'torchtools publish' is deprecated and will be removed in a future release.\n"
        "Use 'torchtools manifest add <torch-path>' instead.",
        err=True,
    )
    ctx.invoke(_manifest_add, torch=torch, yubikey=yubikey, slot=slot, pin=pin)


# ---------------------------------------------------------------------------
# torchtools namespace  (register / show)
# ---------------------------------------------------------------------------

@tools.group("namespace")
def _namespace_group():
    """Namespace registration and inspection commands."""
    pass


@_namespace_group.command("register")
@click.option("--namespace", required=True, help="Namespace to register (e.g. 'us-fda-hfp')")
@click.option("--ipfs-node", default="127.0.0.1", show_default=True)
@click.option("--ipfs-port", default=5001, show_default=True, type=int)
@click.option("--submit-to", multiple=True, metavar="URL",
              help="Log operator URL(s) to submit the genesis block to")
@_yubikey_signer_options
def _namespace_register(namespace, ipfs_node, ipfs_port, submit_to, yubikey, slot, pin):
    """Register a namespace by creating and publishing a genesis block.

    Steps:
      1. Import your Ed25519 key into Kubo's keystore
      2. Create and sign the genesis block
      3. Upload the genesis block to IPFS and pin it
      4. Publish the genesis block CID to IPNS (via your namespace key)
      5. Submit to any --submit-to log operators
      6. Write the IPNS address to ~/.torchbase/config.toml
    """
    import shutil
    import tempfile
    import toml as _toml
    from torchbase.chain import (
        make_genesis_block, block_to_toml, upload_block, pin_cid,
        publish_chain_head, import_key_to_kubo, submit_to_log,
    )

    signer = _build_signer(namespace, yubikey=yubikey, slot=slot, pin=pin)

    # Import key into Kubo keystore (idempotent if already present)
    key_path = Path.home() / ".torchbase" / "keys" / f"{namespace}.key"
    if not yubikey and key_path.exists():
        click.echo(f"Importing key into Kubo keystore as '{namespace}' ...")
        try:
            ipns_address = import_key_to_kubo(namespace, key_path, ipfs_node, ipfs_port)
        except Exception as e:
            raise click.ClickException(f"Key import failed: {e}")
    else:
        # YubiKey or no PEM — derive IPNS address from Kubo if key already there
        # (user must have imported manually); we'll discover it after publishing
        ipns_address = None

    # Build and upload genesis block
    click.echo("Creating genesis block ...")
    genesis = make_genesis_block(namespace, signer)
    genesis_toml = block_to_toml(genesis)

    try:
        genesis_cid = upload_block(genesis_toml, ipfs_node, ipfs_port)
        pin_cid(genesis_cid, ipfs_node, ipfs_port)
    except Exception as e:
        raise click.ClickException(f"IPFS upload failed: {e}")

    # Publish IPNS head
    click.echo(f"Publishing genesis block {genesis_cid} to IPNS ...")
    try:
        publish_chain_head(genesis_cid, namespace, ipfs_node, ipfs_port)
    except Exception as e:
        raise click.ClickException(f"IPNS publish failed: {e}")

    # Discover IPNS address if we didn't get it from import
    if ipns_address is None:
        from torchbase.chain import get_chain_head, _kubo_url
        try:
            resp = requests.post(
                f"{_kubo_url(ipfs_node, ipfs_port)}/key/list",
                timeout=10,
            )
            resp.raise_for_status()
            for entry in resp.json().get("Keys", []):
                if entry.get("Name") == namespace:
                    ipns_address = "/ipns/" + entry["Id"]
                    break
        except Exception:
            pass
        if ipns_address is None:
            ipns_address = f"<run 'ipfs key list' to find IPNS address for '{namespace}'>"

    # Submit to log operators
    all_log_urls = list(submit_to)
    if not all_log_urls:
        from torchbase.config import RegistryConfig
        cfg = RegistryConfig.load()
        all_log_urls = list(cfg.log_operators)
    for log_url in all_log_urls:
        click.echo(f"Submitting to log operator: {log_url}")
        try:
            submit_to_log(genesis_cid, namespace, ipns_address, log_url)
        except Exception as e:
            click.echo(f"  Warning: log submission failed: {e}", err=True)

    # Write IPNS address to user config
    config_path = Path.home() / ".torchbase" / "config.toml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    existing = _toml.load(config_path) if config_path.exists() else {}
    if "namespaces" not in existing:
        existing["namespaces"] = {}
    existing["namespaces"][namespace] = ipns_address
    fd, tmp = tempfile.mkstemp(dir=config_path.parent, suffix=".toml.tmp")
    try:
        with open(fd, "w", encoding="utf-8") as f:
            _toml.dump(existing, f)
        shutil.move(tmp, config_path)
    except Exception:
        Path(tmp).unlink(missing_ok=True)
        raise

    click.echo(f"\nNamespace '{namespace}' registered.")
    click.echo(f"  Genesis CID:  {genesis_cid}")
    click.echo(f"  IPNS address: {ipns_address}")
    click.echo(f"\nConfig snippet for collaborators:")
    click.echo(f'  [namespaces]')
    click.echo(f'  "{namespace}" = "{ipns_address}"')


@_namespace_group.command("show")
@click.argument("namespace", required=True)
@click.option("--ipfs-node", default="127.0.0.1", show_default=True)
@click.option("--ipfs-port", default=5001, show_default=True, type=int)
def _namespace_show(namespace, ipfs_node, ipfs_port):
    """Walk the IPLD chain for a namespace and print its reconstructed manifest."""
    import toml as _toml
    from torchbase.chain import get_chain_head, walk_chain, reconstruct_manifest
    from torchbase.config import RegistryConfig

    config = RegistryConfig.load()
    ipns_address = config.namespaces.get(namespace)
    if not ipns_address:
        raise click.ClickException(
            f"Namespace '{namespace}' not in config. "
            f"Run: torchtools namespace register --namespace {namespace}"
        )

    head_cid = get_chain_head(ipns_address, ipfs_node, ipfs_port)
    if head_cid is None:
        raise click.ClickException(f"Could not resolve IPNS address {ipns_address}")

    click.echo(f"Head CID:     {head_cid}")

    try:
        chain = walk_chain(head_cid, ipfs_node, ipfs_port)
    except ValueError as e:
        raise click.ClickException(f"Chain verification failed: {e}")

    genesis_cid = None
    for block in chain:
        if block.get("type") == "genesis":
            genesis_cid = head_cid if len(chain) == 1 else "..."
            break

    click.echo(f"Block count:  {len(chain)}")
    click.echo(f"Genesis CID:  {chain[0].get('public_key', '')[:16]}... (public key prefix)")
    manifest = reconstruct_manifest(chain)
    if manifest:
        click.echo("\nReconstructed manifest:")
        click.echo(_toml.dumps(manifest))
    else:
        click.echo("\nNo torches published yet.")


# ---------------------------------------------------------------------------
# torchtools manifest  (add / show)
# ---------------------------------------------------------------------------

@tools.group("manifest")
def _manifest_group():
    """Manifest publishing commands (IPLD commit chain)."""
    pass


@_manifest_group.command("add")
@torch
@click.option("--ipfs-node", default="127.0.0.1", show_default=True)
@click.option("--ipfs-port", default=5001, show_default=True, type=int)
@click.option("--submit-to", multiple=True, metavar="URL",
              help="Extra log operator URL(s) to notify")
@_yubikey_signer_options
def _manifest_add(torch, ipfs_node, ipfs_port, submit_to, yubikey, slot, pin):
    """Sign a torch, upload it to IPFS, and publish an update block.

    Steps:
      1. Sign the torch (idempotent — rewrites signature.toml)
      2. Upload the torch directory to IPFS, pin the CID
      3. Sign the CID
      4. Resolve the namespace IPNS head (the previous block CID)
      5. Build and upload an update block referencing the previous block
      6. Publish the new block CID to IPNS
      7. Submit to configured log operators
    """
    import shutil
    import tempfile
    import toml as _toml
    from torchbase.signing import sign_torch, sign_cid
    from torchbase.chain import (
        make_update_block, block_to_toml, upload_block, pin_cid,
        upload_directory, get_chain_head, publish_chain_head,
        walk_chain, submit_to_log,
    )
    from torchbase.config import RegistryConfig

    torch_path = Path(torch)
    if not torch_path.is_dir():
        raise click.ClickException(f"Not a directory: {torch_path}")

    meta = _toml.load(torch_path / "metadata.toml")
    namespace = meta["namespace"]
    name = meta["name"]
    version = str(meta["version"])

    signer = _build_signer(namespace, yubikey=yubikey, slot=slot, pin=pin)

    # Step 1 — sign content
    click.echo("Signing torch ...")
    sig_path = sign_torch(torch_path, signer)
    click.echo(f"  {sig_path}")

    # Step 2 — upload to IPFS
    click.echo("Uploading to IPFS ...")
    try:
        torch_cid = upload_directory(torch_path, ipfs_node, ipfs_port)
        pin_cid(torch_cid, ipfs_node, ipfs_port)
    except Exception as e:
        raise click.ClickException(f"IPFS upload failed: {e}")
    click.echo(f"  torch CID: {torch_cid}")

    # Step 3 — sign the CID
    cid_sig = sign_cid(torch_cid, namespace, version, signer)

    # Step 4 — find previous block (chain head)
    config = RegistryConfig.load()
    ipns_address = config.namespaces.get(namespace)
    if not ipns_address:
        raise click.ClickException(
            f"Namespace '{namespace}' not found in config. "
            f"Run: torchtools namespace register --namespace {namespace}"
        )
    previous_cid = get_chain_head(ipns_address, ipfs_node, ipfs_port)
    if previous_cid is None:
        raise click.ClickException(
            f"Could not resolve chain head for {ipns_address}. "
            f"Run 'torchtools namespace register' first."
        )

    # Step 5 — build update block
    torch_ref = f"{namespace}/{name}"
    entries = {
        torch_ref: {
            version: torch_cid,
            "latest": torch_cid,
            "signatures": {version: cid_sig},
        }
    }
    update_block = make_update_block(namespace, previous_cid, entries, signer)
    block_toml_str = block_to_toml(update_block)

    click.echo("Uploading update block ...")
    try:
        block_cid = upload_block(block_toml_str, ipfs_node, ipfs_port)
        pin_cid(block_cid, ipfs_node, ipfs_port)
    except Exception as e:
        raise click.ClickException(f"Block upload failed: {e}")
    click.echo(f"  block CID: {block_cid}")

    # Step 6 — publish new chain head
    click.echo("Publishing IPNS head ...")
    try:
        publish_chain_head(block_cid, namespace, ipfs_node, ipfs_port)
    except Exception as e:
        raise click.ClickException(f"IPNS publish failed: {e}")

    # Step 7 — submit genesis CID to log operators
    try:
        chain = walk_chain(block_cid, ipfs_node, ipfs_port)
        genesis_block_idx = next(
            (i for i, b in enumerate(chain) if b.get("type") == "genesis"), None
        )
    except Exception:
        chain = []
        genesis_block_idx = None

    all_log_urls = list(submit_to) + list(config.log_operators)
    if all_log_urls and genesis_block_idx is not None:
        # Walk chain to get genesis CID
        from torchbase.chain import fetch_block, _kubo_url
        # Re-walk to find genesis CID (walk_chain returns blocks, not CIDs)
        cur = block_cid
        while True:
            blk = fetch_block(cur, ipfs_node, ipfs_port)
            if blk.get("type") == "genesis":
                genesis_cid_for_log = cur
                break
            cur = blk.get("previous", "")
            if not cur:
                genesis_cid_for_log = None
                break
        if genesis_cid_for_log:
            for log_url in all_log_urls:
                click.echo(f"Submitting to log: {log_url}")
                try:
                    submit_to_log(genesis_cid_for_log, namespace, ipns_address, log_url)
                except Exception as e:
                    click.echo(f"  Warning: log submission failed: {e}", err=True)

    click.echo(f"\nPublished {torch_ref} {version}")
    click.echo(f"  Torch CID: {torch_cid}")
    click.echo(f"  Block CID: {block_cid}")
    click.echo(f"  IPNS:      {ipns_address}")
    click.echo(f"\nManifest entry (for reference):")
    click.echo(f'["{torch_ref}"]')
    click.echo(f'"{version}" = "{torch_cid}"')
    click.echo(f'latest = "{torch_cid}"')
    click.echo(f'["{torch_ref}".signatures]')
    click.echo(f'"{version}" = "{cid_sig}"')


@_manifest_group.command("show")
@click.argument("namespace", required=True)
@click.option("--ipfs-node", default="127.0.0.1", show_default=True)
@click.option("--ipfs-port", default=5001, show_default=True, type=int)
def _manifest_show(namespace, ipfs_node, ipfs_port):
    """Walk the IPLD chain and print the reconstructed manifest as TOML."""
    import toml as _toml
    from torchbase.chain import get_chain_head, walk_chain, reconstruct_manifest
    from torchbase.config import RegistryConfig

    config = RegistryConfig.load()
    ipns_address = config.namespaces.get(namespace)
    if not ipns_address:
        raise click.ClickException(
            f"Namespace '{namespace}' not in config."
        )

    head_cid = get_chain_head(ipns_address, ipfs_node, ipfs_port)
    if head_cid is None:
        raise click.ClickException(f"Could not resolve IPNS for namespace '{namespace}'")

    try:
        chain = walk_chain(head_cid, ipfs_node, ipfs_port)
    except ValueError as e:
        raise click.ClickException(f"Chain verification failed: {e}")

    manifest = reconstruct_manifest(chain)
    if manifest:
        click.echo(_toml.dumps(manifest))
    else:
        click.echo("# No torches published yet")


@cli.command("verify")
@torch
@click.option("--public-key", "public_key_b64", default=None,
              help="Override public key (base64url). Default: resolve from config/embedded.")
@click.option("--require-signature", is_flag=True, default=False,
              help="Fail if no signature is present")
def _verify(torch, public_key_b64, require_signature):
    """Verify the cryptographic signature of a torch."""
    from torchbase.signing import verify_torch, resolve_public_key
    from torchbase.config import RegistryConfig

    torch_path = Path(torch)
    sig_path = torch_path / "signature.toml"

    if not sig_path.exists():
        if require_signature:
            raise click.ClickException("No signature.toml found and --require-signature set")
        click.echo("Warning: no signature.toml found", err=True)
        return

    config = RegistryConfig.load()

    # Resolve public key if not provided
    if public_key_b64 is None:
        import toml as _toml
        meta = _toml.load(torch_path / "metadata.toml")
        namespace = meta["namespace"]
        result = resolve_public_key(namespace, config, torch_path)
        if result is not None:
            public_key_b64, _ = result

    result = verify_torch(torch_path, public_key_b64=public_key_b64)
    if result.valid:
        click.echo(f"OK: {result.message}")
    else:
        raise click.ClickException(f"Invalid signature: {result.message}")


if __name__ == '__main__':
    cli()

