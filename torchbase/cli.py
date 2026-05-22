import csv
import click
import logging
import json
from functools import partial
from pathlib import Path
from xml.etree.ElementTree import ElementTree as xml
from tabulate import tabulate
from subprocess import run

import zstandard as zstd
import zipfile
import gzip
import bz2

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

class ReadsFile(click.Path):
    name = "reads or contigs file"

    def __init__(self):
        super().__init__(exists=True, dir_okay=False, readable=True, resolve_path=True, allow_dash=True, path_type=Path)

    def convert(self, value, param, ctx):
        path = super().convert(value, param, ctx)
        # open the file, try to decompress it, and compress to zstandard
        compressor = zstd.ZstdCompressor()

        def compress_stream(file_obj):
            return compressor.stream_reader(file_obj)

        magic_sigs = (
            (0x1f8b08, gzip.open, compress_stream),
            (0x425a68, bz2.open, compress_stream),
            (0x504b0304, lambda p, m: zipfile.ZipFile(p, m), compress_stream),
            (0x28b52ffd, open, lambda s: s) # zstd doesn't need to be converted
        )

        for signature, method, converter in magic_sigs:
            with open(path, 'rb') as file:
                if file.read(4).startswith(signature):
                    return converter(method(path, 'rb'))
        # otherwise the file is not compressed, compress it
        return compress_stream(open(path, 'rb'))



# We use this a lot

ReadsParam = partial(click.option, 
                     nargs=1, 
                     default="", 
                     type=ReadsFile())

#
# Main running method
#

@cli.command("run", context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option("--cromwell-opts", "cromwell_options", nargs=1, default="", type=click.STRING)
@torch
@click.option("-m" "--method", nargs=1, default="main", type=click.STRING)
@click.option("--workflow", default=None, help="Override workflow torch (namespace/name format)")
@click.option("-o", "--output", default=None, help="Output file for results")
@ReadsParam("-c", "--contigs")
@ReadsParam("-r", "--reads")
@ReadsParam("-pe1", "--paired1", "--pe1")
@ReadsParam("-pe2", "--paired2", "--pe2")
@ReadsParam("-i", "--interlaced")
@ReadsParam("-l", "--longreads")
@click.argument('torch_args', nargs=-1, type=click.UNPROCESSED)
@click.pass_context
def _run(clx, torch, cromwell_options="", method="main", workflow=None, output=None, contigs=None, reads=None, paired1=None, paired2=None, interlaced=None, longreads=None, torch_args=[]):
    "Run the selected torch."
    from torchbase.torchfs import Torch
    from torchbase.registry import RegistryManager
    from torchbase.config import RegistryConfig

    if not (contigs or reads or (paired1 and paired2) or interlaced or longreads):
        if (paired1 and not paired2) or (paired2 and not paired1):
            raise click.Abort("paired-end data requires two files; use -i/--interlaced for single-file paired-end data.")
        raise click.Abort("at least one reads option and file must be given.")
    if len(filter(lambda v: v is not None, (contigs, reads, paired1, interlaced, longreads))) > 1:
        raise click.Abort("provide reads in no more than one layout form.")

    try:
        # Load data torch
        data_torch = Torch.load(torch)

        # Determine workflow to use
        workflow_torch = data_torch

        if workflow:
            # User specified custom workflow
            config = RegistryConfig.load()
            manager = RegistryManager(config)
            try:
                workflow_path = manager.fetch_torch(workflow)
                workflow_torch = Torch.load(workflow_path)
            except Exception as e:
                raise click.ClickException(f"Failed to fetch workflow {workflow}: {str(e)}")
        elif not data_torch.workflow:
            # No workflow in data torch, try default
            try:
                config = RegistryConfig.load()
                manager = RegistryManager(config)
                default_workflow_path = manager.fetch_torch("torchbase/default-workflow")
                workflow_torch = Torch.load(default_workflow_path)
            except Exception as e:
                raise click.ClickException(
                    f"Workflow not found in torch and default workflow fetch failed: {str(e)}"
                )

        # Validate workflow exists and is named main.wdl
        if not workflow_torch.workflow:
            raise click.ClickException("No workflow found (main.wdl) in torch")

        if workflow_torch.workflow.name != "main.wdl":
            raise click.ClickException(
                f"Workflow must be named 'main.wdl', found: {workflow_torch.workflow.name}"
            )

        # Build miniwdl command
        miniwdl_cmd = ['miniwdl', 'run', str(workflow_torch.workflow)]

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

@tools.group("convert")
def convert():
    "Various conversion tools to make torches."
    pass

@convert.command("pubmlst")
@click.argument("scheme", type=click.File())
@click.argument("sequences", type=click.File(), nargs=-1)
def _pubmlst(scheme, sequences=[]):
    "Create a torch from a PubMLST database and schema."
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

