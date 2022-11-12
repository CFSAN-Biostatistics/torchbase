import click
import logging


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

@cli.command("list")
@click.option('-i', '--installed', 'only_installed', default=False)
def _list(only_installed=False):
    "Show available typing frameworks."
    pass

@cli.command("pull")
@torch
@click.option("--force-use-gateway", default=False)
def _pull(torch, force_use_gateway=False):
    "Pull the selected torch via IPFS or an IPFS gateway."
    pass

@cli.command("run", context_settings=dict(ignore_unknown_options=True))
@click.option('-c', "--cromwell-opts", "cromwell_options", nargs=1, default="", type=click.STRING)
@torch
@click.argument('torch_args', nargs=-1, type=click.UNPROCESSED)
def _run(torch, cromwell_options="", torch_args=[]):
    "Run the selected torch."
    pass




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
def _pubcgmlst():
    "Create a torch from a PubMLST cgMLST database and schema."
    pass

@convert.command("chewie-ns")
def _chewie_ns():
    "Create a torch from a Chewie-NS wgMLST database and schema."
    pass





if __name__ == '__main__':
    cli()

