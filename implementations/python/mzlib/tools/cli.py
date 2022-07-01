import logging
import ipdb
import sys
import traceback
import click
import logging

from typing import DefaultDict, List

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.index import MemoryIndex, SQLIndex
from mzlib.backends.text import TextSpectralLibraryWriter
from mzlib.validate import validator
from mzlib.validate.level import RequirementLevel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

logger = logging.getLogger(__name__)


def info(type, value, tb):
    if not sys.stderr.isatty():
        click.secho("Running interactively, not starting debugger", fg='yellow')
        sys.__excepthook__(type, value, tb)
    else:
        traceback.print_exception(type, value, tb)
        ipdb.post_mortem(tb)

sys.excepthook = info
sys.breakpointhook = ipdb.set_trace


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    '''A collection of utilities for inspecting and manipulating
    spectral libraries.
    '''
    format_string = '[%(asctime)s] %(levelname).1s | %(name)s | %(message)s'

    logging.basicConfig(
        level='INFO',
        stream=sys.stderr,
        format=format_string,
        datefmt="%H:%M:%S")


@main.command("describe", short_help=("Produce a minimal textual description"
                                     " of a spectral library"))
@click.argument('path', type=click.Path(exists=True))
@click.option("-d", "--diagnostics", is_flag=True,
              help="Run more diagnostics, greatly increasing runtime but producing additional information")
def describe(path, diagnostics=False):
    '''Produces a minimal textual description of a spectral library.
    '''
    click.echo("Describing \"%s\"" % (path,))
    if SQLIndex.exists(path):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    library = SpectrumLibrary(filename=path, index_type=index_type)
    click.echo(f"Format: {library.format}")
    click.echo(f"Spectrum Count: {library.__len__()}")
    for attr in library.attributes:
        if not attr.group_id:
            click.echo(f"{attr.key}={attr.value}")
        else:
            click.echo(f"[{attr.group_id}]{attr.key}={attr.value}")


@main.command("index", short_help="Build an on-disk index for a spectral library")
@click.argument('inpath', type=click.Path(exists=True))
def build_index(inpath):
    library = SpectrumLibrary(filename=inpath, index_type=SQLIndex)


@main.command("convert", short_help=("Convert a spectral library from one format to another"))
@click.argument('inpath', type=click.Path(exists=True))
@click.argument("outpath", type=click.Path())
@click.option("-f", "--format", type=click.Choice(["text", "json"]), default="text")
def convert(inpath, outpath, format=None):
    '''Convert a spectral library from one format to another. If `outpath` is `-`,
    instead of writing to file, data will instead be sent to STDOUT.
    '''
    if format is None:
        format = "text"
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    click.echo(f"Opening {inpath}")
    library = SpectrumLibrary(filename=inpath, index_type=index_type)
    click.echo(f"Writing to {outpath}")
    fh = click.open_file(outpath, mode='w')
    library.write(fh, format)


@main.command("index", short_help="Build an on-disk index for a spectral library")
@click.argument('inpath', type=click.Path(exists=True))
def build_index(inpath):
    library = SpectrumLibrary(filename=inpath, index_type=SQLIndex)



def progress_logger(iterable, label, increment: int=100):
    n = len(iterable)
    for i, item in enumerate(iterable):
        if i % increment == 0 and i:
            logger.info(f"... {label} {i}/{n} ({i * 100 / n:0.2f}%)")
        yield item


@main.command(short_help="Semantically validate a spectral library")
@click.argument('inpath', type=click.Path(exists=True))
@click.option("-p", "--profile", "profiles", type=click.Choice(
    ["consensus", "single", "silver", "peptide", "gold"],
    case_sensitive=False),
    multiple=True)
def validate(inpath, profiles=None):
    if profiles is None:
        profiles = []
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex

    logger.info(f"Loading library {inpath}...")
    library = SpectrumLibrary(filename=inpath, index_type=index_type)

    logger.info(f"Loading validators...")
    chain = validator.get_validator_for("base")
    chain |= validator.get_object_validator_for("peak_annotations")
    for profile in profiles:
        if profile is None:
            continue
        logger.info(f"... {profile}")
        chain |= validator.get_validator_for(profile)

    logger.info(f"Validating {inpath}...")
    n_spectra = len(library)
    increment = max(min(n_spectra // 10, 5000), 1)
    chain.validate_library(library, progress_logger(library, "Validating spectra", increment))

    by_level: DefaultDict[RequirementLevel, List[validator.ValidationError]] = DefaultDict(list)
    for message in chain.error_log:
        by_level[message.requirement_level].append(message)

    for level, bucket in sorted(by_level.items()):
        logger.info(f"Found {len(bucket)} violations for {level.name.upper()} rules")
        for err in bucket:
            logger.warn(f"... {err.message}")


if __name__ == "__main__":
    main()
