import logging
import sys
import traceback
import click
import logging

from typing import DefaultDict, List

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.index import MemoryIndex, SQLIndex
from mzlib.backends.base import SpectralLibraryBackendBase
from mzlib.validate import validator
from mzlib.validate.level import RequirementLevel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

logger = logging.getLogger(__name__)


def info(type, value, tb):
    if not sys.stderr.isatty():
        click.secho("Running interactively, not starting debugger", fg='yellow')
        sys.__excepthook__(type, value, tb)
    else:
        import pdb
        traceback.print_exception(type, value, tb)
        pdb.post_mortem(tb)


def set_breakpoint_hook():
    try:
        import pdb
        sys.breakpointhook = pdb.set_trace
    except ImportError:
        pass


def _display_tree(tree, indent: int=0):
    if isinstance(tree, dict):
        if not tree:
            return
        for k, v in tree.items():
            if isinstance(v, (list, dict)) and not v:
                continue
            click.echo(f"{' ' * indent}└-{k}", err=True, nl=True)
            _display_tree(v, indent + 2)
    elif isinstance(tree, list):
        for v in tree:
            _display_tree(v, indent + 2)
    else:
        click.echo(f"{' ' * indent}└-{tree}", err=True, nl=True)


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    '''A collection of utilities for inspecting and manipulating
    spectral libraries.
    '''

    format_string = '[%(asctime)s] %(levelname).1s | %(name)s | %(message)s'
    sys.excepthook = info
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
@click.option("-i", "--input-format", type=click.Choice(sorted(SpectralLibraryBackendBase._file_extension_to_implementation)))
@click.option("-f", "--format", type=click.Choice(["text", "json"]), default="text")
@click.option("-k", "--library-attribute", "library_attributes", type=(str, str), multiple=True, help="Specify an attribute to add to the library metadata section. May be repeated.")
def convert(inpath, outpath, format=None, header_file=None, library_attributes=(), input_format=None):
    '''Convert a spectral library from one format to another. If `outpath` is `-`,
    instead of writing to file, data will instead be sent to STDOUT.
    '''
    if format is None:
        format = "text"
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    click.echo(f"Opening {inpath}", err=True)
    library = SpectrumLibrary(filename=inpath, index_type=index_type, format=input_format)
    if library_attributes:
        for k, v in library_attributes:
            library.add_attribute(k, v)
    click.echo(f"Writing to {outpath}", err=True)
    fh = click.open_file(outpath, mode='w')
    library.write(fh, format)

    _display_tree(library.summarize_parsing_errors())


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
