import logging
import sys
import click
import logging
import json
import csv

from typing import DefaultDict, List

from mzlib.spectrum_library import SpectrumLibrary
from mzlib.index import MemoryIndex, SQLIndex
from mzlib.backends.base import FormatInferenceFailure, SpectralLibraryBackendBase
from mzlib.validate import validator
from mzlib.validate.level import RequirementLevel
from mzlib.ontology import ControlledVocabularyResolver

from mzlib.tools.utils import ColoringFormatter

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

logger = logging.getLogger(__name__)


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
    """A collection of utilities for inspecting and manipulating spectral libraries."""
    format_string = '[%(asctime)s] %(levelname).1s | %(name)s | %(message)s'

    logging.basicConfig(
        level='INFO',
        stream=sys.stderr,
        format=format_string,
        datefmt="%H:%M:%S")

    fmtr = ColoringFormatter(format_string, datefmt='%H:%M:%S')

    for handler in logging.getLogger().handlers:
        handler.setFormatter(
            fmtr
        )


@main.command("describe", short_help=("Produce a minimal textual description"
                                      " of a spectral library"))
@click.argument('path', type=click.Path(exists=True))
@click.option("-d", "--diagnostics", is_flag=True,
              help="Run more diagnostics, greatly increasing runtime but producing additional information")
@click.option("-i", "--input-format",
              type=click.Choice(sorted(SpectralLibraryBackendBase._file_extension_to_implementation)),
              help='The file format of the input file. If omitted, will attempt to infer automatically.')
def describe(path, diagnostics=False, input_format=None):
    """Produce a minimal textual description of a spectral library."""
    click.echo("Describing \"%s\"" % (path,))
    if SQLIndex.exists(path):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    try:
        library = SpectrumLibrary(filename=path, index_type=index_type, format=input_format)
    except FormatInferenceFailure as err:
        click.echo(f"{err}", err=True)
        raise click.Abort()
    click.echo(f"Format: {library.format}")
    click.echo(f"Spectrum Count: {library.__len__()}")
    for attr in library.attributes:
        if not attr.group_id:
            click.echo(f"{attr.key}={attr.value}")
        else:
            click.echo(f"[{attr.group_id}]{attr.key}={attr.value}")


@main.command("convert", short_help=("Convert a spectral library from one format to another"))
@click.argument('inpath', type=click.Path(exists=True))
@click.argument("outpath", type=click.Path())
@click.option("-i", "--input-format", type=click.Choice(sorted(SpectralLibraryBackendBase._file_extension_to_implementation)))
@click.option("-f", "--format", type=click.Choice(["text", "json", "msp"]), default="text")
@click.option("-k", "--library-attribute", "library_attributes", type=(str, str), multiple=True,
              help="Specify an attribute to add to the library metadata section. May be repeated.")
@click.option("-K", "--header-file", type=click.Path(readable=True),
              help="Specify a file to read name-value pairs from. May be JSON or TAB-separated")
def convert(inpath, outpath, format=None, header_file=None, library_attributes=(), input_format=None):
    """
    Convert a spectral library from one format to another. If `outpath` is `-`,
    instead of writing to file, data will instead be sent to STDOUT.
    """
    if format is None:
        format = "text"
    if format == 'msp':
        click.secho("MSP conversion is arbitrary and lossy", fg='yellow', err=True)
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex
    click.echo(f"Opening {inpath}", err=True)
    try:
        library = SpectrumLibrary(filename=inpath, index_type=index_type, format=input_format)
    except FormatInferenceFailure as err:
        click.echo(f"{err}", err=True)
        raise click.Abort()
    if header_file:
        library_attributes = list(library_attributes)
        if header_file.endswith(".json"):
            library_attributes.extend(
                json.load(open(header_file, 'rt')).items())
        else:
            library_attributes.extend(
                csv.reader(open(header_file, 'rt'), delimiter='\t'))
    if library_attributes:
        resolver = ControlledVocabularyResolver()
        for k, v in library_attributes:
            k = resolver.attribute_syntax(k)
            library.add_attribute(k, v)
    click.echo(f"Writing to {outpath}", err=True)
    fh = click.open_file(outpath, mode='w')
    try:
        library.write(fh, format)
    except IOError as err:
        if outpath == '-':
            return
        else:
            raise


    _display_tree(library.summarize_parsing_errors())


@main.command("index", short_help="Build an on-disk index for a spectral library")
@click.argument('inpath', type=click.Path(exists=True))
@click.option("-i", "--input-format", type=click.Choice(sorted(SpectralLibraryBackendBase._file_extension_to_implementation)))
def build_index(inpath, input_format=None):
    """Build an external on-disk SQL-based index for the spectral library"""
    try:
        library = SpectrumLibrary(filename=inpath, index_type=SQLIndex, format=input_format)
    except FormatInferenceFailure as err:
        click.echo(f"{err}", err=True)
        raise click.Abort()


def _progress_logger(iterable, label, increment: int=100):
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
@click.option("-i", "--input-format", type=click.Choice(sorted(SpectralLibraryBackendBase._file_extension_to_implementation)))
def validate(inpath, profiles=None, input_format=None):
    """Semantically and structurally validate a spectral library."""
    if profiles is None:
        profiles = []
    if SQLIndex.exists(inpath):
        index_type = SQLIndex
    else:
        index_type = MemoryIndex

    logger.info(f"Loading library {inpath}...")
    try:
        library = SpectrumLibrary(filename=inpath, index_type=index_type, format=input_format)
    except FormatInferenceFailure as err:
        click.echo(f"{err}", err=True)
        raise click.Abort()

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
    chain.validate_library(library, _progress_logger(library, "Validating spectra", increment))

    by_level: DefaultDict[RequirementLevel, List[validator.ValidationError]] = DefaultDict(list)
    for message in chain.error_log:
        by_level[message.requirement_level].append(message)

    for level, bucket in sorted(by_level.items()):
        logger.info(f"Found {len(bucket)} violations for {level.name.upper()} rules")
        for err in bucket:
            logger.warn(f"... {err.message}")


if __name__ == "__main__":
    main()
