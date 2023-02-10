import os
import pathlib
import string


from railroad import (Diagram, Choice, Group, Optional, Terminal,
                      NonTerminal, Sequence, OneOrMore, ZeroOrMore, Stack)
import io

from pyteomics.mass import std_aa_comp

image_dir = pathlib.Path("schema_images/")
os.makedirs(image_dir, exist_ok=True)


UPPER_CASE_LETTER = Choice(0, *string.ascii_uppercase)

LOWER_CASE_LETTER = Choice(0, *string.ascii_lowercase)

SYMBOL = Choice(0, *[(c) for c in string.punctuation])

DIGIT = Choice(0, *string.digits)


ATOM_COUNT = Sequence(
    NonTerminal("UPPER_CASE_LETTER"),
    ZeroOrMore(NonTerminal("LOWER_CASE_LETTER")),
    OneOrMore(NonTerminal("DIGIT"))
)

NUMBER = Sequence(
    OneOrMore(NonTerminal("DIGIT")),
    Optional(Sequence(".", OneOrMore(NonTerminal("DIGIT")))),
    Optional(
        Sequence(
            "e",
            OneOrMore(NonTerminal("DIGIT")),
            Optional(Sequence(".", OneOrMore(NonTerminal("DIGIT"))))
        )
    )
)

ORDINAL = OneOrMore(NonTerminal("DIGIT"))

CHARACTER = Choice(0, NonTerminal("DIGIT"), NonTerminal("UPPER_CASE_LETTER"), NonTerminal("LOWER_CASE_LETTER"))

AMINO_ACID = Choice(0, *list(std_aa_comp)[:-2])

SIGN = Choice(0, "+", "-")

BraceEnclosedContent = Sequence(
    Terminal("["),
    OneOrMore(Choice(0, NonTerminal("CHARACTER"), NonTerminal("SYMBOL"))),
    Terminal("]")
)

IsAuxiliary = Group(Optional(Terminal("&")), "Is Auxiliary")

AnalyteIdentifier = Group(
    Sequence(NonTerminal("ORDINAL"), Terminal("@")), "Analyte Identifier"
)

PeptideIon = Group(
    Sequence(
        Choice(0, *list(map(Terminal, ("a", "b", "c", "x", "y", "z")))),
        NonTerminal("ORDINAL"),
    ),
    "Peptide Ion",
)

ReporterIon = Group(
    Sequence(
        Terminal("r"),
        BraceEnclosedContent
    ),
    "Reporter Ion"
)

InternalIon = Group(
    Sequence(
        Terminal("m"),
        NonTerminal("ORDINAL"),
        Terminal(":"),
        NonTerminal("ORDINAL"),
    ),
    "Internal Peptide Ion",
)

ImmoniumIon = Group(
    Sequence(Terminal("I"), NonTerminal("AMINO_ACID"), Optional(BraceEnclosedContent)),
    "Immonium Ion",
)

PrecursorIon = Group(
    Terminal("p"),
    "Precursor Ion"
)


ChemicalFormula = OneOrMore(NonTerminal('ATOM_COUNT'))


FormulaIon = Group(
    Sequence(
        Terminal("f"),
        Terminal('{'),
        ChemicalFormula,
        Terminal('}')
    ),
    "Formula Ion"
)

ExternalIon = Group(
    Sequence(
        Terminal("_"),
        Terminal('{'),
        OneOrMore(NonTerminal("CHARACTER")),
        Terminal('}'),
    ),
    "External Ion"
)

UnknownIon = Group(
    Sequence(Terminal("?"), Optional(OneOrMore(NonTerminal("DIGIT")))),
    "Unknown Ion"
)

SMILESIon = Group(
    Sequence(
        Terminal("s"),
        Terminal('{'),
        OneOrMore(Terminal("/[^}]/")),
        Terminal('}')
    ),
    "SMILES Ion"
)

IonType = Group(
    Choice(
        0,
        PeptideIon,
        InternalIon,
        ImmoniumIon,
        ReporterIon,
        PrecursorIon,
        FormulaIon,
        ExternalIon,
        SMILESIon,
        UnknownIon,
    ),
    "Ion Type"
)

NeutralLoss = Group(
    Sequence(
        NonTerminal('SIGN'),
        Choice(0, ChemicalFormula, BraceEnclosedContent)
    ),
    "Neutral Loss(es)"
)

Isotope = Group(
    Sequence(
        NonTerminal("SIGN"),
        Optional(NonTerminal("ORDINAL")),
        Terminal("i"),
    ),
    "Isotope"
)

ChargeState = Group(
    Sequence(
        "^",
        NonTerminal("ORDINAL")
    ),
    "Charge State"
)

Adducts = Group(
    Sequence(
        '[',
        'M',
        OneOrMore(
            Sequence(
                NonTerminal('SIGN'),
                ChemicalFormula,
            )
        ),
        ']'
    ),
    "Adducts"
)

MassError = Group(
    Sequence(
        '/',
        NonTerminal("NUMBER"),
        Optional("ppm")
    ),
    "Mass Error"
)

ConfidenceEstimate = Group(
    Sequence(
        "*",
        NonTerminal("NUMBER")
    ),
    "Confidence Estimate"
)


Annotation = (
    Stack(
        IsAuxiliary,
        Optional(
            AnalyteIdentifier
        ),
        IonType,
        ZeroOrMore(NeutralLoss),
        Optional(Isotope),
        Optional(ChargeState),
        Optional(Adducts),
        Optional(MassError),
        Optional(ConfidenceEstimate),
    )
)


def encode_svg(diagram):
    buffer = io.StringIO()
    diagram.writeSvg(buffer.write)
    value: str = buffer.getvalue()
    value = value.replace("\n", " ").replace("'", "&apos;")
    return value


def render_group_to_file(fh, name):
    print("Writing", name)
    tokens = globals()[name]
    pathname: pathlib.Path = (image_dir / name).with_suffix(".svg")
    with pathname.open('wt') as img_fh:
        img_fh.write(encode_svg(Diagram(tokens)))
    fh.write(f"""## {name}\n<img src="{pathname}">\n\n""")


with open("grammar.md", 'wt') as fh:
    fh.write("""# Peak Annotation Grammar\n\n""")
    render_group_to_file(fh, "DIGIT")
    render_group_to_file(fh, "LOWER_CASE_LETTER")
    render_group_to_file(fh, "UPPER_CASE_LETTER")
    render_group_to_file(fh, "SYMBOL")
    render_group_to_file(fh, "ORDINAL")
    render_group_to_file(fh, "NUMBER")
    render_group_to_file(fh, "CHARACTER")
    render_group_to_file(fh, "SIGN")
    render_group_to_file(fh, "ATOM_COUNT")
    render_group_to_file(fh, "AMINO_ACID")
    render_group_to_file(fh, "Annotation")
    fh.write("\n")
