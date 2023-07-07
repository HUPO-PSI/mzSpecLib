import json
from pathlib import Path

import pandas as pd
from jsonschema import validate


def excel_terms_to_dict(excel_file, sheet_name):
    terms_sheet = pd.read_excel(excel_file, sheet_name=sheet_name)
    # Remove trailing spaces from all dataframe values
    terms_sheet = terms_sheet.apply(
        lambda x: x.str.strip() if x.dtype == "object" else x
    )
    terms = terms_sheet.to_dict(orient="records")
    return terms


def parse_terms(terms):
    """Translate terms from excel to metadata rules."""
    for term in terms:
        attribute = {
            "accession": term["accession"],
            "name": term["name"],
            "repeatable": term["repeatable"],
        }

        if isinstance(term["notes"], str):
            attribute["notes"] = term["notes"]

        if term["value"].startswith("Children of"):
            attribute["value"] = {
                "name": "value_is_child_of",
                "accession": term["value"].split("Children of ")[1].strip(),
            }
        elif term["value"].startswith("xsd:"):
            attribute["value"] = {"name": "value_of_type", "value": term["value"]}
        yield {
            "id": "has_" + term["name"].lower().replace(" ", "_"),
            "path": term["path"],
            "attr": [attribute],
            "requirement_level": term["level"],
        }


def main(schema_file, excel_file, sheet_name, json_file):
    # Read Excel
    terms = excel_terms_to_dict(excel_file, sheet_name)

    # Parse terms to validation rules
    rules = {
        "$schema": "./validator-rules-schema.json",
        "name": "all_metadata",
        "rules": [t for t in parse_terms(terms)],
    }

    # Validate
    validate(rules, json.load(open(schema_file)))

    # Write JSON
    with open(json_file, "w") as f:
        json.dump(rules, f, indent=2)


if __name__ == "__main__":
    current_path = Path(__file__).parent
    schema_file = current_path / "../validator-rules-schema.json"
    excel_file = current_path / "terms.xlsx"
    sheet_name = "general"
    json_file = current_path / "../all.json"

    main(schema_file, excel_file, sheet_name, json_file)
