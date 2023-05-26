"""
Generate human-readable documentation from JSON validation rules.
"""

import json
import re
from collections import defaultdict
from functools import partial
from glob import glob

from jsonschema import validate
from tomark import Tomark


def rules_to_markdown(rules):
    rule_dict = defaultdict(lambda: defaultdict(list))
    for rule in rules["rules"]:
        rule_dict[rule["path"]][rule["requirement_level"]].append(rule)

    lines = []
    lines.append(f"# {rules['name'].upper()}\n")
    for path, path_rules in rule_dict.items():
        lines.append(f"## `{path}`\n")
        for level, level_rules in path_rules.items():
            # Combine all rule attributes into single list (#TODO: how to handle combination logic?)
            rule_attrs = []
            for rule in level_rules:
                rule_attrs.extend(rule["attr"])

            fields = []
            for attr in rule_attrs:
                field = dict()

                # Parse name and accession fields
                field["Name"] = f"{attr['name']} ({attr['accession']})"

                # Parse definition field
                info = []
                if "definition" in attr:
                    info.append(attr["definition"])
                if "notes" in attr:
                    info.append(f"**Notes:** {attr['notes']}")
                field["Info"] = "<br /><br />".join(info).replace(
                    "|", "\|"
                )  # Escape markdown table delimiter TODO: Also for other fields?

                # Parse value field
                if "value" in attr:
                    if attr["value"]["name"] == "value_of_type":
                        attr["value"] = attr["value"]["value"]
                    elif attr["value"]["name"] == "value_is_child_of":
                        attr["value"] = f"Children of {attr['value']['accession']}"
                    # TODO add other value types
                    field["Value"] = attr["value"]
                else:
                    field["Value"] = "Undefined"

                # Parse units field
                field["Allowed units"] = (
                    ", ".join(attr["units"]) if "units" in attr else "/"
                )

                # Parse repeatable field
                field["Repeatable"] = attr["repeatable"]

                fields.append(field)

            lines.append(f"### {level}\n")
            # lines.append(f"Combination logic: `{rule['combination_logic']}`\n")
            lines.append("\n")
            lines.append(Tomark.table(fields))
            lines.append("\n")

    # Add link to OLS for all accessions
    add_links = partial(
        re.sub,
        pattern=r"MS:(\d{7})",
        repl=r"[\g<0>](https://www.ebi.ac.uk/ols4/ontologies/ms/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FMS_\g<1>)",
    )
    lines = [add_links(string=line) for line in lines]
    lines.append("\n")
    return lines


def main():
    lines = []
    for json_filename in glob(JSON_FILES_GLOB):
        # Read JSON rules
        with open(json_filename, "rt", encoding="utf-8") as json_file:
            rules = json.load(json_file)

        # Validate
        with open(SCHEMA_FILENAME, "rt", encoding="utf-8") as validator_rules_schema:
            schema = json.load(validator_rules_schema)
        validate(rules, schema)

        # Convert to markdown
        lines.extend(rules_to_markdown(rules))

    # Write to file
    with open(MD_FILENAME, "wt") as md_file:
        for line in lines:
            md_file.write(line)


if __name__ == "__main__":
    SCHEMA_FILENAME = (
        "implementations/python/mzlib/validate/rules/validator-rules-schema.json"
    )
    # JSON_FILES_GLOB = "implementations/python/mzlib/validate/rules/*.json"
    JSON_FILES_GLOB = "implementations/python/mzlib/validate/rules/all.json"
    MD_FILENAME = "docs/metadata-rules.md"

    main()
