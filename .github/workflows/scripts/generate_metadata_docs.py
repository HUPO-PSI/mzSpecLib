"""
Generate human-readable documentation from JSON validation rules.
"""

import json
from glob import glob

from tomark import Tomark

JSON_FILES_GLOB = "implementations/python/mzlib/validate/rules/*.json"
MD_FILENAME = "docs/metadata-rules.md"


def generate_md():
    lines = []
    for json_filename in glob(JSON_FILES_GLOB):
        with open(json_filename, "rt") as json_file:
            rules = json.load(json_file)

        if rules["$schema"] != "./validator-rules-schema.json":
            continue

        lines.append(f"## {rules['name'].upper()}\n")

        for rule in rules["rules"]:
            lines.append(f"### **`{rule['level']}`** {rule['path']}: {rule['id'].replace('_', ' ')}\n")
            lines.append(f"Combination logic: `{rule['combination_logic']}`\n")
            lines.append("\n")
            lines.append(Tomark.table(
                rule["attr"]
            ))
            lines.append("\n")

    with open(MD_FILENAME, "wt") as md_file:
        for line in lines:
            md_file.write(line)


if __name__ == "__main__":
    generate_md()
