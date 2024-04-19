"""Configuration file for the Sphinx documentation builder."""

# Project information
project = "mzSpecLib"
author = "HUPO-PSI"
github_project_url = "https://github.com/HUPO-PSI/mzSpecLib"
github_doc_root = "https://github.com/HUPO-PSI/mzSpecLib/tree/master/docs"

# General configuration
extensions = ["myst_parser"]
source_suffix = [".rst"]
master_doc = "index"
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/HUPO-PSI/mzSpecLib",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        }
    ]
}


def setup(app):
    config = {"enable_eval_rst": True}  # noqa: F841
