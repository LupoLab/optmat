# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Sphinx configuration for the optmat manual.

Markdown throughout (myst-parser), NumPy-style docstrings rendered by numpydoc, an
autosummary API reference and a sphinx-gallery of data and fit figures that executes on
every build, its figures saved as SVG. Warnings are errors and missing references are
warnings (``-W -n``).
"""

import os
import sys
from importlib.metadata import version as _version
from pathlib import Path

# The gallery pages share their plotting through docs/gallery/_helpers.py.
sys.path.insert(0, str(Path(__file__).parent / "gallery"))
# The gallery's image scraper and the source-catalogue extension.
sys.path.insert(0, str(Path(__file__).parent / "_ext"))

project = "optmat"
author = "John C. Travers and contributors"
copyright = "2026, John C. Travers and contributors"
release = _version("optmat")
version = ".".join(release.split(".")[:2])

extensions = [
    "myst_parser",
    "source_catalogue",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_gallery.gen_gallery",
]

templates_path = ["_templates"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
exclude_patterns = [
    "_build",
    "_inventories",
    "Thumbs.db",
    ".DS_Store",
    "gallery/README.rst",
    "gallery/*/README.rst",
]

myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 3

autosummary_generate = True
autosummary_imported_members = False
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False

# Each inventory is read from the vendored copy when the site is unreachable.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", (None, "_inventories/python.inv")),
    "numpy": ("https://numpy.org/doc/stable/", (None, "_inventories/numpy.inv")),
    "scipy": ("https://docs.scipy.org/doc/scipy/", (None, "_inventories/scipy.inv")),
    "matplotlib": (
        "https://matplotlib.org/stable/",
        (None, "_inventories/matplotlib.inv"),
    ),
}
intersphinx_timeout = 10
nitpicky = True
nitpick_ignore = []
nitpick_ignore_regex = [
    ("py:class", r"^jaxtyping\..*"),
    ("py:class", r"^numpy\._typing\..*"),
    ("py:class", r"^ndarray$"),
    ("py:class", r"^'.*'$"),
]

sphinx_gallery_conf = {
    # Every figure as SVG, with a PNG for the thumbnail (docs/_ext/svg_scraper.py).
    "image_scrapers": ("svg_scraper.svg_scraper",),
    "examples_dirs": ["gallery"],
    "gallery_dirs": ["auto_gallery"],
    "subsection_order": [
        "gallery/data",
        "gallery/dispersion",
        "gallery/accuracy",
        "gallery/nonlinear",
        "gallery/ionisation",
        "gallery/solids",
    ],
    "within_subsection_order": "FileNameSortKey",
    "filename_pattern": r"/plot_",
    "ignore_pattern": r"__init__\.py|_helpers\.py",
    "download_all_examples": False,
    "remove_config_comments": True,
    "show_signature": False,
}

html_theme = "pydata_sphinx_theme"
html_title = "optmat"
html_theme_options = {
    "navigation_with_keys": False,
    "show_prev_next": True,
    "header_links_before_dropdown": 5,
    "navigation_depth": 3,
    "show_nav_level": 1,
    "github_url": "https://github.com/LupoLab/optmat",
    "use_edit_page_button": True,
}
html_context = {
    "github_user": "LupoLab",
    "github_repo": "optmat",
    "github_version": "main",
    "doc_path": "docs",
}
# Read the Docs sets the canonical URL of the version it builds.
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")
html_static_path = ["_static"]


def _edit_page_context(_app, pagename, _templatename, context, _doctree):
    # Browse individual records in the alphabetical index; keep the data chapters
    # visible instead of expanding all source records in every page's sidebar.
    if pagename.startswith("generated/sources/"):
        context["theme_navigation_depth"] = 2
    # Generated pages have no source file in GitHub to edit.
    if pagename.startswith(
        ("auto_gallery/", "api/generated/", "generated/")
    ) or pagename.endswith("sg_execution_times"):
        context["theme_use_edit_page_button"] = False
    # These pages include their text from the repository root.
    elif pagename in {"changelog", "register"}:
        context["edit_page_provider_name"] = "GitHub"
        context["edit_page_url_template"] = (
            "https://github.com/LupoLab/optmat/edit/main/"
            "{{ file_name.rsplit('/', 1)[-1].split('.')[0].upper() }}.md"
        )


def setup(app):
    app.connect("html-page-context", _edit_page_context)
