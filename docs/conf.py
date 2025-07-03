# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'Omics Downstream Pipeline'
copyright = '2025, Stanislav Bratchikov'
author = 'Stanislav Bratchikov'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output ------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = 'pipeline_overview.png'

# -- Extension configuration -------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

source_suffix = {
    '.rst': None,
    '.md': 'myst_parser',
}
