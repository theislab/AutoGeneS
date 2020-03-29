import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import sphinx_rtd_theme
import autogenes

# -- Project information -----------------------------------------------------

project = 'AutoGeneS'
copyright = '2020, Hananeh Aliee, Maxim Schmidt'
author = 'Hananeh Aliee, Maxim Schmidt'
version = 'v1.0'

# -- General configuration ---------------------------------------------------

needs_sphinx = '2.0'

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary'
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False
#napoleon_custom_sections = [('Params', 'Parameters')]

templates_path = ['_templates']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
#
html_theme = 'sphinx_rtd_theme'

html_show_sphinx = False

html_logo = '_static/logo.svg'

html_theme_options = {
    'logo_only': True,
    'display_version': True,
    'navigation_depth': 2,
    'prev_next_buttons_location': None
}

html_static_path = ['_static']


def setup(app):
  app.add_stylesheet('custom.css')

# Options

master_doc = 'index'
