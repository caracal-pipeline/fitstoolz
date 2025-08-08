import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'fitstools'
copyright = '2025 Wits Centre for Astrophysics'
author = 'Sphesihle Makhathini, Mika Naidoo, Athanaseus Ramaila'
release = '0.0.1b1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
