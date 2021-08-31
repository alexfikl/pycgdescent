# https://www.sphinx-doc.org/en/master/usage/configuration.html
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
try:
    # python >=3.8 only
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

# {{{ project information

m = metadata.metadata("pycgdescent")
project = m["Name"]
author = m["Author"]
copyright = f"2020 {author}"
version = m["Version"]
release = version

# }}}

# {{{ general configuration

# needed extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]

try:
    import sphinxcontrib.spelling       # noqa: F401
    extensions.append("sphinxcontrib.spelling")
except ImportError:
    pass

# extension for source files
source_suffix = ".rst"
# name of the main (master) document
master_doc = "index"
# min sphinx version (needed for `autodoc_type_aliases`)
needs_sphinx = "3.3"
# files to ignore
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
# highlighting
pygments_style = "sphinx"

# }}}

# {{{ internationalization

language = "en"

# sphinxcontrib.spelling options
spelling_lang = "en_US"
tokenizer_lang = "en_US"
spelling_word_list_filename = "wordlist_en.txt"

# }}

# {{{ output

# html
html_theme = "sphinx_rtd_theme"

# }}}

# {{{ extension settings

autoclass_content = "class"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "show-inheritance": None,
}

autodoc_typehints = "description"
autodoc_type_aliases = {
    "FunType": "pycgdescent.FunType",
    "GradType": "pycgdescent.GradType",
    "FunGradType": "pycgdescent.FunGradType",
    "CallbackType": "pycgdescent.CallbackType"
}

# https://github.com/sphinx-doc/sphinx/issues/8127
nitpick_ignore = [
    ("py:class", "..")
]

# }}}

# {{{ links

intersphinx_mapping = {
    "https://docs.python.org/3": None,
    "https://numpy.org/doc/stable": None,
    "https://docs.scipy.org/doc/scipy/reference": None,
}

# }}}
