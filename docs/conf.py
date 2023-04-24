# https://www.sphinx-doc.org/en/master/usage/configuration.html
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
from importlib import metadata

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
    "autoapi.extension",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
]

try:
    import sphinxcontrib.spelling  # noqa: F401

    extensions.append("sphinxcontrib.spelling")
except ImportError:
    pass

# extension for source files
source_suffix = ".rst"
# name of the main (master) document
master_doc = "index"
# min sphinx version
needs_sphinx = "4.0"
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

autoapi_type = "python"
autoapi_dirs = ["."]
autoapi_add_toctree_entry = False

autoapi_python_class_content = "class"
autoapi_member_order = "bysource"
autoapi_options = [
    "show-inheritance",
]
suppress_warnings = ["autoapi"]

# }}}

# {{{ links

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}

# }}}
