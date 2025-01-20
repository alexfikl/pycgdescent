PYTHON := 'python -X dev'

_default:
    @just --list

# {{{ formatting

alias fmt: format

[doc('Reformat all source code')]
format: isort black pyproject clangfmt justfmt

[doc('Run ruff isort fixes over the source code')]
isort:
    ruff check --fix --select=I src tests examples docs
    ruff check --fix --select=RUF022 src
    @echo -e "\e[1;32mruff isort clean!\e[0m"

[doc('Run ruff format over the source code')]
black:
    ruff format src tests examples docs
    @echo -e "\e[1;32mruff format clean!\e[0m"

[doc('Run pyproject-fmt over the configuration')]
pyproject:
    {{ PYTHON }} -m pyproject_fmt --indent 4 pyproject.toml
    @echo -e "\e[1;32mpyproject clean!\e[0m"

[doc('Run clang-format over wrapper code')]
clangfmt:
    clang-format -i src/wrapper/cg_descent_wrap.cpp
    @echo -e "\e[1;32mclang-format clean!\e[0m"

[doc('Run just --fmt over the justfile')]
justfmt:
    just --unstable --fmt
    @echo -e "\e[1;32mjust --fmt clean!\e[0m"

# }}}
# {{{ linting

[doc('Run all linting checks over the source code')]
lint: typos reuse ruff mypy

[doc('Run typos over the source code and documentation')]
typos:
    typos --sort
    @echo -e "\e[1;32mtypos clean!\e[0m"

[doc('Check REUSE license compliance')]
reuse:
    {{ PYTHON }} -m reuse lint
    @echo -e "\e[1;32mREUSE compliant!\e[0m"

[doc('Run ruff checks over the source code')]
ruff:
    ruff check src tests examples
    @echo -e "\e[1;32mruff clean!\e[0m"

[doc('Run mypy checks over the source code')]
mypy:
    {{ PYTHON }} -m mypy src tests examples
    @echo -e "\e[1;32mmypy clean!\e[0m"

# }}}
# {{{ pin

[private]
requirements_test_txt:
    uv pip compile --upgrade --universal --python-version '3.10' \
        --extra test \
        -o requirements-test.txt pyproject.toml

[private]
requirements_txt:
    uv pip compile --upgrade --universal --python-version '3.10' \
        -o requirements.txt pyproject.toml

[doc('Pin dependency versions to requirements.txt')]
pin: requirements_txt requirements_test_txt

# }}}
# {{{ develop

[doc('Install project in editable mode')]
develop:
    @rm -rf build
    @rm -rf dist
    {{ PYTHON }} -m pip install \
        --verbose \
        --no-build-isolation \
        --config-settings setup-args='-Duse-blas=true' \
        --editable .

[doc("Editable install using pinned dependencies from requirements-test.txt")]
pip-install:
    {{ PYTHON }} -m pip install --upgrade pip pybind11 meson-python ninja poetry
    {{ PYTHON }} -m pip install \
        --verbose \
        --requirement requirements-test.txt \
        --no-build-isolation \
        --config-settings setup-args='-Duse-blas=false' \
        --editable .

[doc("Generate typing stubs for binary module")]
stubgen:
    {{ PYTHON }} -m pybind11_stubgen \
        --numpy-array-use-type-var \
        --output src \
        pycgdescent._cg_descent
    @ruff format --quiet src/pycgdescent/_cg_descent.pyi

[doc("Remove various build artifacts")]
clean:
    rm -rf build dist
    rm -rf docs/_build

[doc("Remove various temporary files")]
purge: clean
    rm -rf .ruff_cache .pytest_cache .mypy_cache tags

[doc("Regenerate ctags")]
ctags:
    ctags --recurse=yes \
        --tag-relative=yes \
        --exclude=.git \
        --exclude=docs \
        --python-kinds=-i \
        --language-force=python

# }}}
# {{{ tests

[doc("Run pytest tests")]
test *PYTEST_ADDOPTS:
    {{ PYTHON }} -m pytest \
        --junit-xml=pytest-results.xml \
        -rswx --durations=25 -v -s \
        {{ PYTEST_ADDOPTS }}

[doc("Run examples with default options")]
examples:
    for ex in `find examples -name "*.py"`; do \
        echo "::group::Running ${ex}"; \
        {{ PYTHON }} "${ex}"; \
        echo "::endgroup::"; \
    done

# }}}