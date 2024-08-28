PYTHON?=python -X dev
PYTEST_ADDOPTS?=
MYPY_ADDOPTS?=

all: help

help: 			## Show this help
	@echo -e "Specify a command. The choices are:\n"
	@grep -E '^[0-9a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[0;36m%-12s\033[m %s\n", $$1, $$2}'
	@echo ""
.PHONY: help

# {{{ linting

format: isort black pyproject clangfmt mesonfmt		## Run all formatting scripts
.PHONY: format

fmt: format
.PHONY: fmt

isort:			## Run ruff isort fixes over the source code
	ruff check --fix --select=I src tests examples docs
	ruff check --fix --select=RUF022 src
	@echo -e "\e[1;32mruff isort clean!\e[0m"
.PHONY: isort

black:			## Run ruff format over the source code
	ruff format src tests examples docs
	@echo -e "\e[1;32mruff format clean!\e[0m"
.PHONY: black

pyproject:		## Run pyproject-fmt over the configuration
	$(PYTHON) -m pyproject_fmt --indent 4 pyproject.toml
	@echo -e "\e[1;32mpyproject clean!\e[0m"
.PHONY: pyproject

clangfmt:		## Format wrapper code
	clang-format -i src/wrapper/cg_descent_wrap.cpp
	@echo -e "\e[1;32mclang-format clean!\e[0m"

mesonfmt: 		## Format meson.build
	meson fmt -i meson.build
	@echo -e "\e[1;32mmeson fmt clean!\e[0m"

lint: typos reuse ruff mypy		## Run linting checks
.PHONY: lint

typos:			## Run typos over the source code and documentation
	typos --sort
	@echo -e "\e[1;32mtypos clean!\e[0m"
.PHONY: typos

reuse:			## Check REUSE license compliance
	$(PYTHON) -m reuse lint
	@echo -e "\e[1;32mREUSE compliant!\e[0m"
.PHONY: reuse

ruff:			## Run ruff checks over the source code
	ruff check src tests examples
	@echo -e "\e[1;32mruff clean!\e[0m"
.PHONY: ruff

mypy:			## Run mypy checks over the source code
	$(PYTHON) -m mypy src tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"
.PHONY: mypy

# }}}

# {{{ testing

REQUIREMENTS=\
	requirements-dev.txt \
	requirements.txt

requirements-dev.txt: pyproject.toml
	uv pip compile --upgrade --universal --python-version '3.10' \
		--extra dev \
		-o $@ $<
.PHONY: requirements-dev.txt

requirements.txt: pyproject.toml
	uv pip compile --upgrade --universal --python-version '3.10' \
		-o $@ $<
.PHONY: requirements.txt

pin: $(REQUIREMENTS)	## Pin dependency versions to requirements.txt
.PHONY: pin

develop:		## Install project in editable mode
	@rm -rf build
	@rm -rf dist
	$(PYTHON) -m pip install \
		--verbose \
		--no-build-isolation \
		--config-settings setup-args='-Duse-blas=true' \
		--editable .

pip-install:	## Install pinned dependencies from requirements.txt
	$(PYTHON) -m pip install --upgrade pip pybind11 meson-python ninja
	$(PYTHON) -m pip install --upgrade poetry
	$(PYTHON) -m pip install \
		--verbose \
		--requirement requirements-dev.txt \
		--no-build-isolation \
		--config-settings setup-args='-Duse-blas=false' \
		--editable .
.PHONY: pip-install

stubgen:		## Generate stubs for binary module
	$(PYTHON) -m pybind11_stubgen \
		--numpy-array-use-type-var \
		--output src \
		pycgdescent._cg_descent

test:			## Run pytest tests
	$(PYTHON) -m pytest \
		--junit-xml=pytest-results.xml \
		-rswx --durations=25 -v -s \
		$(PYTEST_ADDOPTS)
.PHONY: test

run-examples:	## Run examples with default options
	@for ex in $$(find examples -name "*.py"); do \
		echo "::group::Running $${ex}"; \
		$(PYTHON) "$${ex}"; \
		echo "::endgroup::"; \
	done
.PHONY: run-examples

# }}}

ctags:			## Regenerate ctags
	ctags --recurse=yes \
		--tag-relative=yes \
		--exclude=.git \
		--exclude=docs \
		--python-kinds=-i \
		--language-force=python
.PHONY: ctags

clean:			## Remove various build artifacts
	rm -rf build dist
	rm -rf docs/_build
.PHONY: clean

purge: clean	## Remove various temporary files
	rm -rf .ruff_cache .pytest_cache .mypy_cache
.PHONY: purge
