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

fmt: black		## Run all formatting scripts
	$(PYTHON) -m setup_cfg_fmt --include-version-classifiers setup.cfg
	$(PYTHON) -m pyproject_fmt --indent 4 pyproject.toml
	$(PYTHON) -m isort pycgdescent tests examples docs
.PHONY: fmt

black:			## Run black over the source code
	$(PYTHON) -m black \
		--safe --target-version py38 \
		pycgdescent examples tests docs setup.py
.PHONY: black

flake8:			## Run flake8 checks over the source code
	$(PYTHON) -m flake8 pycgdescent examples tests docs setup.py
	@echo -e "\e[1;32mflake8 clean!\e[0m"
.PHONY: flake8

pylint:			## Run pylint checks over the source code
	PYTHONWARNINGS=ignore $(PYTHON) -m pylint pycgdescent tests/*.py examples/*.py
	@echo -e "\e[1;32mpylint clean!\e[0m"
.PHONY: pylint

mypy:			## Run mypy checks over the source code
	$(PYTHON) -m mypy \
		--strict --show-error-codes $(MYPY_ADDOPTS) \
		pycgdescent tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"
.PHONY: mypy

pyright:		## Run pyright checks over the source code
	pyright --stats pycgdescent tests examples
	@echo -e "\e[1;32mpyright clean!\e[0m"
.PHONY: pyright

ruff:			## Run ruff checks over the source code
	ruff pycgdescent tests examples
	@echo -e "\e[1;32mruff clean!\e[0m"
.PHONY: ruff

pytype:			## Run pytype checks over the source code
	$(PYTHON) -m pytype \
		--strict-parameter-checks \
		--strict-primitive-comparisons \
		pycgdescent tests examples
	@echo -e "\e[1;32mpytype clean!\e[0m"
.PHONY: pytype

codespell:		## Run codespell over the source code and documentation
	@codespell --summary \
		--skip _build \
		--ignore-words .codespell-ignore \
		pycgdescent tests examples docs
.PHONY: codespell

reuse:			## Check REUSE license compliance
	@reuse lint
	@echo -e "\e[1;32mREUSE compliant!\e[0m"
.PHONY: reuse

manifest:		## Update MANIFEST.in file
	@check-manifest
	@echo -e "\e[1;32mMANIFEST.in is up to date!\e[0m"
.PHONY: manifest

# }}}

# {{{ testing

REQUIREMENTS=\
	requirements-dev.txt \
	requirements.txt

requirements-dev.txt: setup.cfg
	$(PYTHON) -m piptools compile \
		--resolver=backtracking --upgrade \
		--extra dev \
		-o $@ $<

requirements.txt: setup.cfg
	$(PYTHON) -m piptools compile \
		--resolver=backtracking --upgrade \
		-o $@ $<

pin: $(REQUIREMENTS)	## Pin dependency versions to requirements.txt
.PHONY: pin

pip-install:	## Install pinned dependencies from requirements.txt
	$(PYTHON) -m pip install --upgrade pip wheel setuptools
	$(PYTHON) -m pip install -r requirements-dev.txt -e .
.PHONY: pip-install

test:			## Run pytest tests
	$(PYTHON) -m pytest -rswx --durations=25 -v -s $(PYTEST_ADDOPTS)
.PHONY: test

run-examples:	## Run examples with default options
	@for ex in $$(find examples -name "*.py"); do \
		echo -e "\x1b[1;32m===> \x1b[97mRunning $${ex}\x1b[0m"; \
		$(PYTHON) "$${ex}"; \
		sleep 1; \
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
