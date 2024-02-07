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

format: black isort pyproject					## Run all formatting scripts
.PHONY: format

fmt: format
.PHONY: fmt

pyproject:		## Run pyproject-fmt over the configuration
	$(PYTHON) -m pyproject_fmt --indent 4 pyproject.toml
	@echo -e "\e[1;32mpyproject clean!\e[0m"
.PHONY: pyproject

black:			## Run ruff format over the source code
	ruff format src tests examples docs setup.py
	@echo -e "\e[1;32mruff format clean!\e[0m"
.PHONY: black

isort:			## Run ruff isort fixes over the source code
	ruff check --fix --select=I src tests examples docs setup.py
	ruff check --fix --select=RUF022 src
	@echo -e "\e[1;32mruff isort clean!\e[0m"
.PHONY: isort

lint: ruff mypy codespell reuse					## Run linting checks
.PHONY: lint

ruff:			## Run ruff checks over the source code
	ruff check src tests examples
	@echo -e "\e[1;32mruff clean!\e[0m"
.PHONY: ruff

mypy:			## Run mypy checks over the source code
	$(PYTHON) -m mypy src tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"
.PHONY: mypy

codespell:		## Run codespell over the source code and documentation
	@codespell --summary \
		--skip _build --skip src/wrapper \
		--uri-ignore-words-list '*' \
		--ignore-words .codespell-ignore \
		src tests examples docs
.PHONY: codespell

reuse:			## Check REUSE license compliance
	$(PYTHON) -m reuse lint
	@echo -e "\e[1;32mREUSE compliant!\e[0m"
.PHONY: reuse

# }}}

# {{{ testing

REQUIREMENTS=\
	requirements-dev.txt \
	requirements.txt

requirements-dev.txt: pyproject.toml
	$(PYTHON) -m piptools compile \
		--resolver=backtracking --allow-unsafe \
		--strip-extras --upgrade --extra dev \
		-o $@ $<
.PHONY: requirements-dev.txt

requirements.txt: pyproject.toml
	$(PYTHON) -m piptools compile \
		--resolver=backtracking --allow-unsafe \
		--strip-extras --upgrade \
		-o $@ $<
.PHONY: requirements.txt

pin: $(REQUIREMENTS)	## Pin dependency versions to requirements.txt
.PHONY: pin

pip-install:	## Install pinned dependencies from requirements.txt
	$(PYTHON) -m pip install -r requirements-dev.txt
	$(PYTHON) -m pip install --no-build-isolation -e .
.PHONY: pip-install

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
