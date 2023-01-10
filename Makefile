PYTHON?=python -X dev

all: flake8 pylint mypy

# {{{ linting

black:
	$(PYTHON) -m setup_cfg_fmt --include-version-classifiers setup.cfg
	$(PYTHON) -m pyproject_fmt --indent 4 pyproject.toml
	$(PYTHON) -m isort pycgdescent tests examples docs
	$(PYTHON) -m black --safe --target-version py38 pycgdescent examples tests docs setup.py

flake8:
	$(PYTHON) -m flake8 pycgdescent examples tests docs setup.py
	@echo -e "\e[1;32mflake8 clean!\e[0m"

pylint:
	PYTHONWARNINGS=ignore $(PYTHON) -m pylint pycgdescent tests/*.py examples/*.py
	@echo -e "\e[1;32mpylint clean!\e[0m"

mypy:
	$(PYTHON) -m mypy --strict --show-error-codes pycgdescent tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"

pyright:
	pyright --stats pycgdescent tests examples
	@echo -e "\e[1;32mpyright clean!\e[0m"

ruff:
	ruff pycgdescent tests examples
	@echo -e "\e[1;32mruff clean!\e[0m"

pytype:
	$(PYTHON) -m pytype \
		--strict-parameter-checks \
		--strict-primitive-comparisons \
		pycgdescent tests examples
	@echo -e "\e[1;32mpytype clean!\e[0m"

codespell:
	@codespell --summary \
		--skip _build \
		--ignore-words .codespell-ignore \
		pycgdescent tests examples docs

reuse:
	@reuse lint
	@echo -e "\e[1;32mREUSE compliant!\e[0m"

manifest:
	@check-manifest
	@echo -e "\e[1;32mMANIFEST.in is up to date!\e[0m"

# }}}

# {{{ testing

pin:
	$(PYTHON) -m piptools compile \
		--resolver=backtracking \
		--extra dev --upgrade \
		-o requirements.txt setup.cfg

pip-install:
	$(PYTHON) -m pip install --upgrade pip wheel setuptools
	$(PYTHON) -m pip install -r requirements.txt -e .

test:
	$(PYTHON) -m pytest -rswx --durations=25 -v -s

run-examples:
	@for ex in $$(find examples -name "*.py"); do \
		echo -e "\x1b[1;32m===> \x1b[97mRunning $${ex}\x1b[0m"; \
		$(PYTHON) "$${ex}"; \
		sleep 1; \
	done

# }}}

ctags:
	ctags --recurse=yes \
		--tag-relative=yes \
		--exclude=.git \
		--exclude=docs \
		--python-kinds=-i \
		--language-force=python

.PHONY: all black flake8 mypy pylint test pip-install
