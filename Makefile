PYTHON?=python

all: flake8 pylint mypy

black:
	$(PYTHON) -m black pycgdescent examples tests docs setup.py

flake8:
	$(PYTHON) -m flake8 pycgdescent examples tests docs setup.py
	@echo -e "\e[1;32mflake8 clean!\e[0m"

mypy:
	$(PYTHON) -m mypy --strict --show-error-codes pycgdescent tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"

pylint:
	PYTHONWARNINGS=ignore $(PYTHON) -m pylint pycgdescent tests/*.py examples/*.py
	@echo -e "\e[1;32mpylint clean!\e[0m"

test:
	$(PYTHON) -m pytest -rswx --durations=25 -v -s

test-examples:
	@for ex in $$(find examples -name "*.py"); do \
		echo -e "\x1b[1;32m===> \x1b[97mRunning $${ex}\x1b[0m"; \
		$(PYTHON) "$${ex}"; \
		sleep 1; \
	done

pip-install:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e '.[dev]'

tags:
	ctags -R

.PHONY: all black flake8 mypy pylint test pip-install
