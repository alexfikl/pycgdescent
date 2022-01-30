PYTHON?=python

all: flake8 pylint mypy

black:
	$(PYTHON) -m black pycgdescent examples tests docs setup.py

flake8:
	$(PYTHON) -m flake8 pycgdescent examples tests docs setup.py
	@echo -e "\e[1;32mflake8 clean!\e[0m"

pylint:
	PYTHONWARNINGS=ignore $(PYTHON) -m pylint pycgdescent tests/*.py examples/*.py
	@echo -e "\e[1;32mpylint clean!\e[0m"

mypy:
	$(PYTHON) -m mypy --strict --show-error-codes pycgdescent tests examples
	@echo -e "\e[1;32mmypy clean!\e[0m"

tags:
	ctags -R

.PHONY: all flake8 pylint mypy
