PYTHON ?= python3
export PYTHONPATH = $(CURDIR)/src
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
export VECLIB_MAXIMUM_THREADS = 1
export PYTEST_DISABLE_PLUGIN_AUTOLOAD = 1
export MPLCONFIGDIR = $(CURDIR)/build/matplotlib
.PHONY: test validate reproduce-core figures paper clean

test:
	$(PYTHON) -m pytest -q
validate:
	$(PYTHON) scripts/release.py validate
reproduce-core:
	$(PYTHON) scripts/reproduce_core.py
figures:
	$(PYTHON) scripts/release.py figures
paper: figures
	$(PYTHON) scripts/release.py build
clean:
	$(PYTHON) scripts/release.py clean
