.PHONY: clean data lint test install help
.DEFAULT_GOAL := help

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = renalprog
PYTHON_INTERPRETER = python3
PYTHON_VERSION = 3.9

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Create conda environment with all dependencies (Python + R)
create-env:
	mamba env create -f environment.yml
	@echo ">>> Environment created. Activate with: mamba activate renalprog"
	@echo ">>> Then run: make install"

## Update conda environment from environment.yml
update-env:
	mamba env update -f environment.yml --prune
	@echo ">>> Environment updated."

## Install R packages via conda/mamba
install-r-conda:
	mamba install -c conda-forge r-gprofiler2 r-ggplot2 r-optparse
	@echo ">>> R packages installed via conda."

## Install package and dependencies
install:
	$(PYTHON_INTERPRETER) -m pip install -e .
	@echo ">>> Package installed. Use 'make install-dev' for development dependencies."

## Install package with development dependencies
install-dev:
	$(PYTHON_INTERPRETER) -m pip install -e ".[dev]"
	@echo ">>> Development environment installed."

## Install package with test dependencies only
install-test:
	$(PYTHON_INTERPRETER) -m pip install -e ".[test]"
	@echo ">>> Test dependencies installed."

## Install R dependencies for enrichment analysis (via Rscript)
install-r-cran:
	Rscript scripts/r_analysis/install_r_packages.R
	@echo ">>> R packages installed via CRAN."

## Install R packages (defaults to conda if available, otherwise CRAN)
install-r: install-r-conda
	@echo ">>> R packages installed. Use 'make install-r-cran' for CRAN installation."

## Install all dependencies (Python + R)
install-all: install-dev install-r
	@echo ">>> All dependencies (Python + R) installed."

## Run tests
test:
	pytest tests/ -v

## Run tests with coverage
test-coverage:
	pytest tests/ --cov=renalprog --cov-report=html --cov-report=term
	@echo ">>> Coverage report generated in htmlcov/index.html"

## Run quick tests (exclude slow tests)
test-quick:
	pytest tests/ -v -m "not slow"

## Lint code with flake8
lint:
	flake8 renalprog/ tests/

## Format code with black
format:
	black renalprog/ tests/ --line-length 100

## Sort imports with isort
isort:
	isort renalprog/ tests/ --profile black

## Run all code quality checks
quality: format isort lint
	@echo ">>> Code quality checks complete."

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name ".tox" -exec rm -rf {} +
	rm -rf build/ dist/ htmlcov/ .coverage

## Delete generated data files (keeps raw data)
clean-data:
	rm -rf data/interim/*
	rm -rf data/processed/*

## Delete trained models
clean-models:
	rm -rf models/*

## Full clean (code + data + models)
clean-all: clean clean-data clean-models
	@echo ">>> All generated files removed."

## Build documentation
docs:
	cd docs && make html
	@echo ">>> Documentation built in docs/_build/html/"

## Run preprocessing step
preprocess:
	$(PYTHON_INTERPRETER) scripts/pipeline_steps/01_preprocess.py

## Run train/test split
split:
	$(PYTHON_INTERPRETER) scripts/pipeline_steps/02_train_test_split.py

## Train VAE model
train-vae:
	$(PYTHON_INTERPRETER) scripts/pipeline_steps/03_train_vae.py

## Generate trajectories
generate-trajectories:
	$(PYTHON_INTERPRETER) scripts/pipeline_steps/06_generate_trajectories.py

## Run classification
classify:
	$(PYTHON_INTERPRETER) scripts/pipeline_steps/07_classification.py

## Run gene enrichment analysis (R)
enrichment:
	Rscript scripts/r_analysis/gene_enrichment.R

## Run full pipeline
pipeline:
	@echo ">>> Running full pipeline..."
	$(MAKE) preprocess
	$(MAKE) split
	$(MAKE) train-vae
	$(MAKE) generate-trajectories
	$(MAKE) classify
	@echo ">>> Pipeline complete!"

## Run full pipeline with enrichment analysis
pipeline-full: pipeline enrichment
	@echo ">>> Full pipeline with enrichment analysis complete!"

## Check Python environment
check-env:
	@echo "Python interpreter: $(shell which $(PYTHON_INTERPRETER))"
	@echo "Python version: $(shell $(PYTHON_INTERPRETER) --version)"
	@echo "Pip version: $(shell $(PYTHON_INTERPRETER) -m pip --version)"
	@$(PYTHON_INTERPRETER) -c "import torch; print(f'PyTorch version: {torch.__version__}')" || echo "PyTorch not installed"
	@$(PYTHON_INTERPRETER) -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')" || echo "Cannot check CUDA"

## Initialize project directories
init:
	mkdir -p data/raw data/interim data/processed data/external
	mkdir -p models notebooks reports/figures
	mkdir -p scripts/r_analysis scripts/pipeline_steps
	@echo ">>> Project directories initialized."

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



#################################################################################
# Self Documenting Commands                                                    #
#################################################################################

help:
	@echo "$$HELP_MESSAGE"

define HELP_MESSAGE

Available commands:

Environment Setup:
  make create-env       Create conda environment from environment.yml
  make update-env       Update conda environment

Installation:
  make install          Install package and dependencies
  make install-dev      Install with development dependencies
  make install-test     Install with test dependencies only
  make install-r        Install R packages (via conda, preferred)
  make install-r-conda  Install R packages via conda/mamba
  make install-r-cran   Install R packages via CRAN (Rscript)
  make install-all      Install all dependencies (Python + R)

Testing & Quality:
  make test             Run all tests
  make test-coverage    Run tests with coverage report
  make test-quick       Run quick tests only
  make lint             Lint code with flake8
  make format           Format code with black
  make isort            Sort imports
  make quality          Run all code quality checks

Cleaning:
  make clean            Remove Python artifacts
  make clean-data       Remove intermediate/processed data
  make clean-models     Remove trained models
  make clean-all        Remove all generated files

Pipeline:
  make preprocess       Run preprocessing step
  make split            Create train/test split
  make train-vae        Train VAE model
  make generate-trajectories  Generate synthetic trajectories
  make classify         Run classification
  make enrichment       Run gene enrichment analysis (R)
  make pipeline         Run Python pipeline steps
  make pipeline-full    Run full pipeline including enrichment

Documentation:
  make docs             Build documentation

Utilities:
  make check-env        Check Python environment
  make init             Initialize project directories
  make help             Show this help message

endef
export HELP_MESSAGE

