# Contributing Guidelines

Thank you for your interest in contributing to renalprog! This document provides essential guidelines for contributing to the project.

!!! warning "Respect"
    Be respectful. Any issue or PR that contains disrespectful, unprofessional language, personal attacks, etc. will be disregarded and may be reported.


## Quick Start

1. **Fork** the repository
2. **Clone** your fork: `git clone https://github.com/YOUR_USERNAME/renalprog.git`
3. **Create a branch**: `git checkout -b feature/your-feature-name`
4. **Make changes** following our guidelines
5. **Test** your changes
6. **Commit**: `git commit -m "Description of changes"`
7. **Push**: `git push origin feature/your-feature-name`
8. **Open a Pull Request**

---

## Code Standards

### Python Code

- **Style**: Follow [PEP 8](https://pep8.org/)
- **Formatting**: Use `ruff` for code formatting
- **Type hints**: Add type hints to function signatures
- **Docstrings**: Use [Google-style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).

**Example:**
```python
def process_data(
    data: pd.DataFrame,
    threshold: float = 0.05
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Process input data and return results.

    Args:
        data: Input dataframe with gene expression
        threshold: Significance threshold (default: 0.05)

    Returns:
        Tuple of (processed_data, statistics_dict)

    Examples:
        >>> result, stats = process_data(df, threshold=0.01)
    """
    # Implementation
    pass
```

### R Code

- **Style**: Follow [tidyverse style guide](https://style.tidyverse.org/)
- **Documentation**: Use roxygen2-style comments
- **Naming**: Use snake_case for functions and variables

**Example:**
```r
#' Process Gene Expression Data
#'
#' @param data A data.frame with gene expression values
#' @param alpha Significance threshold (default: 0.05)
#' @return A list with processed data and statistics
#' @export
process_gene_data <- function(data, alpha = 0.05) {
  # Implementation
}
```

---

## Documentation

### Code Documentation

- **All public functions** must have docstrings
- **Include examples** when helpful
- **Document parameters** and return values
- **Note exceptions** that may be raised

### User Documentation

When adding new features, update:

- `README.md` - If changing installation or quick start
- `docs/docs/api/` - API reference for new modules/functions
- `docs/docs/tutorials/` - Tutorial if introducing new workflow
- `CHANGELOG.md` - Document changes

### Building Documentation

```bash
cd docs
mkdocs serve
mkdocs build
```

---

## Commit Messages

Use clear, descriptive commit messages following this format:

```
<type>: <short summary>

<optional detailed description>
```

**Types:**
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `style:` - Code style changes (formatting, no logic change)
- `refactor:` - Code refactoring
- `test:` - Adding or updating tests
- `chore:` - Maintenance tasks

Example:
```
fix: correct control file naming in 1_data_processing.py

Changed from plural (rnaseq_controls.csv) to singular 
(rnaseq_control.csv) to match actual file names in repository.
```

---

## Pull Request Process

### Before Submitting

- [ ] Code follows style guidelines
- [ ] Tests pass locally (`pytest`)
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] Commit messages are clear
- [ ] Branch is up to date with `main`

### PR Description Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Refactoring

## Testing
Describe testing performed

## Checklist
- [ ] Code follows style guidelines
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] All tests pass
```

### Review Process

1. **Automated checks** must pass (tests, linting)
2. **Code review** by maintainer
3. **Changes requested** may need addressing
4. **Approval** and merge by maintainer

---

## Development Setup

### Environment Setup

```bash
# Clone repository
git clone https://github.com/gprolcastelo/renalprog.git
cd renalprog

# Create conda environment
mamba env create -f environment.yml
mamba activate renalprog

# Install in editable mode
pip install -e .

# Install development dependencies
pip install pytest ruff flake8 mypy
```

---

## Code Review Guidelines

### For Contributors

- **Keep PRs focused** - One feature/fix per PR
- **Write clear descriptions** - Explain what and why
- **Respond to feedback** - Address review comments
- **Be patient** 

### For Reviewers

- **Be constructive** - Suggest improvements clearly
- **Explain reasoning** - Help contributors learn
- **Focus on substance** - Not just style
- **Acknowledge good work** - Positive feedback matters


---

## Reporting Issues

### Bug Reports

Include:
- **Description** of the bug
- **Steps to reproduce**
- **Expected behavior**
- **Actual behavior**
- **Environment** (OS, Python version, package versions)
- **Error messages** (full traceback)

**Example:**
```markdown
**Description:** VAE training fails with GPU

**Steps to reproduce:**
1. Set `force_cpu=False`
2. Run `train_vae_with_postprocessing(...)`

**Expected:** Training completes successfully

**Actual:** RuntimeError: CUDA out of memory

**Environment:**
- OS: Ubuntu 20.04
- Python: 3.9.7
- PyTorch: 1.10.0
- GPU: NVIDIA GTX 1080 (8GB)

**Error:**

RuntimeError: CUDA out of memory. Tried to allocate 2.00 GiB...
```

### Feature Requests

Include:
- **Use case** - Why is this needed?
- **Proposed solution** - How should it work?
- **Alternatives** - Other approaches considered?
- **Additional context** - Examples, references

---

## Project Structure

Understanding the project layout:

```
renalprog/
├── renalprog/          # Main package
│   ├── config.py       # Configuration
│   ├── dataset.py      # Data loading
│   ├── features.py     # Feature engineering
│   ├── enrichment.py   # Enrichment analysis
│   ├── plots.py        # Visualization
│   ├── modeling/       # Models (VAE, classification)
│   └── utils/          # Utilities
├── scripts/            # Analysis scripts
│   ├── pipeline_steps/ # Main pipeline
│   └── r_analysis/     # R scripts
├── tests/              # Test suite
├── docs/               # Documentation
├── data/               # Data directories
└── notebooks/          # Jupyter notebooks
```

---

## Questions?

- **Documentation**: Check [docs/](https://gprolcastelo.github.io/renalprog/)
- **Issues**: Search [existing issues](https://github.com/gprolcastelo/renalprog/issues)
- **Discussions**: Start a [discussion](https://github.com/gprolcastelo/renalprog/discussions)
- **Email**: Contact maintainers

---

## Code of Conduct

### Our Standards

- **Be respectful** - Treat everyone with respect
- **Be collaborative** - Work together constructively
- **Be professional** - Maintain professional conduct
- **Be inclusive** - Welcome diverse perspectives

### Unacceptable Behavior

- Harassment or discrimination
- Trolling or insulting comments
- Personal or political attacks
- Publishing others' private information

### Enforcement

Inappropriate behavior may be reported.

---

## License

By contributing, you agree that your contributions will be licensed under the Apache 2.0 License.

---

## Recognition

Contributors will be acknowledged in:
- [Acknowledgments](../acknowledgments.md)
- Release notes for significant contributions
- Documentation where appropriate

Thank you for contributing to renalprog! 🎉
