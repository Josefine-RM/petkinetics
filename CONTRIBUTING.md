# Contributing to petkinetics

Thank you for your interest in contributing. This is a focused toolkit for quantitative medical imaging and kinetic modelling of dynamic PET data — contributions that fit that scope are very welcome.

## Getting started

```bash
git clone https://github.com/Josefine-RM/petkinetics.git
cd petkinetics
pip install -e ".[dev]"
```

## Code style

- Formatting: `black` (line length 88)
- Linting: `ruff`
- Type hints encouraged but not required at this stage

Run before committing:
```bash
black petkinetics/
ruff check petkinetics/
```

## Adding a function

1. Add it to the appropriate submodule (`model_fitting`, `tac_simulation`, etc.)
2. Write a numpy-style docstring with Parameters, Returns, and at least one Example
3. Export it from the submodule's `__init__.py`
4. Add a test in `tests/`

## Tests

```bash
pytest
```

## Opening issues

Please open an issue before submitting a large pull request, so we can discuss the approach first.
