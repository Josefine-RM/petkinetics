# petkinetics

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A Python package for **kinetic modelling of dynamic PET data** — covering graphical analysis, input function handling, TAC simulation, and model selection.

A companion package to [`medimkit`](https://github.com/Josefine-RM/medimkit), which handles the upstream image loading, preprocessing, and ROI extraction needed to produce the time-activity curves this package works with.

Built and maintained by a PhD researcher in parametric imaging and kinetic modelling at Aarhus University Hospital.

---

## Who this is for

- Imaging scientists running kinetic analyses on dynamic PET data
- Researchers building reproducible quantification pipelines for pharmaceutical or clinical trials

---

## What it does

`petkinetics` is organised into five submodules, each targeting a distinct stage of a typical kinetic modelling workflow:

| Submodule | What it handles |
|---|---|
| `petkinetics.input_functions` | Extract image-derived input functions (IDIFs); generate and scale population-based input functions |
| `petkinetics.model_fitting` | Fit linearised kinetic models (Logan, Patlak) to TAC data; generate voxel-wise parametric images |
| `petkinetics.model_selection` | Compare models by Akaike Information Criterion, etc. |
| `petkinetics.tac_simulation` | Simulate time-activity curves from known micro-parameters and a given input function |
| `petkinetics.prototype_bridge` | Call proprietary collaboration executables and DLL-based prototype models (not for distribution; see `bin/README.md`) |

---

## Installation

```bash
# Standard install (once published to PyPI)
pip install petkinetics
```

Or clone and install in editable mode for development:

```bash
git clone https://github.com/Josefine-RM/petkinetics.git
cd petkinetics
pip install -e ".[dev]"
```

`petkinetics` depends on [`medimkit`](https://github.com/Josefine-RM/medimkit) for DICOM type aliases and image I/O utilities. 

**Requirements:** Python 3.9+, NumPy ≥ 1.24, SciPy ≥ 1.10, Matplotlib ≥ 3.7.

---

## Prototype model bridge

The `prototype_bridge` submodule provides Python wrappers for proprietary collaboration executable and DLL-based models that are used in research but cannot be distributed as part of this package. If you have access to these tools, place the binaries in `petkinetics/bin/` and consult `bin/README.md` for setup instructions.

---

## Background

This package grew out of daily analysis work during a PhD on **parametric imaging and kinetic modelling of dynamic PET scans**, conducted at Aarhus University Hospital in collaboration with Siemens Healthineers. The goal is a focused, well-documented set of building blocks for the specific modelling workflows that come up repeatedly in quantitative PET research and imaging studies.

---

## Roadmap

- [ ] `model_fitting`: Logan and Patlak graphical analysis
- [ ] `model_fitting`: voxel-wise parametric image generation
- [ ] `input_functions`: image-derived input function (IDIF) extraction
- [ ] `input_functions`: population-based input function scaling
- [ ] `model_selection`: AIC, BIC, and residual diagnostics
- [ ] `tac_simulation`: forward simulation from micro-parameters
- [ ] `prototype_bridge`: Siemens prototype model wrappers
- [ ] Example notebooks (ROI-level pipeline, voxel-wise parametric imaging)
- [ ] CI with GitHub Actions
- [ ] Publish to PyPI

---

## Contributing

Issues and pull requests are welcome. Please open an issue before submitting a large change. See [CONTRIBUTING.md](CONTRIBUTING.md) for code style and docstring conventions.

---

## Citation

If you use `petkinetics` in a publication, please cite it as:

```
Madsen, J. R. (2025). petkinetics: A Python package for kinetic modelling
of dynamic PET data (v0.1.0). https://github.com/Josefine-RM/petkinetics
```

---

## License

MIT — see [LICENSE](LICENSE) for details.