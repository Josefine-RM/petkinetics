# petkinetics/bin/

This directory holds proprietary collaboration executable and DLL-based tools used
by the `petkinetics.prototype_bridge` submodule.

**This directory is listed in `.gitignore` and must never be committed to
version control.** The files here are not open-source and may not be
redistributed as part of this package.

---

## Why this directory exists

Some kinetic modelling methods used in this project are implemented as
prototype research tools — compiled executables (`.exe`) or dynamic
link libraries (`.dll`) — that are not available as open-source code. This
directory provides a stable, expected location for those binaries so that the
Python wrappers in `petkinetics.prototype_bridge` can find them at runtime via
the `_get_executable()` helper.

---

## What belongs here

| File | Description | How to obtain | Reference |
|---|---|---|---|
| exe/getPetSnake.exe | Executable deriving an input-function from the dynamic PET images. | Siemens Healthineers | Y. Tao, Z. Peng, A. Krishnan and X. S. Zhou, "Robust Learning-Based Parsing and Annotation of Medical Radiographs," in IEEE Transactions on Medical Imaging, vol. 30, no. 2, pp. 338-350, Feb. 2011, doi: 10.1109/TMI.2010.2077740. |
| exe/rawToDICOM.exe | Executable converting a .raw image to a DICOM volume | Siemens Healthineers | - |
| exe/patlak.exe | Executable calculating the net influx (uptake) rate constant Ki (ml/(min*ml)) as slope of the Patlak plot from regional PET time-activity curves. | Turku PET Centre | https://www.turkupetcentre.net/programs/doc/patlak.html |
| exe/logan.exe | Executable calculating the distribution volume (Vt) as slope of the Logan plot from regional PET time-activity curves. | Turku PET Centre | https://www.turkupetcentre.net/programs/doc/logan.html |
| exe/sim_3tcm.exe | Simulation of PET tissue time-radioactivity concentration curves (TACs) from arterial plasma (Ca) and blood (Cb) TACs, based on three-tissue compartmental model, where the compartments are in series | Turku PET Centre | https://www.turkupetcentre.net/tpcclib-doc/v2/sim__3tcm_8c_source.html |
| py/dLogan/dLogan_fits.py | .py containing code for dLogan and additional Logan prototypes | Siemens Healthineers | Madsen JR, Danielsen PB, Dias AH, Gormsen LC, Rodell AB, Panin V, Pigg D, Spottiswoode B, Munk OL. Whole-body parametric PET/CT imaging of the total distribution volume using a new reversible delayed Logan model. EJNMMI Phys. 2026 Mar 12;13(1):40. doi: 10.1186/s40658-026-00853-9. |
| py/Patlak/Patlak_fits.py | .py containing code for Patlak prototypes | Siemens Healthineers | - |
| *(to be documented as tools are integrated)* | | | |


Add a row for each binary as you integrate it into `prototype_bridge/bridge.py`.
Include the tool name, a brief description of what it computes, and where a
collaborator can obtain it.

---

## Setup on a new machine

1. Obtain the required binaries through your relevant research
   collaboration agreement.
2. Place them in this directory (`petkinetics/bin/`).
3. Verify the filenames match exactly what the wrapper functions in
   `prototype_bridge/bridge.py` pass to `_get_executable()`.
4. Test by importing the relevant wrapper and calling it on a small dataset.

If a binary is missing, `_get_executable()` will raise a `FileNotFoundError`
with a message pointing to this file.

---

## `.gitignore` rules

The following patterns in the root `.gitignore` ensure nothing in this
directory is accidentally committed:

```
petkinetics/bin/*.exe
petkinetics/bin/*.dll
```

If you add other binary formats (e.g. `.so`, `.dylib`), extend those rules
accordingly.
