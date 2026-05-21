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
|getPetSnake.exe|Executable deriving an input-function from the dynamic PET images.|Siemens Healthineers|Y. Tao, Z. Peng, A. Krishnan and X. S. Zhou, "Robust Learning-Based Parsing and Annotation of Medical Radiographs," in IEEE Transactions on Medical Imaging, vol. 30, no. 2, pp. 338-350, Feb. 2011, doi: 10.1109/TMI.2010.2077740.|
| *(to be documented as tools are integrated)* | | |

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
