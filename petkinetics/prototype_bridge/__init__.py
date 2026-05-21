"""
petkinetics.prototype_bridge — Prototype Model Bridge
======================================================
Python wrappers for proprietary collaboration executable and DLL-based kinetic
models used in research but not redistributable as part of this package.

How it works
------------
Binaries are located via :func:`_get_executable`, which looks inside the
``petkinetics/bin/`` directory at package root. That directory is listed in
``.gitignore`` and must be populated manually on each machine that needs it.
See ``bin/README.md`` for a list of required files and how to obtain them.

Usage
-----
End users of this submodule will typically call a high-level wrapper function
that handles argument marshalling, subprocess invocation, and output parsing.
The raw :func:`_get_executable` helper is internal and not part of the public
API.

Notes
-----
None of the binaries or DLLs covered by this submodule may be committed to
version control or distributed with the package. This submodule contains only
the Python-side interface code.
"""

from petkinetics.prototype_bridge.bridge import _get_executable

__all__: list[str] = []
