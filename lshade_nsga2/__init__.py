"""
Python reimplementation of the LSHADE + NSGA2 multi-objective docking framework.

This package is being ported incrementally from the original C++ codebase
(`main.cpp`, `parallel_mc.cpp`, `MODEA.cpp`, `lshade.cpp`, `nsga2.cpp`,
`everything.cpp`, etc.).  The goal is to preserve:

- Overall algorithmic structure
- File formats (PDBQT, config files, output PDBQT/log)
- Objective definitions (vina, DLIGAND2, vdw) and multi-objective search

Modules:
- config      : reading Vina-style config files
- types       : basic math and data types (vectors, boxes, enums)
- pdbqt       : PDBQT parsing and model representation (in progress)
- scoring     : scoring terms and energy evaluation (in progress)
- optimizers  : LSHADE, NSGA2, MODEA style optimizers (in progress)
- pipeline    : high-level orchestration equivalent to main_procedure + do_search
"""

from .config import DockingConfig

__all__ = ["DockingConfig"]



