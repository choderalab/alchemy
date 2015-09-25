# TODO list for alchemical tools

## `alchemy.py`

### Code cleanup
* Replace `ONE_4PI_EPS0` with import from `simtk.openmm.constants` once OpenMM internal constants are available there

### Alchemical protocol specification
* Allow alchemical protocols to be specified by YAML files (for YANK)
* Eliminate default alchemical protocols in code and replace them with a library of YAML files for different kinds of protocols

### Tests
* Add tests for alchemically softening loops (for loop model refinement) and proteins (for model refinement) using no-cutoff methods
* Add thorough tests for overlap with original (unmodified) states
* Test to make sure energies are never NaN

### Forcefield support

#### PME
* Can we get around the neglect of reciprocal-space components of electrostatics?

#### Implicit solvent models
* Make support for implicit solvent models generic by rewriting force terms in a general way

#### AMOEBA
* Finish support for AMOEBA

#### Other electrstatics models
* Add support for Wolf-like electrostatics?