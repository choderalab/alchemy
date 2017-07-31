[![DOI](https://zenodo.org/badge/42878787.svg)](https://zenodo.org/badge/latestdoi/42878787)

# Deprecated and Moved to OpenMMTools

This software has been deprecated and its functionality has been moved to [OpenMMTools](https://github.com/choderalab/openmmtools)

This page is kept alive for legacy purposes, but the package is no longer developed here. 


## Alchemical tools for OpenMM 

This package contains [factories](https://en.wikipedia.org/wiki/Factory_(object-oriented_programming)) for creating [alchemically-modified](http://www.alchemistry.org/wiki/Best_Practices#Guideline_1:_Always_use_soft-core_potentials_while_decoupling_or_annihilating_Lennard-Jones_interactions) versions of OpenMM systems.

### Installation

The easiest way to install is through the [conda](http://conda.pydata.org/) package manager through the [omnia](http://omnia.md) [anaconda cloud repository](https://anaconda.org/omnia/alchemy):
```bash
conda install -c omnia alchemy
```

### Examples

Create alchemical intermediates for an absolute binding free energy calculation for p-xylene binding to T4 lysozyme L99A in implicit solvent:

```python
# Create a reference system for T4 lysozyme L99A
from openmmtools import testsystems
complex = testsystems.LysozymeImplicit()
[reference_system, positions] = [complex.system, complex.positions]
receptor_atoms = range(0,2603) # T4 lysozyme L99A
ligand_atoms = range(2603,2621) # p-xylene

# Create a factory to produce alchemical intermediates.
factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)

# Get the default protocol for 'denihilating' in complex in implicit solvent.
protocol = factory.defaultComplexProtocolImplicit()

# Create the perturbed systems using this protocol.
systems = factory.createPerturbedSystems(protocol)
```

### Features

#### Absolute alchemical factory (`AbsoluteAlchemicalFactory`)

Features:

* Absolute alchemical intermediates for a single specified region (can be a separate molecule, or part of a larger molecule)
* Support for vacuum, explicit solvent (currently `PME` and `Ewald`), and implicit solvent (currently only `GBSAOBCForce`)
* Exposure of alchemical parameters for softcore Lennard-Jones (`alpha, a, b, c`) and softcore electrostatics (`beta, d, e, f`).
* Separate alchemical parameter for Lennard-Jones (`lambda_sterics`), electrostatics (`lambda_electrostatics`), torsions (`lambda_torsions`), angles (`lambda_angles`), and bonds (`lambda_bonds`)

Known issues:

* PME and Ewald support do not currently include the reciprocal-space contributions from the alchemically-modified region
* Reaction field (`CutoffPeriodic`) currently subtracts a `lambda_electrostatics`-dependent value to ensure the energy is identically zero at the cutoff, which causes problems with alchemical free energy calculations utilizing reaction-field electrostatics
* Alternative electrostatics models based on `CustomNonbondedForce` are not yet supported
* Support for alchemically-modified implicit solvent models based on `CustomGBForce` has not yet been added
* Poor efficiency for softcore electrostatics has been noted

### Changelog

#### 1.2.3 - Bugfix release to update defaults
* `softcore_beta` now defaults to 0.0 because we have not yet found a good alchemcial path with `softcore_beta > 0`
* `lambda_restraints` now defaults to 1.0 to avoid the need to explicitly specify restraints in implicit solvent free energy calculations

#### 1.2.2 - Minor bugfix release
* `alpha_ewald` is now properly computed if error tolerance is zero
* better PME tests
* atom lists containing `numpy.int64` integers are now properly handled
* adds barostat to overlap tests
* adds benchmarking script

#### 1.2.1 - Critical bugfix release and overhaul of exceptions and testing
This version overhauls the way exceptions are handled, fixing a variety of issues present in earlier versions.
More thorough testing of fully-interacting and non-interacting systems is implemented, though not all tests run on travis.

#### 1.2 - Expose softcore parameters as context parameters
Alchemical softcore parameters are now exposed as context parameters, and can be tweaked on the fly.
The default selection of softcore c was also changed.

#### 1.1 - Allow bonds, angles, and torsions to be alchemically softened
This release allows specified bonds, angles, and torsions to be alchemically softened.
1.1.1 was a critical bugfix release.

#### 1.0 - Initial release
This release breaks out `alchemy.py` from the [`yank`](http://github.com/choderalab/yank) project
