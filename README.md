[![Build Status](https://travis-ci.org/choderalab/alchemy.svg?branch=master)](https://travis-ci.org/choderalab/alchemy)
[![Build Status](https://jenkins.choderalab.org/buildStatus/icon?job=test-alchemy-linux)](https://jenkins.choderalab.org/job/test-alchemy-linux/)
[![Anaconda Badge](https://anaconda.org/omnia/alchemy/badges/version.svg)](https://anaconda.org/omnia/alchemy)

# Alchemical tools for OpenMM

This package contains [factories](https://en.wikipedia.org/wiki/Factory_(object-oriented_programming)) for creating [alchemically-modified](http://www.alchemistry.org/wiki/Best_Practices#Guideline_1:_Always_use_soft-core_potentials_while_decoupling_or_annihilating_Lennard-Jones_interactions) versions of OpenMM systems.

## Installation

The easiest way to install is through the [conda](http://conda.pydata.org/) package manager through the [omnia](http://omnia.md) [binstar repository](https://binstar.org/omnia/alchemy):
```bash
conda install -c omnia alchemy
```

## Examples

Create alchemical intermediates for p-xylene binding to T4 lysozyme L99A in implicit solvent.

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

## Changelog

### 1.2.1 - Critical bugfix release and overhaul of exceptions and testing
This version overhauls the way exceptions are handled, fixing a variety of issues present in earlier versions.
More thorough testing of fully-interacting and non-interacting systems is implemented, though not all tests run on travis.

### 1.2 - Expose softcore parameters as context parameters
Alchemical softcore parameters are now exposed as context parameters, and can be tweaked on the fly.
The default selection of softcore c was also changed.

### 1.1 - Allow bonds, angles, and torsions to be alchemically softened
This release allows specified bonds, angles, and torsions to be alchemically softened.
1.1.1 was a critical bugfix release.

### 1.0 - Initial release
This release breaks out `alchemy.py` from the [`yank`](http://github.com/choderalab/yank) project
