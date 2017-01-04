#!/usr/bin/python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Alchemical factory for free energy calculations that operates directly on OpenMM System objects.

DESCRIPTION

This module contains enumerative factories for generating alchemically-modified System objects
usable for the calculation of free energy differences of hydration or ligand binding.

* `AbsoluteAlchemicalFactory` uses fused elecrostatic and steric alchemical modifications.

TODO

* Remove default protocol class methods, since these are no longer needed.
* Generalize treatment of nonbonded sterics/electrostatics intra-alchemical forces to support arbitrary mixing rules.
  Can we eliminate decoupling to something simpler?
* Add support for other GBSA models.
* Add functions for the automatic optimization of alchemical states?
* Can we store serialized form of Force objects so that we can save time in reconstituting
  Force objects when we make copies?  We can even manipulate the XML representation directly.
* Allow protocols to automatically be resized to arbitrary number of states, to
  allow number of states to be enlarged to be an integral multiple of number of GPUs.
* Finish AMOEBA support.
* Can alchemically-modified System objects share unmodified Force objects to avoid overhead
  of duplicating Forces that are not modified?

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import copy
import time
import itertools

import simtk.openmm as openmm
import simtk.unit as unit

import logging
logger = logging.getLogger(__name__)

#=============================================================================================
# PARAMETERS
#=============================================================================================

ONE_4PI_EPS0 = 138.935456 # OpenMM constant for Coulomb interactions (openmm/platforms/reference/include/SimTKOpenMMRealType.h) in OpenMM units
                          # TODO: Replace this with an import from simtk.openmm.constants once these constants are available there

#=============================================================================================
# MODULE UTILITIES
#=============================================================================================

def _is_periodic(system):
    """
    Determine whether the specified system is periodic.

    Parameters
    ----------
    system : simtk.openmm.System
         The system to check.

    Returns
    -------
    is_periodic : bool
        If True, the system uses a nonbonded method recognized as being periodic.

    Examples
    --------

    >>> from openmmtools import testsystems

    Periodic water box.

    >>> waterbox = testsystems.WaterBox()
    >>> print(_is_periodic(waterbox.system))
    True

    Non-periodic Lennard-Jones cluster.

    >>> cluster = testsystems.LennardJonesCluster()
    >>> print(_is_periodic(cluster.system))
    False

    Notes
    -----
    This method simply checks to see if any nonbonded methods recognized to be periodic are in use.
    Addition of new nonbonded methods or force types may require adding to recognized_periodic_methods.

    """
    forces = { system.getForce(index).__class__.__name__ for index in range(system.getNumForces()) }

    # Only the following nonbonded methods are recognized as being periodic.
    recognized_periodic_methods = {
        'NonbondedForce' : [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.Ewald, openmm.NonbondedForce.PME],
        'CustomNonbondedForce' : [openmm.CustomNonbondedForce.CutoffPeriodic],
        'AmoebaVdwForce' : [openmm.AmoebaVdwForce.CutoffPeriodic],
        }

    for force_index in range(system.getNumForces()):
        force = system.getForce(force_index)
        force_name = force.__class__.__name__
        if force_name in recognized_periodic_methods:
            method = force.getNonbondedMethod()
            if method in recognized_periodic_methods[force_name]:
                return True

    # Nothing we recognize was found, so assume system was not periodic.
    return False

#=============================================================================================
# AlchemicalState
#=============================================================================================

class AlchemicalState(dict):
    """
    Alchemical state description.

    These parameters describe the parameters that affect computation of the energy.

    Attributes
    ----------
    lambda_restraints : float
        Scaling factor for remaining receptor-ligand relative restraint terms (to help keep ligand near protein).
    lambda_electrostatics : float
        Scaling factor for ligand charges, intrinsic Born radii, and surface area term.
    lambda_sterics : float
        Scaling factor for ligand sterics (Lennard-Jones and Halgren) interactions.
    labmda_torsions : float
        Scaling factor for alchemically-softened torsions.
    labmda_angles : float
        Scaling factor for alchemically-softened angles.
    labmda_bonds : float
        Scaling factor for alchemically-softened bonds.

    """
    def __init__(self, **kwargs):
        self['lambda_restraints'] = 1.0
        self['lambda_electrostatics'] = 1.0
        self['lambda_sterics'] = 1.0
        self['lambda_torsions'] = 1.0
        self['lambda_angles'] = 1.0
        self['lambda_bonds'] = 1.0

        for key in kwargs.keys():
            # Raise an exception if we don't know how to handle a specified parameter.
            if key not in self:
                raise Exception("AlchemicalState parameter '%s' unknown" % key)

            self[key] = kwargs[key]

#=============================================================================================
# AbsoluteAlchemicalFactory
#=============================================================================================

class AbsoluteAlchemicalFactory(object):
    """
    Factory for generating OpenMM System objects that have been alchemically perturbed for absolute binding free energy calculation.

    The context parameters created are:
    * softcore_alpha - factor controlling softcore lengthscale for Lennard-Jones
    * softcore_beta - factor controlling softcore lengthscale for Coulomb
    * softcore_a - softcore Lennard-Jones parameter from Eq. 13 of Ref [1]
    * softcore_b - softcore Lennard-Jones parameter from Eq. 13 of Ref [1]
    * softcore_c - softcore Lennard-Jones parameter from Eq. 13 of Ref [1]
    * softcore_d - softcore electrostatics parameter
    * softcore_e - softcore electrostatics parameter
    * softcore_f - softcore electrostatics parameter

    Examples
    --------

    Create alchemical intermediates for default alchemical protocol for p-xylene in T4 lysozyme L99A in GBSA.

    >>> # Create a reference system.
    >>> from openmmtools import testsystems
    >>> complex = testsystems.LysozymeImplicit()
    >>> [reference_system, positions] = [complex.system, complex.positions]
    >>> # Create a factory to produce alchemical intermediates.
    >>> receptor_atoms = range(0,2603) # T4 lysozyme L99A
    >>> ligand_atoms = range(2603,2621) # p-xylene
    >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)
    >>> # Get the default protocol for 'denihilating' in complex in explicit solvent.
    >>> protocol = factory.defaultComplexProtocolImplicit()
    >>> # Create the perturbed systems using this protocol.
    >>> systems = factory.createPerturbedSystems(protocol)

    Create alchemical intermediates for default alchemical protocol for one water in a water box.

    >>> # Create a reference system.
    >>> from openmmtools import testsystems
    >>> waterbox = testsystems.WaterBox()
    >>> [reference_system, positions] = [waterbox.system, waterbox.positions]
    >>> # Create a factory to produce alchemical intermediates.
    >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])
    >>> # Get the default protocol for 'denihilating' in solvent.
    >>> protocol = factory.defaultSolventProtocolExplicit()
    >>> # Create the perturbed systems using this protocol.
    >>> systems = factory.createPerturbedSystems(protocol)

    Alchemically modify some angles and torsions in alanine dipeptide

    >>> # Create an alchemically-perturbed test system
    >>> from openmmtools import testsystems
    >>> testsystem = testsystems.AlanineDipeptideVacuum()
    >>> from alchemy import AbsoluteAlchemicalFactory
    >>> factory = AbsoluteAlchemicalFactory(testsystem.system, ligand_atoms=[0], alchemical_torsions=[0,1,2], alchemical_angles=[0,1,2], annihilate_sterics=True, annihilate_electrostatics=True)
    >>> # Create an alchemically-perturbed system.
    >>> alchemical_system = factory.createPerturbedSystem()
    >>> # Create a Context to make sure this works.
    >>> integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    >>> context = openmm.Context(alchemical_system, integrator)
    >>> del context

    Alchemically modify a bond, angles, and torsions in toluene by automatically selecting bonds involving alchemical atoms.

    >>> # Create an alchemically-perturbed test system.
    >>> from openmmtools import testsystems
    >>> testsystem = testsystems.TolueneImplicit()
    >>> from alchemy import AbsoluteAlchemicalFactory
    >>> factory = AbsoluteAlchemicalFactory(testsystem.system, ligand_atoms=[0,1], alchemical_torsions=True, alchemical_angles=True, annihilate_sterics=True, annihilate_electrostatics=True)
    >>> # Create an alchemically-perturbed system.
    >>> alchemical_system = factory.createPerturbedSystem()
    >>> # Create a Context to make sure this works.
    >>> integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    >>> context = openmm.Context(alchemical_system, integrator)
    >>> del context

    References
    ----------
    [1] Pham TT and Shirts MR. Identifying low variance pathways for free energy calculations of molecular transformations in solution phase.
    JCP 135:034114, 2011. http://dx.doi.org/10.1063/1.3607597

    """

    # Factory initialization.
    def __init__(self, reference_system, ligand_atoms=None, receptor_atoms=None,
                 alchemical_torsions=None, alchemical_angles=None, alchemical_bonds=None,
                 annihilate_electrostatics=True, annihilate_sterics=False,
                 softcore_alpha=0.5, softcore_beta=0.0, softcore_a=1, softcore_b=1, softcore_c=6, softcore_d=1, softcore_e=1, softcore_f=2,
                 alchemical_functions=None,
                 test_positions=None, platform=None, consistent_exceptions=False):
        """
        Initialize absolute alchemical intermediate factory with reference system.

        The reference system will not be modified when alchemical intermediates are generated.

        Parmeters
        ---------
        reference_system : simtk.openmm.System
            The reference system that is to be alchemically modified.
        ligand_atoms : list of int, optional, default = None
            List of atoms to be designated as 'ligand' for alchemical modification; everything else in system is considered the 'environment'.
        receptor_atoms : list of int, optional, default = None
            List of atoms to be designated as 'receptor' for alchemical modification.
            UNUSED; WILL BE DEPRECATED
        alchemical_torsions : list of int, optional, default = None
            If a list of torsion indices are specified, these PeriodicTorsionForce entries are softened with 'lambda_torsions'.
            If set to True, this list is autogenerated to include al proper torsions involving any alchemical atoms; improper torsions are not softened.
        alchemical_angles : list of int, optional, default = None
            If a list of angle indices are specified, these HarmonicAngleForce entries are softened with 'lambda_angles'.
            If set to True, this list is autogenerated to include all angles involving any alchemical atoms.
        alchemical_bonds : list of int, optional, default = None
            If a list of bond indices are specified, these HarmonicBondForce entries are softened with 'lambda_bonds'.
            If set to True, this list is autogenerated to include all bonds involving any alchemical atoms.
        annihilateElectrostatics : bool
            If True, electrostatics should be annihilated, rather than decoupled.
        annihilateSterics : bool
            If True, sterics (Lennard-Jones or Halgren potential) will be annihilated, rather than decoupled.
        softcore_alpha : float, optional, default = 0.5
            Alchemical softcore parameter for Lennard-Jones.
        softcore_beta : float, optional, default = 0.0
            Alchemical softcore parameter for electrostatics.
            Set this to zero to recover standard electrostatic scaling.
        softcore_a, softcore_b, softcore_c : float, optional, default=1
            Parameters modifying softcore Lennard-Jones form.
            Introduced in Eq. 13 of Ref. [1]
        softcore_d, softcore_e, softcore_f : float, optional, default=1
            Parameters modifying softcore electrostatics form.
            r_eff = sigma*((softcore_beta*(lambda_electrostatics-1)^softcore_e + (r/sigma)^softcore_f))^(1/softcore_f)
        alchemical_functions : dict, optional, default=None
            If not None, this dict specifies a mapping from context parameters to one or more globally-controlled parameters.
            This allows groups of alchemical parameters to be slaved to one or more global context parameters.
            To slave everything to a single lambda:
              alchemical_functions = { 'lambda_sterics' : 'lambda', 'lambda_electrostatics' : 'lambda', 'lambda_bonds' : 'lambda', 'lambda_angles' : 'lambda', 'lambda_torsions' : 'lambda' }
            For a two-stage function:
              alchemical_functions = { 'lambda_sterics' : '2*lambda * step(0.5 - lambda)', 'lambda_electrostatics' : '2*(lambda - 0.5) * step(lambda - 0.5)' }
        test_positions : simtk.unit.Quantity of dimension (natoms,3) with units compatible with nanometers, optional, default=None
            If provided, these coordinates will be used to test alchemically-modified system to ensure the potential energy is finite.
            If the potential energy is NaN, the energy for each force component will be computed for the Reference platform to aid in debugging.
        platform : simtk.openmm.Platform, optionl default=None
            If provided, this Platform will be used to check energies are finite.
        consistent_exceptions : bool, optional, default = False
            If True, the same functional form of the System's nonbonded method will be use to determine
            the electrostatics contribution to the potential energy of 1,4 exceptions instead of the
            classical q1*q2/(4*epsilon*epsilon0*pi*r).

        TODO:
        * Can we use a Topology object to simplify this?
        * Can we replace ligand_atoms and receptor_atoms with just alchemical_atoms?
        * Can we specify multiple groups of alchemically-modified atoms that have different alchemical parameters associated with them, like `lambda_sterics_[groupname]`?
        * Can we collect related parameters (e.g. softcore parameters) into a dict?

        """
        # If no ligand atom set is specified, create an empty list.
        if ligand_atoms is None:
            ligand_atoms = list()

        # Ligand atoms must be list of numpy int type
        ligand_atoms = [ int(index) for index in ligand_atoms ]

        # Store annihilation/decoupling information.
        self.annihilate_electrostatics = annihilate_electrostatics
        self.annihilate_sterics = annihilate_sterics
        self.softcore_alpha = softcore_alpha
        self.softcore_beta = softcore_beta
        self.softcore_a = softcore_a
        self.softcore_b = softcore_b
        self.softcore_c = softcore_c
        self.softcore_d = softcore_d
        self.softcore_e = softcore_e
        self.softcore_f = softcore_f
        self.alchemical_functions = alchemical_functions
        if self.alchemical_functions == None:
            self.alchemical_functions = dict()
        self.consistent_exceptions = consistent_exceptions

        # Store serialized form of reference system.
        self.reference_system = copy.deepcopy(reference_system)

        # Store reference forces.
        self.reference_forces = { self.reference_system.getForce(index).__class__.__name__ : self.reference_system.getForce(index) for index in range(self.reference_system.getNumForces()) }

        # Store copy of atom sets.
        all_particles_set = set(range(reference_system.getNumParticles()))
        if not (set(ligand_atoms).issubset(all_particles_set)):
            msg  = 'Some specified ligand atom indices >= number of particles (%d)\n' % reference_system.getNumParticles()
            msg += 'These specified atoms are not in the system: %s\n' % str(set(ligand_atoms).difference(all_particles_set))
            raise Exception(msg)
        self.ligand_atoms = copy.deepcopy(ligand_atoms)

        # Store atom sets
        self.ligand_atomset = set(self.ligand_atoms)

        # Store specified lists of alchemical bonds, angles, and torsions to soften (or None).
        self.alchemical_bonds = alchemical_bonds
        self.alchemical_angles = alchemical_angles
        self.alchemical_torsions = alchemical_torsions

        # If True was specified, build lists of bonds, angles, or torsions involving alchemical atoms.
        if self.alchemical_bonds is True:
            self.alchemical_bonds = self._buildAlchemicalBondList(self.ligand_atomset)
        if self.alchemical_angles is True:
            self.alchemical_angles = self._buildAlchemicalAngleList(self.ligand_atomset)
        if self.alchemical_torsions is True:
            self.alchemical_torsions = self._buildAlchemicalTorsionList(self.ligand_atomset)

        # Create an alchemically-modified system to cache
        [self.alchemically_modified_system, self.force_labels] = self._createAlchemicallyModifiedSystem(self.reference_system)

        # Build a list of all alchemical parameters available in this system.
        self.alchemical_parameters = self._getSystemGlobalParameters(self.alchemically_modified_system, prefix='lambda_')

        # Store information for use in aiding debugging of alchemical factory
        self.test_positions = test_positions
        self.platform = platform
        if self.test_positions is not None:
            self._checkEnergyIsFinite(self.alchemically_modified_system, test_positions, platform=platform)

        # DEBUG: Write XML
        debug_write_xml = False
        if debug_write_xml:
            import os, os.path
            def write_file(filename, contents):
                with open(filename, 'w') as outfile:
                    outfile.write(contents)

            logger.info("Serializing to XML...")
            system_filename = os.path.join('setup', 'alchemical-system.xml')
            write_file(system_filename, openmm.XmlSerializer.serialize(self.alchemically_modified_system))

        return

    @classmethod
    def _getSystemGlobalParameters(cls, system, prefix='lambda_'):
        """
        Get a list of available alchemical parameters beginning with 'prefix'

        Parameters
        ----------
        prefix : str, optional, default='lambda_'
            The prefix that all Context parameters found must begin with.

        """
        platform = openmm.Platform.getPlatformByName('Reference')
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        context = openmm.Context(system, integrator, platform)
        alchemical_parameters = [ parameter for parameter in context.getState(getParameters=True).getParameters().keys() if parameter.startswith(prefix) ]
        del context, integrator
        return alchemical_parameters

    def NoninteractingAlchemicalState(self):
        """
        Return a noninteracting alchemical state where all alchemical parmeters are set to zero.
        """
        kwargs = { parameter : 0 for parameter in self.alchemical_parameters }
        return AlchemicalState(**kwargs)

    def FullyInteractingAlchemicalState(self):
        """
        Return a noninteracting alchemical state where all alchemical parmeters are set to unity.
        """
        kwargs = { parameter : 1 for parameter in self.alchemical_parameters }
        return AlchemicalState(**kwargs)

    @classmethod
    def _tabulateBonds(cls, system):
        """
        Tabulate bonds for the specified system.

        Parameters
        ----------
        system : simtk.openmm.System
            The system for which bonds are to be tabulated.

        Returns
        -------
        bonds : list of set
            bonds[i] is the set of bonds to atom i

        TODO:
        * Could we use a Topology object to simplify this?

        """
        bonds = [ set() for particle_index in range(system.getNumParticles()) ]

        forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }

        # Process HarmonicBondForce
        bond_force = forces['HarmonicBondForce']
        for bond_index in range(bond_force.getNumBonds()):
            [particle1, particle2, r, K] = bond_force.getBondParameters(bond_index)
            bonds[particle1].add(particle2)
            bonds[particle2].add(particle1)
        # Process constraints.
        for constraint_index in range(system.getNumConstraints()):
            [particle1, particle2, r] = system.getConstraintParameters(constraint_index)
            bonds[particle1].add(particle2)
            bonds[particle2].add(particle1)

        # TODO: Process CustomBondForce?

        return bonds

    def _buildAlchemicalTorsionList(self, alchemical_atomset):
        """
        Build a list of proper torsion indices that involve any alchemical atom.

        Parameters
        ----------
        alchemical_atomset : set of int
            The set of alchemically modified atoms

        Returns
        -------
        torsion_list : list of int
            The list of torsion indices that should be alchemically softened

        """

        # Tabulate all bonds
        bonds = self._tabulateBonds(self.reference_system)
        def is_bonded(i,j):
            if j in bonds[i]:
                return True
            return False
        def is_proper_torsion(i,j,k,l):
            if is_bonded(i,j) and is_bonded(j,k) and is_bonded(k,l):
                return True
            return False

        # Create a list of proper torsions that involve any alchemical atom.
        torsion_list = list()
        force = self.reference_forces['PeriodicTorsionForce']
        for torsion_index in range(force.getNumTorsions()):
            [particle1, particle2, particle3, particle4, periodicity, phase, k] = force.getTorsionParameters(torsion_index)
            if set([particle1,particle2,particle3,particle4]).intersection(alchemical_atomset):
                if is_proper_torsion(particle1,particle2,particle3,particle4):
                    torsion_list.append(torsion_index)

        return torsion_list

    def _buildAlchemicalAngleList(self, alchemical_atomset):
        """
        Build a list of angle indices that involve any alchemical atom.

        Parameters
        ----------
        alchemical_atomset : set of int
            The set of alchemically modified atoms

        Returns
        -------
        angle_list : list of int
            The list of angle indices that should be alchemically softened

        """
        angle_list = list()
        force = self.reference_forces['HarmonicAngleForce']
        for angle_index in range(force.getNumAngles()):
            [particle1, particle2, particle3, theta0, K] = force.getAngleParameters(angle_index)
            if set([particle1,particle2,particle3]).intersection(alchemical_atomset):
                angle_list.append(angle_index)

        return angle_list

    def _buildAlchemicalBondList(self, alchemical_atomset):
        """
        Build a list of bond indices that involve any alchemical atom, allowing a list of bonds to override.

        Parameters
        ----------
        alchemical_atomset : set of int
            The set of alchemically modified atoms

        Returns
        -------
        bond_list : list of int
            The list of bond indices that should be alchemically softened

        """
        bond_list = list()
        force = self.reference_forces['HarmonicBondForce']
        for bond_index in range(force.getNumBonds()):
            [particle1, particle2, r, K] = force.getBondParameters(bond_index)
            if set([particle1,particle2]).intersection(alchemical_atomset):
                bond_list.append(bond_index)

        return bond_list

    def _checkEnergyIsFinite(self, system, positions, platform=None):
        """
        Test that the potential energy is finite for the provided positions.
        If energy is NaN, compute energy from each force component on Reference platform to aid in debugging.

        Parameters
        ----------
        system : simtk.openmm.System
             The alchemical system to check.
        positions : simtk.unit.Quantity of dimension (natoms,3) with units compatible with nanometers
            Coordinates to use for energy test.
        platform : simtk.openmm.Platform, optional, default=None
            If specified, this platform will be used to compute test energy (but not by force component).

        """
        logger.debug("Checking alchemical system produces finite energies.")

        def compute_potential_energy(system, positions, platform=None, groups=-1):
            # Compute potential energy of reference system.
            integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
            if platform is not None:
                # Compute potential energy with requested platform.
                context = openmm.Context(system, integrator, platform)
            else:
                # Use fastest available platform.
                context = openmm.Context(system, integrator)
            context.setPositions(positions)
            potential = context.getState(getEnergy=True,groups=groups).getPotentialEnergy()
            del context, integrator
            return potential

        # Compute potential energy error between reference and alchemical system.
        reference_potential = compute_potential_energy(self.reference_system, positions, platform)
        alchemical_potential = compute_potential_energy(system, positions, platform)
        energy_error = alchemical_potential - reference_potential

        # If alchemical potential energy is NaN, compute energy by component on Reference platform.
        if np.isnan(alchemical_potential / unit.kilocalories_per_mole):
            logger.debug("Energy for alchemically modified system is NaN.")
            logger.debug("Decomposing energy by force component on Reference platform:")

            # Assign all forces to separate components.
            system = self.alchemically_modified_system
            force_groups = list()
            for (force_index, force) in enumerate(system.getForces()):
                force_groups.append(force.getForceGroup())
                force.setForceGroup(force_index)

            # Compute potential energy for each force using Reference platform.
            reference_platform = openmm.Platform.getPlatformByName('Reference')
            for force_index in range(system.getNumForces()):
                groups = 1 << force_index # group index bitwise selector
                potential = compute_potential_energy(system, positions, reference_platform, groups)
                force_classname = system.getForce(force_index).__class__.__name__
                logger.debug("Force %5d / %5d [%24s] %12.3f kcal/mol" % (force_index, system.getNumForces(), force_classname, potential / unit.kilocalories_per_mole))

            for (force, force_group) in zip(system.getForces(), force_groups):
                force.setForceGroup(force_index)

            # Clean up
            del context, integrator

            raise Exception("Energy for alchemically modified system is NaN.")

        # Return the energy error.
        logger.debug("Difference between alchemical and reference potential energy is %8.3f kcal/mol" % (energy_error / unit.kilocalories_per_mole))
        return energy_error

    @classmethod
    def defaultComplexProtocolImplicit(cls):
        """
        Return the default protocol for 'denihilating' a ligand in complex with a protein in implicit solvent.

        Returns
        -------
        alchemical_states : list of AlchemicalState
            The list of alchemical states defining the protocol.

        Notes
        -----
        The unrestrained, fully interacting system is always listed first.

        """

        alchemical_states = list()
        lambda_values = [1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.88, 0.86, 0.84, 0.81,
                         0.78, 0.74, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025, 0.01, 0.00]

        for lambda_value in lambda_values:
            alchemical_state = AlchemicalState()
            alchemical_state['lambda_restraints'] = 1.0
            alchemical_state['lambda_electrostatics'] = lambda_value
            alchemical_state['lambda_sterics'] = lambda_value
            alchemical_states.append(alchemical_state)

        return alchemical_states

    @classmethod
    def defaultComplexProtocolExplicit(cls):
        """
        Return the default protocol for 'denihilating' a ligand in complex with a protein in explicit solvent.

        Returns
        -------
        alchemical_states : list of AlchemicalState
            The list of alchemical states defining the protocol.

        Notes
        -----
        The unrestrained, fully interacting system is always listed first.

        TODO
        ----
        * Update this with optimized set of alchemical states.

        """

        alchemical_states = list()

        lambda_values = [1.0, 0.97, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025, 0.01, 0.0]

        for lambda_value in lambda_values:
            alchemical_state = AlchemicalState()
            alchemical_state['lambda_electrostatics'] = lambda_value
            alchemical_state['lambda_sterics'] = lambda_value
            alchemical_states.append(alchemical_state)

        return alchemical_states

    @classmethod
    def defaultSolventProtocolImplicit(cls):
        """
        Return the default protocol for 'denihilating' a ligand in imlicit solvent.

        Returns
        -------
        alchemical_states : list of AlchemicalState
            The list of alchemical states defining the protocol.

        Notes
        -----
        The unrestrained, fully interacting system is always listed first.

        TODO
        ----
        * Update this with optimized set of alchemical states.

        """

        alchemical_states = list()

        lambda_values = [1.0, 0.97, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025, 0.01, 0.0]

        for lambda_value in lambda_values:
            alchemical_state = AlchemicalState()
            alchemical_state['lambda_electrostatics'] = lambda_value
            alchemical_state['lambda_sterics'] = lambda_value
            alchemical_states.append(alchemical_state)

        return alchemical_states

    @classmethod
    def defaultVacuumProtocol(cls):
        """
        Return the default protocol for 'denihilating' a ligand in vacuum.

        Returns
        -------
        alchemical_states : list of AlchemicalState
            The list of alchemical states defining the protocol.

        Notes
        -----
        The unrestrained, fully interacting system is always listed first.

        """

        alchemical_states = list()

        lambda_values = [1.0, 0.0]

        for lambda_value in lambda_values:
            alchemical_state = AlchemicalState()
            alchemical_state['lambda_electrostatics'] = lambda_value
            alchemical_state['lambda_sterics'] = lambda_value
            alchemical_states.append(alchemical_state)

        return alchemical_states

    @classmethod
    def defaultSolventProtocolExplicit(cls):
        """
        Return the default protocol for 'denihilating' a ligand in explicit solvent.

        Returns
        -------
        alchemical_states : list of AlchemicalState
            The list of alchemical states defining the protocol.

        Notes
        -----
        The unrestrained, fully interacting system is always listed first.

        TODO
        ----
        * Update this with optimized set of alchemical states.

        """

        alchemical_states = list()

        lambda_values = [1.0, 0.97, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025, 0.01, 0.00]

        for lambda_value in lambda_values:
            alchemical_state = AlchemicalState()
            alchemical_state['lambda_electrostatics'] = lambda_value
            alchemical_state['lambda_sterics'] = lambda_value
            alchemical_states.append(alchemical_state)

        return alchemical_states

    def _alchemicallyModifyPeriodicTorsionForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of PeriodicTorsionForce.

        Parameters
        ----------
        system : simtk.openmm.System
            The new alchemically-modified System object being built.  This object will be modified.
        reference_force : simtk.openmm.PeriodicTorsionForce
            The reference copy of the PeriodicTorsionForce to be alchemically-modified.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        """

        # Create PeriodicTorsionForce to handle unmodified torsions.
        force = openmm.PeriodicTorsionForce()

        # Create CustomTorsionForce to handle alchemically modified torsions.
        energy_function = "lambda_torsions*k*(1+cos(periodicity*theta-phase))"
        custom_force = openmm.CustomTorsionForce(energy_function)
        custom_force.addGlobalParameter('lambda_torsions', 1.0)
        custom_force.addPerTorsionParameter('periodicity')
        custom_force.addPerTorsionParameter('phase')
        custom_force.addPerTorsionParameter('k')
        # Process reference torsions.
        for torsion_index in range(reference_force.getNumTorsions()):
            # Retrieve parameters.
            [particle1, particle2, particle3, particle4, periodicity, phase, k] = reference_force.getTorsionParameters(torsion_index)
            # Create torsions.
            if set([particle1,particle2,particle3,particle4]).issubset(self.ligand_atomset):
                # Alchemically modified torsion.
                custom_force.addTorsion(particle1, particle2, particle3, particle4, [periodicity, phase, k])
            else:
                # Standard torsion.
                force.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)

        # Add newly-populated forces to system.
        force_index = system.addForce(force)
        force_labels['unmodified PeriodicTorsionForce'] = force_index

        force_index = system.addForce(custom_force)
        force_labels['alchemically modified PeriodicTorsionForce'] = force_index

    def _alchemicallyModifyHarmonicAngleForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of HarmonicAngleForce

        Parameters
        ----------
        system : simtk.openmm.System
            The new alchemically-modified System object being built.  This object will be modified.
        reference_force : simtk.openmm.HarmonicAngleForec
            The reference copy of the HarmonicAngleForce to be alchemically-modified.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        """

        # Create standard HarmonicAngleForce to handle unmodified angles.
        force = openmm.HarmonicAngleForce()

        # Create CustomAngleForce to handle alchemically modified angles.
        energy_function = "lambda_angles*(K/2)*(theta-theta0)^2;"
        custom_force = openmm.CustomAngleForce(energy_function)
        custom_force.addGlobalParameter('lambda_angles', 1.0)
        custom_force.addPerAngleParameter('theta0')
        custom_force.addPerAngleParameter('K')
        # Process reference torsions.
        for angle_index in range(reference_force.getNumAngles()):
            # Retrieve parameters.
            [particle1, particle2, particle3, theta0, K] = reference_force.getAngleParameters(angle_index)
            if angle_index in self.alchemical_angles:
                # Alchemically modified angle.
                custom_force.addAngle(particle1, particle2, particle3, [theta0, K])
            else:
                # Standard torsion.
                force.addAngle(particle1, particle2, particle3, theta0, K)

        # Add newly-populated forces to system.
        force_index = system.addForce(force)
        force_labels['unmodified HarmonicAngleForce'] = force_index

        force_index = system.addForce(custom_force)
        force_labels['alchemically modified HarmonicAngleForce'] = force_index

    def _alchemicallyModifyHarmonicBondForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of HarmonicBondForce

        Parameters
        ----------
        system : simtk.openmm.System
            The new alchemically-modified System object being built.  This object will be modified.
        reference_force : simtk.openmm.HarmonicBondForec
            The reference copy of the HarmonicBondForce to be alchemically-modified.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        """

        # Create standard HarmonicBondForce to handle unmodified bonds.
        force = openmm.HarmonicBondForce()

        # Create CustomBondForce to handle alchemically modified bonds.
        energy_function = "lambda_bonds*(K/2)*(r-r0)^2;"
        custom_force = openmm.CustomBondForce(energy_function)
        custom_force.addGlobalParameter('lambda_bonds', 1.0)
        custom_force.addPerBondParameter('r0')
        custom_force.addPerBondParameter('K')
        # Process reference torsions.
        for bond_index in range(reference_force.getNumBonds()):
            # Retrieve parameters.
            [particle1, particle2, theta0, K] = reference_force.getBondParameters(bond_index)
            if bond_index in self.alchemical_bonds:
                # Alchemically modified torsion.
                custom_force.addBond(particle1, particle2, [theta0, K])
            else:
                # Standard torsion.
                force.addBond(particle1, particle2, theta0, K)

        # Add newly-populated forces to system.
        force_index = system.addForce(force)
        force_labels['unmodified HarmonicBondForce'] = force_index

        force_index = system.addForce(custom_force)
        force_labels['alchemically modified HarmonicBondForce'] = force_index

    def _alchemicallyModifyNonbondedForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of NonbondedForce.

        Parameters
        ----------
        system : simtk.openmm.System
            Alchemically-modified system being built.  This object will be modified.
        nonbonded_force : simtk.openmm.NonbondedForce
            The NonbondedForce used as a template.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        TODO
        ----
        Change softcore_beta to a dimensionless scalar to multiply some intrinsic length-scale, like Lennard-Jones alpha.
        Try using a single, common "reff" effective softcore distance for both Lennard-Jones and Coulomb.

        References
        ----------
        [1] Pham TT and Shirts MR. Identifying low variance pathways for free energy calculations of molecular transformations in solution phase.
        JCP 135:034114, 2011. http://dx.doi.org/10.1063/1.3607597

        """

        nonbonded_method = reference_force.getNonbondedMethod()

        # --------------------------------------------------
        # Determine energy expression for all custom forces
        # --------------------------------------------------

        # Form energy expression for slaved context parameters.
        alchemical_function_expression = ""
        for variable in self.alchemical_functions:
            expression = self.alchemical_functions[variable]
            alchemical_function_expression = " %s = %s;" % (variable, expression)

        # Soft-core Lennard-Jones
        sterics_energy_expression = "U_sterics = (lambda_sterics^softcore_a)*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;"

        # Electrostatics energy expression for NoCutoff nonbonded method
        nocutoff_electrostatics_energy_expression = "U_electrostatics = (lambda_electrostatics^softcore_d)*ONE_4PI_EPS0*chargeprod/reff_electrostatics;"

        # Electrostatics energy expression will change according to the nonbonded method.
        electrostatics_energy_expression = ""

        # Energy expression for 1,4 electrostatics exceptions can be either electrostatics_energy_expression
        # or nocutoff_electrostatics_energy_expression if self.consistent_exceptions is True or False respectively
        exceptions_electrostatics_energy_expression = ""

        # Select electrostatics functional form based on nonbonded method.
        if nonbonded_method in [openmm.NonbondedForce.NoCutoff]:
            # soft-core Coulomb
            electrostatics_energy_expression += nocutoff_electrostatics_energy_expression
        elif nonbonded_method in [openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic]:
            # reaction-field electrostatics
            epsilon_solvent = reference_force.getReactionFieldDielectric()
            r_cutoff = reference_force.getCutoffDistance()
            electrostatics_energy_expression += "U_electrostatics = (lambda_electrostatics^softcore_d)*ONE_4PI_EPS0*chargeprod*(reff_electrostatics^(-1) + k_rf*reff_electrostatics^2 - c_rf);"
            k_rf = r_cutoff**(-3) * ((epsilon_solvent - 1) / (2*epsilon_solvent + 1))
            c_rf = r_cutoff**(-1) * ((3*epsilon_solvent) / (2*epsilon_solvent + 1))
            electrostatics_energy_expression += "k_rf = %f;" % (k_rf.value_in_unit_system(unit.md_unit_system))
            electrostatics_energy_expression += "c_rf = %f;" % (c_rf.value_in_unit_system(unit.md_unit_system))
        elif nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
            # Ewald direct-space electrostatics
            [alpha_ewald, nx, ny, nz] = reference_force.getPMEParameters()
            if (alpha_ewald/alpha_ewald.unit) == 0.0:
                # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance.
                tol = reference_force.getEwaldErrorTolerance()
                alpha_ewald = (1.0/reference_force.getCutoffDistance()) * np.sqrt(-np.log(2.0*tol))
            electrostatics_energy_expression += "U_electrostatics = (lambda_electrostatics^softcore_d)*ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*reff_electrostatics)/reff_electrostatics;"
            electrostatics_energy_expression += "alpha_ewald = %f;" % (alpha_ewald.value_in_unit_system(unit.md_unit_system))
            # TODO: Handle reciprocal-space electrostatics for alchemically-modified particles.  These are otherwise neglected.
            # NOTE: There is currently no way to do this in OpenMM.
        else:
            raise Exception("Nonbonded method %s not supported yet." % str(nonbonded_method))

        # Add additional definitions common to all methods.
        sterics_energy_expression += "reff_sterics = sigma*((softcore_alpha*(1.0-lambda_sterics)^softcore_b + (r/sigma)^softcore_c))^(1/softcore_c);" # effective softcore distance for sterics
        reff_electrostatics_expression = "reff_electrostatics = sigma*((softcore_beta*(1.0-lambda_electrostatics)^softcore_e + (r/sigma)^softcore_f))^(1/softcore_f);" # effective softcore distance for electrostatics
        reff_electrostatics_expression += "ONE_4PI_EPS0 = %f;" % ONE_4PI_EPS0 # already in OpenMM units
        electrostatics_energy_expression += reff_electrostatics_expression

        # Define functional form of 1,4 electrostatic exceptions
        if self.consistent_exceptions:
            exceptions_electrostatics_energy_expression += electrostatics_energy_expression
        else:
            exceptions_electrostatics_energy_expression += nocutoff_electrostatics_energy_expression
            exceptions_electrostatics_energy_expression += reff_electrostatics_expression

        # Define mixing rules.
        sterics_mixing_rules = ""
        sterics_mixing_rules += "epsilon = sqrt(epsilon1*epsilon2);" # mixing rule for epsilon
        sterics_mixing_rules += "sigma = 0.5*(sigma1 + sigma2);" # mixing rule for sigma
        electrostatics_mixing_rules = ""
        electrostatics_mixing_rules += "chargeprod = charge1*charge2;" # mixing rule for charges
        electrostatics_mixing_rules += "sigma = 0.5*(sigma1 + sigma2);" # mixing rule for sigma

        # ------------------------------------------------------------
        # Create and configure all forces to add to alchemical system
        # ------------------------------------------------------------

        # Interactions and exceptions will be distributed according to the following table

        # --------------------------------------------------------------------------------------------------
        # FORCE                                    | INTERACTION GROUP                                     |
        # --------------------------------------------------------------------------------------------------
        # nonbonded_force (unmodified)             | all interactions nonalchemical/nonalchemical          |
        #                                          | all exceptions nonalchemical/nonalchemical            |
        # --------------------------------------------------------------------------------------------------
        # aa_sterics_custom_nonbonded_force        | sterics interactions alchemical/alchemical            |
        # --------------------------------------------------------------------------------------------------
        # aa_electrostatics_custom_nonbonded_force | electrostatics interactions alchemical/alchemical     |
        # --------------------------------------------------------------------------------------------------
        # na_sterics_custom_nonbonded_force        | sterics interactions non-alchemical/alchemical        |
        # --------------------------------------------------------------------------------------------------
        # na_electrostatics_custom_nonbonded_force | electrostatics interactions non-alchemical/alchemical |
        # --------------------------------------------------------------------------------------------------
        # aa_sterics_custom_bond_force             | sterics exceptions alchemical/alchemical              |
        # --------------------------------------------------------------------------------------------------
        # aa_electrostatics_custom_bond_force      | electrostatics exceptions alchemical/alchemical       |
        # --------------------------------------------------------------------------------------------------
        # na_sterics_custom_bond_force             | sterics exceptions non-alchemical/alchemical          |
        # --------------------------------------------------------------------------------------------------
        # na_electrostatics_custom_bond_force      | electrostatics exceptions non-alchemical/alchemical   |
        # --------------------------------------------------------------------------------------------------

        # Create a copy of the NonbondedForce to handle particle interactions and
        # 1,4 exceptions between non-alchemical/non-alchemical atoms (nn).
        nonbonded_force = copy.deepcopy(reference_force)
        force_index = system.addForce(nonbonded_force)
        force_labels['unmodified NonbondedForce'] = force_index

        # Create CustomNonbondedForces to handle sterics particle interactions between
        # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
        # to 1.0 for decoupled interactions
        basic_sterics_expression = "U_sterics;" + sterics_energy_expression +\
                                   sterics_mixing_rules + alchemical_function_expression

        na_sterics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_sterics_expression)
        na_sterics_custom_nonbonded_force.addGlobalParameter("lambda_sterics", 1.0)
        if self.annihilate_sterics:
            aa_sterics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_sterics_expression)
            aa_sterics_custom_nonbonded_force.addGlobalParameter("lambda_sterics", 1.0)
        else:  # for decoupling fix lambda_sterics to 1.0
            aa_sterics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_sterics_expression +
                                                                            'lambda_sterics=1.0;')

        # Add parameters and configure CustomNonbondedForces to match reference force
        is_method_periodic = nonbonded_method in [openmm.NonbondedForce.Ewald, openmm.NonbondedForce.PME,
                                                  openmm.NonbondedForce.CutoffPeriodic]

        for force in [na_sterics_custom_nonbonded_force, aa_sterics_custom_nonbonded_force]:
            force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
            force.addPerParticleParameter("epsilon")  # Lennard-Jones epsilon
            force.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction())
            force.setCutoffDistance(nonbonded_force.getCutoffDistance())
            force.setSwitchingDistance(nonbonded_force.getSwitchingDistance())
            force.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())

            if is_method_periodic:
                force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            else:
                force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

        # Create CustomNonbondedForces to handle electrostatics particle interactions between
        # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
        # to 1.0 for decoupled interactions
        basic_electrostatics_expression = "U_electrostatics;" + electrostatics_energy_expression +\
                                          electrostatics_mixing_rules + alchemical_function_expression

        na_electrostatics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_electrostatics_expression)
        na_electrostatics_custom_nonbonded_force.addGlobalParameter("lambda_electrostatics", 1.0)
        if self.annihilate_electrostatics:
            aa_electrostatics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_electrostatics_expression)
            aa_electrostatics_custom_nonbonded_force.addGlobalParameter("lambda_electrostatics", 1.0)
        else:  # for decoupling fix lambda_electrostatics to 1.0
            aa_electrostatics_custom_nonbonded_force = openmm.CustomNonbondedForce(basic_electrostatics_expression +
                                                                                   'lambda_electrostatics=1.0;')

        # Add parameters and configure CustomNonbondedForces to match reference force
        for force in [na_electrostatics_custom_nonbonded_force, aa_electrostatics_custom_nonbonded_force]:
            force.addPerParticleParameter("charge")  # partial charge
            force.addPerParticleParameter("sigma")  # Lennard-Jones sigma
            force.setUseSwitchingFunction(False)  # no switch for electrostatics, since NonbondedForce doesn't use it
            force.setCutoffDistance(nonbonded_force.getCutoffDistance())
            force.setUseLongRangeCorrection(False)  # long-range dispersion correction is meaningless for electrostatics

            if is_method_periodic:
                force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            else:
                force.setNonbondedMethod(nonbonded_force.getNonbondedMethod())

        # Add CustomNonbondedForces modelling particle interactions to alchemical system.
        force_index = system.addForce(na_sterics_custom_nonbonded_force)
        force_labels['alchemically modified NonbondedForce for non-alchemical/alchemical sterics'] = force_index

        force_index = system.addForce(aa_sterics_custom_nonbonded_force)
        force_labels['alchemically modified NonbondedForce for alchemical/alchemical sterics'] = force_index

        force_index = system.addForce(na_electrostatics_custom_nonbonded_force)
        force_labels['alchemically modified NonbondedForce for non-alchemical/alchemical electrostatics'] = force_index

        force_index = system.addForce(aa_electrostatics_custom_nonbonded_force)
        force_labels['alchemically modified NonbondedForce for alchemical/alchemical electrostatics'] = force_index

        # Create CustomBondForces to handle sterics 1,4 exceptions interactions between
        # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
        # to 1.0 for decoupled interactions
        basic_sterics_expression = "U_sterics;" + sterics_energy_expression + alchemical_function_expression

        na_sterics_custom_bond_force = openmm.CustomBondForce(basic_sterics_expression)
        na_sterics_custom_bond_force.addGlobalParameter("lambda_sterics", 1.0)
        if self.annihilate_sterics:
            aa_sterics_custom_bond_force = openmm.CustomBondForce(basic_sterics_expression)
            aa_sterics_custom_bond_force.addGlobalParameter("lambda_sterics", 1.0)
        else:  # for decoupling fix lambda_sterics to 1.0
            aa_sterics_custom_bond_force = openmm.CustomBondForce(basic_sterics_expression + 'lambda_sterics=1.0;')

        for force in [na_sterics_custom_bond_force, aa_sterics_custom_bond_force]:
            force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma
            force.addPerBondParameter("epsilon")  # Lennard-Jones effective epsilon

        # Create CustomBondForces to handle electrostatics 1,4 exceptions interactions between
        # non-alchemical/alchemical atoms (na) and alchemical/alchemical atoms (aa). Fix lambda
        # to 1.0 for decoupled interactions
        if self.consistent_exceptions:
            basic_electrostatics_expression = "U_electrostatics;" + electrostatics_energy_expression +\
                                              alchemical_function_expression
        else:
            basic_electrostatics_expression = "U_electrostatics;" + exceptions_electrostatics_energy_expression +\
                                              alchemical_function_expression

        na_electrostatics_custom_bond_force = openmm.CustomBondForce(basic_electrostatics_expression)
        na_electrostatics_custom_bond_force.addGlobalParameter("lambda_electrostatics", 1.0)
        if self.annihilate_electrostatics:
            aa_electrostatics_custom_bond_force = openmm.CustomBondForce(basic_electrostatics_expression)
            aa_electrostatics_custom_bond_force.addGlobalParameter("lambda_electrostatics", 1.0)
        else:  # for decoupling fix lambda_electrostatics to 1.0
            aa_electrostatics_custom_bond_force = openmm.CustomBondForce(basic_electrostatics_expression +
                                                                         'lambda_electrostatics=1.0;')

        # Create CustomBondForce to handle exceptions for electrostatics
        for force in [na_electrostatics_custom_bond_force, aa_electrostatics_custom_bond_force]:
            force.addPerBondParameter("chargeprod")  # charge product
            force.addPerBondParameter("sigma")  # Lennard-Jones effective sigma

        # Add CustomNonbondedForces modelling particle interactions to alchemical system.
        force_index = system.addForce(na_sterics_custom_bond_force)
        force_labels['alchemically modified BondForce for non-alchemical/alchemical sterics exceptions'] = force_index

        force_index = system.addForce(aa_sterics_custom_bond_force)
        force_labels['alchemically modified BondForce for alchemical/alchemical sterics exceptions'] = force_index

        force_index = system.addForce(na_electrostatics_custom_bond_force)
        force_labels['alchemically modified BondForce for non-alchemical/alchemical electrostatics exceptions'] = force_index

        force_index = system.addForce(aa_electrostatics_custom_bond_force)
        force_labels['alchemically modified BondForce for alchemical/alchemical electrostatics exceptions'] = force_index

        # -------------------------------------------------------------------------------
        # Distribute particle interactions contributions in appropriate nonbonded forces
        # -------------------------------------------------------------------------------

        # Fix any NonbondedForce issues with Lennard-Jones sigma = 0 (epsilon = 0), which should have sigma > 0.
        for particle_index in range(nonbonded_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
            # Check particle sigma is not zero.
            if (sigma == 0.0 * unit.angstrom):
                logger.warning("particle %d has Lennard-Jones sigma = 0 (charge=%s, sigma=%s, epsilon=%s); setting sigma=1A" % (particle_index, str(charge), str(sigma), str(epsilon)))
                sigma = 1.0 * unit.angstrom
                # Fix it.
                nonbonded_force.setParticleParameters(particle_index, charge, sigma, epsilon)
        for exception_index in range(nonbonded_force.getNumExceptions()):
            # Retrieve parameters.
            [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
            # Check particle sigma is not zero.
            if (sigma == 0.0 * unit.angstrom):
                logger.warning("exception %d has Lennard-Jones sigma = 0 (iatom=%d, jatom=%d, chargeprod=%s, sigma=%s, epsilon=%s); setting sigma=1A" % (exception_index, iatom, jatom, str(chargeprod), str(sigma), str(epsilon)))
                sigma = 1.0 * unit.angstrom
                # Fix it.
                nonbonded_force.setExceptionParameters(exception_index, iatom, jatom, chargeprod, sigma, epsilon)

        # Copy NonbondedForce particle terms for alchemically-modified particles to CustomNonbondedForces.
        # On CUDA, for efficiency reasons, all nonbonded forces (custom and not) must have the same particles.
        for particle_index in range(nonbonded_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
            if (sigma / unit.angstroms) == 0.0:
                raise Exception('sigma is %s for particle %d; sigma must be positive' % (str(sigma), particle_index))
            na_sterics_custom_nonbonded_force.addParticle([sigma, epsilon])
            aa_sterics_custom_nonbonded_force.addParticle([sigma, epsilon])
            na_electrostatics_custom_nonbonded_force.addParticle([charge, sigma])
            aa_electrostatics_custom_nonbonded_force.addParticle([charge, sigma])

        # Create atom groups.
        alchemical_atomset = self.ligand_atomset.copy()
        all_atomset = set(range(system.getNumParticles()))  # all atoms, including alchemical region
        nonalchemical_atomset = all_atomset.difference(alchemical_atomset)

        # Turn off interactions contribution from alchemically-modified particles in unmodified
        # NonbondedForce that will be handled by all other forces
        for particle_index in range(nonbonded_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
            if particle_index in alchemical_atomset:
                nonbonded_force.setParticleParameters(particle_index, abs(0*charge), sigma, abs(0*epsilon))

        # Restrict interaction evaluation of CustomNonbondedForces to their respective atom groups.
        na_sterics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomset)
        na_electrostatics_custom_nonbonded_force.addInteractionGroup(nonalchemical_atomset, alchemical_atomset)
        aa_sterics_custom_nonbonded_force.addInteractionGroup(alchemical_atomset, alchemical_atomset)
        aa_electrostatics_custom_nonbonded_force.addInteractionGroup(alchemical_atomset, alchemical_atomset)

        # ---------------------------------------------------------------
        # Distribute exceptions contributions in appropriate bond forces
        # ---------------------------------------------------------------

        # Move all NonbondedForce exception terms for alchemically-modified particles to CustomBondForces.
        for exception_index in range(nonbonded_force.getNumExceptions()):
            # Retrieve parameters.
            [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)

            # Exclude this atom pair in CustomNonbondedForces. All nonbonded forces
            # must have the same number of exceptions/exclusions on CUDA platform.
            na_sterics_custom_nonbonded_force.addExclusion(iatom, jatom)
            aa_sterics_custom_nonbonded_force.addExclusion(iatom, jatom)
            na_electrostatics_custom_nonbonded_force.addExclusion(iatom, jatom)
            aa_electrostatics_custom_nonbonded_force.addExclusion(iatom, jatom)

            # Check how many alchemical atoms we have
            both_alchemical = iatom in alchemical_atomset and jatom in alchemical_atomset
            only_one_alchemical = (iatom in alchemical_atomset) != (jatom in alchemical_atomset)

            # Check if this is an exception or an exclusion
            is_exception_epsilon = abs(epsilon.value_in_unit_system(unit.md_unit_system)) > 0.0
            is_exception_chargeprod = abs(chargeprod.value_in_unit_system(unit.md_unit_system)) > 0.0

            # If exception (and not exclusion), add special CustomBondForce terms to
            # handle alchemically-modified Lennard-Jones and electrostatics exceptions
            if both_alchemical:
                if is_exception_epsilon:
                    aa_sterics_custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])
                if is_exception_chargeprod:
                    aa_electrostatics_custom_bond_force.addBond(iatom, jatom, [chargeprod, sigma])
            elif only_one_alchemical:
                if is_exception_epsilon:
                    na_sterics_custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])
                if is_exception_chargeprod:
                    na_electrostatics_custom_bond_force.addBond(iatom, jatom, [chargeprod, sigma])
            # else: both particles are non-alchemical, leave them in the unmodified NonbondedForce

        # Turn off all exception contributions from alchemical atoms in the NonbondedForce
        # modelling non-alchemical atoms only
        for exception_index in range(nonbonded_force.getNumExceptions()):
            [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
            if iatom in alchemical_atomset or jatom in alchemical_atomset:
                nonbonded_force.setExceptionParameters(exception_index, iatom, jatom,
                                                       abs(0.0*chargeprod), sigma, abs(0.0*epsilon))

        # Add global parameters to forces.
        def add_global_parameters(force):
            force.addGlobalParameter('softcore_alpha', self.softcore_alpha)
            force.addGlobalParameter('softcore_beta', self.softcore_beta)
            force.addGlobalParameter('softcore_a', self.softcore_a)
            force.addGlobalParameter('softcore_b', self.softcore_b)
            force.addGlobalParameter('softcore_c', self.softcore_c)
            force.addGlobalParameter('softcore_d', self.softcore_d)
            force.addGlobalParameter('softcore_e', self.softcore_e)
            force.addGlobalParameter('softcore_f', self.softcore_f)

            # Add control variables.
            control_variables = set(self.alchemical_functions.values())
            for variable in control_variables:
                force.addGlobalParameter(variable, 1.0)

        for force in [na_sterics_custom_nonbonded_force, na_electrostatics_custom_nonbonded_force,
                      aa_sterics_custom_nonbonded_force, aa_electrostatics_custom_nonbonded_force,
                      na_sterics_custom_bond_force, na_electrostatics_custom_bond_force,
                      aa_sterics_custom_bond_force, aa_electrostatics_custom_bond_force]:
            add_global_parameters(force)

        return

    def _alchemicallyModifyAmoebaMultipoleForce(self, system, reference_force):
        raise Exception("Not implemented; needs CustomMultipleForce")
        alchemical_atom_indices = self.ligand_atoms


    def _alchemicallyModifyAmoebaVdwForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of AmoebaVdwForce

        Parameters
        ----------
        system : simtk.openmm.System
            Alchemically-modified system being built.  This object will be modified.
        nonbonded_force : simtk.openmm.NonbondedForce
            The NonbondedForce used as a template.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        TODO
        ----
        * Supported periodic boundary conditions need to be handled correctly.
        * Exceptions/exclusions need to be dealt with.

        """
        # This feature is incompletely implemented, so raise an exception.
        raise Exception("Not implemented")

        # Softcore Halgren potential from Eq. 3 of
        # Shi, Y., Jiao, D., Schnieders, M.J., and Ren, P. (2009). Trypsin-ligand binding free energy calculation with AMOEBA. Conf Proc IEEE Eng Med Biol Soc 2009, 2328-2331.
        energy_expression = 'lambda^5 * epsilon * (1.07^7 / (0.7*(1-lambda)^2+(rho+0.07)^7)) * (1.12 / (0.7*(1-lambda)^2 + rho^7 + 0.12) - 2);'
        energy_expression += 'epsilon = 4*epsilon1*epsilon2 / (sqrt(epsilon1) + sqrt(epsilon2))^2;'
        energy_expression += 'rho = r / R0;'
        energy_expression += 'R0 = (R01^3 + R02^3) / (R01^2 + R02^2);'
        energy_expression += 'lambda = vdw_lambda * (ligand1*(1-ligand2) + ligand2*(1-ligand1)) + ligand1*ligand2;'
        energy_expression += 'vdw_lambda = %f;' % vdw_lambda

        softcore_force = openmm.CustomNonbondedForce(energy_expression)
        softcore_force.addPerParticleParameter('epsilon')
        softcore_force.addPerParticleParameter('R0')
        softcore_force.addPerParticleParameter('ligand')

        for particle_index in range(system.getNumParticles()):
            # Retrieve parameters from vdW force.
            [parentIndex, sigma, epsilon, reductionFactor] = force.getParticleParameters(particle_index)
            # Add parameters to CustomNonbondedForce.
            if particle_index in ligand_atoms:
                softcore_force.addParticle([epsilon, sigma, 1])
            else:
                softcore_force.addParticle([epsilon, sigma, 0])

            # Deal with exclusions.
            excluded_atoms = force.getParticleExclusions(particle_index)
            for jatom in excluded_atoms:
                if (particle_index < jatom):
                    softcore_force.addExclusion(particle_index, jatom)

        # Make sure periodic boundary conditions are treated the same way.
        # TODO: Handle PBC correctly.
        softcore_force.setNonbondedMethod( openmm.CustomNonbondedForce.CutoffPeriodic )
        softcore_force.setCutoffDistance( force.getCutoff() )

        # Add the softcore force.
        force_index = system.addForce( softcore_force )
        force_labels['alchemically modified AmoebaVdwForce'] = force_index

        # Turn off vdW interactions for alchemically-modified atoms.
        for particle_index in ligand_atoms:
            # Retrieve parameters.
            [parentIndex, sigma, epsilon, reductionFactor] = force.getParticleParameters(particle_index)
            epsilon = 1.0e-6 * epsilon # TODO: For some reason, we cannot set epsilon to 0.
            force.setParticleParameters(particle_index, parentIndex, sigma, epsilon, reductionFactor)

        # Deal with exceptions here.
        # TODO

        return system

    def _alchemicallyModifyGBSAOBCForce(self, system, reference_force, force_labels, sasa_model='ACE'):
        """
        Create alchemically-modified version of GBSAOBCForce.

        Parameters
        ----------
        system : simtk.openmm.System
            Alchemically-modified System object being built.  This object will be modified.
        reference_force : simtk.openmm.GBSAOBCForce
            Reference force to use for template.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`
        sasa_model : str, optional, default='ACE'
            Solvent accessible surface area model.

        """

        custom_force = openmm.CustomGBForce()

        # Add per-particle parameters.
        custom_force.addGlobalParameter("lambda_electrostatics", 1.0);
        custom_force.addPerParticleParameter("charge");
        custom_force.addPerParticleParameter("radius");
        custom_force.addPerParticleParameter("scale");
        custom_force.addPerParticleParameter("alchemical");

        # Set nonbonded method.
        custom_force.setNonbondedMethod(reference_force.getNonbondedMethod())
        custom_force.setCutoffDistance(reference_force.getCutoffDistance())

        # Add global parameters.
        custom_force.addGlobalParameter("solventDielectric", reference_force.getSolventDielectric())
        custom_force.addGlobalParameter("soluteDielectric", reference_force.getSoluteDielectric())
        custom_force.addGlobalParameter("offset", 0.009)

        custom_force.addComputedValue("I",  "(lambda_electrostatics*alchemical2 + (1-alchemical2))*step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
                                "U=r+sr2;"
                                "L=max(or1, D);"
                                "D=abs(r-sr2);"
                                "sr2 = scale2*or2;"
                                "or1 = radius1-offset; or2 = radius2-offset", openmm.CustomGBForce.ParticlePairNoExclusions)

        custom_force.addComputedValue("B", "1/(1/or-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);"
                                  "psi=I*or; or=radius-offset", openmm.CustomGBForce.SingleParticle)

        custom_force.addEnergyTerm("-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*(lambda_electrostatics*alchemical+(1-alchemical))*charge^2/B", openmm.CustomGBForce.SingleParticle)
        if sasa_model == 'ACE':
            custom_force.addEnergyTerm("(lambda_electrostatics*alchemical+(1-alchemical))*28.3919551*(radius+0.14)^2*(radius/B)^6", openmm.CustomGBForce.SingleParticle)

        custom_force.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*(lambda_electrostatics*alchemical1+(1-alchemical1))*charge1*(lambda_electrostatics*alchemical2+(1-alchemical2))*charge2/f;"
                             "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", openmm.CustomGBForce.ParticlePairNoExclusions);

        # Add particle parameters.
        for particle_index in range(reference_force.getNumParticles()):
            # Retrieve parameters.
            [charge, radius, scaling_factor] = reference_force.getParticleParameters(particle_index)
            # Set particle parameters.
            if particle_index in self.ligand_atoms:
                parameters = [charge, radius, scaling_factor, 1.0]
            else:
                parameters = [charge, radius, scaling_factor, 0.0]
            custom_force.addParticle(parameters)

        # Add alchemically-modified GBSAOBCForce to system.
        force_index = system.addForce(custom_force)
        force_labels['alchemically modified GBSAOBCForce'] = force_index

    def _alchemicallyModifyCustomGBForce(self, system, reference_force, force_labels):
        """
        Create alchemically-modified version of CustomGBForce by metaprogramming GB functions.

        The following rules are applied:
        - 'lambda_electrostatics' is added as a global parameter.
        - 'alchemical' is added as a per-particle parameter.
           All atoms in the alchemical group have this parameter set to 1; otherwise 0.
        - Any single-particle energy term (`CustomGBForce.SingleParticle`) is scaled by `(lambda_electrostatics*alchemical+(1-alchemical))`
        - Any two-particle energy term (`CustomGBForce.ParticlePairNoExclusions`) has charge 1 (`charge1`) replaced by `(lambda_electrostatics*alchemical1+(1-alchemical1))*charge1` and charge 2 (`charge2`) replaced by `(lambda_electrostatics*alchemical2+(1-alchemical2))*charge2`.
        - Any single-particle computed value (`CustomGBForce.SingleParticle`) remains unmodified
        - Any two-particle computed value (`CustomGBForce.ParticlePairNoExclusions`) is scaled by `(lambda_electrostatics*alchemical2 + (1-alchemical2))`

        Scaling of a term should always prepend and capture the value with an intermediate variable.
        For example, prepending `scaling * unscaled; unscaled =` will capture the value of the expression as `unscaled` and multiple by `scaled`.
        This avoids the need to identify the head expression and add parentheses.

        WARNING: This may not work correctly for all GB models.

        Parameters
        ----------
        system : simtk.openmm.System
            Alchemically-modified System object being built.  This object will be modified.
        reference_force : simtk.openmm.GBSAOBCForce
            Reference force to use for template.
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        """

        custom_force = openmm.CustomGBForce()

        # Add global parameters
        for index in range(reference_force.getNumGlobalParameters()):
            name = reference_force.getGlobalParameterName(index)
            default_value = reference_force.getGlobalParameterDefaultValue(index)
            custom_force.addGlobalParameter(name, default_value)
        custom_force.addGlobalParameter("lambda_electrostatics", 1.0);

        # Add per-particle parameters.
        for index in range(reference_force.getNumPerParticleParameters()):
            name = reference_force.getPerParticleParameterName(index)
            custom_force.addPerParticleParameter(name)
        custom_force.addPerParticleParameter("alchemical");

        # Set nonbonded methods.
        custom_force.setNonbondedMethod(reference_force.getNonbondedMethod())
        custom_force.setCutoffDistance(reference_force.getCutoffDistance())

        # Add computations.
        for index in range(reference_force.getNumComputedValues()):
            [name, expression, computation_type] = reference_force.getComputedValueParameters(index)

            # Alter expression for particle pair terms only.
            if not (computation_type == openmm.CustomGBForce.SingleParticle):
                prepend = 'alchemical_scaling*unscaled; alchemical_scaling = (lambda_electrostatics*alchemical2 + (1-alchemical2)); unscaled = '
                expression = prepend + expression

            custom_force.addComputedValue(name, expression, computation_type)

        # Add energy terms.
        for index in range(reference_force.getNumEnergyTerms()):
            [expression, computation_type] = reference_force.getEnergyTermParameters(index)

            # Alter expressions
            if (computation_type == openmm.CustomGBForce.SingleParticle):
                prepend = 'alchemical_scaling*unscaled; alchemical_scaling = (lambda_electrostatics*alchemical + (1-alchemical)); unscaled = '
                expression = prepend + expression
            else:
                expression.replace('charge1', 'alchemically_scaled_charge1')
                expression.replace('charge2', 'alchemically_scaled_charge2')
                expression += ' ; alchemically_scaled_charge1 = (lambda_electrostatics*alchemical1+(1-alchemical1)) * charge1;'
                expression += ' ; alchemically_scaled_charge2 = (lambda_electrostatics*alchemical2+(1-alchemical2)) * charge2;'

            custom_force.addEnergyTerm(expression, computation_type)

        # Add particle parameters
        for particle_index in range(reference_force.getNumParticles()):
            parameters = reference_force.getParticleParameters(particle_index)
            # Append alchemical parameter
            parameters = list(parameters)
            if particle_index in self.ligand_atoms:
                parameters.append(1.0)
            else:
                parameters.append(0.0)
            custom_force.addParticle(parameters)

        # Add tabulated functions
        for function_index in range(reference_force.getNumTabulatedFunctions()):
            name = reference_force.getTabulatedFunctionName(function_index)
            function = reference_force.getTabulatedFunction(function_index)
            function_copy = copy.deepcopy(function)
            custom_force.addTabulatedFunction(name, function_copy)

        # Add alchemically-modified CustomGBForce to system
        force_index = system.addForce(custom_force)
        force_labels['alchemically modified CustomGBForce'] = force_index

    def _createAlchemicallyModifiedSystem(self, mm=None):
        """
        Create an alchemically modified version of the reference system with global parameters encoding alchemical parameters.

        Returns
        -------
        system : simtk.openmm.System
            Alchemically-modified version of self.reference_system
        force_labels : dict of int : str
            force_labels[name] is the force index in the alchemically modified system of the modified force `name`

        TODO
        ----
        * This could be streamlined if it was possible to modify System or Force objects.
        * isinstance(mm.NonbondedForce) and related expressions won't work if reference system was created with a different OpenMM implementation.
          Use class names instead.

        """
        # TODO: Check that the provided system hasn't already been alchemically modified.
        # We don't want to allow a System to be modified twice!

        # Record timing statistics.
        initial_time = time.time()
        logger.debug("Creating alchemically modified system...")

        reference_system = self.reference_system

        # Create new deep copy reference system to modify.
        system = openmm.System()

        # Set periodic box vectors.
        [a,b,c] = reference_system.getDefaultPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(a,b,c)

        # Add atoms.
        for atom_index in range(reference_system.getNumParticles()):
            mass = reference_system.getParticleMass(atom_index)
            system.addParticle(mass)

        # Add constraints
        for constraint_index in range(reference_system.getNumConstraints()):
            [iatom, jatom, r0] = reference_system.getConstraintParameters(constraint_index)
            system.addConstraint(iatom, jatom, r0)

        # Keep track of which force indices of alchemically-modifed system arise from which contributions.
        # e.g. force_labels['alchemically modified HarmonicBondForce'] is force index of alchemically-modified HarmonicBondForce
        from collections import OrderedDict
        force_labels = OrderedDict()

        # Modify forces as appropriate, copying other forces without modification.
        # TODO: Use introspection to automatically dispatch registered modifiers?
        nforces = reference_system.getNumForces()
        for force_index in range(nforces):
            reference_force = reference_system.getForce(force_index)
            if isinstance(reference_force, openmm.PeriodicTorsionForce) and (self.alchemical_torsions is not None):
                self._alchemicallyModifyPeriodicTorsionForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.HarmonicAngleForce) and (self.alchemical_angles is not None):
                self._alchemicallyModifyHarmonicAngleForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.HarmonicBondForce) and (self.alchemical_bonds is not None):
                self._alchemicallyModifyHarmonicBondForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.NonbondedForce):
                self._alchemicallyModifyNonbondedForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.CustomGBForce):
                self._alchemicallyModifyCustomGBForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.GBSAOBCForce):
                self._alchemicallyModifyGBSAOBCForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.AmoebaMultipoleForce):
                self._alchemicallyModifyAmoebaMultipoleForce(system, reference_force, force_labels)
            elif isinstance(reference_force, openmm.AmoebaVdwForce):
                self._alchemicallyModifyAmoebaVdwForce(system, reference_force, force_labels)
            else:
                # Copy force without modification.
                force = copy.deepcopy(reference_force)
                force_index = system.addForce(force)
                force_name = force.__class__.__name__
                force_labels['unmodified %s' % force_name] = force_index

        # Record timing statistics.
        final_time = time.time()
        elapsed_time = final_time - initial_time
        logger.debug("Elapsed time %.3f s." % (elapsed_time))

        return [system, force_labels]

    @classmethod
    def perturbSystem(cls, system, alchemical_state):
        """
        Perturb the specified system (previously generated by this alchemical factory) by setting default global parameters.

        Parameters
        ----------
        system : simtk.openmm.System created by AlchemicalFactory
            The alchemically-modified system whose default global parameters are to be perturbed.
        alchemical_state : AlchemicalState
            The alchemical state to create from the reference system.

        Examples
        --------

        Create alchemical intermediates for 'denihilating' one water in a water box.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> waterbox = testsystems.WaterBox()
        >>> [reference_system, positions] = [waterbox.system, waterbox.positions]
        >>> # Create a factory to produce alchemical intermediates.
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])
        >>> # Create an alchemically-perturbed state corresponding to fully-interacting.
        >>> alchemical_state = AlchemicalState(lambda_sterics=1.00, lambda_electrostatics=1.00)
        >>> # Create the perturbed system.
        >>> alchemical_system = factory.createPerturbedSystem(alchemical_state)
        >>> # Perturb this system.
        >>> alchemical_state = AlchemicalState(lambda_sterics=0.90, lambda_electrostatics=0.90)
        >>> factory.perturbSystem(alchemical_system, alchemical_state)

        """

        for force_index in range(system.getNumForces()):
            force = system.getForce(force_index)
            if hasattr(force, 'getNumGlobalParameters'):
                for parameter_index in range(force.getNumGlobalParameters()):
                    parameter_name = force.getGlobalParameterName(parameter_index)
                    if parameter_name in alchemical_state:
                        force.setGlobalParameterDefaultValue(parameter_index, alchemical_state[parameter_name])

        return

    def getEnergyComponents(self, alchemical_state, positions, box_vectors=None, use_all_parameters=True):
        """
        Compute potential energy by Force component for the corresponding alchemically-modified system.

        Parameters
        ----------
        alchemical_state : AlchemicalState
            The alchemical state to est the Context to.
        use_all_parameters : bool, optional, default=False
            If True, will ensure that all parameters are used or raise an Exception if not.

        """
        if not hasattr(self, 'alchemically_modified_system_with_force_groups'):
            # Create deep copy of alchemical system.
            system = copy.deepcopy(self.alchemically_modified_system)
            # Separate all forces into separate force groups.
            assert system.getNumForces() <= 32, "self.alchemically_modified_system has more than 32 force groups; " \
                                                "can't compute individual force component energies."
            for (force_index, force) in enumerate(system.getForces()):
                force.setForceGroup(force_index)
            # Store system
            self.alchemically_modified_system_with_force_groups = system

        # Retrieve system with distinct force groups
        system = self.alchemically_modified_system_with_force_groups

        # Create a Context
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        context = openmm.Context(system, integrator)
        # Set alchemical state
        AbsoluteAlchemicalFactory.perturbContext(context, alchemical_state, use_all_parameters=use_all_parameters)
        # Set positions and box vectors
        if box_vectors is not None:
            context.setDefaultPeriodicBoxVectors(*box_vectors)
        context.setPositions(positions)
        # Get energy component  s
        from collections import OrderedDict
        energy_components = OrderedDict()
        for (force_label, force_index) in self.force_labels.items():
            energy_components[force_label] = context.getState(getEnergy=True,groups=2**force_index).getPotentialEnergy()
        # Clean up
        del context, integrator
        # Return energy components
        return energy_components

    @classmethod
    def perturbContext(cls, context, alchemical_state, use_all_parameters=False):
        """
        Perturb the specified context (using a system previously generated by this alchemical factory) by setting context parameters.

        Parameters
        ----------
        context : simtk.openmm.Context
            The Context object associated with a System that has been alchemically modified.
        alchemical_state : AlchemicalState
            The alchemical state to est the Context to.
        use_all_parameters : bool, optional, default=False
            If True, will ensure that all parameters are used or raise an Exception if not.

        Examples
        --------

        Create an alchemically-modified water box and set alchemical parameters appropriately.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> waterbox = testsystems.WaterBox()
        >>> [reference_system, positions] = [waterbox.system, waterbox.positions]
        >>> # Create a factory to produce alchemical intermediates.
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])
        >>> # Create an alchemically-perturbed state corresponding to fully-interacting.
        >>> alchemical_state = AlchemicalState()
        >>> # Create the perturbed system.
        >>> alchemical_system = factory.createPerturbedSystem(alchemical_state)
        >>> # Create a Context.
        >>> integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
        >>> context = openmm.Context(alchemical_system, integrator)
        >>> # Set alchemical state.
        >>> alchemical_state = AlchemicalState(lambda_electrostatics=0.5, lambda_sterics=0.5)
        >>> factory.perturbContext(context, alchemical_state)

        """

        # Set parameters in Context.
        for parameter in alchemical_state.keys():
            try:
                context.setParameter(parameter, alchemical_state[parameter])
            except:
                if use_all_parameters:
                    # Print some useful output about what parameters were available.
                    available_parameters = context.getState(getParameters=True).getParameters().keys()
                    raise Exception("Parameter '%s' not available in Context; available parameters are: %s" % (parameter, str(available_parameters)))
                else:
                    # It's OK if some parameters are not found.
                    pass
        return

    def createPerturbedSystem(self, alchemical_state=None, mm=None):
        """
        Create a perturbed copy of the system given the specified alchemical state.

        Parameters
        ----------
        alchemical_state : AlchemicalState, optional, default=None
            The alchemical state to create from the reference system; if None, will create a fully-interacting alchemically-modified version.

        TODO
        ----
        * This could be streamlined if it was possible to modify System or Force objects.
        * isinstance(mm.NonbondedForce) and related expressions won't work if reference system was created with a different OpenMM implementation.
          Use class names instead.

        Examples
        --------

        Create alchemical intermediates for 'denihilating' one water in a water box.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> waterbox = testsystems.WaterBox()
        >>> [reference_system, positions] = [waterbox.system, waterbox.positions]
        >>> # Create a factory to produce alchemical intermediates.
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])
        >>> # Create an alchemically-perturbed state corresponding to fully-interacting.
        >>> alchemical_state = AlchemicalState()
        >>> # Create the perturbed system.
        >>> alchemical_system = factory.createPerturbedSystem(alchemical_state)

        Create alchemical intermediates for 'denihilating' p-xylene in T4 lysozyme L99A in GBSA.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> complex = testsystems.LysozymeImplicit()
        >>> [reference_system, positions] = [complex.system, complex.positions]
        >>> # Create a factory to produce alchemical intermediates.
        >>> receptor_atoms = range(0,2603) # T4 lysozyme L99A
        >>> ligand_atoms = range(2603,2621) # p-xylene
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)
        >>> # Create an alchemically-perturbed state corresponding to fully-interacting.
        >>> alchemical_state = AlchemicalState()
        >>> # Create the perturbed systems using this protocol.
        >>> alchemical_system = factory.createPerturbedSystem(alchemical_state)

        """

        if alchemical_state == None:
            # TODO: Also set any other alchemical parameters defined for this system to be fully interacting.
            alchemical_state = AlchemicalState()

        # Record timing statistics.
        initial_time = time.time()
        logger.debug("Creating alchemically modified intermediate...")

        # Return an alchemically modified copy.
        system = copy.deepcopy(self.alchemically_modified_system)

        # Perturb the default global parameters for this system according to the alchemical parameters.
        self.perturbSystem(system, alchemical_state)

        # Test the system energy if requested.
        if self.test_positions is not None:
            self._checkEnergyIsFinite(system, self.test_positions, self.platform)

        # Record timing statistics.
        final_time = time.time()
        elapsed_time = final_time - initial_time
        logger.debug("Elapsed time %.3f s." % (elapsed_time))

        return system

    def createPerturbedSystems(self, alchemical_states):
        """
        Create a list of perturbed copies of the system given a specified set of alchemical states.

        Parameters
        ----------
        states : list of AlchemicalState
            List of alchemical states to generate.

        Returns
        -------
        systems : list of simtk.openmm.System
            List of alchemically-modified System objects.  The cached reference system will be unmodified.

        Examples
        --------

        Create alchemical intermediates for 'denihilating' p-xylene in T4 lysozyme L99A in GBSA.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> complex = testsystems.LysozymeImplicit()
        >>> [reference_system, positions] = [complex.system, complex.positions]
        >>> # Create a factory to produce alchemical intermediates.
        >>> receptor_atoms = range(0,2603) # T4 lysozyme L99A
        >>> ligand_atoms = range(2603,2621) # p-xylene
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=ligand_atoms)
        >>> # Get the default protocol for 'denihilating' in complex in explicit solvent.
        >>> protocol = factory.defaultComplexProtocolImplicit()
        >>> # Create the perturbed systems using this protocol.
        >>> systems = factory.createPerturbedSystems(protocol)

        """

        initial_time = time.time()

        systems = list()
        for (state_index, alchemical_state) in enumerate(alchemical_states):
            logger.debug("Creating alchemical system %d / %d..." % (state_index, len(alchemical_states)))
            system = self.createPerturbedSystem(alchemical_state)
            systems.append(system)

        # Report timing.
        final_time = time.time()
        elapsed_time = final_time - initial_time
        logger.debug("createPerturbedSystems: Elapsed time %.3f s." % elapsed_time)

        return systems

    def _is_restraint(self, valence_atoms):
        """
        Determine whether specified valence term connects the ligand with its environment.

        Parameters
        ----------
        valence_atoms : list of int
            Atom indices involved in valence term (bond, angle or torsion).

        Returns
        -------
        is_restraint : bool
            True if the set of atoms includes at least one ligand atom and at least one non-ligand atom; False otherwise

        Examples
        --------

        Various tests for a simple system.

        >>> # Create a reference system.
        >>> from openmmtools import testsystems
        >>> alanine_dipeptide = testsystems.AlanineDipeptideImplicit()
        >>> [reference_system, positions] = [alanine_dipeptide.system, alanine_dipeptide.positions]
        >>> # Create a factory.
        >>> factory = AbsoluteAlchemicalFactory(reference_system, ligand_atoms=[0, 1, 2])
        >>> factory._is_restraint([0,1,2])
        False
        >>> factory._is_restraint([1,2,3])
        True
        >>> factory._is_restraint([3,4])
        False
        >>> factory._is_restraint([2,3,4,5])
        True

        """

        valence_atomset = set(valence_atoms)
        intersection = set.intersection(valence_atomset, self.ligand_atomset)
        if (len(intersection) >= 1) and (len(intersection) < len(valence_atomset)):
            return True

        return False
