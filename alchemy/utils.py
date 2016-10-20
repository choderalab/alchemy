#!/usr/bin/python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Utility functions

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import copy
import time
import itertools

from simtk import openmm, unit
import openmmtools
import argparse

import logging
logger = logging.getLogger(__name__)

#=============================================================================================
# PARAMETERS
#=============================================================================================

#=============================================================================================
# MODULE UTILITIES
#=============================================================================================

def benchmark(reference_system, positions, platform_name=None, nsteps=500, timestep=1.0*unit.femtoseconds, factory_args=None):
    """
    Benchmark performance of alchemically modified system relative to original system.

    Parameters
    ----------
    reference_system : simtk.openmm.System
       The reference System object to compare with
    positions : simtk.unit.Quantity with units compatible with nanometers
       The positions to assess energetics for.
    platform_name : str, optional, default=None
       The name of the platform to use for benchmarking.
    nsteps : int, optional, default=500
       Number of molecular dynamics steps to use for benchmarking.
    timestep : simtk.unit.Quantity with units compatible with femtoseconds, optional, default=1*femtoseconds
       Timestep to use for benchmarking.
    factory_args : dict(), optional, default=None
       Arguments passed to AbsoluteAlchemicalFactory.

    """
    from alchemy import AbsoluteAlchemicalFactory, AlchemicalState

    # Create a factory to produce alchemical intermediates.
    print("Creating alchemical factory...")
    initial_time = time.time()
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    print("AbsoluteAlchemicalFactory initialization took %.3f s" % elapsed_time)

    # Create an alchemically-perturbed state corresponding to nearly fully-interacting.
    # NOTE: We use a lambda slightly smaller than 1.0 because the AlchemicalFactory does not use Custom*Force softcore versions if lambda = 1.0 identically.
    lambda_value = 1.0 - 1.0e-6
    alchemical_state = AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value)

    platform = None
    if platform_name:
        platform = openmm.Platform.getPlatformByName(platform_name)

    # Create the perturbed system.
    print("Creating alchemically-modified state...")
    initial_time = time.time()
    alchemical_system = factory.createPerturbedSystem(alchemical_state)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    # Compare energies.
    print("Computing reference energies...")
    reference_integrator = openmm.VerletIntegrator(timestep)
    if platform:
        reference_context = openmm.Context(reference_system, reference_integrator, platform)
    else:
        reference_context = openmm.Context(reference_system, reference_integrator)
    reference_context.setPositions(positions)
    reference_state = reference_context.getState(getEnergy=True)
    reference_potential = reference_state.getPotentialEnergy()
    print("Computing alchemical energies...")
    alchemical_integrator = openmm.VerletIntegrator(timestep)
    if platform:
        alchemical_context = openmm.Context(alchemical_system, alchemical_integrator, platform)
    else:
        alchemical_context = openmm.Context(alchemical_system, alchemical_integrator)
    alchemical_context.setPositions(positions)
    alchemical_state = alchemical_context.getState(getEnergy=True)
    alchemical_potential = alchemical_state.getPotentialEnergy()
    delta = alchemical_potential - reference_potential

    # Make sure all kernels are compiled.
    reference_integrator.step(2)
    alchemical_integrator.step(2)

    # Time simulations.
    print("Simulating reference system...")
    initial_time = time.time()
    reference_integrator.step(nsteps)
    reference_state = reference_context.getState(getEnergy=True)
    reference_potential = reference_state.getPotentialEnergy()
    final_time = time.time()
    reference_time = final_time - initial_time
    print("Simulating alchemical system...")
    initial_time = time.time()
    alchemical_integrator.step(nsteps)
    alchemical_state = alchemical_context.getState(getEnergy=True)
    alchemical_potential = alchemical_state.getPotentialEnergy()
    final_time = time.time()
    alchemical_time = final_time - initial_time

    seconds_per_day = (1.*unit.day)/(1.*unit.seconds)

    print("TIMINGS")
    print("reference system       : %12.3f s for %8d steps (%12.3f ms/step; %12.3f ns/day)" % (reference_time, nsteps, reference_time/nsteps*1000, nsteps*timestep*(seconds_per_day/reference_time)/unit.nanoseconds))
    print("alchemical system      : %12.3f s for %8d steps (%12.3f ms/step; %12.3f ns/day)" % (alchemical_time, nsteps, alchemical_time/nsteps*1000, nsteps*timestep*(seconds_per_day/alchemical_time)/unit.nanoseconds))
    print("alchemical simulation is %12.3f x slower than unperturbed system" % (alchemical_time / reference_time))

    return delta

def run_benchmark():
    """
    Benchmark fully interacting system versus alchemically-modified system.
    """
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark alchemically modified system against unmodified system.')
    parser.add_argument('--platform', dest='platform_name', action='store', default=None, help='platform name to benchmark (default: None)')
    options = parser.parse_args()

    from sams.tests import testsystems
    for testsystem_name in ['AblImatinibExplicitAlchemical']:
        cls = getattr(testsystems, testsystem_name)
        testsystem = cls()
        factory_args = { 'ligand_atoms' : testsystem.alchemical_atoms, 'receptor_atoms' : range(0,4266) }
        benchmark(testsystem.system, testsystem.positions, platform_name=options.platform_name, nsteps=5000, timestep=2.0*unit.femtoseconds, factory_args=factory_args)
