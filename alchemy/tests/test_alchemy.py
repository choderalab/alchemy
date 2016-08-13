#!/usr/bin/python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for alchemical factory in `alchemy.py`.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os, os.path
import numpy as np
import copy
import time
from functools import partial

from simtk import unit, openmm
from simtk.openmm import app

from nose.plugins.attrib import attr
import pymbar

import logging
logger = logging.getLogger(__name__)

from openmmtools import testsystems

from alchemy import AlchemicalState, AbsoluteAlchemicalFactory

from nose.plugins.skip import Skip, SkipTest

#=============================================================================================
# CONSTANTS
#=============================================================================================

kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA # Boltzmann constant
temperature = 300.0 * unit.kelvin # reference temperature
#MAX_DELTA = 0.01 * kB * temperature # maximum allowable deviation
MAX_DELTA = 1.0 * kB * temperature # maximum allowable deviation

#=============================================================================================
# SUBROUTINES FOR TESTING
#=============================================================================================

def config_root_logger(verbose, log_file_path=None, mpicomm=None):
    """Setup the the root logger's configuration.
     The log messages are printed in the terminal and saved in the file specified
     by log_file_path (if not None) and printed. Note that logging use sys.stdout
     to print logging.INFO messages, and stderr for the others. The root logger's
     configuration is inherited by the loggers created by logging.getLogger(name).
     Different formats are used to display messages on the terminal and on the log
     file. For example, in the log file every entry has a timestamp which does not
     appear in the terminal. Moreover, the log file always shows the module that
     generate the message, while in the terminal this happens only for messages
     of level WARNING and higher.
    Parameters
    ----------
    verbose : bool
        Control the verbosity of the messages printed in the terminal. The logger
        displays messages of level logging.INFO and higher when verbose=False.
        Otherwise those of level logging.DEBUG and higher are printed.
    log_file_path : str, optional, default = None
        If not None, this is the path where all the logger's messages of level
        logging.DEBUG or higher are saved.
    mpicomm : mpi4py.MPI.COMM communicator, optional, default=None
        If specified, this communicator will be used to determine node rank.
    """

    class TerminalFormatter(logging.Formatter):
        """
        Simplified format for INFO and DEBUG level log messages.
        This allows to keep the logging.info() and debug() format separated from
        the other levels where more information may be needed. For example, for
        warning and error messages it is convenient to know also the module that
        generates them.
        """

        # This is the cleanest way I found to make the code compatible with both
        # Python 2 and Python 3
        simple_fmt = logging.Formatter('%(asctime)-15s: %(message)s')
        default_fmt = logging.Formatter('%(asctime)-15s: %(levelname)s - %(name)s - %(message)s')

        def format(self, record):
            if record.levelno <= logging.INFO:
                return self.simple_fmt.format(record)
            else:
                return self.default_fmt.format(record)

    # Check if root logger is already configured
    n_handlers = len(logging.root.handlers)
    if n_handlers > 0:
        root_logger = logging.root
        for i in xrange(n_handlers):
            root_logger.removeHandler(root_logger.handlers[0])

    # If this is a worker node, don't save any log file
    if mpicomm:
        rank = mpicomm.rank
    else:
        rank = 0

    if rank != 0:
        log_file_path = None

    # Add handler for stdout and stderr messages
    terminal_handler = logging.StreamHandler()
    terminal_handler.setFormatter(TerminalFormatter())
    if rank != 0:
        terminal_handler.setLevel(logging.WARNING)
    elif verbose:
        terminal_handler.setLevel(logging.DEBUG)
    else:
        terminal_handler.setLevel(logging.INFO)
    logging.root.addHandler(terminal_handler)

    # Add file handler to root logger
    if log_file_path is not None:
        file_format = '%(asctime)s - %(levelname)s - %(name)s - %(message)s'
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(file_format))
        logging.root.addHandler(file_handler)

    # Do not handle logging.DEBUG at all if unnecessary
    if log_file_path is not None:
        logging.root.setLevel(logging.DEBUG)
    else:
        logging.root.setLevel(terminal_handler.level)

def dump_xml(system=None, integrator=None, state=None):
    """
    Dump system, integrator, and state to XML for debugging.
    """
    from simtk.openmm import XmlSerializer
    def write_file(filename, contents):
        outfile = open(filename, 'w')
        outfile.write(contents)
        outfile.close()
    if system: write_file('system.xml', XmlSerializer.serialize(system))
    if integrator: write_file('integrator.xml', XmlSerializer.serialize(integrator))
    if state: write_file('state.xml', XmlSerializer.serialize(state))
    return

def compute_energy(system, positions, platform=None, precision=None):
    timestep = 1.0 * unit.femtoseconds
    integrator = openmm.VerletIntegrator(timestep)
    if platform:
        context = openmm.Context(system, integrator, platform)
    else:
        context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    potential = state.getPotentialEnergy()
    del context, integrator, state
    return potential

def check_waterbox(platform=None, precision=None, nonbondedMethod=openmm.NonbondedForce.CutoffPeriodic):
    """Compare annihilated states in vacuum and a large box.
    """
    platform_name = platform.getName()
    from openmmtools import testsystems
    testsystem = testsystems.WaterBox()
    system = testsystem.system
    positions = testsystem.positions

    # Use reaction field
    for force in system.getForces():
        if force.__class__.__name__ == 'NonbondedForce':
            force.setNonbondedMethod(nonbondedMethod)

    factory_args = {'ligand_atoms' : [], 'receptor_atoms' : [],
        'annihilate_sterics' : False, 'annihilate_electrostatics' : True }

    # Create alchemically-modified system
    factory = AbsoluteAlchemicalFactory(system, **factory_args)
    alchemical_system = factory.createPerturbedSystem()

    # Compare energies
    system_energy = compute_energy(system, positions, platform=platform, precision=precision)
    alchemical_1_energy = compute_energy(alchemical_system, positions, platform=platform, precision=precision)

    # Set lambda = 0
    lambda_value = 0.0
    alchemical_state = AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value)
    AbsoluteAlchemicalFactory.perturbSystem(alchemical_system, alchemical_state)
    alchemical_0_energy = compute_energy(alchemical_system, positions, platform=platform, precision=precision)

    # Check deviation.
    logger.info("========")
    logger.info("Platform %s" % platform_name)
    logger.info("Alchemically-modified WaterBox with no alchemical atoms")
    logger.info('real system : %8.3f kcal/mol' % (system_energy / unit.kilocalories_per_mole))
    logger.info('lambda = 1  : %8.3f kcal/mol' % (alchemical_1_energy / unit.kilocalories_per_mole))
    logger.info('lambda = 0  : %8.3f kcal/mol' % (alchemical_0_energy / unit.kilocalories_per_mole))
    delta = alchemical_1_energy - alchemical_0_energy
    logger.info("ERROR       : %8.3f kcal/mol" % (delta / unit.kilocalories_per_mole))
    if (abs(delta) > MAX_DELTA):
        raise Exception("Maximum allowable deviation on platform %s exceeded (was %.8f kcal/mol; allowed %.8f kcal/mol); test failed." % (platform_name, delta / unit.kilocalories_per_mole, MAX_DELTA / unit.kilocalories_per_mole))

def test_waterbox():
    for platform_index in range(openmm.Platform.getNumPlatforms()):
        for nonbondedMethod in [openmm.NonbondedForce.PME, openmm.NonbondedForce.CutoffPeriodic]:
            platform = openmm.Platform.getPlatform(platform_index)
            f = partial(check_waterbox, platform=platform, nonbondedMethod=nonbondedMethod)
            yield f

def compare_platforms(system, positions, factory_args=dict()):
    # Create annihilated version of vacuum system.
    factory = AbsoluteAlchemicalFactory(system, **factory_args)
    alchemical_system = factory.createPerturbedSystem()

    def set_lambda(alchemical_system, lambda_value):
        alchemical_state = AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value)
        AbsoluteAlchemicalFactory.perturbSystem(alchemical_system, alchemical_state)

    # Compare energies
    energies = dict()
    platform_names = list()
    for platform_index in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(platform_index)
        platform_name = platform.getName()
        if platform_name != 'Reference':
            platform_names.append(platform_name)
        energies[platform_name] = dict()
        energies[platform_name]['full'] = compute_energy(system, positions, platform=platform)
        set_lambda(alchemical_system, 1.0)
        energies[platform_name]['lambda = 1'] = compute_energy(alchemical_system, positions, platform=platform)
        set_lambda(alchemical_system, 0.0)
        energies[platform_name]['lambda = 0'] = compute_energy(alchemical_system, positions, platform=platform)

    # Check deviations.
    for platform_name in platform_names:
        for energy_name in ['full', 'lambda = 1', 'lambda = 0']:
            delta = energies[platform_name][energy_name] - energies['Reference'][energy_name]
            if (abs(delta) > MAX_DELTA):
                raise Exception("Maximum allowable deviation on platform %s exceeded (was %.8f kcal/mol; allowed %.8f kcal/mol); test failed." % (platform_name, delta / unit.kilocalories_per_mole, MAX_DELTA / unit.kilocalories_per_mole))

def test_annihilated_states(platform_name=None, precision=None):
    """Compare annihilated states in vacuum and a large box.
    """
    from openmmtools import testsystems
    testsystem = testsystems.TolueneVacuum()
    vacuum_system = testsystem.system
    positions = testsystem.positions

    factory_args = {'ligand_atoms' : range(0,15), 'receptor_atoms' : [],
        'annihilate_sterics' : False, 'annihilate_electrostatics' : True }

    # Create annihilated version of vacuum system.
    factory = AbsoluteAlchemicalFactory(vacuum_system, **factory_args)
    vacuum_alchemical_system = factory.createPerturbedSystem()

    # Make copy of system that has periodic boundaries and uses reaction field.
    periodic_system = copy.deepcopy(vacuum_system)
    box_edge = 18.5 * unit.angstroms
    from simtk.openmm import Vec3
    periodic_system.setDefaultPeriodicBoxVectors(Vec3(box_edge,0,0), Vec3(0,box_edge,0), Vec3(0,0,box_edge))
    for force in periodic_system.getForces():
        if force.__class__.__name__ == 'NonbondedForce':
            force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(9.0 * unit.angstroms)
            force.setUseDispersionCorrection(False)
            force.setReactionFieldDielectric(1.0)
    factory = AbsoluteAlchemicalFactory(periodic_system, **factory_args)
    periodic_alchemical_system = factory.createPerturbedSystem()

    # Compare energies
    platform = None
    if platform_name:
        platform = openmm.Platform.getPlatformByName(platform_name)

    vacuum_alchemical_1_energy = compute_energy(vacuum_alchemical_system, positions, platform=platform, precision=precision)
    periodic_alchemical_1_energy = compute_energy(periodic_alchemical_system, positions, platform=platform, precision=precision)

    #compareSystemEnergies(positions, [vacuum_alchemical_system, periodic_alchemical_system], ['vacuum (fully interacting)', 'periodic (fully interacting)'], platform=platform, precision=precision)

    # Set lambda = 0
    lambda_value = 0.0
    alchemical_state = AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value)
    AbsoluteAlchemicalFactory.perturbSystem(vacuum_alchemical_system, alchemical_state)
    AbsoluteAlchemicalFactory.perturbSystem(periodic_alchemical_system, alchemical_state)

    #compareSystemEnergies(positions, [vacuum_alchemical_system, periodic_alchemical_system], ['vacuum (noninteracting)', 'periodic (noninteracting)'], platform=platform, precision=precision)

    vacuum_alchemical_0_energy = compute_energy(vacuum_alchemical_system, positions, platform=platform, precision=precision)
    periodic_alchemical_0_energy = compute_energy(periodic_alchemical_system, positions, platform=platform, precision=precision)

    logger.info('vacuum   lambda = 1 : %8.3f kcal/mol' % (vacuum_alchemical_1_energy / unit.kilocalories_per_mole))
    logger.info('vacuum   lambda = 0 : %8.3f kcal/mol' % (vacuum_alchemical_0_energy / unit.kilocalories_per_mole))
    logger.info('difference          : %8.3f kcal/mol' % ((vacuum_alchemical_1_energy - vacuum_alchemical_0_energy) / unit.kilocalories_per_mole))

    logger.info('periodic lambda = 1 : %8.3f kcal/mol' % (periodic_alchemical_1_energy / unit.kilocalories_per_mole))
    logger.info('periodic lambda = 0 : %8.3f kcal/mol' % (periodic_alchemical_0_energy / unit.kilocalories_per_mole))
    logger.info('difference          : %8.3f kcal/mol' % ((periodic_alchemical_1_energy - periodic_alchemical_0_energy) / unit.kilocalories_per_mole))


def compareSystemEnergies(positions, systems, descriptions, platform=None, precision=None):
    # Compare energies.
    timestep = 1.0 * unit.femtosecond

    if platform:
        platform_name = platform.getName()
        if precision:
            if platform_name == 'CUDA':
                platform.setDefaultPropertyValue('CudaPrecision', precision)
            elif platform_name == 'OpenCL':
                platform.setDefaultPropertyValue('OpenCLPrecision', precision)

    potentials = list()
    states = list()
    for system in systems:
        #dump_xml(system=system)
        integrator = openmm.VerletIntegrator(timestep)
        #dump_xml(integrator=integrator)
        if platform:
            context = openmm.Context(system, integrator, platform)
        else:
            context = openmm.Context(system, integrator)
        context.setPositions(positions)
        state = context.getState(getEnergy=True, getPositions=True)
        #dump_xml(system=system, integrator=integrator, state=state)
        potential = state.getPotentialEnergy()
        potentials.append(potential)
        states.append(state)
        del context, integrator, state

    logger.info("========")
    for i in range(len(systems)):
        logger.info("%32s : %24.8f kcal/mol" % (descriptions[i], potentials[i] / unit.kilocalories_per_mole))
        if (i > 0):
            delta = potentials[i] - potentials[0]
            logger.info("%32s : %24.8f kcal/mol" % ('ERROR', delta / unit.kilocalories_per_mole))
            if (abs(delta) > MAX_DELTA):
                raise Exception("Maximum allowable deviation exceeded (was %.8f kcal/mol; allowed %.8f kcal/mol); test failed." % (delta / unit.kilocalories_per_mole, MAX_DELTA / unit.kilocalories_per_mole))

    return potentials

def alchemical_factory_check(reference_system, positions, platform_name=None, precision=None, factory_args=None):
    """
    Compare energies of reference system and fully-interacting alchemically modified system.

    ARGUMENTS

    reference_system : simtk.openmm.System
       The reference System object to compare with
    positions : simtk.unit.Quantity of dimentsion [natoms,3] with units compatible with angstroms
       The positions to assess energetics for
    precision : str, optional, default=None
       Precision model, or default if not specified. ('single', 'double', 'mixed')
    factory_args : dict(), optional, default=None
       Arguments passed to AbsoluteAlchemicalFactory.

    """

    # Create a factory to produce alchemical intermediates.
    logger.info("Creating alchemical factory...")
    initial_time = time.time()
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    logger.info("AbsoluteAlchemicalFactory initialization took %.3f s" % elapsed_time)
    platform = None
    if platform_name:
        platform = openmm.Platform.getPlatformByName(platform_name)
    alchemical_system = factory.createPerturbedSystem()
    compareSystemEnergies(positions, [reference_system, alchemical_system], ['reference', 'alchemical'], platform=platform, precision=precision)
    return

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

    # Create a factory to produce alchemical intermediates.
    logger.info("Creating alchemical factory...")
    initial_time = time.time()
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    logger.info("AbsoluteAlchemicalFactory initialization took %.3f s" % elapsed_time)

    # Create an alchemically-perturbed state corresponding to nearly fully-interacting.
    # NOTE: We use a lambda slightly smaller than 1.0 because the AlchemicalFactory does not use Custom*Force softcore versions if lambda = 1.0 identically.
    lambda_value = 1.0 - 1.0e-6
    alchemical_state = AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value)

    platform = None
    if platform_name:
        platform = openmm.Platform.getPlatformByName(platform_name)

    # Create the perturbed system.
    logger.info("Creating alchemically-modified state...")
    initial_time = time.time()
    alchemical_system = factory.createPerturbedSystem(alchemical_state)
    final_time = time.time()
    elapsed_time = final_time - initial_time
    # Compare energies.
    logger.info("Computing reference energies...")
    reference_integrator = openmm.VerletIntegrator(timestep)
    if platform:
        reference_context = openmm.Context(reference_system, reference_integrator, platform)
    else:
        reference_context = openmm.Context(reference_system, reference_integrator)
    reference_context.setPositions(positions)
    reference_state = reference_context.getState(getEnergy=True)
    reference_potential = reference_state.getPotentialEnergy()
    logger.info("Computing alchemical energies...")
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
    reference_integrator.step(1)
    alchemical_integrator.step(1)

    # Time simulations.
    logger.info("Simulating reference system...")
    initial_time = time.time()
    reference_integrator.step(nsteps)
    reference_state = reference_context.getState(getEnergy=True)
    reference_potential = reference_state.getPotentialEnergy()
    final_time = time.time()
    reference_time = final_time - initial_time
    logger.info("Simulating alchemical system...")
    initial_time = time.time()
    alchemical_integrator.step(nsteps)
    alchemical_state = alchemical_context.getState(getEnergy=True)
    alchemical_potential = alchemical_state.getPotentialEnergy()
    final_time = time.time()
    alchemical_time = final_time - initial_time

    logger.info("TIMINGS")
    logger.info("reference system       : %12.3f s for %8d steps (%12.3f ms/step)" % (reference_time, nsteps, reference_time/nsteps*1000))
    logger.info("alchemical system      : %12.3f s for %8d steps (%12.3f ms/step)" % (alchemical_time, nsteps, alchemical_time/nsteps*1000))
    logger.info("alchemical simulation is %12.3f x slower than unperturbed system" % (alchemical_time / reference_time))

    return delta

def overlap_check(reference_system, positions, platform_name=None, precision=None, nsteps=50, nsamples=200, factory_args=None, cached_trajectory_filename=None):
    """
    Test overlap between reference system and alchemical system by running a short simulation.

    Parameters
    ----------
    reference_system : simtk.openmm.System
       The reference System object to compare with
    positions : simtk.unit.Quantity with units compatible with nanometers
       The positions to assess energetics for.
    platform_name : str, optional, default=None
       The name of the platform to use for benchmarking.
    nsteps : int, optional, default=50
       Number of molecular dynamics steps between samples.
    nsamples : int, optional, default=100
       Number of samples to collect.
    factory_args : dict(), optional, default=None
       Arguments passed to AbsoluteAlchemicalFactory.
    cached_trajectory_filename : str, optional, default=None
       If specified, attempt to cache (or reuse) trajectory.

    """

    # Create a fully-interacting alchemical state.
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    alchemical_state = AlchemicalState()
    alchemical_system = factory.createPerturbedSystem(alchemical_state)

    temperature = 300.0 * unit.kelvin
    collision_rate = 5.0 / unit.picoseconds
    timestep = 2.0 * unit.femtoseconds
    kT = (kB * temperature)

    # Select platform.
    platform = None
    if platform_name:
        platform = openmm.Platform.getPlatformByName(platform_name)

    # Create integrators.
    reference_integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    alchemical_integrator = openmm.VerletIntegrator(timestep)

    # Create contexts.
    if platform:
        reference_context = openmm.Context(reference_system, reference_integrator, platform)
        alchemical_context = openmm.Context(alchemical_system, alchemical_integrator, platform)
    else:
        reference_context = openmm.Context(reference_system, reference_integrator)
        alchemical_context = openmm.Context(alchemical_system, alchemical_integrator)

    ncfile = None
    if cached_trajectory_filename:
        cache_mode = 'write'

        # Try reading from cache
        from netCDF4 import Dataset
        if os.path.exists(cached_trajectory_filename):
            try:
                ncfile = Dataset(cached_trajectory_filename, 'r')
                if (ncfile.variables['positions'].shape == (nsamples, reference_system.getNumParticles(), 3)):
                    # Read the cache if everything matches
                    cache_mode = 'read'
            except:
                pass

        if cache_mode == 'write':
            # If anything went wrong, create a new cache.
            try:
                (pathname, filename) = os.path.split(cached_trajectory_filename)
                if not os.path.exists(pathname): os.makedirs(pathname)
                ncfile = Dataset(cached_trajectory_filename, 'w', format='NETCDF4')
                ncfile.createDimension('samples', 0)
                ncfile.createDimension('atoms', reference_system.getNumParticles())
                ncfile.createDimension('spatial', 3)
                ncfile.createVariable('positions', 'f4', ('samples', 'atoms', 'spatial'))
            except Exception as e:
                logger.info(str(e))
                logger.info('Could not create a trajectory cache (%s).' % cached_trajectory_filename)
                ncfile = None

    # Collect simulation data.
    reference_context.setPositions(positions)
    du_n = np.zeros([nsamples], np.float64) # du_n[n] is the
    print()
    import click
    with click.progressbar(range(nsamples)) as bar:
        for sample in bar:
            if cached_trajectory_filename and (cache_mode == 'read'):
                # Load cached frames.
                positions = unit.Quantity(ncfile.variables['positions'][sample,:,:], unit.nanometers)
                reference_context.setPositions(positions)
            else:
                # Run dynamics.
                reference_integrator.step(nsteps)

            # Get reference energies.
            reference_state = reference_context.getState(getEnergy=True, getPositions=True)
            reference_potential = reference_state.getPotentialEnergy()
            if np.isnan(reference_potential/kT):
                raise Exception("Reference potential is NaN")

            # Get alchemical energies.
            alchemical_context.setPositions(reference_state.getPositions(asNumpy=True))
            alchemical_state = alchemical_context.getState(getEnergy=True)
            alchemical_potential = alchemical_state.getPotentialEnergy()
            if np.isnan(alchemical_potential/kT):
                raise Exception("Alchemical potential is NaN")

            du_n[sample] = (alchemical_potential - reference_potential) / kT

            if cached_trajectory_filename and (cache_mode == 'write') and (ncfile is not None):
                ncfile.variables['positions'][sample,:,:] = reference_state.getPositions(asNumpy=True) / unit.nanometers

    # Clean up.
    del reference_context, alchemical_context
    if cached_trajectory_filename and (ncfile is not None):
        ncfile.close()

    # Discard data to equilibration and subsample.
    from pymbar import timeseries
    [t0, g, Neff] = timeseries.detectEquilibration(du_n)
    indices = timeseries.subsampleCorrelatedData(du_n, g=g)
    du_n = du_n[indices]

    # Compute statistics.
    from pymbar import EXP
    [DeltaF, dDeltaF] = EXP(du_n)

    # Raise an exception if the error is larger than 3kT.
    MAX_DEVIATION = 3.0 # kT
    if (dDeltaF > MAX_DEVIATION):
        report = "DeltaF = %12.3f +- %12.3f kT (%5d samples, g = %6.1f)" % (DeltaF, dDeltaF, Neff, g)
        raise Exception(report)

    return

def rstyle(ax):
    '''Styles x,y axes to appear like ggplot2
    Must be called after all plot and axis manipulation operations have been
    carried out (needs to know final tick spacing)

    From:
    http://nbviewer.ipython.org/github/wrobstory/climatic/blob/master/examples/ggplot_styling_for_matplotlib.ipynb
    '''
    import pylab
    import matplotlib
    import matplotlib.pyplot as plt

    #Set the style of the major and minor grid lines, filled blocks
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
    ax.grid(True, 'minor', color='0.99', linestyle='-', linewidth=0.7)
    ax.patch.set_facecolor('0.90')
    ax.set_axisbelow(True)

    #Set minor tick spacing to 1/2 of the major ticks
    ax.xaxis.set_minor_locator((pylab.MultipleLocator((plt.xticks()[0][1]
                                -plt.xticks()[0][0]) / 2.0 )))
    ax.yaxis.set_minor_locator((pylab.MultipleLocator((plt.yticks()[0][1]
                                -plt.yticks()[0][0]) / 2.0 )))

    #Remove axis border
    for child in ax.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_alpha(0)

    #Restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)

    #Remove the minor tick lines
    for line in (ax.xaxis.get_ticklines(minor=True) +
                 ax.yaxis.get_ticklines(minor=True)):
        line.set_markersize(0)

    #Only show bottom left ticks, pointing out of axis
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def lambda_trace(reference_system, positions, platform_name=None, precision=None, nsteps=100, factory_args=None):
    """
    Compute potential energy as a function of lambda.

    """

    # Create a factory to produce alchemical intermediates.
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)

    platform = None
    if platform_name:
        # Get platform.
        platform = openmm.Platform.getPlatformByName(platform_name)

    if precision:
        if platform_name == 'CUDA':
            platform.setDefaultPropertyValue('CudaPrecision', precision)
        elif platform_name == 'OpenCL':
            platform.setDefaultPropertyValue('OpenCLPrecision', precision)

    # Take equally-sized steps.
    delta = 1.0 / nsteps

    def compute_potential(system, positions, platform=None):
        timestep = 1.0 * unit.femtoseconds
        integrator = openmm.VerletIntegrator(timestep)
        if platform:
            context = openmm.Context(system, integrator, platform)
        else:
            context = openmm.Context(system, integrator)
        context.setPositions(positions)
        state = context.getState(getEnergy=True)
        potential = state.getPotentialEnergy()
        del integrator, context
        return potential

    # Compute unmodified energy.
    u_original = compute_potential(reference_system, positions, platform)

    # Scan through lambda values.
    lambda_i = np.zeros([nsteps+1], np.float64) # lambda values for u_i
    u_i = unit.Quantity(np.zeros([nsteps+1], np.float64), unit.kilocalories_per_mole) # u_i[i] is the potential energy for lambda_i[i]
    for i in range(nsteps+1):
        lambda_value = 1.0-i*delta # compute lambda value for this step
        alchemical_system = factory.createPerturbedSystem(AlchemicalState(lambda_electrostatics=lambda_value, lambda_sterics=lambda_value, lambda_torsions=lambda_value))
        lambda_i[i] = lambda_value
        u_i[i] = compute_potential(alchemical_system, positions, platform)
        logger.info("%12.9f %24.8f kcal/mol" % (lambda_i[i], u_i[i] / unit.kilocalories_per_mole))

    # Write figure as PDF.
    import pylab
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    with PdfPages('lambda-trace.pdf') as pdf:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        plt.plot(1, u_original / unit.kilocalories_per_mole, 'ro', label='unmodified')
        plt.plot(lambda_i, u_i / unit.kilocalories_per_mole, 'k.', label='alchemical')
        plt.title('T4 lysozyme L99A + p-xylene : AMBER96 + OBC GBSA')
        plt.ylabel('potential (kcal/mol)')
        plt.xlabel('lambda')
        ax.legend()
        rstyle(ax)
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

    return

def generate_trace(test_system):
    lambda_trace(test_system['test'].system, test_system['test'].positions, test_system['receptor_atoms'], test_system['ligand_atoms'])
    return

#=============================================================================================
# TEST SYSTEM DEFINITIONS
#=============================================================================================

test_systems = dict()
test_systems['Lennard-Jones cluster'] = {
    'test' : testsystems.LennardJonesCluster(),
    'factory_args' : {'ligand_atoms' : range(0,1), 'receptor_atoms' : range(1,2) }}
test_systems['Lennard-Jones cluster with modified softcore parameters'] = {
    'test' : testsystems.LennardJonesCluster(),
    'factory_args' : {'ligand_atoms' : range(0,1), 'receptor_atoms' : range(1,2), 'softcore_alpha' : 1, 'softcore_beta' : 1, 'softcore_a' : 2, 'softcore_b' : 2, 'softcore_c' : 2, 'softcore_d' : 2, 'softcore_e' : 2, 'softcore_f' : 2 }}
test_systems['Lennard-Jones fluid without dispersion correction'] = {
    'test' : testsystems.LennardJonesFluid(dispersion_correction=False),
    'factory_args' : {'ligand_atoms' : range(0,1), 'receptor_atoms' : range(1,2) }}
test_systems['Lennard-Jones fluid with dispersion correction'] = {
    'test' : testsystems.LennardJonesFluid(dispersion_correction=True),
    'factory_args' : {'ligand_atoms' : range(0,1), 'receptor_atoms' : range(1,2) }}
test_systems['TIP3P with reaction field, no charges, no switch, no dispersion correction'] = {
    'test' : testsystems.DischargedWaterBox(dispersion_correction=False, switch=False, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6) }}
test_systems['TIP3P with reaction field, switch, no dispersion correction'] = {
    'test' : testsystems.WaterBox(dispersion_correction=False, switch=True, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6) }}
test_systems['TIP3P with reaction field, no switch, dispersion correction'] = {
    'test' : testsystems.WaterBox(dispersion_correction=True, switch=False, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6) }}
test_systems['TIP3P with reaction field, no switch, dispersion correction, no alchemical atoms'] = {
    'test' : testsystems.WaterBox(dispersion_correction=True, switch=False, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : [], 'receptor_atoms' : [] }}
test_systems['TIP3P with reaction field, switch, dispersion correction'] = {
    'test' : testsystems.WaterBox(dispersion_correction=True, switch=True, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6) }}
test_systems['TIP3P with reaction field, switch, dispersion correction, no alchemical atoms'] = {
    'test' : testsystems.WaterBox(dispersion_correction=True, switch=True, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : [], 'receptor_atoms' : [] }}
test_systems['TIP3P with reaction field, switch, dispersion correctionm, electrostatics scaling followed by softcore Lennard-Jones'] = {
    'test' : testsystems.WaterBox(dispersion_correction=True, switch=True, nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6), 'softcore_beta' : 0.0, 'alchemical_functions' : { 'lambda_sterics' : '2*lambda * step(0.5 - lambda)', 'lambda_electrostatics' : '2*(lambda - 0.5) * step(lambda - 0.5)' }}}
test_systems['alanine dipeptide in vacuum'] = {
    'test' : testsystems.AlanineDipeptideVacuum(),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22) }}
test_systems['alanine dipeptide in vacuum with annihilated bonds, angles, and torsions'] = {
    'test' : testsystems.AlanineDipeptideVacuum(),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22),
    'alchemical_torsions' : True, 'alchemical_angles' : True, 'alchemical_bonds' : True }}
test_systems['alanine dipeptide in vacuum with annihilated sterics'] = {
    'test' : testsystems.AlanineDipeptideVacuum(),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22),
    'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}
test_systems['alanine dipeptide in OBC GBSA'] = {
    'test' : testsystems.AlanineDipeptideImplicit(),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22) }}
test_systems['alanine dipeptide in OBC GBSA, with sterics annihilated'] = {
    'test' : testsystems.AlanineDipeptideImplicit(),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22),
    'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}
test_systems['alanine dipeptide in TIP3P with reaction field'] = {
    'test' : testsystems.AlanineDipeptideExplicit(nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,22), 'receptor_atoms' : range(22,22) }}
test_systems['T4 lysozyme L99A with p-xylene in OBC GBSA'] = {
    'test' : testsystems.LysozymeImplicit(),
    'factory_args' : {'ligand_atoms' : range(2603,2621), 'receptor_atoms' : range(0,2603) }}
test_systems['DHFR in explicit solvent with reaction field, annihilated'] = {
    'test' : testsystems.DHFRExplicit(nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,2849), 'receptor_atoms' : [],
    'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}
test_systems['Src in TIP3P with reaction field, with Src sterics annihilated'] = {
    'test' : testsystems.SrcExplicit(nonbondedMethod=app.CutoffPeriodic),
    'factory_args' : {'ligand_atoms' : range(0,4428), 'receptor_atoms' : [],
    'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}
test_systems['Src in GBSA'] = {
    'test' : testsystems.SrcImplicit(),
    'factory_args' : {'ligand_atoms' : range(0,4427), 'receptor_atoms' : [],
    'annihilate_sterics' : False, 'annihilate_electrostatics' : False }}
test_systems['Src in GBSA, with Src sterics annihilated'] = {
    'test' : testsystems.SrcImplicit(),
    'factory_args' : {'ligand_atoms' : range(0,4427), 'receptor_atoms' : [],
    'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}

# Problematic tests: PME is not fully implemented yet
test_systems['TIP3P with PME, no switch, no dispersion correction'] = {
    'test' : testsystems.WaterBox(dispersion_correction=False, switch=False, nonbondedMethod=app.PME),
    'factory_args' : {'ligand_atoms' : range(0,3), 'receptor_atoms' : range(3,6) }}
test_systems['TIP3P with PME, no switch, no dispersion correction, no alchemical atoms'] = {
    'test' : testsystems.WaterBox(dispersion_correction=False, switch=False, nonbondedMethod=app.PME),
    'factory_args' : {'ligand_atoms' : [], 'receptor_atoms' : [] }}

test_systems['toluene in implicit solvent'] = {
    'test' : testsystems.TolueneImplicit(),
    'factory_args' : {'ligand_atoms' : [0,1], 'receptor_atoms' : list(),
    'alchemical_torsions' : True, 'alchemical_angles' : True, 'annihilate_sterics' : True, 'annihilate_electrostatics' : True }}

# Slow tests
#test_systems['Src in OBC GBSA'] = {
#    'test' : testsystems.SrcImplicit(),
#    'ligand_atoms' : range(0,21), 'receptor_atoms' : range(21,7208) }
#test_systems['Src in TIP3P with reaction field'] = {
#    'test' : testsystems.SrcExplicit(nonbondedMethod=app.CutoffPeriodic),
#    'ligand_atoms' : range(0,21), 'receptor_atoms' : range(21,4091) }

accuracy_testsystem_names = [
    'Lennard-Jones cluster',
    'Lennard-Jones fluid without dispersion correction',
    'Lennard-Jones fluid with dispersion correction',
    'TIP3P with reaction field, no charges, no switch, no dispersion correction',
    'TIP3P with reaction field, switch, no dispersion correction',
    'TIP3P with reaction field, switch, dispersion correction',
    'alanine dipeptide in vacuum with annihilated sterics',
    'toluene in implicit solvent',
]

overlap_testsystem_names = [
    'Lennard-Jones cluster',
    'Lennard-Jones fluid without dispersion correction',
    'Lennard-Jones fluid with dispersion correction',
    'TIP3P with reaction field, no charges, no switch, no dispersion correction',
    'TIP3P with reaction field, switch, no dispersion correction',
    'TIP3P with reaction field, switch, dispersion correction',
    'alanine dipeptide in vacuum with annihilated sterics',
    'TIP3P with PME, no switch, no dispersion correction', # PME still lacks reciprocal space component; known energy comparison failure
    'toluene in implicit solvent',
]

#=============================================================================================
# Test various options to AbsoluteAlchemicalFactory
#=============================================================================================

def test_alchemical_functions():
    """
    Testing alchemical slave functions
    """
    alchemical_functions = { 'lambda_sterics' : 'lambda', 'lambda_electrostatics' : 'lambda', 'lambda_bonds' : 'lambda', 'lambda_angles' : 'lambda', 'lambda_torsions' : 'lambda' }
    name = 'Lennard-Jones fluid with dispersion correction'
    test_system = copy.deepcopy(test_systems[name])
    reference_system = test_system['test'].system
    positions = test_system['test'].positions
    factory_args = test_system['factory_args']
    factory_args['alchemical_functions'] = alchemical_functions
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    alchemical_system = factory.createPerturbedSystem()
    compareSystemEnergies(positions, [reference_system, alchemical_system], ['reference', 'alchemical'])

def test_softcore_parameters():
    """
    Testing alchemical slave functions
    """
    alchemical_functions = { 'lambda_sterics' : 'lambda', 'lambda_electrostatics' : 'lambda', 'lambda_bonds' : 'lambda', 'lambda_angles' : 'lambda', 'lambda_torsions' : 'lambda' }
    name = 'Lennard-Jones fluid with dispersion correction'
    test_system = copy.deepcopy(test_systems[name])
    reference_system = test_system['test'].system
    positions = test_system['test'].positions
    factory_args = test_system['factory_args']
    factory_args.update({ 'softcore_alpha' : 1.0, 'softcore_beta' : 1.0, 'softcore_a' : 1.0, 'softcore_b' : 1.0, 'softcore_c' : 1.0, 'softcore_d' : 1.0, 'softcore_e' : 1.0, 'softcore_f' : 1.0 })
    factory = AbsoluteAlchemicalFactory(reference_system, **factory_args)
    alchemical_system = factory.createPerturbedSystem()
    compareSystemEnergies(positions, [reference_system, alchemical_system], ['reference', 'alchemical'])

#=============================================================================================
# NOSETEST GENERATORS
#=============================================================================================

@attr('slow')
def test_overlap():
    """
    Generate nose tests for overlap for all alchemical test systems.
    """
    for name in overlap_testsystem_names:
        test_system = test_systems[name]
        reference_system = test_system['test'].system
        positions = test_system['test'].positions
        factory_args = test_system['factory_args']
        cached_trajectory_filename = os.path.join(os.environ['HOME'], '.cache', 'alchemy', 'tests', name + '.nc')
        f = partial(overlap_check, reference_system, positions, factory_args=factory_args, cached_trajectory_filename=cached_trajectory_filename)
        f.description = "Testing reference/alchemical overlap for %s..." % name
        yield f

def test_alchemical_accuracy():
    """
    Generate nose tests for overlap for all alchemical test systems.
    """
    for name in accuracy_testsystem_names:
        test_system = test_systems[name]
        reference_system = test_system['test'].system
        positions = test_system['test'].positions
        factory_args = test_system['factory_args']
        f = partial(alchemical_factory_check, reference_system, positions, factory_args=factory_args)
        f.description = "Testing alchemical fidelity of %s..." % name
        yield f

def test_alchemical_accuracy():
    """
    Generate nose tests for overlap for all alchemical test systems.
    """
    for name in accuracy_testsystem_names:
        test_system = test_systems[name]
        reference_system = test_system['test'].system
        positions = test_system['test'].positions
        factory_args = test_system['factory_args']
        f = partial(compare_platforms, reference_system, positions, factory_args=factory_args)
        f.description = "Comparing platforms for alchemically-modified forms of %s..." % name
        yield f

#=============================================================================================
# MAIN FOR MANUAL DEBUGGING
#=============================================================================================

if __name__ == "__main__":
    #generate_trace(test_systems['TIP3P with reaction field, switch, dispersion correction'])
    config_root_logger(True)

    logging.basicConfig(level=logging.INFO)
    #test_waterbox()
    test_annihilated_states()


    #name = 'Lennard-Jones fluid with dispersion correction'
    #name = 'Src in GBSA, with Src sterics annihilated'
    #name = 'Src in GBSA'
    #name = 'alanine dipeptide in OBC GBSA, with sterics annihilated'
    #name = 'alanine dipeptide in OBC GBSA'
    name = 'Src in TIP3P with reaction field, with Src sterics annihilated'
    test_system = test_systems[name]
    reference_system = test_system['test'].system
    positions = test_system['test'].positions
    ligand_atoms = test_system['ligand_atoms']
    receptor_atoms = test_system['receptor_atoms']
    alchemical_factory_check(reference_system, positions, receptor_atoms, ligand_atoms)
