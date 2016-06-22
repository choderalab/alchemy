from alchemy import relative
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit
import copy
from pkg_resources import resource_filename
import numpy as np
import os
try:
    from urllib.request import urlopen
    from io import StringIO
except:
    from urllib2 import urlopen
    from cStringIO import StringIO

temperature = 300*unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
# Compute kT and inverse temperature.
kT = kB * temperature
beta = 1.0 / kT

def load_pdbid_to_openmm(pdbid):
    """
    create openmm topology without pdb file
    lifted from pandegroup/pdbfixer
    """
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
    file = urlopen(url)
    contents = file.read().decode('utf-8')
    file.close()
    file = StringIO(contents)

    if _guessFileFormat(file, url) == 'pdbx':
        pdbx = app.PDBxFile(contents)
        topology = pdbx.topology
        positions = pdbx.positions
    else:
        pdb = app.PDBFile(file)
        topology = pdb.topology
        positions = pdb.positions

    return topology, positions

def get_data_filename(relative_path):
    """Get the full path to one of the reference files shipped for testing
    In the source distribution, these files are in ``perses/data/*/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the openmoltools folder).
    """

    fn = resource_filename('perses', relative_path)

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn

def extractPositionsFromOEMOL(molecule):
    positions = unit.Quantity(np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
    coords = molecule.GetCoords()
    for index in range(molecule.NumAtoms()):
        positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
    return positions

def generate_initial_molecule(mol_smiles):
    """
    Generate an oemol with a geometry
    """
    import openeye.oechem as oechem
    import openeye.oeomega as oeomega
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, mol_smiles)
    mol.SetTitle("MOL")
    oechem.OEAddExplicitHydrogens(mol)
    oechem.OETriposAtomNames(mol)
    oechem.OETriposBondTypeNames(mol)
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega(mol)
    return mol

def oemol_to_omm_ff(oemol, molecule_name):
    from perses.rjmc import topology_proposal
    from openmoltools import forcefield_generators
    gaff_xml_filename = get_data_filename('data/gaff.xml')
    system_generator = topology_proposal.SystemGenerator([gaff_xml_filename])
    topology = forcefield_generators.generateTopologyFromOEMol(oemol)
    system = system_generator.build_system(topology)
    positions = extractPositionsFromOEMOL(oemol)
    return system, positions, topology

def _guessFileFormat(file, filename):
    """
    Guess whether a file is PDB or PDBx/mmCIF based on its filename and contents.
    authored by pandegroup
    """
    filename = filename.lower()
    if '.pdbx' in filename or '.cif' in filename:
        return 'pdbx'
    if '.pdb' in filename:
        return 'pdb'
    for line in file:
        if line.startswith('data_') or line.startswith('loop_'):
            file.seek(0)
            return 'pdbx'
        if line.startswith('HEADER') or line.startswith('REMARK') or line.startswith('TITLE '):
            file.seek(0)
            return 'pdb'
    file.seek(0)
    return 'pdb'

def test_run_point_mutation_propose():
    """
    Propose a random mutation in insulin
    """
    import perses.rjmc.topology_proposal as topology_proposal

    pdbid = "2HIU"
    topology, positions = load_pdbid_to_openmm(pdbid)
    modeller = app.Modeller(topology, positions)
    for chain in modeller.topology.chains():
        pass

    modeller.delete([chain])

    ff_filename = "amber99sbildn.xml"
    max_point_mutants = 1

    ff = app.ForceField(ff_filename)
    system = ff.createSystem(modeller.topology)
    chain_id = 'A'

    system_generator = topology_proposal.SystemGenerator([ff_filename])

    pm_top_engine = topology_proposal.PointMutationEngine(modeller.topology, system_generator, chain_id, max_point_mutants=max_point_mutants)
    return pm_top_engine.propose(system, modeller.topology)

def test_small_molecule_proposals():
    """
    Make sure the small molecule proposal engine generates molecules
    """
    from perses.rjmc import topology_proposal
    from openmoltools import forcefield_generators
    import openeye.oechem as oechem
    list_of_smiles = ['CCCC','CCCCC','CCCCCC']
    gaff_xml_filename = get_data_filename('data/gaff.xml')
    stats_dict = {smiles : 0 for smiles in list_of_smiles}
    system_generator = topology_proposal.SystemGenerator([gaff_xml_filename])
    proposal_engine = topology_proposal.SmallMoleculeSetProposalEngine(list_of_smiles, system_generator)
    initial_molecule = generate_initial_molecule('CCCC')
    initial_system, initial_positions, initial_topology = oemol_to_omm_ff(initial_molecule, "MOL")
    return proposal_engine.propose(initial_system, initial_topology)

def _relative_factory(top_proposal):
    system1 = top_proposal.old_system
    system2 = top_proposal.new_system
    topology1 = top_proposal.old_topology
    topology2 = top_proposal.new_topology
    positions1 = np.zeros((topology1.getNumAtoms(), 3))
    positions1 = unit.Quantity(value=positions1, unit=unit.nanometer)
    positions2 = np.zeros((topology2.getNumAtoms(), 3))
    positions2 = unit.Quantity(value=positions2, unit=unit.nanometer)
    atom_mapping_1to2 = top_proposal.old_to_new_atom_map

    atom_list = list(topology1.atoms())
    atom_list2 = list(topology2.atoms())
    printed_resmap = False
    for atom in atom_list:
        try:
            atom2 = atom_mapping_1to2[atom.index]
            core = True
        except:
            core = False
        if core:
            atom2 = atom_list2[atom2]
            if not atom.residue.name == atom2.residue.name:
                if not printed_resmap:
                    print(atom.residue.name, atom2.residue.name)
                    for resatom in atom.residue.atoms():
                        try:
                            print(resatom.name, atom_list2[atom_mapping_1to2[resatom.index]].name)
                        except:
                            pass
                    printed_resmap = True

    print("New topology chemical state key:")
    print(top_proposal.new_chemical_state_key)
    print("Old topology chemical state key:")
    print(top_proposal.old_chemical_state_key)

    if topology1 == topology2:
        return

    alchemical_factory = relative.HybridTopologyFactory(system1, system2, topology1, topology2, positions1, positions2, atom_mapping_1to2)
    hybrid_system, hybrid_topology, hybrid_positions = alchemical_factory.createPerturbedSystem()

def test_relative_factory_point_mutation():
    top_proposal = test_run_point_mutation_propose()
    _relative_factory(top_proposal)

def test_relative_factory_small_molecule():
    top_proposal = test_small_molecule_proposals()
    _relative_factory(top_proposal)


if __name__ == "__main__":
    for i in range(50):
        test_relative_factory_point_mutation()
        test_relative_factory_small_molecule()





