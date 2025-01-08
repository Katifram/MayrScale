from pyscf.tools import cubegen

from utils import draw_isosurface, parse_cube
from pyscf import gto, lo

import pyscf
from pyscf import gto
from pyscf import scf
from pyscf.geomopt.geometric_solver import optimize

from rdkit import Chem
from rdkit.Chem import AllChem

def create_xyz_string(name, smiles):  # Add 'self' as the first parameter
        """Generates a 3D structure and returns it in XYZ format as a string without atom count or molecule name."""

        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule object
        mol = Chem.AddHs(mol)  # Add hydrogen atoms

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        # Prepare XYZ string format without the number of atoms or molecule name
        conf = mol.GetConformer()  # Get the conformation of the molecule
        xyz_string = ""
        for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                # Format each line with increased spacing and align values for a cleaner look
                xyz_string += f"{atom.GetSymbol():<2}    {pos.x: .10f}    {pos.y: .10f}    {pos.z: .10f}\n"

        return xyz_string.strip()  # Remove any trailing newline

def init_calculation(xyz, charge):

    atom = xyz

    # Build the molecule without GPU-specific settings
    mol = gto.M(atom=atom, basis='def2-tzvpp', charge=charge)

    # Perform SCF calculation without GPU acceleration
    mf = scf.RHF(mol)

    # optimize geometry
    mol_opt = optimize(mf, maxsteps=1)
    mf = scf.RHF(mol_opt)

    mf.kernel()

    return mol, mf

def calc_ESP(mf, mol, cube_filename):
    # No GPU processing, directly use the density matrix on CPU
    rdm1_cpu = mf.make_rdm1()  # Directly using the density matrix on the CPU
               
    cubegen.density(mol, cube_filename, rdm1_cpu)
    parsed_cube = parse_cube(cube_filename)

    return parsed_cube

def calc_LMO(mf, mol, cube_filename, orbital_number, method):
    mo_coeff_cpu = mf.mo_coeff

    # Select the localization method
    if method == "Boys":
        localizer = lo.Boys(mol, mo_coeff_cpu)
    elif method == "Edmiston-Ruedenberg":
        localizer = lo.EdmistonRuedenberg(mol, mo_coeff_cpu)
    else:
        raise ValueError("Invalid localization method. Choose 'Boys' or 'Edmiston-Ruedenberg'.")

    # Localize orbitals
    localized_mos = localizer.kernel()

    # Analyze localized orbitals to identify lone pairs
    cubegen.orbital(mol, cube_filename, localized_mos[:, orbital_number], nx=80)
    parsed_cube = parse_cube(cube_filename)

    return parsed_cube
