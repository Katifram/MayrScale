import os
import tkinter as tk
from rdkit import Chem
from rdkit.Chem import AllChem

import numpy as np

import pyscf
from pyscf import gto
from pyscf.tools import cubegen
from pyscf.geomopt.geometric_solver import optimize

from gpu4pyscf import scf as gpu_scf
from gpu4pyscf.dft import rks, numint, gen_grid

from utils import parse_cube
from utils import draw_isosurface

from cube_tools import cube


class CalculationWindow:
    def __init__(self, root, name, smiles, n_params, sn_params, charge):
        self.root = root
        self.root.title("Calculation Window")

        # Create the XYZ file for the molecule
        self.xyz = self.create_xyz_string(name, smiles)

        # Store charge as an instance attribute
        self.charge = charge

        # Display the molecule's information
        name_label = tk.Label(root, text=f"Name: {name}", font=('Arial', 14))
        name_label.pack(pady=10)

        n_params_label = tk.Label(root, text=f"N Params: {n_params}", font=('Arial', 12))
        n_params_label.pack(pady=5)

        sn_params_label = tk.Label(root, text=f"sN Params: {sn_params}", font=('Arial', 12))
        sn_params_label.pack(pady=5)

        # Add buttons or fields for calculation logic here
        calc_button = tk.Button(root, text="Start Calculation", command=self.start_calculation)
        calc_button.pack(pady=10)

    def create_xyz_string(self, name, smiles):  # Add 'self' as the first parameter
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

    def start_calculation(self):

        atom = self.xyz

        print(atom)
        # Build the molecule using GPU-enabled settings
        mol = gto.M(atom=atom, basis='def2-tzvpp', charge=self.charge)

        # Perform SCF calculation using GPU-accelerated SCF
        mf = gpu_scf.RHF(mol)

        mol_opt = optimize(mf, maxsteps=100)
        mf = gpu_scf.RHF(mol_opt)

        mf.kernel()

        # Convert GPU-accelerated density matrix to CPU-compatible format
        rdm1_cpu = mf.make_rdm1().get()  # Use .get() to convert CuPy to NumPy  # Convert to a NumPy array on CPU

        cube_filename = "mol.cube"

        cubegen.density(mol, cube_filename, rdm1_cpu)
        draw_isosurface(cube_filename)

        

        

