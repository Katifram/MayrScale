import tkinter as tk
from rdkit import Chem
from rdkit.Chem import AllChem


class CalculationWindow:
    def __init__(self, root, name, smiles, n_params, sn_params):
        self.root = root
        self.root.title("Calculation Window")

        # Create the XYZ file for the molecule
        self.create_xyz_file(name, smiles)

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


    def create_xyz_file(self, name, smiles):
        """Generates a 3D structure and writes it to an XYZ file."""
        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule object
        mol = Chem.AddHs(mol)  # Add hydrogen atoms

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

        # Write to XYZ format
        with open(f"xyz/{name}.xyz", "w") as f:
            conf = mol.GetConformer()  # Get the conformation of the molecule
            num_atoms = mol.GetNumAtoms()
            f.write(f"{num_atoms}\n")  # First line: number of atoms
            f.write(f"{name}\n")  # Second line: molecule name or comment
            for i in range(num_atoms):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

    def start_calculation(self):
        # Placeholder for calculation logic
        result_label = tk.Label(self.root, text="Calculations started...")
        result_label.pack(pady=10)
