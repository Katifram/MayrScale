from rdkit import Chem
from rdkit.Chem import AllChem


# Function to generate 3D coordinates and save as .xyz file
def save_as_xyz(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return

    mol = Chem.AddHs(mol)  # Add hydrogens
    success = AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    if success == -1:
        print(f"Failed to generate 3D coordinates for: {smiles}")
        return

    AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry using UFF force field
    conformer = mol.GetConformer()

    # Writing the XYZ file
    with open(filename, 'w') as xyz_file:
        xyz_file.write(f"{mol.GetNumAtoms()}\n")  # Atom count
        xyz_file.write(f"{smiles}\n")  # Comment line (SMILES or molecule name)

        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            xyz_file.write(f"{symbol} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

    print(f"XYZ file saved as: {filename}")