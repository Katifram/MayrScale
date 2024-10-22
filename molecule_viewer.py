import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageTk
import tkinter as tk

from utils import save_as_xyz  # Import save_as_xyz function from utils

class MoleculeViewer:
    def __init__(self, root):
        self.root = root
        self.current_index = 0
        self.extracted_data = pd.DataFrame()  # Empty DataFrame

        # GUI Elements
        self.img_label = tk.Label(root)
        self.img_label.pack(pady=20)

        self.name_label = tk.Label(root, font=('Arial', 14))
        self.name_label.pack(pady=10)

        self.n_params_label = tk.Label(root, font=('Arial', 12))
        self.n_params_label.pack(pady=5)

        self.sn_params_label = tk.Label(root, font=('Arial', 12))
        self.sn_params_label.pack(pady=5)

        # Buttons for switching files
        tk.Button(root, text="N-Nucleophiles", command=lambda: self.switch_file('data/N_Nucleophiles/Nucleophile.csv')).pack(side='top', padx=20, pady=5)
        tk.Button(root, text="C-Nucleophiles", command=lambda: self.switch_file('data/C_Nucleophiles/Nucleophile.csv')).pack(side='top', padx=20, pady=5)
        tk.Button(root, text="O-Nucleophiles", command=lambda: self.switch_file('data/O_Nucleophiles/Nucleophile.csv')).pack(side='top', padx=20, pady=5)

        # Navigation buttons
        tk.Button(root, text="Vorheriges Molekül", command=self.previous_molecule).pack(side='left', padx=20)
        tk.Button(root, text="Nächstes Molekül", command=self.next_molecule).pack(side='right', padx=20)

        # Initially show an empty molecule
        self.show_molecule(self.current_index)

    def load_data(self, file_path):
        df = pd.read_csv(file_path, sep=';', engine='python', quotechar='"')
        df['Name'] = df['Name'].str.split(r' \(', n=1).str[0]
        df = df.dropna(subset=['Smiles'])
        columns_of_interest = ['Name', 'Smiles', 'N Params', 'sN Params']
        self.extracted_data = df[columns_of_interest]
        self.show_molecule(self.current_index)

    def show_molecule(self, index):
        if self.extracted_data.empty:
            self.name_label.config(text="No data loaded")
            self.n_params_label.config(text="")
            self.sn_params_label.config(text="")
            self.img_label.config(image='')  # Clear the image label
            return

        # Ensure the index is within bounds
        index = index % len(self.extracted_data)
        smiles = self.extracted_data.iloc[index]['Smiles']
        name = self.extracted_data.iloc[index]['Name']
        n_params = self.extracted_data.iloc[index]['N Params']
        sn_params = self.extracted_data.iloc[index]['sN Params']

        # Create molecule image
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            img.save('resources/current_molecule.png')  # Save the image temporarily
            self.display_image()

        # Display molecule data
        self.name_label.config(text=name)
        self.n_params_label.config(text=f"N Params: {n_params}")
        self.sn_params_label.config(text=f"sN Params: {sn_params}")

    def display_image(self):
        img = Image.open('resources/current_molecule.png')
        img = img.resize((300, 300))  # Resize
        img_tk = ImageTk.PhotoImage(img)
        self.img_label.config(image=img_tk)
        self.img_label.image = img_tk  # Keep a reference

    def next_molecule(self):
        self.current_index = (self.current_index + 1) % len(self.extracted_data)
        self.show_molecule(self.current_index)

    def previous_molecule(self):
        self.current_index = (self.current_index - 1) % len(self.extracted_data)
        self.show_molecule(self.current_index)

    def switch_file(self, file_path):
        self.current_index = 0  # Reset index
        self.load_data(file_path)