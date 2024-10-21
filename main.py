import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import tkinter as tk
from PIL import Image, ImageTk

# Funktion zum Laden der CSV-Datei
def load_data(file_path):
    global extracted_data
    df = pd.read_csv(file_path, sep=';', engine='python', quotechar='"')
    df['Name'] = df['Name'].str.split(r' \(', n=1).str[0]
    df = df.dropna(subset=['Smiles'])
    columns_of_interest = ['Name', 'Smiles', 'N Params', 'sN Params']
    extracted_data = df[columns_of_interest]
    show_molecule(current_index)

# Funktion zum Anzeigen des Moleküls
def show_molecule(index):
    # Check if there is any data in extracted_data
    if extracted_data.empty:
        name_label.config(text="No data loaded")
        n_params_label.config(text="")
        sn_params_label.config(text="")
        img_label.config(image='')  # Clear the image label
        return

    # Ensure the index is within bounds
    index = index % len(extracted_data)

    smiles = extracted_data.iloc[index]['Smiles']
    name = extracted_data.iloc[index]['Name']
    n_params = extracted_data.iloc[index]['N Params']
    sn_params = extracted_data.iloc[index]['sN Params']

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        img.save('current_molecule.png')  # Speichere das Bild zwischen
        display_image()

    name_label.config(text=name)
    n_params_label.config(text=f"N Params: {n_params}")
    sn_params_label.config(text=f"sN Params: {sn_params}")


# Funktion zum Anzeigen des Bildes
def display_image():
    img = Image.open('current_molecule.png')
    img = img.resize((300, 300))  # Größe anpassen
    img_tk = ImageTk.PhotoImage(img)
    img_label.config(image=img_tk)
    img_label.image = img_tk  # Referenz halten

# Button-Funktion für nächste und vorherige Moleküle
def next_molecule():
    global current_index
    current_index = (current_index + 1) % len(extracted_data)
    show_molecule(current_index)

def previous_molecule():
    global current_index
    current_index = (current_index - 1) % len(extracted_data)
    show_molecule(current_index)

# Button-Funktion zum Wechseln der Dateien
def switch_file(file_path):
    global current_index
    current_index = 0  # Setze den Index zurück
    load_data(file_path)

# Hauptfenster erstellen
root = tk.Tk()
root.title("Molecule Viewer")

# Initialer Index
current_index = 0
extracted_data = pd.DataFrame()  # Leeres DataFrame initialisieren

# Bildlabel
img_label = tk.Label(root)
img_label.pack(pady=20)

# Namenlabel
name_label = tk.Label(root, font=('Arial', 14))
name_label.pack(pady=10)

# N Params Label
n_params_label = tk.Label(root, font=('Arial', 12))
n_params_label.pack(pady=5)

# sN Params Label
sn_params_label = tk.Label(root, font=('Arial', 12))
sn_params_label.pack(pady=5)

# Buttons zum Wechseln der Dateien
n_button = tk.Button(root, text="N-Nucleophiles", command=lambda: switch_file('N_Nucleophiles/Nucleophile.csv'))
n_button.pack(side='top', padx=20, pady=5)

c_button = tk.Button(root, text="C-Nucleophiles", command=lambda: switch_file('C_Nucleophiles/Nucleophile.csv'))
c_button.pack(side='top', padx=20, pady=5)

p_button = tk.Button(root, text="O-Nucleophiles", command=lambda: switch_file('O_Nucleophiles/Nucleophile.csv'))
p_button.pack(side='top', padx=20, pady=5)

# Button für vorheriges Molekül
prev_button = tk.Button(root, text="Vorheriges Molekül", command=previous_molecule)
prev_button.pack(side='left', padx=20)

# Button für nächstes Molekül
next_button = tk.Button(root, text="Nächstes Molekül", command=next_molecule)
next_button.pack(side='right', padx=20)

# Zeige das erste Molekül (anfänglich leer)
show_molecule(current_index)

# Hauptloop starten
root.mainloop()
