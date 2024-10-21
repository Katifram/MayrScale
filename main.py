import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import tkinter as tk
from PIL import Image, ImageTk

# Lade die CSV-Datei mit Semikolon als Trennzeichen
file_path = 'Nucleophile.csv'
df = pd.read_csv(file_path, sep=';', engine='python', quotechar='"')

# Bereinige die Spalte "Name", indem alles nach " (" entfernt wird
df['Name'] = df['Name'].str.split(' \(', n=1).str[0]

# Verwerfe Zeilen ohne SMILES
df = df.dropna(subset=['Smiles'])

# Extrahiere die gewünschten Spalten
columns_of_interest = ['Name', 'Smiles', 'N Params', 'sN Params']
extracted_data = df[columns_of_interest]

# Funktion zum Anzeigen des Moleküls
def show_molecule(index):
    smiles = extracted_data.iloc[index]['Smiles']
    name = extracted_data.iloc[index]['Name']
    n_params = extracted_data.iloc[index]['N Params']
    sn_params = extracted_data.iloc[index]['sN Params']

    # Erzeuge das Molekül-Objekt aus dem SMILES-String
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Erzeuge das Molekülbild
        img = Draw.MolToImage(mol)
        img.save('current_molecule.png')  # Speichere das Bild zwischen
        display_image()

    # Update den Namen und die Parameter
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
    current_index = (current_index + 1) % len(extracted_data)  # Zyklisch
    show_molecule(current_index)

def previous_molecule():
    global current_index
    current_index = (current_index - 1) % len(extracted_data)  # Zyklisch
    show_molecule(current_index)

# Hauptfenster erstellen
root = tk.Tk()
root.title("Molecule Viewer")

# Initialer Index
current_index = 0

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

# Button für vorheriges Molekül
prev_button = tk.Button(root, text="Vorheriges Molekül", command=previous_molecule)
prev_button.pack(side='left', padx=20)

# Button für nächstes Molekül
next_button = tk.Button(root, text="Nächstes Molekül", command=next_molecule)
next_button.pack(side='right', padx=20)

# Zeige das erste Molekül
show_molecule(current_index)

# Hauptloop starten
root.mainloop()
