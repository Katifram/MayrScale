import tkinter as tk
from molecule_viewer import MoleculeViewer

# Hauptfenster erstellen
root = tk.Tk()
root.title("Molecule Viewer")

# MoleculeViewer Instanz erzeugen
viewer = MoleculeViewer(root)

# Hauptloop starten
root.mainloop()
