import os
import tkinter as tk
from tkinter import ttk

import numpy as np

from DFT import init_calculation, calc_ESP, calc_LMO, create_xyz_string
from utils import draw_isosurface

from cube_tools import cube


class CalculationWindow:
    def __init__(self, root, name, smiles, n_params, sn_params, charge):
        self.root = root
        self.root.title("Calculation Window")

        # Create the XYZ file for the molecule
        self.xyz = create_xyz_string(name, smiles)

        # Store charge as an instance attribute
        self.charge = charge

        # Initial DFT calculation
        self.mol, self.mf = init_calculation(self.xyz, self.charge)

        # Display the molecule's information
        name_label = tk.Label(root, text=f"Name: {name}", font=('Arial', 14))
        name_label.pack(pady=10)

        n_params_label = tk.Label(root, text=f"N Params: {n_params}", font=('Arial', 12))
        n_params_label.pack(pady=5)

        sn_params_label = tk.Label(root, text=f"sN Params: {sn_params}", font=('Arial', 12))
        sn_params_label.pack(pady=5)

        # Dropdown for selecting calculation type
        self.calc_type = tk.StringVar()
        calc_dropdown_label = tk.Label(root, text="Calculation Type:", font=('Arial', 12))
        calc_dropdown_label.pack(pady=5)

        calc_dropdown = ttk.Combobox(root, textvariable=self.calc_type, font=('Arial', 12))
        calc_dropdown['values'] = ['ESP', 'LMO']
        calc_dropdown.set('Select Calculation')
        calc_dropdown.pack(pady=5)
        calc_dropdown.bind("<<ComboboxSelected>>", self.update_inputs)  # Update inputs on selection

        # Frame to hold dynamic inputs
        self.dynamic_frame = tk.Frame(root)
        self.dynamic_frame.pack(pady=10)

        # "iso_value" input (permanent)
        iso_value_label = tk.Label(root, text="Iso Value:", font=('Arial', 12))
        iso_value_label.pack(pady=5)
        self.iso_value = tk.Entry(root)
        self.iso_value.insert(0, "0.3")  # Default iso value
        self.iso_value.pack(pady=5)

        # Start calculation button
        calc_button = tk.Button(root, text="Start Calculation", command=self.start_calculation)
        calc_button.pack(pady=10)

    def update_inputs(self, event=None):
        # Clear previous inputs
        for widget in self.dynamic_frame.winfo_children():
            widget.destroy()

        calc_type = self.calc_type.get()

        if calc_type == "ESP":
            # Inputs for ESP
            esp_label = tk.Label(self.dynamic_frame, text="ESP Grid Spacing:", font=('Arial', 12))
            esp_label.pack(pady=5)
            self.esp_spacing = tk.Entry(self.dynamic_frame)
            self.esp_spacing.pack(pady=5)

        elif calc_type == "LMO":
            # Inputs for LMO
            lmo_label = tk.Label(self.dynamic_frame, text="Localization Method:", font=('Arial', 12))
            lmo_label.pack(pady=5)

            # Set default method to "Boys"
            self.lmo_method = ttk.Combobox(self.dynamic_frame, values=["Boys", "Edmiston-Ruedenberg"])
            self.lmo_method.set("Boys")
            self.lmo_method.pack(pady=5)

            # Add integer input for orbital number
            lmo_int_label = tk.Label(self.dynamic_frame, text="Orbital Number:", font=('Arial', 12))
            lmo_int_label.pack(pady=5)

            # Create the Entry widget for the orbital number with default value 0
            self.lmo_orbital_number = tk.Entry(self.dynamic_frame)
            self.lmo_orbital_number.insert(0, "0")
            self.lmo_orbital_number.pack(pady=5)


    def start_calculation(self):    
        calc_type = self.calc_type.get()
        cube_filename = "mol.cube"

        # Get iso_value from the entry widget
        try:
            iso_value = float(self.iso_value.get())  # Convert the iso value to float
        except ValueError:
            print("Please enter a valid number for the Iso Value.")
            return  # Stop the calculation if iso_value is invalid

        if calc_type == "ESP":
            esp_spacing = self.esp_spacing.get()
            print(f"Running ESP calculation with grid spacing: {esp_spacing}")
            parsed_cube = calc_ESP(self.mf, self.mol, cube_filename)

            draw_isosurface(parsed_cube, iso_value)

        elif calc_type == "LMO":
            lmo_method = self.lmo_method.get()
            print(f"Running LMO calculation using method: {lmo_method}")
            
            # Retrieve the orbital number, convert to integer
            try:
                orbital_number = int(self.lmo_orbital_number.get())  # Convert string input to integer
                parsed_cube = calc_LMO(self.mf, self.mol, cube_filename, orbital_number, lmo_method)

                draw_isosurface(parsed_cube, iso_value)
            except ValueError:
                print("Please enter a valid integer for the orbital number.")
