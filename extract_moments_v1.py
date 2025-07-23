# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 10:43:16 2025

extract_moments_v1.py

Version 1

User gives desired reference isotope and target isotope and scan.
Program finds reference data and scan/target data, then calculates and outputs target nuclear moments.

Author: Shay Rowlinson
"""

import os
import pandas as pd
from tkinter import ttk, StringVar, messagebox, Entry

class Moments(ttk.Frame):
    """
    Tab for extracting nuclear moments (μ, Q) from fitted hyperfine parameters, using a reference isotope.
    User loads a reference isotope and a scan, and moments are calculated automatically.
    """
    def __init__(self, parent):
        """
        Initialise the moments extraction tab and widgets.
        """
        super().__init__(parent)
        self._build_gui()
        self.ref_data = None
        self.scan_data = None
        self.I_scan = None  # Nuclear spin for scan isotope
        

    def _build_gui(self):
        """
        Build the layout for the reference/scan input and calculated moments output.
        """
        # Reference Isotope Frame 
        ref_frame = ttk.LabelFrame(self, text="Reference Isotope", padding=(10,6))
        ref_frame.grid(row=0, column=0, sticky="ew", padx=2, pady=(10,0))
        
        # Top row: Element, Mass, Load
        ttk.Label(ref_frame, text="Element:").grid(row=0, column=0)
        self.ref_elem_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_elem_var, width=8).grid(row=0, column=1)
        ttk.Label(ref_frame, text="Mass #:").grid(row=0, column=2)
        self.ref_mass_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_mass_var, width=8).grid(row=0, column=3)
        ttk.Button(ref_frame, text="Load Reference", command=self.load_reference).grid(row=0, column=4, padx=12)
        
        # Second row:  μ_ref, Q_ref I_ref
        ttk.Label(ref_frame, text="μ_ref / μN:").grid(row=1, column=0)
        self.ref_mu_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_mu_var, width=11, state="readonly").grid(row=1, column=1)
        ttk.Label(ref_frame, text="Q_ref / b:").grid(row=1, column=2)
        self.ref_Q_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_Q_var, width=11, state="readonly").grid(row=1, column=3)
        ttk.Label(ref_frame, text="I_ref:").grid(row=1, column=4)
        self.ref_I_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_I_var, width=11, state="readonly").grid(row=1, column=5)
        
        # Third row: A_refs and B_refs 
        ttk.Label(ref_frame, text="A_u ref:").grid(row=2, column=0)
        self.ref_Au_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_Au_var, width=11, state="readonly").grid(row=2, column=1)
        ttk.Label(ref_frame, text="A_l ref:").grid(row=2, column=2)
        self.ref_Al_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_Al_var, width=11, state="readonly").grid(row=2, column=3)
        ttk.Label(ref_frame, text="B_u ref:").grid(row=2, column=4)
        self.ref_Bu_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_Bu_var, width=11, state="readonly").grid(row=2, column=5)
        ttk.Label(ref_frame, text="B_l ref:").grid(row=2, column=6)
        self.ref_Bl_var = StringVar()
        Entry(ref_frame, textvariable=self.ref_Bl_var, width=11, state="readonly").grid(row=2, column=7)
    
        # Target Scan Frame 
        scan_frame = ttk.LabelFrame(self, text="Target Scan", padding=(10,6))
        scan_frame.grid(row=1, column=0, sticky="ew", padx=2)
        
        # Top row: Element, Mass, Scan, Load
        ttk.Label(scan_frame, text="Element:").grid(row=0, column=0)
        self.elem_var = StringVar()
        Entry(scan_frame, textvariable=self.elem_var, width=8).grid(row=0, column=1)
        ttk.Label(scan_frame, text="Mass #:").grid(row=0, column=2)
        self.mass_var = StringVar()
        Entry(scan_frame, textvariable=self.mass_var, width=8).grid(row=0, column=3)
        ttk.Label(scan_frame, text="Scan #:").grid(row=0, column=4)
        self.scan_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_var, width=8).grid(row=0, column=5)
        ttk.Button(scan_frame, text="Load Scan", command=self.load_scan).grid(row=0, column=6, padx=12)
    
        # Second row: I, A_u, A_l
        ttk.Label(scan_frame, text="I:").grid(row=1, column=0)
        self.scan_I_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_I_var, width=11, state="readonly").grid(row=1, column=1)
        ttk.Label(scan_frame, text="A_u:").grid(row=1, column=2)
        self.scan_Au_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_Au_var, width=11, state="readonly").grid(row=1, column=3)
        ttk.Label(scan_frame, text="A_l:").grid(row=1, column=4)
        self.scan_Al_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_Al_var, width=11, state="readonly").grid(row=1, column=5)
    
        # Third row: B_u, B_l
        ttk.Label(scan_frame, text="B_u:").grid(row=2, column=0)
        self.scan_Bu_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_Bu_var, width=11, state="readonly").grid(row=2, column=1)
        ttk.Label(scan_frame, text="B_l:").grid(row=2, column=2)
        self.scan_Bl_var = StringVar()
        Entry(scan_frame, textvariable=self.scan_Bl_var, width=11, state="readonly").grid(row=2, column=3)
    
        # Calculated Moments 
        out_frame = ttk.LabelFrame(self, text="Calculated Moments", padding=(10,6))
        out_frame.grid(row=2, column=0, sticky="ew", padx=2, pady=10)
        ttk.Label(out_frame, text="μ_l / μN:").grid(row=0, column=0, sticky='e')
        self.mu_l_var = StringVar()
        Entry(out_frame, textvariable=self.mu_l_var, width=14, state="readonly").grid(row=0, column=1)
        ttk.Label(out_frame, text="μ_u / μN:").grid(row=0, column=2, sticky='e')
        self.mu_u_var = StringVar()
        Entry(out_frame, textvariable=self.mu_u_var, width=14, state="readonly").grid(row=0, column=3)
        ttk.Label(out_frame, text="Q_l / b:").grid(row=1, column=0, sticky='e')
        self.Q_l_var = StringVar()
        Entry(out_frame, textvariable=self.Q_l_var, width=14, state="readonly").grid(row=1, column=1)
        ttk.Label(out_frame, text="Q_u / b:").grid(row=1, column=2, sticky='e')
        self.Q_u_var = StringVar()
        Entry(out_frame, textvariable=self.Q_u_var, width=14, state="readonly").grid(row=1, column=3)


    def load_reference(self):
        """
        Load reference isotope parameters (μ_ref, Q_ref, I_ref, A_ref, B_ref) from a CSV, and fill GUI fields.
        """
        elem = self.ref_elem_var.get().strip()
        mass = self.ref_mass_var.get().strip()
        if not elem or not mass:
            messagebox.showerror("Missing Input", "Enter both element and mass number for reference isotope.")
            return
        path = os.path.join(".", "references", f"{elem}_{mass}.csv")
        try:
            df = pd.read_csv(path)
            row = df.iloc[0]
            self.ref_I_var.set(str(row.get("I", "")))
            self.ref_mu_var.set(str(row.get("mu", "")))
            self.ref_Q_var.set(str(row.get("Q", "")))
            self.ref_Au_var.set(str(row.get("A_u", "")))
            self.ref_Al_var.set(str(row.get("A_l", "")))
            self.ref_Bu_var.set(str(row.get("B_u", "")))
            self.ref_Bl_var.set(str(row.get("B_l", "")))
            
            # Store floats for calculation:
            self.ref_data = {
                "I_ref": self._parse_fraction(row.get("I", "")),
                "mu_ref": float(row.get("mu", 0)),
                "Q_ref": float(row.get("Q", 0)),
                "A_l_ref": float(row.get("A_l", 0)),
                "A_u_ref": float(row.get("A_u", 0)),
                "B_l_ref": float(row.get("B_l", 0)),
                "B_u_ref": float(row.get("B_u", 0)),
            }
            self.try_calculate()
        except Exception as e:
            messagebox.showerror("Reference Load Error", f"Could not load reference:\n{e}")

    def load_scan(self):
        """
        Load target scan fit parameters A, B and nuclear spin from CSV, and fill GUI fields.
        """
        elem = self.elem_var.get().strip()
        mass = self.mass_var.get().strip()
        scan = self.scan_var.get().strip()
        if not elem or not mass or not scan:
            messagebox.showerror("Missing Input", "Enter element, mass, and scan number for target scan.")
            return
        # Load scan fit parameters
        path = os.path.join(".", f"{elem.capitalize()}_Results", str(mass), "Saved_Parameters.csv")
        try:
            df = pd.read_csv(path)
            match = df[df['scan'].astype(str) == str(scan)]
            if match.empty:
                raise ValueError(f"Scan {scan} not found in {path}")
            row = match.iloc[0]
            self.scan_Au_var.set(str(row.get("A_u", "")))
            self.scan_Al_var.set(str(row.get("A_l", "")))
            self.scan_Bu_var.set(str(row.get("B_u", "")))
            self.scan_Bl_var.set(str(row.get("B_l", "")))
            
            # Now extract I from Elements CSV (as before)
            I_scan = self._get_I_from_elements(elem, mass)
            self.scan_I_var.set(str(I_scan))
            self.I_scan = I_scan
            
            # Store floats for calculation:
            self.scan_data = {
                "A_l": float(row.get("A_l", 0)),
                "A_u": float(row.get("A_u", 0)),
                "B_l": float(row.get("B_l", 0)),
                "B_u": float(row.get("B_u", 0)),
            }
            self.try_calculate()
        except Exception as e:
            messagebox.showerror("Scan Load Error", f"Could not load scan parameters:\n{e}")

    def try_calculate(self):
        """
        If all reference and scan parameters are loaded, perform nuclear moments calculation.
        """
        if self.ref_data is not None and self.scan_data is not None and self.I_scan is not None:
            self.calculate_moments()

    def calculate_moments(self):
        """
        Calculate nuclear moments μ and Q for the target scan using the reference isotope.
        """
        try:
            # Reference (already parsed to float)
            mu_ref = self.ref_data["mu_ref"]
            Q_ref = self.ref_data["Q_ref"]
            I_ref = self.ref_data["I_ref"]
            A_l_ref = self.ref_data["A_l_ref"]
            A_u_ref = self.ref_data["A_u_ref"]
            B_l_ref = self.ref_data["B_l_ref"]
            B_u_ref = self.ref_data["B_u_ref"]
            # Scan
            A_l = self.scan_data["A_l"]
            A_u = self.scan_data["A_u"]
            B_l = self.scan_data["B_l"]
            B_u = self.scan_data["B_u"]
            I = self.I_scan

            # μ = (A * I * μ_ref) / (A_ref * I_ref)
            mu_l = (A_l * I * mu_ref) / (A_l_ref * I_ref) if (A_l_ref and I_ref) else "n/a"
            mu_u = (A_u * I * mu_ref) / (A_u_ref * I_ref) if (A_u_ref and I_ref) else "n/a"
            # Q = (B * Q_ref) / B_ref
            Q_l = (B_l * Q_ref) / B_l_ref if B_l_ref else "n/a"
            Q_u = (B_u * Q_ref) / B_u_ref if B_u_ref else "n/a"

            self.mu_l_var.set(f"{mu_l:.5g}" if isinstance(mu_l, float) else mu_l)
            self.mu_u_var.set(f"{mu_u:.5g}" if isinstance(mu_u, float) else mu_u)
            self.Q_l_var.set(f"{Q_l:.5g}" if isinstance(Q_l, float) else Q_l)
            self.Q_u_var.set(f"{Q_u:.5g}" if isinstance(Q_u, float) else Q_u)
        except Exception as e:
            self.mu_l_var.set("Err")
            self.mu_u_var.set("Err")
            self.Q_l_var.set("Err")
            self.Q_u_var.set("Err")
            messagebox.showerror("Calculation Error", f"Could not calculate moments:\n{e}")

    def _get_I_from_elements(self, elem, mass):
        """
        Extract nuclear spin I for the isotope from Elements/{Element}.csv for given mass number.
        """
        # Load from ./Elements/{Element}.csv, get I for mass number
        path = os.path.join(".", "Elements", f"{elem.capitalize()}.csv")
        df = pd.read_csv(path)
        match = df[df['Mass'] == int(mass)]
        if match.empty:
            raise ValueError(f"Mass {mass} not found in {path}")
        I_val = match.iloc[0].get("I", "")
        return self._parse_fraction(I_val)

    def _parse_fraction(self, val):
        """
        Convert fractional or decimal string to float for nuclear spin values.
        """
        try:
            if isinstance(val, (int, float)):
                return float(val)
            if "/" in str(val):
                num, denom = val.split("/")
                return float(num) / float(denom)
            return float(val)
        except Exception:
            return None