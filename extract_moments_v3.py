# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 10:43:16 2025

extract_moments_v3.py

Version 3

User gives desired reference isotope and target isotope and scan.
Program finds reference data and scan/target data, then calculates and outputs target nuclear moments.
Now with moment vs isotope plots for Q and mu.

Author: Shay Rowlinson
"""

import os
import pandas as pd
from tkinter import ttk, StringVar, messagebox, Entry
from uncertainties import ufloat
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

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
        
        self.plot_frame = ttk.LabelFrame(self, text="Moments vs Isotope", padding=(8,4))
        self.plot_frame.grid(row=0, column=1, rowspan=4, sticky="ns", padx=(12,4), pady=(10,2))
        
        self.elem_plot_var = StringVar()
        ttk.Label(self.plot_frame, text="Element:").grid(row=0, column=0, sticky='e', pady=3)
        ttk.Entry(self.plot_frame, textvariable=self.elem_plot_var, width=6).grid(row=0, column=1, sticky='w', pady=3)
        ttk.Button(self.plot_frame, text="Plot", command=self.plot_moments_vs_isotope).grid(row=0, column=2, padx=8)
        
        self.moments_fig = Figure(figsize=(4.7,4))
        self.moments_ax = self.moments_fig.subplots()
        self.moments_canvas = FigureCanvasTkAgg(self.moments_fig, master=self.plot_frame)
        self.moments_canvas.get_tk_widget().grid(row=1, column=0, columnspan=3, sticky="nsew")
        

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
    
        # Target Scan  
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
  
        # Uncertainty rows
        ttk.Label(out_frame, text="σ(μ_l):").grid(row=0, column=4, sticky='e')
        self.mu_l_err_var = StringVar()
        Entry(out_frame, textvariable=self.mu_l_err_var, width=12, state="readonly").grid(row=0, column=5)
        
        ttk.Label(out_frame, text="σ(μ_u):").grid(row=0, column=6, sticky='e')
        self.mu_u_err_var = StringVar()
        Entry(out_frame, textvariable=self.mu_u_err_var, width=12, state="readonly").grid(row=0, column=7)
        
        ttk.Label(out_frame, text="σ(Q_l):").grid(row=1, column=4, sticky='e')
        self.Q_l_err_var = StringVar()
        Entry(out_frame, textvariable=self.Q_l_err_var, width=12, state="readonly").grid(row=1, column=5)
        
        ttk.Label(out_frame, text="σ(Q_u):").grid(row=1, column=6, sticky='e')
        self.Q_u_err_var = StringVar()
        Entry(out_frame, textvariable=self.Q_u_err_var, width=12, state="readonly").grid(row=1, column=7)
        
        # Export button
        btn_frame = ttk.Frame(self)
        btn_frame.grid(row=3, column=0, sticky="w", padx=2, pady=(0,10))
        ttk.Button(btn_frame, text="Export moments", command=self.export_moments).grid(row=0, column=0)
        
    def _uf(self, val, err=0.0):
        try:
            v = float(val)
            e = float(err) if err not in (None, "", "nan") else 0.0
            if e < 0 or not (e == e):  # NaN/<0 guard
                e = 0.0
            return ufloat(v, e)
        except Exception:
            return ufloat(0.0, 0.0)


    def load_reference(self):
        """
        Load reference isotope parameters (μ_ref, Q_ref, I_ref, A_ref, B_ref) from a CSV, and fill GUI fields.
        """
        elem = self.ref_elem_var.get().strip()
        mass = self.ref_mass_var.get().strip()
        if not elem or not mass:
            messagebox.showerror("Missing Input", "Enter both element and mass number for reference isotope.")
            return
        path = os.path.join(".", "References", f"{elem}_{mass}.csv")
        try:
            df = pd.read_csv(path)
            row = df.iloc[0]
            required_cols = ["I", "mu", "Q", "A_l", "A_u", "B_l", "B_u"]
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                messagebox.showerror("Reference CSV Error", f"Missing columns in {path}:\n{', '.join(missing)}")
                return
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
                "mu_ref": self._uf(row.get("mu", 0), row.get("mu_err", 0)),
                "Q_ref":  self._uf(row.get("Q", 0),  row.get("Q_err", 0)),
                "A_l_ref": self._uf(row.get("A_l", 0), row.get("A_l_err", 0)),
                "A_u_ref": self._uf(row.get("A_u", 0), row.get("A_u_err", 0)),
                "B_l_ref": self._uf(row.get("B_l", 0), row.get("B_l_err", 0)),
                "B_u_ref": self._uf(row.get("B_u", 0), row.get("B_u_err", 0)),
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
            required_cols = ["scan", "A_l", "A_u", "B_l", "B_u"]
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                messagebox.showerror("Scan CSV Error", f"Missing columns in {path}:\n{', '.join(missing)}")
                return
            match = df[df['scan'].astype(str) == str(scan)]

            if match.empty:
                raise ValueError(f"Scan {scan} not found in {path}")
            self.current_scan_number = str(scan)
            row = match.iloc[0]
            self.scan_Au_var.set(str(row.get("A_u", "")))
            self.scan_Al_var.set(str(row.get("A_l", "")))
            self.scan_Bu_var.set(str(row.get("B_u", "")))
            self.scan_Bl_var.set(str(row.get("B_l", "")))
            
            I_scan = None
            if "I" in row.index or "I" in df.columns:
                try:
                    I_scan_val = row.get("I", "")
                    if I_scan_val not in ("", None):
                        I_scan = self._parse_fraction(I_scan_val)
                except Exception:
                    I_scan = None
            
            # If not found in Saved_Parameters, fall back to Elements CSV
            if I_scan is None:
                I_scan = self._get_I_from_elements(elem, mass)
            
            self.scan_I_var.set(str(I_scan) if I_scan is not None else "")
            self.I_scan = I_scan
            
            # Store for calculation:
            self.scan_data = {
                "A_l": self._uf(row.get("A_l", 0), row.get("A_l_err", 0)),
                "A_u": self._uf(row.get("A_u", 0), row.get("A_u_err", 0)),
                "B_l": self._uf(row.get("B_l", 0), row.get("B_l_err", 0)),
                "B_u": self._uf(row.get("B_u", 0), row.get("B_u_err", 0)),
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
            mu_ref  = self.ref_data["mu_ref"]
            Q_ref   = self.ref_data["Q_ref"]
            I_ref   = float(self.ref_data["I_ref"]) if self.ref_data["I_ref"] is not None else None
            A_l_ref = self.ref_data["A_l_ref"]
            A_u_ref = self.ref_data["A_u_ref"]
            B_l_ref = self.ref_data["B_l_ref"]
            B_u_ref = self.ref_data["B_u_ref"]
        
            A_l = self.scan_data["A_l"]
            A_u = self.scan_data["A_u"]
            B_l = self.scan_data["B_l"]
            B_u = self.scan_data["B_u"]
            I   = float(self.I_scan) if self.I_scan is not None else None
        
            # Guards
            if None in (I, I_ref) or A_l_ref.n == 0 or A_u_ref.n == 0 or B_l_ref.n == 0 or B_u_ref.n == 0:
                raise ValueError("Missing spin or zero reference constants prevent moment calculation.")
        
            mu_l = (A_l * I * mu_ref) / (A_l_ref * I_ref)
            mu_u = (A_u * I * mu_ref) / (A_u_ref * I_ref)
            Q_l  = (B_l * Q_ref) / B_l_ref
            Q_u  = (B_u * Q_ref) / B_u_ref
        
            self.mu_l_var.set(f"{mu_l.n:.6g}"); self.mu_l_err_var.set(f"{mu_l.s:.3g}")
            self.mu_u_var.set(f"{mu_u.n:.6g}"); self.mu_u_err_var.set(f"{mu_u.s:.3g}")
            self.Q_l_var.set(f"{Q_l.n:.6g}");   self.Q_l_err_var.set(f"{Q_l.s:.3g}")
            self.Q_u_var.set(f"{Q_u.n:.6g}");   self.Q_u_err_var.set(f"{Q_u.s:.3g}")
        
            # Stash for exporter
            self._last_moments = dict(mu_l=mu_l, mu_u=mu_u, Q_l=Q_l, Q_u=Q_u)
    
        except Exception as e:
            for v in (self.mu_l_var, self.mu_u_var, self.Q_l_var, self.Q_u_var,
                      self.mu_l_err_var, self.mu_u_err_var, self.Q_l_err_var, self.Q_u_err_var):
                v.set("Err")
            self._last_moments = None
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
        
    def export_moments(self):
        """
        Append calculated moments with 1 stdev uncertainties to
        ./<Element>_Results/<Mass>/Saved_moments.csv
        Create CSV if not found
        """
        try:
            elem = self.elem_var.get().strip()
            mass = self.mass_var.get().strip()
            if not elem or not mass:
                messagebox.showerror("Missing fields", "Element and Mass must be set before exporting.")
                return
            if not hasattr(self, "current_scan_number"):
                messagebox.showerror("Missing scan", "Load a scan first so we know the scan number.")
                return
            if self._last_moments is None:
                messagebox.showerror("No moments", "Calculate moments before exporting.")
                return
    
            root_dir = os.path.dirname(os.path.abspath(__file__))
            save_dir = os.path.join(root_dir, f"{elem.capitalize()}_Results", str(mass))
            os.makedirs(save_dir, exist_ok=True)
            save_path = os.path.join(save_dir, "Saved_moments.csv")
    
            # Prepare row
            row = {
                "scan": self.current_scan_number,
                "I": str(self.I_scan) if self.I_scan is not None else "",
                "mu_l": f"{self._last_moments['mu_l'].n:.12g}",
                "mu_l_err": f"{self._last_moments['mu_l'].s:.12g}",
                "mu_u": f"{self._last_moments['mu_u'].n:.12g}",
                "mu_u_err": f"{self._last_moments['mu_u'].s:.12g}",
                "Q_l": f"{self._last_moments['Q_l'].n:.12g}",
                "Q_l_err": f"{self._last_moments['Q_l'].s:.12g}",
                "Q_u": f"{self._last_moments['Q_u'].n:.12g}",
                "Q_u_err": f"{self._last_moments['Q_u'].s:.12g}",
            }
            cols = ["scan","I","mu_l","mu_l_err","mu_u","mu_u_err","Q_l","Q_l_err","Q_u","Q_u_err"]
    
            import pandas as pd
            if os.path.exists(save_path):
                df = pd.read_csv(save_path, dtype=str)
                # Replace row if scan already present
                df = df[df["scan"] != str(self.current_scan_number)]
                df = pd.concat([df, pd.DataFrame([row], columns=cols)], ignore_index=True)
                df.to_csv(save_path, index=False)
            else:
                pd.DataFrame([row], columns=cols).to_csv(save_path, index=False)
    
            messagebox.showinfo("Exported", f"Moments saved to:\n{save_path}")
        except Exception as e:
            messagebox.showerror("Export error", f"Could not export moments:\n{e}")
    
    def plot_moments_vs_isotope(self):
        def weighted_mean_and_err(vals, errs):
            vals = np.array(vals, dtype=float)
            errs = np.array(errs, dtype=float)
            mask = (errs > 0) & np.isfinite(vals) & np.isfinite(errs)
            if not np.any(mask):
                return np.nan, np.nan
            weights = 1.0 / (errs[mask] ** 2)
            mean = np.sum(vals[mask] * weights) / np.sum(weights)
            err = np.sqrt(1.0 / np.sum(weights))
            return mean, err
        element = self.elem_plot_var.get().strip()
        if not element:
            messagebox.showerror("Element?", "Enter element symbol")
            return
    
        # Find ./<Element>_Results/<mass>/Saved_moments.csv
        root_dir = os.path.dirname(os.path.abspath(__file__))
        elem_dir = os.path.join(root_dir, f"{element.capitalize()}_Results")
        if not os.path.exists(elem_dir):
            messagebox.showerror("Not found", f"No results directory for {element}")
            return
    
        mass_list = []
        mu_l_means, mu_l_errs = [], []
        mu_u_means, mu_u_errs = [], []
        Q_l_means, Q_l_errs = [], []
        Q_u_means, Q_u_errs = [], []
    
        for sub in sorted(os.listdir(elem_dir), key=lambda x: int(x) if x.isdigit() else 999999):
            subdir = os.path.join(elem_dir, sub)
            if not os.path.isdir(subdir):
                continue
            csv_path = os.path.join(subdir, "Saved_moments.csv")
            if not os.path.exists(csv_path):
                continue
            try:
                df = pd.read_csv(csv_path)
                # For each moment type, collect all values & errors, then weighted mean
                for col, errcol, collect, collect_err in [
                    ("mu_l", "mu_l_err", mu_l_means, mu_l_errs),
                    ("mu_u", "mu_u_err", mu_u_means, mu_u_errs),
                    ("Q_l",  "Q_l_err",  Q_l_means,  Q_l_errs),
                    ("Q_u",  "Q_u_err",  Q_u_means,  Q_u_errs)
                ]:
                    vals = pd.to_numeric(df[col], errors='coerce')
                    errs = pd.to_numeric(df[errcol], errors='coerce')
                    mean, err = weighted_mean_and_err(vals, errs)
                    collect.append(mean)
                    collect_err.append(err)
                mass_list.append(int(sub))
            except Exception as e:
                print(f"Error loading {csv_path}: {e}")
    
        if not mass_list:
            messagebox.showinfo("No Data", f"No isotopes with Saved_moments.csv found for {element}.")
            return
    
        # Plotting
        self.moments_fig.clear()
        ax_mu = self.moments_fig.add_subplot(1, 2, 1)
        ax_q = self.moments_fig.add_subplot(1, 2, 2)
        
        # Plot mu moments
        ax_mu.errorbar(mass_list, mu_l_means, yerr=mu_l_errs, fmt='o', ms=6, capsize=3, color='red', label='$\mu_l$')
        ax_mu.errorbar(mass_list, mu_u_means, yerr=mu_u_errs, fmt='o', ms=6, capsize=3, color='orange', label='$\mu_u$')
        ax_mu.set_xlabel("Isotope mass number")
        ax_mu.set_ylabel("Magnetic moment / $\mu_N$")
        ax_mu.set_title(f"{element} μ vs isotope")
        ax_mu.grid(True, linestyle=':', alpha=0.32)
        ax_mu.legend()
        
        # Plot Q moments
        ax_q.errorbar(mass_list, Q_l_means, yerr=Q_l_errs, fmt='s', ms=6, capsize=3, color='blue', label='$Q_l$')
        ax_q.errorbar(mass_list, Q_u_means, yerr=Q_u_errs, fmt='s', ms=6, capsize=3, color='purple', label='$Q_u$')
        ax_q.set_xlabel("Isotope mass number")
        ax_q.set_ylabel("Quadrupole moment / b")
        ax_q.set_title(f"{element} Q vs isotope")
        ax_q.grid(True, linestyle=':', alpha=0.32)
        ax_q.legend()
        
        self.moments_fig.tight_layout()
        self.moments_canvas.draw()
        
