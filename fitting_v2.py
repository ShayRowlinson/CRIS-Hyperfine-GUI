# -*- coding: utf-8 -*-
"""
Created on Mon Jul 7 12:21:57 2025

fitting_v2.py

Version 2

Performs hyperfine structure fitting capabilities of the interactive GUI.
User imports scan data in format exported from binning tab, and inputs starting parameters.
Program exports saved results to /{Elemental Symbol}_Results/{Mass Number}
Now with ismoeric fitting capability. Exports parameters of the currently selected fit.
Author: Shay Rowlinson
"""

import os
from tkinter import ttk, filedialog, messagebox, DoubleVar, BooleanVar, StringVar
import tkinter as tk
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import satlas2
from matplotlib.figure import Figure

class Fitting(ttk.Frame):
    """
    Fitting tab for hyperfine structure scan fitting in CRIS data, with isomeric (two-fit) support.
    User can import scan, set parameters for Fit 1 and Fit 2, and view/fit both overlays on same plot.
    """
    def __init__(self, parent):
        """
        Initialise fitting frame and widgets.
        """
        super().__init__(parent)
        self.filepath = None
        self.data = None

        # Fit tracking: each fit stores its state
        self.selected_fit = 1  # Which fit is selected for editing 
        self.fits = {
            1: {
                "active": False,
                "estimate_visible": False,
                "fit_curve": None,
                "estimate_curve": None,
                "fit_x": None,
                "last_fit_hfs": None,
                "fit_result": None
            },
            2: {
                "active": False,
                "estimate_visible": False,
                "fit_curve": None,
                "estimate_curve": None,
                "fit_x": None,
                "last_fit_hfs": None,
                "fit_result": None
            }
        }

        # Variables for parameters (per fit)
        self.I_vars = {1: DoubleVar(), 2: DoubleVar()}
        self.A_l_vars = {1: DoubleVar(value=0), 2: DoubleVar(value=0)}
        self.A_u_vars = {1: DoubleVar(value=0), 2: DoubleVar(value=0)}
        self.B_l_vars = {1: DoubleVar(value=0), 2: DoubleVar(value=0)}
        self.B_u_vars = {1: DoubleVar(value=0), 2: DoubleVar(value=0)}
        self.A_l_fix = {1: BooleanVar(), 2: BooleanVar()}
        self.A_u_fix = {1: BooleanVar(), 2: BooleanVar()}
        self.B_l_fix = {1: BooleanVar(), 2: BooleanVar()}
        self.B_u_fix = {1: BooleanVar(), 2: BooleanVar()}
        self.Bu_Bl_fix = {1: BooleanVar(), 2: BooleanVar()}
        self.df_vars = {1: DoubleVar(value=0), 2: DoubleVar(value=0)}
        self.scale_vars = {1: DoubleVar(value=10), 2: DoubleVar(value=10)}
        self.FWHMg_vars = {1: DoubleVar(value=100), 2: DoubleVar(value=100)}
        self.FWHMl_vars = {1: DoubleVar(value=100), 2: DoubleVar(value=100)}
        self.intensity_vars = {1: DoubleVar(value=1.0), 2: DoubleVar(value=1.0)}
        self.racah_int = {1: BooleanVar(), 2: BooleanVar()}
        self.background_vars = {1: DoubleVar(value=0.0005), 2: DoubleVar(value=0.0005)}
        self.A_ratio = {1: None, 2: None}
        self.B_ratio = {1: None, 2: None}
        self._updating_A = False
        self._updating_B = False

        # Shared variables
        self.element_var = StringVar()
        self.massnumber_var = StringVar()
        self.Jl_var = DoubleVar()
        self.Ju_var = DoubleVar()
        self.fit1_active = BooleanVar(value=False)
        self.fit2_active = BooleanVar(value=False)

        self._build_gui()
        self._add_param_traces()

    def _build_gui(self):
        """
        Build the full GUI: fit selection, plot, parameter frame, buttons, etc.
        """
        # --- Top fit selection controls ---
        control_frame = ttk.Frame(self)
        control_frame.grid(row=0, column=0, sticky="ew", padx=2, pady=(10, 0), columnspan=2)

        # Fit radiobuttons
        self.fit_select_var = tk.IntVar(value=1)
        ttk.Radiobutton(control_frame, text="Fit 1", variable=self.fit_select_var, value=1, command=lambda: self._on_fit_toggle(1)).grid(row=0, column=0, padx=(0,6))
        ttk.Radiobutton(control_frame, text="Fit 2", variable=self.fit_select_var, value=2, command=lambda: self._on_fit_toggle(2)).grid(row=0, column=1, padx=(0,10))
        # Active checkboxes
        ttk.Checkbutton(control_frame, text="Active", variable=self.fit1_active, command=lambda: self._on_active_toggle(1)).grid(row=0, column=2, padx=(0,12))
        ttk.Checkbutton(control_frame, text="Active", variable=self.fit2_active, command=lambda: self._on_active_toggle(2)).grid(row=0, column=3, padx=(0,12))
        # Show Estimate
        ttk.Button(control_frame, text="Show Estimate (Fit 1)", command=lambda: self._on_show_estimate(1)).grid(row=0, column=4, padx=(6,0))
        ttk.Button(control_frame, text="Show Estimate (Fit 2)", command=lambda: self._on_show_estimate(2)).grid(row=0, column=5, padx=(6,0))

        # --- Parameter inputs ---
        param_frame = ttk.LabelFrame(self, text="Spin and isotope parameters", padding=(10, 6))
        param_frame.grid(row=1, column=0, sticky="ew", padx=2)
        ttk.Label(param_frame, text="Element symbol:").grid(row=0, column=0, sticky='e', pady=3)
        ttk.Entry(param_frame, textvariable=self.element_var, width=6).grid(row=0, column=1, sticky='w', pady=3)
        ttk.Label(param_frame, text="Mass number:").grid(row=0, column=2, sticky='e', pady=3)
        ttk.Entry(param_frame, textvariable=self.massnumber_var, width=7).grid(row=0, column=3, sticky='w', pady=3)
        ttk.Label(param_frame, text="J_l").grid(row=1, column=0, sticky="e")
        ttk.Entry(param_frame, textvariable=self.Jl_var, width=5).grid(row=1, column=1, sticky="w")
        ttk.Label(param_frame, text="J_u").grid(row=1, column=2, sticky="e")
        ttk.Entry(param_frame, textvariable=self.Ju_var, width=5).grid(row=1, column=3, sticky="w")
        ttk.Label(param_frame, text="I (selected fit)").grid(row=2, column=0, sticky="e")
        self.I_entry = ttk.Entry(param_frame, textvariable=self.I_vars[self.selected_fit], width=5)
        self.I_entry.grid(row=2, column=1, sticky="w")

        # --- Main plost ---
        plot_frame = ttk.LabelFrame(self, text="Scan Fit", padding=(6,4))
        plot_frame.grid(row=2, column=0, sticky="nsew", pady=8)
        self.fig = Figure(figsize=(6,2.5))
        self.ax = self.fig.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")

        # --- Fit parameter control ---
        self.fitparam_frame = ttk.LabelFrame(self, text="Fit parameters (selected fit)", padding=(10, 6))
        self.fitparam_frame.grid(row=1, column=1, rowspan=2, sticky="ns", padx=(8,2), pady=(0,8))
        self._render_param_widgets()

        # --- Fitted parameter display area ---
        self.fit_results_label = ttk.Label(self.fitparam_frame, text="Fit results (selected fit):", font=('Segoe UI', 10, 'bold'))
        self.fit_results_label.grid(row=12, column=0, columnspan=4, sticky="w", pady=(10,0))
        self.fit_results_box = tk.Text(self.fitparam_frame, height=10, width=44, font=("Consolas", 9), state='disabled', background="#f7f7ff")
        self.fit_results_box.grid(row=13, column=0, columnspan=5, sticky="ew", pady=(0,6))

        # --- Buttons ---
        button_frame = ttk.Frame(self)
        button_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(8,2))
        ttk.Button(button_frame, text="Import scan(s)", command=self._import_scan).grid(row=0, column=0, padx=5)
        ttk.Button(button_frame, text="Fit Peaks", command=self._fit_peaks).grid(row=0, column=1, padx=5)
        ttk.Button(button_frame, text="Import parameters", command=self._import_parameters).grid(row=0, column=2, padx=10)
        ttk.Button(button_frame, text="Save fit parameters", command=self._save_fit_parameters).grid(row=0, column=3, padx=5)
        ttk.Button(button_frame, text="Expand Plot", command=self._expand_plot).grid(row=0, column=4, padx=5)

    def _render_param_widgets(self):
        """
        Clears and re-renders the parameter controls for the currently selected fit.
        Called on startup and every time the selected fit changes.
        """
        # Remove old widgets
        for child in self.fitparam_frame.winfo_children():
            info = child.grid_info()
            if info.get("row") not in (12, 13):
                child.destroy()
        fit = self.selected_fit  # 1 or 2
        self._add_param_slider(self.fitparam_frame, 0, "A_l", self.A_l_vars[fit], self.A_l_fix[fit])
        self._add_param_slider(self.fitparam_frame, 1, "A_u", self.A_u_vars[fit], self.A_u_fix[fit])
        ttk.Checkbutton(self.fitparam_frame, text="A_u/A_l Fix", variable=self.A_u_fix[fit]).grid(row=1, column=4)
        self._add_param_slider(self.fitparam_frame, 2, "B_l", self.B_l_vars[fit], self.B_l_fix[fit])
        self._add_param_slider(self.fitparam_frame, 3, "B_u", self.B_u_vars[fit], self.B_u_fix[fit])
        ttk.Checkbutton(self.fitparam_frame, text="B_u/B_l Fix", variable=self.Bu_Bl_fix[fit]).grid(row=3, column=4)
        self._add_param_slider(self.fitparam_frame, 4, "df (MHz)", self.df_vars[fit])
        self._add_param_slider(self.fitparam_frame, 5, "scale", self.scale_vars[fit], from_=0, to=20)
        self._add_param_slider(self.fitparam_frame, 6, "FWHM_g", self.FWHMg_vars[fit], from_=0, to=5000)
        self._add_param_slider(self.fitparam_frame, 7, "FWHM_l", self.FWHMl_vars[fit], from_=0, to=5000)
        self._add_param_slider(self.fitparam_frame, 8, "Intensity", self.intensity_vars[fit], from_=0, to=5)
        self._add_param_slider(self.fitparam_frame, 9, "Background", self.background_vars[fit], from_=0, to=0.00001)
        ttk.Checkbutton(self.fitparam_frame, text="Racah int.", variable=self.racah_int[fit]).grid(row=10, column=3)

    def _add_param_slider(self, parent, row, label, var, fix_var=None, from_=-2000, to=2000):
        """Helper to add a label, slider, entry and (optionally) fix checkbox for a fit parameter."""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="e")
        slider = ttk.Scale(parent, from_=from_, to=to, variable=var, orient="horizontal", length=80)
        slider.grid(row=row, column=1, sticky="ew")
        entry = ttk.Entry(parent, textvariable=var, width=7)
        entry.grid(row=row, column=2, sticky="w")
        if fix_var is not None:
            ttk.Checkbutton(parent, text="Fix", variable=fix_var).grid(row=row, column=3, padx=2)

    def _add_param_traces(self):
        """
        Sets up traces for parameter sliders and fix ratio logic for both fits.
        """
        for fit in (1,2):
            self.A_l_vars[fit].trace_add("write", lambda *args, fit=fit: self._on_A_var_changed("A_l", fit))
            self.A_u_vars[fit].trace_add("write", lambda *args, fit=fit: self._on_A_var_changed("A_u", fit))
            self.A_u_fix[fit].trace_add("write", lambda *args, fit=fit: self._on_A_ratio_fix_toggled(fit))
            self.B_l_vars[fit].trace_add("write", lambda *args, fit=fit: self._on_B_var_changed("B_l", fit))
            self.B_u_vars[fit].trace_add("write", lambda *args, fit=fit: self._on_B_var_changed("B_u", fit))
            self.Bu_Bl_fix[fit].trace_add("write", lambda *args, fit=fit: self._on_B_ratio_fix_toggled(fit))
            for var in [
                self.df_vars[fit], self.scale_vars[fit], self.FWHMg_vars[fit], self.FWHMl_vars[fit],
                self.intensity_vars[fit], self.background_vars[fit]
            ]:
                var.trace_add("write", lambda *args, fit=fit: self._slider_update_estimate(fit))

    def _on_fit_toggle(self, fit_num):
        """Called when Fit 1 or Fit 2 is selected for parameter editing."""
        self.selected_fit = fit_num
        self.I_entry.config(textvariable=self.I_vars[fit_num])
        self._render_param_widgets()
        self._update_fit_results_box()

    def _on_active_toggle(self, fit_num):
        """Toggle the 'active' state for a fit (enables/disables for fitting)."""
        self.fits[fit_num]["active"] = not self.fits[fit_num]["active"]

    def _on_show_estimate(self, fit_num):
        """Show/hide estimate for given fit; overlays estimates from both fits."""
        self.fits[fit_num]["estimate_visible"] = not self.fits[fit_num]["estimate_visible"]
        self._update_estimate(fit_num)

    def _update_estimate(self, fit_num):
        """Updates and displays SATLAS2 estimate for current parameter values of fit_num, overlays on main plot."""
        if self.data is None or not self.fits[fit_num]["estimate_visible"]:
            self.fits[fit_num]["estimate_curve"] = None
            self._plot_data()
            return
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
        params = self._collect_fit_params(fit_num)
        hfs = satlas2.HFS(
            I=params["I"],
            J=[params["J_l"], params["J_u"]],
            A=[params["A_l"], params["A_u"]],
            B=[params["B_l"], params["B_u"]],
            C=[0,0],
            df=params["df"],
            scale=params["scale"],
            racah=params["racah"],
            fwhmg=params["FWHMg"],
            fwhml=params["FWHMl"],
            name='hfs',
        )
        background = params["background"]
        bkg = satlas2.Polynomial([background], name='bkg')
        source = satlas2.Source(x=x, y=y, yerr=yerr, name="source")
        source.addModel(hfs)
        source.addModel(bkg)
        fit_x = np.linspace(x.min(), x.max(), 500)
        estimate_curve = source.evaluate(fit_x)
        self.fits[fit_num]["estimate_curve"] = estimate_curve
        self.fits[fit_num]["fit_x"] = fit_x
        self._plot_data()

    def _slider_update_estimate(self, fit_num):
        """Only update estimate if estimate is visible for given fit."""
        if self.fits[fit_num]["estimate_visible"]:
            self._update_estimate(fit_num)

    def _import_scan(self):
        """Opens a dialogue to import a binned scan CSV and plots the data."""
        path = filedialog.askopenfilename(title="Select binned scan CSV", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        self.filepath = path
        self.data = pd.read_csv(self.filepath, sep=";")
        for fit_num in (1,2):
            self.fits[fit_num]["fit_curve"] = None
            self.fits[fit_num]["estimate_curve"] = None
            self.fits[fit_num]["fit_x"] = None
        self._plot_data()
        # Show all estimates that are toggled on
        for fit_num in (1,2):
            if self.fits[fit_num]["estimate_visible"]:
                self._update_estimate(fit_num)

    def _plot_data(self):
        """Plots the imported scan data and overlays the fits/estimates for both fits."""
        self.fig.clf()
        self.ax = self.fig.subplots()
        if self.data is None:
            self.ax.set_title("No scan loaded.")
            self.canvas.draw()
            return
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
        self.ax.errorbar(x, y, yerr=yerr, fmt='o', color='red', markersize=2, ecolor='k', capsize=2, label='Data')
        # Overlay estimate/fit for both fits if present
        colors = {1: 'green', 2: 'purple'}
        fit_labels = {1: 'Fit 1', 2: 'Fit 2'}
        estimate_labels = {1: 'Estimate 1', 2: 'Estimate 2'}
        for fit_num in (1,2):
            if self.fits[fit_num]["estimate_visible"] and self.fits[fit_num]["estimate_curve"] is not None and self.fits[fit_num]["fit_x"] is not None:
                self.ax.plot(self.fits[fit_num]["fit_x"], self.fits[fit_num]["estimate_curve"], color=colors[fit_num], linestyle='--', label=estimate_labels[fit_num])
            if self.fits[fit_num]["fit_curve"] is not None and self.fits[fit_num]["fit_x"] is not None:
                self.ax.plot(self.fits[fit_num]["fit_x"], self.fits[fit_num]["fit_curve"], color=colors[fit_num], label=fit_labels[fit_num])
        self.ax.set_ylabel('Counts per bunch')
        self.ax.legend()
        self.ax.set_title(os.path.basename(self.filepath) if self.filepath else "Scan")
        self.ax.set_xlabel("Frequency offset / MHz")
        self.ax.grid(True)
        self.fig.tight_layout()
        self.canvas.draw()

    def _fit_peaks(self):
        """
        Runs SATLAS2 fitting for all fits currently set as active.
        Overlays both fits on the same plot if both are active.
        """
        if self.data is None:
            messagebox.showerror("No scan loaded", "Please import a scan CSV first.")
            return
        fits_to_run = []
        if self.fit1_active.get():
            fits_to_run.append(1)
        if self.fit2_active.get():
            fits_to_run.append(2)
        if not fits_to_run:
            messagebox.showinfo("No fit active", "Please activate at least one fit.")
            return
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
        for fit_num in fits_to_run:
            params = self._collect_fit_params(fit_num)
            hfs = satlas2.HFS(
                I=params["I"],
                J=[params["J_l"], params["J_u"]],
                A=[params["A_l"], params["A_u"]],
                B=[params["B_l"], params["B_u"]],
                C=[0,0],
                df=params["df"],
                scale=params["scale"],
                racah=params["racah"],
                fwhmg=params["FWHMg"],
                fwhml=params["FWHMl"],
                name='hfs',
            )
            hfs.params['scale'].vary = True
            hfs.params['Al'].vary = not self.A_l_fix[fit_num].get()
            hfs.params['Au'].vary = not self.A_u_fix[fit_num].get()
            hfs.params['Bl'].vary = not self.B_l_fix[fit_num].get()
            hfs.params['Bu'].vary = not self.B_u_fix[fit_num].get()
            hfs.params['Cl'].vary = False
            hfs.params['Cu'].vary = False
            background = params["background"]
            bkg = satlas2.Polynomial([background], name='bkg')
            source = satlas2.Source(x=x, y=y, yerr=yerr, name="source")
            source.addModel(hfs)
            source.addModel(bkg)
            f = satlas2.Fitter()
            f.addSource(source)
            f.fit()
            report = f.reportFit()
            self.fits[fit_num]["last_fit_hfs"] = hfs  # Store fitted HFS object for later use
            if "uncertainties could not be estimated" in report:
                messagebox.showwarning(
                    "Uncertainties Not Estimated",
                    f"Warning: uncertainties could not be estimated for fit {fit_num}.\nCheck parameter values and try again."
                )
            print(report)
            fit_x = np.linspace(x.min(), x.max(), 500)
            fit_y = source.evaluate(fit_x)
            self.fits[fit_num]["fit_x"] = fit_x
            self.fits[fit_num]["fit_curve"] = fit_y
        for fit_num in fits_to_run:
            self.fits[fit_num]["estimate_visible"] = False
            self.fits[fit_num]["estimate_curve"] = None
        self._plot_data()
        self._update_fit_results_box()
        messagebox.showinfo("Fit finished", "Fit(s) finished! Select 'Expand Plot' for residuals.")

    def _collect_fit_params(self, fit_num):
        """Collects all fit parameter values from GUI for fit_num (1 or 2)."""
        return {
            "A_l": self.A_l_vars[fit_num].get(),
            "A_u": self.A_u_vars[fit_num].get(),
            "B_l": self.B_l_vars[fit_num].get(),
            "B_u": self.B_u_vars[fit_num].get(),
            "df": self.df_vars[fit_num].get(),
            "scale": self.scale_vars[fit_num].get(),
            "FWHMg": self.FWHMg_vars[fit_num].get(),
            "FWHMl": self.FWHMl_vars[fit_num].get(),
            "I": self.I_vars[fit_num].get(),
            "J_l": self.Jl_var.get(),
            "J_u": self.Ju_var.get(),
            "racah": self.racah_int[fit_num].get(),
            "intensity": self.intensity_vars[fit_num].get(),
            "background": self.background_vars[fit_num].get(),
        }

    def _import_parameters(self):
        """Extracts A_l, A_u, B_l, B_u, and I from element CSV using given element symbol and mass number,
        and fills in the fit parameter fields for the selected fit."""
        element = self.element_var.get().strip()
        mass_number = self.massnumber_var.get().strip()
        if not element or not mass_number:
            messagebox.showerror("Error", "Please provide both element symbol and mass number.")
            return
        elements_folder = os.path.join(os.path.dirname(__file__), "Elements")
        csv_path = os.path.join(elements_folder, f"{element.capitalize()}.csv")
        if not os.path.exists(csv_path):
            messagebox.showerror("Error", f"No data file found for {element} at {csv_path}")
            return
        try:
            df = pd.read_csv(csv_path)
            match = df[df['Mass'] == int(mass_number)]
            if match.empty:
                messagebox.showerror("Error", f"No isotope with mass {mass_number} found for {element}")
                return
            fit = self.selected_fit
            self.A_l_vars[fit].set(float(match['A_l'].values[0]))
            self.A_u_vars[fit].set(float(match['A_u'].values[0]))
            self.B_l_vars[fit].set(float(match['B_l'].values[0]))
            self.B_u_vars[fit].set(float(match['B_u'].values[0]))
            if 'I' in match:
                self.I_vars[fit].set(float(match['I'].values[0]))
            if 'centroid' in match:
                self.df_vars[fit].set(float(match['centroid'].values[0]))
            messagebox.showinfo("Parameters Imported",
                                f"Extracted I, A_l, A_u, B_l, B_u for {element}-{mass_number} (fit {fit}).")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to import parameters: {e}")

    def _save_fit_parameters(self):
        """Saves fit parameters and uncertainties for the currently selected fit as a new row in
        (current directory)/<Element>_Results/<mass>/Saved_Parameters.csv"""
        import csv
        fit = self.selected_fit
        if not self.fits[fit]["last_fit_hfs"]:
            messagebox.showerror("No fit", "No fit has been performed yet for this fit. Fit data before saving parameters.")
            return
        element = self.element_var.get().strip()
        mass = self.massnumber_var.get().strip()
        if not element or not mass:
            messagebox.showerror("Error", "Element and mass number required in parameter fields.")
            return
        scan_number = tk.simpledialog.askstring("Scan Number", "Enter scan number for this fit (e.g. 1234):")
        if not scan_number:
            messagebox.showerror("Cancelled", "Scan number is required to save parameters.")
            return
        try:
            int(scan_number)
        except ValueError:
            messagebox.showerror("Invalid input", "Scan number must be an integer.")
            return
        root_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(root_dir, f"{element.capitalize()}_Results", str(mass))
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, "Saved_Parameters.csv")
        fit_params = self.fits[fit]["last_fit_hfs"].params
        centroid_param = fit_params["centroid"] if "centroid" in fit_params else fit_params["df"]
        param_dict = {
            "scan": scan_number,
            "A_l": fit_params["Al"].value,
            "A_l_err": fit_params["Al"].unc if getattr(fit_params["Al"], 'unc', None) is not None else "",
            "A_u": fit_params["Au"].value,
            "A_u_err": fit_params["Au"].unc if getattr(fit_params["Au"], 'unc', None) is not None else "",
            "B_l": fit_params["Bl"].value,
            "B_l_err": fit_params["Bl"].unc if getattr(fit_params["Bl"], 'unc', None) is not None else "",
            "B_u": fit_params["Bu"].value,
            "B_u_err": fit_params["Bu"].unc if getattr(fit_params["Bu"], 'unc', None) is not None else "",
            "centroid": centroid_param.value,
            "centroid_err": centroid_param.unc if getattr(centroid_param, 'unc', None) is not None else "",
        }
        columns = [
            "scan",
            "A_l", "A_l_err",
            "A_u", "A_u_err",
            "B_l", "B_l_err",
            "B_u", "B_u_err",
            "centroid", "centroid_err"
        ]
        new_row = [param_dict[c] for c in columns]
        file_exists = os.path.exists(save_path)
        if file_exists:
            df = pd.read_csv(save_path, dtype=str)
            if set(columns).issubset(df.columns):
                existing = df[df['scan'] == str(scan_number)]
                if not existing.empty:
                    if not messagebox.askyesno("Overwrite?", f"Scan {scan_number} already saved. Overwrite?"):
                        messagebox.showinfo("Aborted", "Save cancelled. Existing scan not overwritten.")
                        return
                    df = df[df['scan'] != str(scan_number)]
                df = pd.concat([df, pd.DataFrame([dict(zip(columns, new_row))])], ignore_index=True)
                df.to_csv(save_path, index=False)
                messagebox.showinfo("Success", f"Fit parameters (with uncertainties) saved in:\n{save_path}")
                return
        with open(save_path, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if not file_exists:
                writer.writerow(columns)
            writer.writerow(new_row)
        messagebox.showinfo("Success", f"Fit parameters (with uncertainties) saved to:\n{save_path}")

    def _expand_plot(self):
        """Opens a new top-level window with a larger version of the current plot, including all fits/estimates."""
        if self.data is None:
            messagebox.showerror("No scan loaded", "Import a scan first.")
            return
        top = tk.Toplevel(self)
        top.title("Expanded Plot")
        fig = Figure(figsize=(8,4.25))
        any_fit = any(self.fits[fit_num]["fit_curve"] is not None for fit_num in (1,2))
        if any_fit:
            axs = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
            ax_main, ax_resid = axs
        else:
            ax_main = fig.subplots()
            ax_resid = None
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
        ax_main.errorbar(x, y, yerr=yerr, fmt='o', color='red', markersize=2, ecolor='k', capsize=2, label='Data')
        colors = {1: 'green', 2: 'purple'}
        fit_labels = {1: 'Fit 1', 2: 'Fit 2'}
        estimate_labels = {1: 'Estimate 1', 2: 'Estimate 2'}
        for fit_num in (1,2):
            if self.fits[fit_num]["estimate_visible"] and self.fits[fit_num]["estimate_curve"] is not None and self.fits[fit_num]["fit_x"] is not None:
                ax_main.plot(self.fits[fit_num]["fit_x"], self.fits[fit_num]["estimate_curve"], color=colors[fit_num], linestyle='--', label=estimate_labels[fit_num])
            if self.fits[fit_num]["fit_curve"] is not None and self.fits[fit_num]["fit_x"] is not None:
                ax_main.plot(self.fits[fit_num]["fit_x"], self.fits[fit_num]["fit_curve"], color=colors[fit_num], label=fit_labels[fit_num])
        ax_main.set_ylabel('Counts per bunch')
        ax_main.legend()
        ax_main.set_title(os.path.basename(self.filepath) if self.filepath else "Scan")
        ax_main.grid(True)
        if ax_resid is not None and self.fits[self.selected_fit]["fit_curve"] is not None:
            from scipy.interpolate import interp1d
            fit_interp = interp1d(self.fits[self.selected_fit]["fit_x"], self.fits[self.selected_fit]["fit_curve"], kind='linear', fill_value="extrapolate")
            residuals = y - fit_interp(x)
            ax_resid.axhline(0, color='grey', lw=1, linestyle='--')
            ax_resid.errorbar(x, residuals, yerr=yerr, fmt='o', markersize=1, color='red', ecolor='black', capsize=2)
            ax_resid.set_ylabel("Residuals")
            ax_resid.set_xlabel("Frequency offset / MHz")
            ax_resid.grid(True)
        else:
            ax_main.set_xlabel("Frequency offset / MHz")
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=False)
        toolbar = NavigationToolbar2Tk(canvas, top)
        toolbar.update()
        toolbar.pack(side="top", fill="x")

    def _on_A_ratio_fix_toggled(self, fit_num):
        """Stores ratio for A_u/A_l for fit_num when ratio fixing is enabled/disabled."""
        if self.A_u_fix[fit_num].get():
            Al = self.A_l_vars[fit_num].get()
            Au = self.A_u_vars[fit_num].get()
            self.A_ratio[fit_num] = Au / Al if Al != 0 else None
        else:
            self.A_ratio[fit_num] = None

    def _on_B_ratio_fix_toggled(self, fit_num):
        """Stores ratio for B_u/B_l for fit_num when ratio fixing is enabled/disabled."""
        if self.Bu_Bl_fix[fit_num].get():
            Bl = self.B_l_vars[fit_num].get()
            Bu = self.B_u_vars[fit_num].get()
            self.B_ratio[fit_num] = Bu / Bl if Bl != 0 else None
        else:
            self.B_ratio[fit_num] = None

    def _on_A_var_changed(self, changed, fit_num):
        """fix for A_l/A_u for fit_num when enabled."""
        if self.A_u_fix[fit_num].get() and self.A_ratio[fit_num] is not None and not self._updating_A:
            try:
                self._updating_A = True
                Al = self.A_l_vars[fit_num].get()
                Au = self.A_u_vars[fit_num].get()
                if changed == "A_l":
                    self.A_u_vars[fit_num].set(Al * self.A_ratio[fit_num])
                elif changed == "A_u":
                    if self.A_ratio[fit_num] != 0:
                        self.A_l_vars[fit_num].set(Au / self.A_ratio[fit_num])
            finally:
                self._updating_A = False
        self._slider_update_estimate(fit_num)

    def _on_B_var_changed(self, changed, fit_num):
        """fix for B_l/B_u for fit_num when enabled."""
        if self.Bu_Bl_fix[fit_num].get() and self.B_ratio[fit_num] is not None and not self._updating_B:
            try:
                self._updating_B = True
                Bl = self.B_l_vars[fit_num].get()
                Bu = self.B_u_vars[fit_num].get()
                if changed == "B_l":
                    self.B_u_vars[fit_num].set(Bl * self.B_ratio[fit_num])
                elif changed == "B_u":
                    if self.B_ratio[fit_num] != 0:
                        self.B_l_vars[fit_num].set(Bu / self.B_ratio[fit_num])
            finally:
                self._updating_B = False
        self._slider_update_estimate(fit_num)

    def _update_fit_results_box(self):
        """Updates the fit results box with the latest results for the currently selected fit."""
        fit = self.selected_fit
        try:
            hfs = self.fits[fit]["last_fit_hfs"]
            if hfs is None:
                self.fit_results_box.config(state='normal')
                self.fit_results_box.delete(1.0, tk.END)
                self.fit_results_box.insert(tk.END, "No fit result for selected fit yet.")
                self.fit_results_box.config(state='disabled')
                return
            fit_params = hfs.params
            self.fit_results_box.config(state='normal')
            self.fit_results_box.delete(1.0, tk.END)
            self.fit_results_box.insert(tk.END, "---- Fit Results (fit {}) ----\n".format(fit))
            for k in fit_params:
                v = fit_params[k].value
                self.fit_results_box.insert(tk.END, f"{k:12s}: {v:.5g}\n")
            self.fit_results_box.config(state='disabled')
        except Exception as e:
            self.fit_results_box.config(state='normal')
            self.fit_results_box.delete(1.0, tk.END)
            self.fit_results_box.insert(tk.END, f"Could not display results: {e}")
            self.fit_results_box.config(state='disabled')