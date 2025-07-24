# -*- coding: utf-8 -*-
"""
Created on Mon Jul 7 12:21:57 2025

fitting_v1.py

Version 1

Performs hyperfine structure fitting capabilities of the interactive GUI.
User imports scan data in format exported from binning tab, and inputs starting parameters.

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
    Fitting tab for hyperfine structure scan fitting in CRIS data.
    User can import scan, set HFS parameters, view/fit data, and export results.
    """
    def __init__(self, parent):
        """
        Initialise fitting frame and widgets.
        """
        super().__init__(parent)
        self.filepath = None
        self.data = None
        self.fit_result = None
        self.estimate_visible = False  # Estimate toggling
        self.estimate_curve = None
        self.fit_curve = None
        self.fit_x = None
        self.A_ratio = None
        self.B_ratio = None
        self._updating_A = False
        self._updating_B = False
        self._build_gui()
        self._add_param_traces()

    def _build_gui(self):
        """
        Build layout for the Fitting tab of GUI
        """
        # --- Parameter inputs (I,J, element, mass) ---
        param_frame = ttk.LabelFrame(self, text="Spin and isotope parameters", padding=(10, 6))
        param_frame.grid(row=0, column=0, sticky="ew", padx=2, pady=(10, 0))

        # Element Symbol and Mass Number
        ttk.Label(param_frame, text="Element symbol:").grid(row=0, column=0, sticky='e', pady=3)
        self.element_var = StringVar()
        element_entry = ttk.Entry(param_frame, textvariable=self.element_var, width=6)
        element_entry.grid(row=0, column=1, sticky='w', pady=3)
        ttk.Label(param_frame, text="Mass number:").grid(row=0, column=2, sticky='e', pady=3)
        self.massnumber_var = StringVar()
        massnumber_entry = ttk.Entry(param_frame, textvariable=self.massnumber_var, width=7)
        massnumber_entry.grid(row=0, column=3, sticky='w', pady=3)

        # Jl, Ju, I 
        ttk.Label(param_frame, text="J_l").grid(row=1, column=0, sticky="e")
        self.Jl_var = DoubleVar()
        ttk.Entry(param_frame, textvariable=self.Jl_var, width=5).grid(row=1, column=1, sticky="w")
        ttk.Label(param_frame, text="J_u").grid(row=1, column=2, sticky="e")
        self.Ju_var = DoubleVar()
        ttk.Entry(param_frame, textvariable=self.Ju_var, width=5).grid(row=1, column=3, sticky="w")
        ttk.Label(param_frame, text="I").grid(row=2, column=0, sticky="e")
        self.I_var = DoubleVar()
        ttk.Entry(param_frame, textvariable=self.I_var, width=5).grid(row=2, column=1, sticky="w")

        # --- Main plot ---
        plot_frame = ttk.LabelFrame(self, text="Scan Fit", padding=(6,4))
        plot_frame.grid(row=1, column=0, sticky="nsew", pady=8)
        self.fig = Figure(figsize=(6,2.5))
        self.ax = self.fig.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")

        # --- Fit parameter controls ---
        fitparam_frame = ttk.LabelFrame(self, text="Fit parameters", padding=(10, 6))
        fitparam_frame.grid(row=0, column=1, rowspan=2, sticky="ns", padx=(8,2), pady=(0,8))
        self.A_l_var = DoubleVar(value=-0); self.A_l_fix = BooleanVar()
        self.A_u_var = DoubleVar(value=0); self.A_u_fix = BooleanVar()
        self.B_l_var = DoubleVar(value=0); self.B_l_fix = BooleanVar()
        self.B_u_var = DoubleVar(value=-0); self.B_u_fix = BooleanVar()
        self.Bu_Bl_fix = BooleanVar()
        self.df_var = DoubleVar(value=0)
        self.scale_var = DoubleVar(value=10)
        self.FWHMg_var = DoubleVar(value=100)
        self.FWHMl_var = DoubleVar(value=100)
        self.intensity_var = DoubleVar(value=1.0)
        self.racah_int = BooleanVar()
        self.background_var = DoubleVar(value=0.0005)
        self._add_param_slider(fitparam_frame, 0, "A_l", self.A_l_var, self.A_l_fix)
        self._add_param_slider(fitparam_frame, 1, "A_u", self.A_u_var, self.A_u_fix)
        ttk.Checkbutton(fitparam_frame, text="A_u/A_l Fix", variable=self.A_u_fix).grid(row=1, column=4)
        self._add_param_slider(fitparam_frame, 2, "B_l", self.B_l_var, self.B_l_fix)
        self._add_param_slider(fitparam_frame, 3, "B_u", self.B_u_var, self.B_u_fix)
        ttk.Checkbutton(fitparam_frame, text="B_u/B_l Fix", variable=self.Bu_Bl_fix).grid(row=3, column=4)
        self._add_param_slider(fitparam_frame, 4, "df (MHz)", self.df_var)
        self._add_param_slider(fitparam_frame, 5, "scale", self.scale_var, from_=0, to=20)
        self._add_param_slider(fitparam_frame, 6, "FWHM_g", self.FWHMg_var, from_=0, to=5000)
        self._add_param_slider(fitparam_frame, 7, "FWHM_l", self.FWHMl_var, from_=0, to=5000)
        self._add_param_slider(fitparam_frame, 8, "Intensity", self.intensity_var, from_=0, to=5)
        self._add_param_slider(fitparam_frame, 9, "Background", self.background_var,from_=0, to=0.00001)
        ttk.Checkbutton(fitparam_frame, text="Racah int.", variable=self.racah_int).grid(row=10, column=3)

        # --- Fitted parameter display area ---
        self.fit_results_label = ttk.Label(fitparam_frame, text="Fit results:", font=('Segoe UI', 10, 'bold'))
        self.fit_results_label.grid(row=12, column=0, columnspan=4, sticky="w", pady=(10,0))
        self.fit_results_box = tk.Text(fitparam_frame, height=10, width=44, font=("Consolas", 9), state='disabled', background="#f7f7ff")
        self.fit_results_box.grid(row=13, column=0, columnspan=5, sticky="ew", pady=(0,6))

        # --- Buttons ---
        button_frame = ttk.Frame(self)
        button_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8,2))
        ttk.Button(button_frame, text="Import scan(s)", command=self._import_scan).grid(row=0, column=0, padx=5)
        ttk.Button(button_frame, text="Show Estimate", command=self._show_estimate).grid(row=0, column=1, padx=5)
        ttk.Button(button_frame, text="Fit Peaks", command=self._fit_peaks).grid(row=0, column=2, padx=5)
        ttk.Button(button_frame, text="Import parameters", command=self._import_parameters).grid(row=0, column=3, padx=10)
        ttk.Button(button_frame, text="Save fit parameters", command=self._save_fit_parameters).grid(row=0, column=4, padx=5)
        ttk.Button(button_frame, text="Expand Plot", command=self._expand_plot).grid(row=0, column=5, padx=5)

    def _add_param_slider(self, parent, row, label, var, fix_var=None, from_=-2000, to=2000):
        """
        Helper to add a label, slider, entry and fix checkbox for a fit parameter.
        """
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="e")
        slider = ttk.Scale(parent, from_=from_, to=to, variable=var, orient="horizontal", length=80)
        slider.grid(row=row, column=1, sticky="ew")
        entry = ttk.Entry(parent, textvariable=var, width=7)
        entry.grid(row=row, column=2, sticky="w")
        if fix_var is not None:
            ttk.Checkbutton(parent, text="Fix", variable=fix_var).grid(row=row, column=3, padx=2)

    def _add_param_traces(self):
        """
        Sets up traces for parameter sliders and fix ratio logic.
        """
        # A ratios
        self.A_l_var.trace_add("write", lambda *args: self._on_A_var_changed("A_l"))
        self.A_u_var.trace_add("write", lambda *args: self._on_A_var_changed("A_u"))
        self.A_u_fix.trace_add("write", self._on_A_ratio_fix_toggled)
        # B ratios
        self.B_l_var.trace_add("write", lambda *args: self._on_B_var_changed("B_l"))
        self.B_u_var.trace_add("write", lambda *args: self._on_B_var_changed("B_u"))
        self.Bu_Bl_fix.trace_add("write", self._on_B_ratio_fix_toggled)
        # All other parameter traces
        for var in [
            self.df_var, self.scale_var, self.FWHMg_var, self.FWHMl_var,
            self.intensity_var, self.background_var
        ]:
            var.trace_add("write", lambda *args: self._slider_update_estimate())

    def _on_A_ratio_fix_toggled(self, *args):
        """
        Stores ratio for A_u/A_l when ratio fixing is enabled/disabled.
        """
        if self.A_u_fix.get():
            Al = self.A_l_var.get()
            Au = self.A_u_var.get()
            self.A_ratio = Au / Al if Al != 0 else None
        else:
            self.A_ratio = None

    def _on_B_ratio_fix_toggled(self, *args):
        """
        Stores ratio for B_u/B_l when ratio fixing is enabled/disabled.
        """
        if self.Bu_Bl_fix.get():
            Bl = self.B_l_var.get()
            Bu = self.B_u_var.get()
            self.B_ratio = Bu / Bl if Bl != 0 else None
        else:
            self.B_ratio = None

    def _on_A_var_changed(self, changed):
        """
        fix for A_l/A_u when enabled.
        """
        if self.A_u_fix.get() and self.A_ratio is not None and not self._updating_A:
            try:
                self._updating_A = True
                Al = self.A_l_var.get()
                Au = self.A_u_var.get()
                if changed == "A_l":
                    self.A_u_var.set(Al * self.A_ratio)
                elif changed == "A_u":
                    if self.A_ratio != 0:
                        self.A_l_var.set(Au / self.A_ratio)
            finally:
                self._updating_A = False
        self._slider_update_estimate()

    def _on_B_var_changed(self, changed):
        """
        fix for B_l/B_u when enabled.
        """
        if self.Bu_Bl_fix.get() and self.B_ratio is not None and not self._updating_B:
            try:
                self._updating_B = True
                Bl = self.B_l_var.get()
                Bu = self.B_u_var.get()
                if changed == "B_l":
                    self.B_u_var.set(Bl * self.B_ratio)
                elif changed == "B_u":
                    if self.B_ratio != 0:
                        self.B_l_var.set(Bu / self.B_ratio)
            finally:
                self._updating_B = False
        self._slider_update_estimate()

    def _slider_update_estimate(self):
        """
        Only update estimate if estimate is visible.
        """
        if self.estimate_visible:
            self._update_estimate()

    def _import_scan(self):
        """
        Opens a dialogue to import a binned scan CSV and plots the data.
        """
        path = filedialog.askopenfilename(title="Select binned scan CSV", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        self.filepath = path
        self.data = pd.read_csv(self.filepath, sep=";")
        self.fit_curve = None
        self.estimate_curve = None
        self.fit_x = None
        self._plot_data()
        if self.estimate_visible:
            self._update_estimate()

    def _plot_data(self, fit_curve=None, estimate_curve=None, fit_x=None):
        """
        Plots the imported scan data and the fit and/or estimate.
        """
        self.fig.clf()
        self.ax = self.fig.subplots()  # <--- Always create a new axes after clf()
        if self.data is None:
            self.ax.set_title("No scan loaded.")
            self.canvas.draw()
            return
    
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
    
        # Main spectrum plot
        self.ax.errorbar(x, y, yerr=yerr, fmt='o', color='red', markersize=2, ecolor='k', capsize=2, label='Data')
        if estimate_curve is not None and fit_x is not None:
            self.ax.plot(fit_x, estimate_curve, color='green', label='Estimate')
        if fit_curve is not None and fit_x is not None:
            self.ax.plot(fit_x, fit_curve, color='blue', label='SATLAS2 Fit')
    
        self.ax.set_ylabel('Countrate per bunch')
        self.ax.legend()
        self.ax.set_title(os.path.basename(self.filepath) if self.filepath else "Scan")
        self.ax.set_xlabel("Frequency offset / MHz")
        self.ax.grid(True)
    
        self.fig.tight_layout()
        self.canvas.draw()

    def _fit_peaks(self):
        """
        SATLAS2 fitting to the imported scan data.
        """
        if self.data is None:
            messagebox.showerror("No scan loaded", "Please import a scan CSV first.")
            return
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()

        # Prepare fit parameters
        params = self._collect_fit_params()
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
        hfs.params['Al'].vary = not self.A_l_fix.get()
        hfs.params['Au'].vary = not self.A_u_fix.get()
        hfs.params['Bl'].vary = not self.B_l_fix.get()
        hfs.params['Bu'].vary = not self.B_u_fix.get()
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
        self.last_fit_hfs = hfs  #Store the fitted HFS object for later use
        print(report)

        fit_x = np.linspace(x.min(), x.max(), 500)
        fit_y = source.evaluate(fit_x)

        # Store and plot fit
        self.fit_x = fit_x
        self.fit_curve = fit_y
        
        # Hide estimate after fit if present
        if self.estimate_visible:
            self.estimate_visible = False
            self.estimate_curve = None
       
        # Display the fit curve
        self._plot_data(fit_curve=self.fit_curve,
                        estimate_curve=self.estimate_curve if self.estimate_visible else None,
                        fit_x=self.fit_x)

        # Display fit results panel 
        try:
            fit_params = hfs.params
            self.fit_results_box.config(state='normal')
            self.fit_results_box.delete(1.0, tk.END)
            self.fit_results_box.insert(tk.END, "---- Fit Results ----\n")
            for k in fit_params:
                v = fit_params[k].value
                self.fit_results_box.insert(tk.END, f"{k:12s}: {v:.5g}\n")
            self.fit_results_box.config(state='disabled')
        except Exception as e:
            print("Error writing to results box:", e)
            self.fit_results_box.config(state='normal')
            self.fit_results_box.delete(1.0, tk.END)
            self.fit_results_box.insert(tk.END, f"Could not display results: {e}")
            self.fit_results_box.config(state='disabled')

        messagebox.showinfo("Fit finished", "Fit finished! Expand plot for more detail.")

    def _show_estimate(self):
        """
        Toggles estimate on or off. If on, update and allow live slider updates. If off, hide and stop updating.
        """
        if not self.estimate_visible:
            # Turn ON
            self.estimate_visible = True
            self._update_estimate()
        else:
            # Turn OFF
            self.estimate_visible = False
            self.estimate_curve = None
            self._plot_data(fit_curve=self.fit_curve, estimate_curve=None, fit_x=self.fit_x)

    def _update_estimate(self):
        """
        Always updates and displays the SATLAS2 spectrum for current parameter values (without fitting).
        Only visible if self.estimate_visible is True.
        """
        if self.data is None or not self.estimate_visible:
            return
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
        params = self._collect_fit_params()
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

        self.estimate_curve = estimate_curve
        self.fit_x = fit_x

        self._plot_data(
            fit_curve=self.fit_curve,
            estimate_curve=self.estimate_curve if self.estimate_visible else None,
            fit_x=self.fit_x
        )

    def _collect_fit_params(self):
        """
        Collects all fit parameter values from GUI.
        Returns parameter values.
        """
        return {
            "A_l": self.A_l_var.get(),
            "A_u": self.A_u_var.get(),
            "B_l": self.B_l_var.get(),
            "B_u": self.B_u_var.get(),
            "df": self.df_var.get(),
            "scale": self.scale_var.get(),
            "FWHMg": self.FWHMg_var.get(),
            "FWHMl": self.FWHMl_var.get(),
            "I": self.I_var.get(),
            "J_l": self.Jl_var.get(),
            "J_u": self.Ju_var.get(),
            "racah": self.racah_int.get(),
            "intensity": self.intensity_var.get(),
            "background": self.background_var.get(),
        }

    def _import_parameters(self):
        """
        Extracts A_l, A_u, B_l, B_u, and I from element CSV using given element symbol and mass number,
        and fills in/overwrites those fields.
        """
        element = self.element_var.get().strip()
        mass_number = self.massnumber_var.get().strip()
        if not element or not mass_number:
            messagebox.showerror("Error", "Please provide both element symbol and mass number.")
            return
        # Path to the Elements CSV folder
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
            # Pull values
            self.A_l_var.set(float(match['A_l'].values[0]))
            self.A_u_var.set(float(match['A_u'].values[0]))
            self.B_l_var.set(float(match['B_l'].values[0]))
            self.B_u_var.set(float(match['B_u'].values[0]))
            # Extract I (nuclear spin) and centroud
            if 'I' in match:
                self.I_var.set(float(match['I'].values[0]))
            if 'centroid' in match:
                self.df_var.set(float(match['centroid'].values[0]))
            messagebox.showinfo("Parameters Imported",
                                f"Extracted I, A_l, A_u, B_l, B_u for {element}-{mass_number}.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to import parameters: {e}")

    def _save_fit_parameters(self):
        """
        Saves fit parameters AND uncertainties as a new row in
        (current directory)/<Element>_Results/<mass>/Saved_Parameters.csv
        Columns: scan, A_l, A_l_err, A_u, A_u_err, B_l, B_l_err, B_u, B_u_err, centroid, centroid_err
        If scan exists: prompt to overwrite.
        Creates folder/file if missing.
        """
        import pandas as pd
        import csv

        if not hasattr(self, "last_fit_hfs") or self.last_fit_hfs is None:
            messagebox.showerror("No fit", "No fit has been performed yet. Fit data before saving parameters.")
            return

        element = self.element_var.get().strip()
        mass = self.massnumber_var.get().strip()
        if not element or not mass:
            messagebox.showerror("Error", "Element and mass number required in parameter fields.")
            return

        scan_number = ""
        scan_number = tk.simpledialog.askstring("Scan Number", "Enter scan number for this fit (e.g. 1234):")
        if not scan_number:
            messagebox.showerror("Cancelled", "Scan number is required to save parameters.")
            return
        # Check scan_number is an integer
        try:
            int(scan_number)
        except ValueError:
            messagebox.showerror("Invalid input", "Scan number must be an integer.")
            return

        root_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(root_dir, f"{element.capitalize()}_Results", str(mass))
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, "Saved_Parameters.csv")

        fit_params = self.last_fit_hfs.params
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

        # Overwrite option if scan exists already
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

        # Otherwise: append as new row
        with open(save_path, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            if not file_exists:
                writer.writerow(columns)
            writer.writerow(new_row)
        messagebox.showinfo("Success", f"Fit parameters (with uncertainties) saved to:\n{save_path}")
        
    def _expand_plot(self):
        """
        Opens a new top-level window with a larger version of the current plot,
        including residuals if a fit is present.
        """
        if self.data is None:
            messagebox.showerror("No scan loaded", "Import a scan first.")
            return
        
        top = tk.Toplevel(self)
        top.title("Expanded Plot")
        fig = Figure(figsize=(10, 6))  
    
        # If fit present, show main and residuals, else just main
        if self.fit_curve is not None and self.fit_x is not None:
            axs = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
            ax_main, ax_resid = axs
        else:
            ax_main = fig.subplots()
            ax_resid = None
    
        x = self.data['x'].to_numpy()
        y = (self.data['y'] / self.data['bunches']).to_numpy()
        yerr = (self.data['yerr'] / self.data['bunches']).to_numpy()
    
        # Plot main data
        ax_main.errorbar(x, y, yerr=yerr, fmt='o', color='red', markersize=2, ecolor='k', capsize=2, label='Data')
    
        # Plot estimate if desired
        if self.estimate_visible and self.estimate_curve is not None and self.fit_x is not None:
            ax_main.plot(self.fit_x, self.estimate_curve, color='green', label='Estimate')
    
        # Plot fit if desired
        if self.fit_curve is not None and self.fit_x is not None:
            ax_main.plot(self.fit_x, self.fit_curve, color='blue', label='SATLAS2 Fit')
    
        ax_main.set_ylabel('Countrate per bunch')
        ax_main.legend()
        ax_main.set_title(os.path.basename(self.filepath) if self.filepath else "Scan")
        ax_main.grid(True)
    
        # Plot residuals when fit present
        if ax_resid is not None:
            # Interpolate fit to data x-points for residuals
            from scipy.interpolate import interp1d
            fit_interp = interp1d(self.fit_x, self.fit_curve, kind='linear', fill_value="extrapolate")
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
        canvas.get_tk_widget().pack(fill='both', expand=True)
        toolbar = NavigationToolbar2Tk(canvas, top)
        toolbar.update()
        toolbar.pack(side="top", fill="x")
