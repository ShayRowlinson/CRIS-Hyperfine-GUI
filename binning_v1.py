# -*- coding: utf-8 -*-
"""
Created on Mon Jul 7 12:22:47 2025

binning_v1.py

Version 1

Performs binning and file combining capabilities (and associated plots) of the interactive GUI.
User can export data single scans or combined scans formatted ready for fitting and analysis.

Author: Shay Rowlinson
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, StringVar
import threading
import pandas as pd
import matplotlib
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from Dopplershift_analysis import Dopplershift_data


class Binning(ttk.Frame):
    """
    Dopplershifts and bins CRIS scan data.
    User selects scan folder and parameters, bins and plots data, and can export binned as a CSV.
    Now supports combining binned and raw data from multiple scans for a unified combined export and plot.
    """
    def __init__(self, parent):
        """
        Initialises the binning frame and widgets.
        """
        super().__init__(parent)
        self.grid_columnconfigure(0, weight=1)
        self._build_gui()
        self.binned_data = None
        self.combined_data_list = []
        self.combined_rawdata_list = []
        self.combined_scan_labels = []

    def _build_gui(self):
        """
        Builds the layout for the Binning tab.
        """
        header = ttk.Label(
            self, text="CRIS scan binning for hyperfine analysis ",
            font=("Segoe UI", 18, "bold")
        )
        header.grid(row=0, column=0, pady=(5, 8), sticky="ew")
        self._build_param_frame()
        self.status_var = tk.StringVar(value="Ready.")
        status = ttk.Label(self, textvariable=self.status_var, font=("Segoe UI", 10, "italic"), anchor='w', foreground="#555")
        status.grid(row=2, column=0, sticky="ew", pady=(2, 8))
        self._build_plot_frame()

        btn_frame = ttk.Frame(self)
        btn_frame.grid(row=4, column=0, sticky="e", pady=(8,2))
        self.export_button = ttk.Button(
            btn_frame, text="Export data as CSV",
            command=self.export_csv, state="disabled"
        )
        self.export_button.grid(row=0, column=0, padx=6)
    
    def _build_param_frame(self):
        """
        Creates and places all input widgets for user parameters in a labeled frame.
        """
        param_frame = ttk.LabelFrame(self, text="Parameters", padding=(16,10))
        param_frame.grid(row=1, column=0, sticky="ew", padx=6, pady=4)
        for i in range(8): param_frame.columnconfigure(i, weight=1)
       
        # Scan Folder
        ttk.Label(param_frame, text="Scan folder:").grid(row=0, column=0, sticky='e', pady=3)
        self.folder_var = StringVar()
        folder_entry = ttk.Entry(param_frame, textvariable=self.folder_var, width=20, state='readonly')
        folder_entry.grid(row=0, column=1, columnspan=5, sticky='ew', padx=(0,6), pady=3)
        browse_btn = ttk.Button(param_frame, text="Browse...", command=self.browse_scan_folder)
        browse_btn.grid(row=0, column=6, padx=(0,6))
       
        # Element Symbol 
        ttk.Label(param_frame, text="Element symbol:").grid(row=1, column=0, sticky='e', pady=3)
        self.element_var = StringVar()
        element_entry = ttk.Entry(param_frame, textvariable=self.element_var, width=6)
        element_entry.grid(row=1, column=1, sticky='w', pady=3)
        
        # Mass Number 
        ttk.Label(param_frame, text="Mass number:").grid(row=1, column=2, sticky='e', pady=3)
        self.massnumber_var = StringVar()
        massnumber_entry = ttk.Entry(param_frame, textvariable=self.massnumber_var, width=7)
        massnumber_entry.grid(row=1, column=3, sticky='w', pady=3)
       
        # Harmonic
        ttk.Label(param_frame, text="Harmonic:").grid(row=2, column=2, sticky='e', pady=3)
        self.harmonic_var = StringVar(value="4")
        harmonic_entry = ttk.Entry(param_frame, textvariable=self.harmonic_var, width=8)
        harmonic_entry.grid(row=2, column=3, sticky='w', pady=3)
       
        # TOF Gate
        ttk.Label(param_frame, text="TOF gate (μs):").grid(row=2, column=4, sticky='e', pady=3)
        self.tof_lower_var = StringVar(value="0")
        self.tof_upper_var = StringVar(value="30")
        tof_lower_entry = ttk.Entry(param_frame, textvariable=self.tof_lower_var, width=7)
        tof_upper_entry = ttk.Entry(param_frame, textvariable=self.tof_upper_var, width=7)
        tof_lower_entry.grid(row=2, column=5, sticky='w', padx=(0,1), pady=3)
        ttk.Label(param_frame, text="–").grid(row=2, column=6, sticky='w')
        tof_upper_entry.grid(row=2, column=7, sticky='w', pady=3)
      
        # Bin size
        ttk.Label(param_frame, text="Bin size (MHz):").grid(row=3, column=0, sticky='e', pady=3)
        self.bin_size_var = StringVar(value="30")
        bin_size_entry = ttk.Entry(param_frame, textvariable=self.bin_size_var, width=10)
        bin_size_entry.grid(row=3, column=1, sticky='w', pady=3)
        
        # Rebin Button
        self.rebin_button = ttk.Button(param_frame, text="Rebin", command=self.run_binning_thread)
        self.rebin_button.grid(row=3, column=6, sticky='e', pady=6)

    def _build_plot_frame(self):
        """
        Creates and places the matplotlib Figure and canvas inside a labeled frame.
        """
        plot_frame = ttk.LabelFrame(self, text="Results", padding=(6,4))
        plot_frame.grid(row=3, column=0, sticky="nsew", padx=6, pady=2)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        self.fig = Figure(figsize=(3.5,4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")

    def browse_scan_folder(self):
        """
        Dialogue for the user to select the scan folder.
        """
        path = filedialog.askdirectory(title="Select scan folder")
        if path:
            self.folder_var.set(path)

    def get_exact_mass(self, element, mass_number):
        # Path to the Elements CSV folder
        elements_folder = os.path.join(os.path.dirname(__file__), "Elements")
        csv_path = os.path.join(elements_folder, f"{element.capitalize()}.csv")
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"No data file found for {element} at {csv_path}")
        df = pd.read_csv(csv_path)
        # Try to find row with this mass number
        match = df[df['Mass'] == int(mass_number)]
        if match.empty:
            raise ValueError(f"No isotope with mass {mass_number} found for {element}")
        return float(match['ExactMass'].values[0])

    def export_csv(self):
        """
        Opens save dialogue and writes current binned data to CSV.
        Uses pop-up dialogs for success and error messages.
        """
        import tkinter.messagebox as mb
    
        if self.binned_data is not None:
            save_path = filedialog.asksaveasfilename(
                defaultextension=".csv", filetypes=[("CSV", "*.csv")]
            )
            if save_path:
                try:
                    self.binned_data.to_csv(save_path, sep=';')
                    mb.showinfo("Export Successful", f"Binned data saved to:\n{save_path}")
                except Exception as e:
                    mb.showerror("Save Failed", f"Failed to save file:\n{e}")
        else:
            mb.showwarning("No Data", "No binned data to export.")

    def set_status(self, msg):
        """
        Updates the status message label.
        """
        self.status_var.set(msg)
        self.update_idletasks()

    def set_running_state(self, running=True):
        """
        Enables or disables buttons while running binning.
        """
        state = 'disabled' if running else 'normal'
        self.rebin_button.config(state=state)
        self.export_button.config(state=state if running else ('normal' if self.binned_data is not None else 'disabled'))
        
    def run_binning_thread(self):
        """
        Launches binning processing in a background thread so GUI isn't blocked.
        """
        self.set_running_state(True)
        self.set_status("Running binning...")
        threading.Thread(target=self.run_binning, daemon=True).start()

    def run_binning(self):
        """
        Loads scan folder and parameters from user input, extracts exact mass from Elements CSV using element symbol and mass number, 
        runs Dopplershift_data binning, plots results.
        """
        try:
            scan_folder = self.folder_var.get()
            if not scan_folder:
                raise ValueError("No scan folder selected.")
    
            element = self.element_var.get().strip()
            mass_number = self.massnumber_var.get().strip()
            if not element or not mass_number:
                raise ValueError("Element symbol and mass number must be provided.")
    
            # Get exact mass from Elements CSV
            exact_mass = self.get_exact_mass(element, mass_number)
    
            bin_size_MHz = float(self.bin_size_var.get())
            tof_lower = float(self.tof_lower_var.get())
            tof_upper = float(self.tof_upper_var.get())
            harmonic = int(self.harmonic_var.get())
        except Exception:
            self.set_status("Parameter error.")
            self.set_running_state(False)
            return
    
        save_data_path = scan_folder
        save_fig_path = scan_folder
    
        try:
            D = Dopplershift_data(
                mass=int(mass_number), scan=1, voltage_scanning=False,
                wn_channel='wavenumber_1', exact_mass=exact_mass,
                data_folder=scan_folder, save_data_path=save_data_path, save_fig_path=save_fig_path
            )
            D._exact_mass = exact_mass
            D._PATH = scan_folder if scan_folder.endswith(os.sep) else scan_folder + os.sep
            D.bin_size_voltage = 4.0
            D.bin_size_MHz = bin_size_MHz
    
            self.set_status("Extracting and processing data...")
            data = D.extract_raw_data(devices_to_read=D._devices, path=D._PATH)
            data = D.AdvCutNoise(data, threshold=0.15)
            data = D.gate_tof(data, manual=[tof_lower, tof_upper])
            data = D.filter_scatter(data, filename='iscool2', method='avg', ISCOOL_voltage_multiplier=10001.645)
            data = D.calibrate_CRIS_voltagemeter(data, calibration_factor=1.005030)
            #data = D.apply_wavenumber_correction(data)
            data = D.gate_wavenumber(data, wavenumber=D._wn_channel)
            self.binned_data = D.bin_wm_data(data, freq_multiplier=harmonic)
    
            self.set_status("Plotting...")
            self.after(0, self.plot_results, data, D)
            self.set_status("Binning complete. Data may now be exported.")
            self.set_running_state(False)
        except Exception:
            self.set_status("Error during binning.")
            import tkinter.messagebox as mb
            mb.showerror(
                "Binning Error",
                "Error: Either the mass/scan combination is wrong or the path is wrong."
            )
            self.set_running_state(False)
        self.after(0, self.plot_results, data)

    def plot_results(self, data, D):
        """
        Plots as:
        - Top: Doppler-shifted spectrum (binned, Counts/bunch vs freq)
        - Bottom left: TOF histogram (Counts vs TOF)
        - Bottom right: 2D histogram (freq vs TOF), sharing axes.
        Includes print debugging to diagnose frequency axis alignment.
        """
        import matplotlib.pyplot as plt
        import numpy as np
    
        self.fig.clf()
        gs = GridSpec(
            4, 5,  # rows, cols
            figure=self.fig,
            width_ratios=[1, 4, 4, 4, 0.25],  # last col = colorbar
            hspace=0.25,
            wspace=0.25
        )
    
        ax_main = self.fig.add_subplot(gs[1:, 1:4])
        ax_top = self.fig.add_subplot(gs[0, 1:4], sharex=ax_main)
        ax_left = self.fig.add_subplot(gs[1:, 0], sharey=ax_main)
        cax = self.fig.add_subplot(gs[1:, 4])
    
        # Hide x tick labels on ax_top and y tick labels on ax_left
        plt.setp(ax_top.get_xticklabels(), visible=False)
        plt.setp(ax_left.get_yticklabels(), visible=False)
    
        # Data extraction
        if 'delta_t' in data.columns and 'wavenumber_1' in data.columns:
            valid = (data['delta_t'] > 1) & (data['wavenumber_1'].notna())
            freq = data.loc[valid, 'wavenumber_1'] * 29979.2458
            transition_frequency = D._transition_wavenumber * 29979.2458
            rel_freq = freq - transition_frequency
            time = data.loc[valid, 'delta_t']
        else:
            rel_freq = np.array([])
            time = np.array([])
    
        #  Main 2D histogram 
        h, xedges, yedges, img = ax_main.hist2d(
            rel_freq, time, bins=[75, 100], cmap='Blues', cmin=1
        )
        cb = self.fig.colorbar(img, cax=cax)
        cb.set_label("Counts", fontsize=11)
        cb.ax.tick_params(labelsize=9)
    
        ax_main.set_xlabel("Relative Frequency / MHz", fontsize=12)
        ax_main.set_ylabel("TOF (μs)", fontsize=13)
    
        #  Spectrum above (Counts/bunch vs freq) 
        if self.binned_data is not None and not self.binned_data.empty:
            ax_top.errorbar(
                self.binned_data['x'],
                self.binned_data['y'] / self.binned_data['bunches'],
                xerr=self.binned_data['xerr'],
                yerr=self.binned_data['yerr'] / self.binned_data['bunches'],
                fmt='.', markersize=3, capsize=2,
                color='k', ecolor='k'
            )
        ax_top.set_ylabel("Counts per bunch", fontsize=8)
        ax_top.set_title("Hyperfine Spectrum", fontsize=13)
        ax_top.grid(True, linestyle='--', alpha=0.25)
        ax_top.tick_params(axis='x', labelbottom=False, which='both')
        ax_top.tick_params(axis='y', labelsize=9)
    
        #  TOF histogram 
        if time is not None and len(time) > 0:
            ax_left.hist(
                time[time > 0],
                bins=50,
                orientation='horizontal',
                color='cornflowerblue',
                edgecolor='navy',
                alpha=0.85,
                linewidth=1.0
            )
        ax_left.set_xlabel("Counts", fontsize=11)
        ax_left.set_ylabel("Time / μs", fontsize=11)
        ax_left.invert_xaxis()  # To match your style: counts increasing leftward
        ax_left.tick_params(axis='y', labelleft=False, which='both')
        ax_left.tick_params(axis='x', labelsize=9)
    
        # Remove overlapping ticks for clarity
        ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax_left.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        ax_main.tick_params(axis='both', labelsize=10)
    
        self.fig.tight_layout()
        self.canvas.draw()
