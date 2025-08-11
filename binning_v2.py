# -*- coding: utf-8 -*-
"""
Created on Mon Jul 7 12:22:47 2025

binning_v2.py

Version 2

Performs binning and file combining capabilities (and associated plots) of the interactive GUI.
User can export data single scans or combined scans formatted ready for fitting and analysis.
User adds queue list of scans which are combined as raw data and processed together.
Toggle to show TOF histogram or just the hyperfine spectra.

Author: Shay Rowlinson
"""


import os
import tkinter as tk
from tkinter import ttk, filedialog, StringVar, messagebox
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
    User can add scans to a queue, then bin all as one dataset.
    """
    def __init__(self, parent):
        """
        Initialises the binning frame and widgets.
        """
        super().__init__(parent)
        self.grid_columnconfigure(0, weight=1)
        self.binned_data = None
        self.scan_queue = []
        self._build_gui()

    def _build_gui(self):
        """
        Builds the layout for the Binning tab.
        """
        header = ttk.Label(
            self, text="CRIS scan binning for hyperfine analysis",
            font=("Segoe UI", 18, "bold")
        )
        header.grid(row=0, column=0, pady=(5, 8), sticky="ew")
        self._build_param_frame()
        self.status_var = tk.StringVar(value="Ready.")
        status = ttk.Label(self, textvariable=self.status_var, font=("Segoe UI", 10, "italic"), anchor='w', foreground="#555")
        status.grid(row=2, column=0, sticky="ew", pady=(2, 8))
        self._build_plot_frame()
        self._build_bottom_controls()

    def _build_param_frame(self):
        """
        Creates and places all input widgets for user parameters in a labeled frame.
        """
        param_frame = ttk.LabelFrame(self, text="Parameters", padding=(16,10))
        param_frame.grid(row=1, column=0, sticky="ew", padx=6, pady=4)
        for i in range(8): param_frame.columnconfigure(i, weight=1)

        ttk.Label(param_frame, text="Scan folder:").grid(row=0, column=0, sticky='e', pady=3)
        self.folder_var = StringVar()
        folder_entry = ttk.Entry(param_frame, textvariable=self.folder_var, width=20, state='readonly')
        folder_entry.grid(row=0, column=1, columnspan=5, sticky='ew', padx=(0,6), pady=3)
        browse_btn = ttk.Button(param_frame, text="Browse...", command=self.browse_scan_folder)
        browse_btn.grid(row=0, column=6, padx=(0,6))

        ttk.Label(param_frame, text="Element symbol:").grid(row=1, column=0, sticky='e', pady=3)
        self.element_var = StringVar()
        element_entry = ttk.Entry(param_frame, textvariable=self.element_var, width=6)
        element_entry.grid(row=1, column=1, sticky='w', pady=3)

        ttk.Label(param_frame, text="Mass number:").grid(row=1, column=2, sticky='e', pady=3)
        self.massnumber_var = StringVar()
        massnumber_entry = ttk.Entry(param_frame, textvariable=self.massnumber_var, width=7)
        massnumber_entry.grid(row=1, column=3, sticky='w', pady=3)

        ttk.Label(param_frame, text="Harmonic:").grid(row=2, column=2, sticky='e', pady=3)
        self.harmonic_var = StringVar(value="4")
        harmonic_entry = ttk.Entry(param_frame, textvariable=self.harmonic_var, width=8)
        harmonic_entry.grid(row=2, column=3, sticky='w', pady=3)

        ttk.Label(param_frame, text="TOF gate (μs):").grid(row=2, column=4, sticky='e', pady=3)
        self.tof_lower_var = StringVar(value="0")
        self.tof_upper_var = StringVar(value="30")
        tof_lower_entry = ttk.Entry(param_frame, textvariable=self.tof_lower_var, width=7)
        tof_upper_entry = ttk.Entry(param_frame, textvariable=self.tof_upper_var, width=7)
        tof_lower_entry.grid(row=2, column=5, sticky='w', padx=(0,1), pady=3)
        ttk.Label(param_frame, text="–").grid(row=2, column=6, sticky='w')
        tof_upper_entry.grid(row=2, column=7, sticky='w', pady=3)

        ttk.Label(param_frame, text="Bin size (MHz):").grid(row=3, column=0, sticky='e', pady=3)
        self.bin_size_var = StringVar(value="30")
        bin_size_entry = ttk.Entry(param_frame, textvariable=self.bin_size_var, width=10)
        bin_size_entry.grid(row=3, column=1, sticky='w', pady=3)
        
        bin_selected_btn = ttk.Button(param_frame, text="Bin Current Scan", command=self.bin_current_scan)
        bin_selected_btn.grid(row=0, column=7, padx=(16,0), sticky="ew")
        
        self.show_spectra_only = tk.BooleanVar(value=False)
        spectra_mode_check = ttk.Checkbutton(
            param_frame, text="Show Only Spectra", variable=self.show_spectra_only,
            command=self._redraw_plot_mode
        )
        spectra_mode_check.grid(row=3, column=7, sticky='e', padx=(6, 4), pady=(2, 0))

    def _build_plot_frame(self):
        """
        Creates canvas for both figure set ups.
        """
        plot_frame = ttk.LabelFrame(self, text="Results", padding=(6,4))
        plot_frame.grid(row=3, column=0, sticky="nsew", padx=6, pady=2)
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        self.fig = Figure(figsize=(3.45,4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")

    def _build_bottom_controls(self):
        """
        Controls for multiscan/exports
        """
        bottom_frame = ttk.Frame(self)
        bottom_frame.grid(row=4, column=0, sticky="ew", padx=8, pady=6)
        bottom_frame.columnconfigure(1, weight=1)

        # Add scans to queue
        add_btn = ttk.Button(bottom_frame, text="Add Scan", command=self.add_scan_to_queue)
        add_btn.grid(row=0, column=0, sticky="w", padx=(0, 8))

        # Listbox for queued scan
        lb_frame = ttk.Frame(bottom_frame)
        lb_frame.grid(row=0, column=1, padx=(0, 8), sticky="ew")
        self.scan_listbox = tk.Listbox(lb_frame, height=4, width=60)
        self.scan_listbox.pack(side="left", fill="both", expand=True)
        lb_scroll = ttk.Scrollbar(lb_frame, orient="vertical", command=self.scan_listbox.yview)
        lb_scroll.pack(side="right", fill="y")
        self.scan_listbox.config(yscrollcommand=lb_scroll.set)

        # Remove scans
        remove_btn = ttk.Button(bottom_frame, text="Remove Selected", command=self.remove_selected_scan)
        remove_btn.grid(row=0, column=2, sticky="w", padx=(0, 8))

        # Bin all scan
        binall_btn = ttk.Button(bottom_frame, text="Bin All", command=self.run_binning_thread)
        binall_btn.grid(row=0, column=3, sticky="w", padx=(0, 8))

        # Export CSVs
        self.export_button = ttk.Button(
            bottom_frame, text="Export data as CSV",
            command=self.export_csv, state="disabled"
        )
        self.export_button.grid(row=0, column=4, padx=(0,0))

    def browse_scan_folder(self):
        """
        Select scan location
        """
        path = filedialog.askdirectory(title="Select scan folder")
        if path:
            self.folder_var.set(path)

    def get_exact_mass(self, element, mass_number):
        """
        Find the exact mass from elements CSV for dopplershifting.
        """
        elements_folder = os.path.join(os.path.dirname(__file__), "Elements")
        csv_path = os.path.join(elements_folder, f"{element.capitalize()}.csv")
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"No data file found for {element} at {csv_path}")
        df = pd.read_csv(csv_path)
        match = df[df['Mass'] == int(mass_number)]
        if match.empty:
            raise ValueError(f"No isotope with mass {mass_number} found for {element}")
        return float(match['ExactMass'].values[0])

    def add_scan_to_queue(self):
        try:
            scan_folder = self.folder_var.get()
            if not scan_folder:
                raise ValueError("No scan folder selected.")
            element = self.element_var.get().strip()
            mass_number = self.massnumber_var.get().strip()
            if not element or not mass_number:
                raise ValueError("Element symbol and mass number must be provided.")
            exact_mass = self.get_exact_mass(element, mass_number)
            bin_size_MHz = float(self.bin_size_var.get())
            tof_lower = float(self.tof_lower_var.get())
            tof_upper = float(self.tof_upper_var.get())
            harmonic = int(self.harmonic_var.get())
        except Exception as e:
            messagebox.showerror("Parameter error", str(e))
            return

        # Save parameters for this scan
        scan_info = dict(
            scan_folder=scan_folder,
            element=element,
            mass_number=mass_number,
            exact_mass=exact_mass,
            bin_size_MHz=bin_size_MHz,
            tof_lower=tof_lower,
            tof_upper=tof_upper,
            harmonic=harmonic
        )
        self.scan_queue.append(scan_info)
        display = f"{os.path.basename(scan_folder)} | {element}-{mass_number} | TOF: {tof_lower}-{tof_upper}μs"
        self.scan_listbox.insert(tk.END, display)
        self.set_status(f"Added scan: {display}")

    def remove_selected_scan(self):
        selection = self.scan_listbox.curselection()
        if not selection:
            return
        idx = selection[0]
        self.scan_listbox.delete(idx)
        del self.scan_queue[idx]
        self.set_status("Removed scan from queue.")

    def export_csv(self):
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
        Update the GUI’s status message
        """
        self.status_var.set(msg)
        self.update_idletasks()

    def set_running_state(self, running=True):
        """
        Enable or disable the export button depending on whether a process is running.
        """
        state = 'disabled' if running else 'normal'
        self.export_button.config(state=state if running else ('normal' if self.binned_data is not None else 'disabled'))

    def run_binning_thread(self):
        """
        Start binning in a background thread so the GUI doesn’t freeze.
        """
        self.set_running_state(True)
        self.set_status("Running binning on all queued scans...")
        threading.Thread(target=self.run_binning, daemon=True).start()

    def run_binning(self):
        """
        Run binning.
        """
        if not self.scan_queue:
            self.set_status("No scans in queue to bin.")
            self.set_running_state(False)
            return

        try:
            all_data = []
            D_template = None
            for idx, scan in enumerate(self.scan_queue):
                D = Dopplershift_data(
                    mass=int(scan["mass_number"]), scan=idx+1, voltage_scanning=False,
                    wn_channel='wavenumber_1', exact_mass=scan["exact_mass"],
                    data_folder=scan["scan_folder"], save_data_path=scan["scan_folder"], save_fig_path=scan["scan_folder"]
                )
                D._exact_mass = scan["exact_mass"]
                D._PATH = scan["scan_folder"] if scan["scan_folder"].endswith(os.sep) else scan["scan_folder"] + os.sep
                D.bin_size_voltage = 4.0
                D.bin_size_MHz = scan["bin_size_MHz"]

                data = D.extract_raw_data(devices_to_read=D._devices, path=D._PATH)
                data = D.AdvCutNoise(data, threshold=0.15)
                data = D.gate_tof(data, manual=[scan["tof_lower"], scan["tof_upper"]])
                data = D.filter_scatter(data, filename='iscool2', method='avg', ISCOOL_voltage_multiplier=10001.645)
                data = D.calibrate_CRIS_voltagemeter(data, calibration_factor=1.005030)
                data = D.gate_wavenumber(data, wavenumber=D._wn_channel)
                all_data.append(data)
                if D_template is None:
                    D_template = D  # use first scan's D for binning

            combined_data = pd.concat(all_data, ignore_index=True)
            self.set_status("Binning combined data...")
            self.binned_data = D_template.bin_wm_data(combined_data, freq_multiplier=scan["harmonic"])
            self.set_status("Plotting...")
            self.after(0, self.plot_results, combined_data, D_template)
            self.set_status("Binning complete. Data may now be exported.")
            self.set_running_state(False)
        except Exception as e:
            self.set_status(f"Error during binning: {e}")
            messagebox.showerror(
                "Binning Error",
                f"Error: {e}\nEither the mass/scan combination is wrong or the path is wrong."
            )
            self.set_running_state(False)
            
    def bin_current_scan(self):
        """
        Bin only the currently selected scan, not any in the queue.
        """
        try:
            scan_folder = self.folder_var.get()
            if not scan_folder:
                raise ValueError("No scan folder selected.")
            element = self.element_var.get().strip()
            mass_number = self.massnumber_var.get().strip()
            if not element or not mass_number:
                raise ValueError("Element symbol and mass number must be provided.")
            exact_mass = self.get_exact_mass(element, mass_number)
            bin_size_MHz = float(self.bin_size_var.get())
            tof_lower = float(self.tof_lower_var.get())
            tof_upper = float(self.tof_upper_var.get())
            harmonic = int(self.harmonic_var.get())
        except Exception as e:
            messagebox.showerror("Parameter error", str(e))
            return
    
        try:
            D = Dopplershift_data(
                mass=int(mass_number), scan=1, voltage_scanning=False,
                wn_channel='wavenumber_1', exact_mass=exact_mass,
                data_folder=scan_folder, save_data_path=scan_folder, save_fig_path=scan_folder
            )
            D._exact_mass = exact_mass
            D._PATH = scan_folder if scan_folder.endswith(os.sep) else scan_folder + os.sep
            D.bin_size_voltage = 4.0
            D.bin_size_MHz = bin_size_MHz
    
            data = D.extract_raw_data(devices_to_read=D._devices, path=D._PATH)
            data = D.AdvCutNoise(data, threshold=0.15)
            data = D.gate_tof(data, manual=[tof_lower, tof_upper])
            data = D.filter_scatter(data, filename='iscool2', method='avg', ISCOOL_voltage_multiplier=10001.645)
            data = D.calibrate_CRIS_voltagemeter(data, calibration_factor=1.005030)
            data = D.gate_wavenumber(data, wavenumber=D._wn_channel)
            binned = D.bin_wm_data(data, freq_multiplier=harmonic)
            self.binned_data = binned
            self.set_status(f"Binned single scan: {element}-{mass_number} from {os.path.basename(scan_folder)}")
            self.after(0, self.plot_results, data, D)
            self.set_running_state(False)
        except Exception as e:
            self.set_status(f"Error: {e}")
            messagebox.showerror("Error", f"Could not bin current scan:\n{e}")

    def plot_results(self, data, D):
        """
        Plot spectra OR full 3D layout
        """
        import matplotlib.pyplot as plt
        import numpy as np
    
        # Store for redrawing
        self.last_data = data
        self.last_D = D
    
        self.fig.clf()
    
        if hasattr(self, "show_spectra_only") and self.show_spectra_only.get():
            # --- Only plot the top hyperfine spectrum ---
            ax_top = self.fig.subplots()
    
            if self.binned_data is not None and not self.binned_data.empty:
                ax_top.errorbar(
                    self.binned_data['x'],
                    self.binned_data['y'] / self.binned_data['bunches'],
                    xerr=self.binned_data['xerr'],
                    yerr=self.binned_data['yerr'] / self.binned_data['bunches'],
                    fmt='o', markersize=2, capsize=2,
                    color='red', ecolor='k'
                )
    
            ax_top.set_ylabel("Counts per bunch", fontsize=8)
            ax_top.set_title("Hyperfine Spectrum", fontsize=13)
            ax_top.grid(True, linestyle='--', alpha=0.25)
            ax_top.set_xlabel("Relative Frequency / MHz", fontsize=12)
            ax_top.tick_params(axis='both', labelsize=9)
    
            self.fig.tight_layout()
            self.canvas.draw()
            return
    
        # --- Defalt: full 3 plot layout---
        gs = GridSpec(
            4, 5,
            figure=self.fig,
            width_ratios=[1, 4, 4, 4, 0.25],
            hspace=0.25,
            wspace=0.25
        )
        ax_main = self.fig.add_subplot(gs[1:, 1:4])
        ax_top = self.fig.add_subplot(gs[0, 1:4], sharex=ax_main)
        ax_left = self.fig.add_subplot(gs[1:, 0], sharey=ax_main)
        cax = self.fig.add_subplot(gs[1:, 4])
    
        plt.setp(ax_top.get_xticklabels(), visible=False)
        plt.setp(ax_left.get_yticklabels(), visible=False)
    
        if 'delta_t' in data.columns and 'wavenumber_1' in data.columns:
            valid = (data['delta_t'] > 1) & (data['wavenumber_1'].notna())
            freq = data.loc[valid, 'wavenumber_1'] * 29979.2458
            transition_frequency = D._transition_wavenumber * 29979.2458
            rel_freq = freq - transition_frequency
            time = data.loc[valid, 'delta_t']
        else:
            rel_freq = np.array([])
            time = np.array([])
    
        h, xedges, yedges, img = ax_main.hist2d(
            rel_freq, time, bins=[75, 100], cmap='Blues', cmin=1
        )
        cb = self.fig.colorbar(img, cax=cax)
        cb.set_label("Counts", fontsize=11)
        cb.ax.tick_params(labelsize=9)
    
        ax_main.set_xlabel("Relative Frequency / MHz", fontsize=12)
        ax_main.set_ylabel("TOF (μs)", fontsize=13)
    
        if self.binned_data is not None and not self.binned_data.empty:
            ax_top.errorbar(
                self.binned_data['x'],
                self.binned_data['y'] / self.binned_data['bunches'],
                xerr=self.binned_data['xerr'],
                yerr=self.binned_data['yerr'] / self.binned_data['bunches'],
                fmt='.', markersize=2, capsize=2,
                color='k', ecolor='white'
            )
        ax_top.set_ylabel("Counts per bunch", fontsize=8)
        ax_top.set_title("Hyperfine Spectrum", fontsize=13)
        ax_top.grid(True, linestyle='--', alpha=0.25)
        ax_top.tick_params(axis='x', labelbottom=False, which='both')
        ax_top.tick_params(axis='y', labelsize=9)
    
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
        ax_left.invert_xaxis()
        ax_left.tick_params(axis='y', labelleft=False, which='both')
        ax_left.tick_params(axis='x', labelsize=9)
    
        ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        ax_left.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        ax_main.tick_params(axis='both', labelsize=10)
    
        self.fig.tight_layout()
        self.canvas.draw()

        
    def _redraw_plot_mode(self):
        """
        Redraw plot if mode changed.
        """
        if hasattr(self, "last_data") and hasattr(self, "last_D"):
            self.plot_results(self.last_data, self.last_D)
