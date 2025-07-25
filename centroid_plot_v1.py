# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 10:49:03 2025

Centroid_plot_v1.py

Version 1

User enters desired isotope and program extracts centroid values from saved 
fitting parameter data bank, then plots centroid value vs scan number

Author: Shay Rowlinson
"""

import os
import numpy as np
from tkinter import ttk, messagebox, StringVar
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class CentroidPlot(ttk.Frame):
    """
    Tab for plotting centroid frequency vs scan number for chosen isotope.
    """
    def __init__(self, parent):
        """
        Initialise centroid plot tab and widgets.
        """
        super().__init__(parent)
        self._build_gui()
        self.data = None

    def _build_gui(self):
        """
        Build the layout for the centroid plotting tab.
        """
        # --- Input area ---
        param_frame = ttk.LabelFrame(self, text="Isotope selection", padding=(10, 6))
        param_frame.grid(row=0, column=0, sticky="ew", padx=2, pady=(10, 0))

        ttk.Label(param_frame, text="Element symbol:").grid(row=0, column=0, sticky='e', pady=3)
        self.element_var = StringVar()
        ttk.Entry(param_frame, textvariable=self.element_var, width=6).grid(row=0, column=1, sticky='w', pady=3)

        ttk.Label(param_frame, text="Mass number:").grid(row=0, column=2, sticky='e', pady=3)
        self.massnumber_var = StringVar()
        ttk.Entry(param_frame, textvariable=self.massnumber_var, width=7).grid(row=0, column=3, sticky='w', pady=3)

        ttk.Button(param_frame, text="Load & Plot", command=self._load_and_plot).grid(row=0, column=4, padx=12)

        # --- Plot area ---
        plot_frame = ttk.LabelFrame(self, text="Centroid Frequency Plot", padding=(6,4))
        plot_frame.grid(row=1, column=0, sticky="nsew", padx=8, pady=8)
        self.fig = Figure(figsize=(8,3))
        self.ax = self.fig.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")

    def _load_and_plot(self):
        """
        Loads centroid data from the correct CSV based on user input and plots centroid vs scan number.
        Calculates and plots the mean and ±1σ lines.
        """
        # Find correct CSV based on element and mass
        element = self.element_var.get().strip()
        mass = self.massnumber_var.get().strip()
        if not element or not mass:
            messagebox.showerror("Error", "Element and mass number required.")
            return
    
        root_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(root_dir, f"{element.capitalize()}_Results", str(mass))
        save_path = os.path.join(save_dir, "Saved_Parameters.csv")
    
        if not os.path.exists(save_path):
            messagebox.showerror("Not found", f"No saved parameters found for {element}-{mass} at:\n{save_path}")
            self.ax.clear()
            self.ax.set_title("No data loaded")
            self.canvas.draw()
            return
    
        try:
            df = pd.read_csv(save_path, dtype={'scan': str})
            if 'scan' not in df.columns or 'centroid' not in df.columns:
                raise ValueError("File missing required columns.")
            df = df.dropna(subset=['scan', 'centroid'])
            df['scan_int'] = df['scan'].astype(int)
            df = df.sort_values('scan_int')
    
            scan_x = df['scan_int'].to_numpy()
            centroid_y = df['centroid'].astype(float).to_numpy()
            centroid_err = None
            if 'centroid_err' in df.columns:
                try:
                    centroid_err = pd.to_numeric(df['centroid_err'], errors='coerce').to_numpy()
                except:
                    centroid_err = None
    
            # --- Plot ---
            self.ax.clear()
    
            # Plot error bars
            if centroid_err is not None and not np.all(np.isnan(centroid_err)):
                self.ax.errorbar(scan_x, centroid_y, yerr=centroid_err, fmt='none', ecolor='k', capsize=4)
    
            # Plot red circles for data
            self.ax.plot(scan_x, centroid_y, 'o', color='red', markersize=3, label='Centroid frequency')
    
            # mean and std
            if centroid_err is not None and not np.all(np.isnan(centroid_err)):
                valid = ~np.isnan(centroid_y) & ~np.isnan(centroid_err) & (centroid_err > 0)
                if np.any(valid):
                    weights = 1 / centroid_err[valid]**2
                    mean = np.average(centroid_y[valid], weights=weights)
                    var = np.average((centroid_y[valid] - mean) ** 2, weights=weights)
                    std = np.sqrt(var)
                else:
                    mean = np.mean(centroid_y)
                    std = np.std(centroid_y, ddof=1)
            else:
                mean = np.mean(centroid_y)
                std = np.std(centroid_y, ddof=1)
    
            # Plot mean and ±1σ lines
            self.ax.axhline(mean, color='green', linestyle='-', label=f"Mean: {mean:.3g} / MHz")
            self.ax.axhline(mean + 0.5 * std, color='green', linestyle=':', label="Width between lines = σ")
            self.ax.axhline(mean - 0.5 * std, color='green', linestyle=':')
    
            self.ax.set_xlabel("Scan number")
            self.ax.set_ylabel("Centroid frequency / MHz")
            self.ax.set_title(f"{element}-{mass}: Centroid frequency vs scan")
            self.ax.grid(True)
            self.ax.legend()
            self.fig.tight_layout()
            self.canvas.draw()
    
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load or plot data:\n{e}")
            self.ax.clear()
            self.ax.set_title("Error plotting")
            self.canvas.draw()
