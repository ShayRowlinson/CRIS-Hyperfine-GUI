# -*- coding: utf-8 -*-
"""
Created on Mon Jul 7 12:20:31 2025

main_v1.py

Version 1

launcher of the interactive GUI for CRIS scan binning and hyperfine analysis.
Requires binning_v1.py, fitting_v1.py, centroid_plot_v1.py, extract_moments_v1.py and Dopplershift_analysis.py, interpolate_ISCOOL.py

Author: Shay Rowlinson
"""

import tkinter as tk
from tkinter import ttk
from binning_v1 import Binning
from fitting_v1 import Fitting
from centroid_plot_v1 import CentroidPlot
from extract_moments_v1 import Moments

class MainGUI(tk.Tk):
    """
    Main GUI for CRIS Analysis with tabbed interface.
    Tabs:
        Binning, Fitting, Centroid plot and Nuclear moment extraction.
    """
    def __init__(self):
        super().__init__()
        self.title("Graphical interface for CRIS scan binning and hyperfine analysis")
        self.geometry("1920x1080")

        # Canvas and scrollbars
        self.canvas = tk.Canvas(self, borderwidth=0)
        self.vscroll = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.hscroll = ttk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=self.vscroll.set, xscrollcommand=self.hscroll.set)
        self.vscroll.pack(side="right", fill="y")
        self.hscroll.pack(side="bottom", fill="x")
        self.canvas.pack(side="left", fill="both", expand=True)

        # Add frame inside canvas
        self.frame = ttk.Frame(self.canvas)
        self.frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.create_window((0,0), window=self.frame, anchor="nw")

        # Notebook inside frame
        self.notebook = ttk.Notebook(self.frame)
        self.notebook.pack(fill='both', expand=True)
        self.binning_tab = Binning(self.notebook)
        self.fitting_tab = Fitting(self.notebook)
        self.centroid_tab = CentroidPlot(self.notebook)
        self.moments_tab = Moments(self.notebook)
        self.notebook.add(self.binning_tab, text="Binning")
        self.notebook.add(self.fitting_tab, text="Fitting")
        self.notebook.add(self.centroid_tab, text="Centroid Plot")
        self.notebook.add(self.moments_tab, text="Extract Moments")

        # Bind mousewheel to scroll
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")


MainGUI().mainloop()