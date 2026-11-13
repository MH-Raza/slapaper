"""Fig. S5 - closure against pressure for all 25 geometries (single cycles and their mean).
Data: valve logs (or results/valve_closure_per_cycle.csv).
Run from anywhere: python FigS5_closure_vs_pressure.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import valve_closure_analysis as VAL

v = L.valves()
fig = VAL.figSV_cycles(v["D"])
# save_tiff(fig, os.path.join(C.FIGURES, 'FigS5.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
