"""Fig. 4 - valve closure vs geometry and pressure; minimum sealing pressure vs as-fabricated aspect ratio.
Data: valve logs (or results/valve_closure_per_run.csv) and the profilometry result table (as-fabricated widths/depths).
Run from anywhere: python Fig4_valve_closure.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import valve_closure_analysis as VAL

v = L.valves()
fig = VAL.fig4(v["G"], v["T"])
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig4.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
