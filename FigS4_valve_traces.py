"""Fig. S4 - normalised flow of all 125 valve runs with the open/actuated analysis windows.
Data: valve logs (raw data needed).
Run from anywhere: python FigS4_valve_traces.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import valve_closure_analysis as VAL

if not L._raw('valves'): L._need(C.VALVE_DIR, 'Fig. S4')
fig = VAL.figSV_traces(C.VALVE_DIR, C.VALVE_EXTRA)
# save_tiff(fig, os.path.join(C.FIGURES, 'FigS4.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
