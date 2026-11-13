"""Fig. S6 - closure against design and as-fabricated aspect ratio; separation for four width definitions.
Data: valve result table + profilometry result table (raw data optional).
Run from anywhere: python FigS6_aspect_ratio.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import valve_closure_analysis as VAL

v = L.valves()
fig = VAL.figSV_ar(v["G"], v["T"], v["prof_path"])
# save_tiff(fig, os.path.join(C.FIGURES, 'FigS6.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
