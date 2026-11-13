"""Fig. S2 - every flow-layer mold and PDMS replica scan (F1-F5).
Data: profilometry scans (raw data needed).
Run from anywhere: python FigS2_profiles_flow.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import profilometry_analysis as PRO

d = L.profilometry(need_raw=True)
fig = PRO.figS_individual(d["R"], d["store"], "F")
# save_tiff(fig, os.path.join(C.FIGURES, 'FigS2.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
