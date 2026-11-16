"""Fig. 2 - push-up valve regimes (schematic) and representative mold/replica profiles of C1 and F3.
Data: profilometry scans (raw data needed: panels b-c draw the registered scan lines).
Run from anywhere: python Fig2_regimes_profiles.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import fig2_regimes_profiles as F2

d = L.profilometry(need_raw=True)
fig = F2.fig2(d["R"], d["F"], d["store"])
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig2.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
