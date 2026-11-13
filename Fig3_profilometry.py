"""Fig. 3 - dimensional fidelity: footprint offsets of control and flow features, height transfer mold -> PDMS.
Data: profilometry (raw scans, or results/profilometry_per_file.csv).
Run from anywhere: python Fig3_profilometry.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import profilometry_analysis as PRO

d = L.profilometry()
fig = PRO.fig3(d["R"], d["F"], d["PO"], d["FIT"], d["T_D"], d["T_F"], d["T_L"])
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig3.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
