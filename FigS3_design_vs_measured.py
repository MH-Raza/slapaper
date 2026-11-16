"""Fig. S3 - design silhouettes against the measured mean profiles of molds and replicas.
Data: profilometry scans (raw data needed).
Run from anywhere: python FigS3_design_vs_measured.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff

import profilometry_analysis as PRO

d = L.profilometry(need_raw=True)
fig = PRO.figS_cad(d["F"], d["store"])
# save_tiff(fig, os.path.join(C.FIGURES, 'FigS3.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
