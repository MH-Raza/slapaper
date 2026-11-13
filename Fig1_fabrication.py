"""Fig. 1 - fabrication workflow for multilayer push-up valve chips from laser-SLA molds.
Schematic only, no data. The drawing code is in lib/fig1_fabrication.py.
Run from anywhere: python Fig1_fabrication.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
from figio import save_tiff

import fig1_fabrication as F1

fig = F1.make_figure()
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig1.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
