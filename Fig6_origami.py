"""Fig. 6 - DNA nanodevice activation on the chip: (a) operating principle of the FRET nanodevice, (b) layout and photograph of
the origami chip, (c) FRET efficiency on chip and in tube as log2 fold change relative to origami alone (mean ± SD, n = 3).
Data: data/FRET_Chip_and_tube.xlsx (included), or results_reference/fret_replicates.csv.
Images: assets/fig6a_nanodevice_*.png (config.FIG6_RENDERS; placeholders if missing), assets/fig6b_chip_photo.jpg (config.FIG6_PHOTO).
Run from anywhere: python Fig6_origami.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff, report_placeholders

import fig6_panels as F6

df = L.fret()
fig = F6.figure6(df, photo=C.FIG6_PHOTO, renders=C.FIG6_RENDERS, c_units=C.FIG6C_UNITS)
report_placeholders(fig)
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig6.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
