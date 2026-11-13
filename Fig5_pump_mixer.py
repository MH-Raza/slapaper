"""Fig. 5 - peristaltic pump (a flow rate, b volume per cycle), mixer schematic (c), video frames of the mixer (d, placeholders),
mixing index (e).
Data: pump gravimetry (data/0_5_mm_spacing.xlsx, or the result tables); mixing videos (data/Mixing data/) for the full panel e
(line + band), otherwise the mixing curves from the result tables (line only). Panel d: placeholders, see config.FIG5_VIDEO_FRAMES.
Run from anywhere: python Fig5_pump_mixer.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib.pyplot as plt
import config as C
import loaders as L
from figio import save_tiff, report_placeholders

import fig5_panels as F5

S = L.pump()["S"]
res = L.mixing()
fig = F5.composite(res, S, frames=C.FIG5_VIDEO_FRAMES)
report_placeholders(fig)
# save_tiff(fig, os.path.join(C.FIGURES, 'Fig5.tif'))   # save figure (uncomment to export: 600 dpi RGB TIFF)
plt.show()
