"""Shared figure style for all figures (Arial, bold axis labels,
gold = control layer, teal = flow layer). Import this module in every plotting script."""
import numpy as np
import matplotlib as mpl
from matplotlib import font_manager

# ---- colours (identical to Figure_1_Rev_V2.ipynb / Revised.ipynb)
C_CTRL, C_FLOW, C_TEXT = '#C9A84C', '#2A8A8A', '#1A2535'
C_MOLD, C_PDMS_BG, C_BRIDGE, C_NORMAL = '#7F7F7F', '#E2E2E2', '#C0392B', '#1A1A1A'
C_GOLD_D, C_TEAL_D = '#8C6D1F', '#1E6868'
PRESSURE_COLORS = {500: '#D4B96A', 750: '#C9A84C', 1000: '#89A874', 1250: '#3A9A96', 1500: '#1E6868'}
PRESSURE_MARKERS = {500: 'o', 750: 's', 1000: '^', 1250: 'D', 1500: 'v'}
PRESSURES = [500, 750, 1000, 1250, 1500]

BASE_FS = 7
MM = 1 / 25.4                                   # mm -> inch
_avail = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ('Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans') if f in _avail), 'DejaVu Sans')

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': [FONT, 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': BASE_FS, 'axes.titlesize': BASE_FS, 'axes.labelsize': BASE_FS, 'axes.labelweight': 'bold',
    'xtick.labelsize': BASE_FS - 1, 'ytick.labelsize': BASE_FS - 1, 'legend.fontsize': BASE_FS - 1.2,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5,
    'ytick.major.size': 2.5, 'xtick.direction': 'out', 'ytick.direction': 'out', 'xtick.top': False, 'ytick.right': False,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'text.color': C_TEXT, 'axes.labelcolor': C_TEXT, 'xtick.color': C_TEXT, 'ytick.color': C_TEXT, 'axes.edgecolor': C_TEXT,
    # maths in the same face as the text (the earlier notebooks rendered maths in the default serif face)
    'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold',
    'mathtext.default': 'it', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
})


def panel_label(fig_or_ax, text, x=None, y=None, ax=None):
    """Bold panel letter. With an Axes: axes coordinates (default just above the top-left corner)."""
    a = fig_or_ax
    a.text(-0.16 if x is None else x, 1.04 if y is None else y, text, transform=a.transAxes, fontsize=10, fontweight='bold',
           va='bottom', ha='left', color=C_TEXT)


def strip_axes(fig, spec):
    """Axes reserved for a panel title + legend + statistics line, placed above the plotting axes (no data in it)."""
    ax = fig.add_subplot(spec); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    return ax
