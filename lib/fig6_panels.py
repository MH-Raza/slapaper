"""Fig. 6 - DNA nanodevice activation on the chip.
 a  operating principle of the FRET nanodevice: inactive (high FRET) -> + sRNA mimic -> active (low FRET). 3D renders of the
    DNA origami nanodevice (assets/fig6a_nanodevice_*.png, config.FIG6_RENDERS), adapted from Zadegan & Hughes (ACS, 2017);
    labels are drawn here. Without the render files the code draws placeholders of the same size.
 b  chip layout (top) and photograph of the origami chip (bottom; assets/fig6b_chip_photo.jpg, config.FIG6_PHOTO)
 c  FRET efficiency on chip and in tube, log2 fold change relative to origami alone (config.FIG6C_UNITS)
Layout (17.1 cm wide): three columns a | b | c. The top edges of the pictures of a and b and of the plot area of c line up,
and so do their bottom edges.
Style: lib/ch2_style.py (same as Figs 1-5)."""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
from PIL import Image
from ch2_style import *
from figio import placeholder

ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets')
LAYOUT = os.path.join(ASSETS, 'fig6b_layout.png')                    # chip layout as in the manuscript
PHOTO_ASPECT = 600 / 434                                             # width/height of the photograph box (as in the manuscript)

COL_CHIP, COL_TUBE = '#298c8c', '#d9c077'                            # flow-layer teal, control-layer gold
CATS = ['Origami + sRNA mimic', 'Origami', 'Origami + random sequence']
CAT_LABELS = ['Origami +\nsRNA mimic', 'Origami', 'Origami +\nrandom sequence']
SHOW_REPLICATES = True                                               # the three replicates over each bar (n = 3)
CM = 1 / 2.54
FS_ANN = BASE_FS - 1                                                 # annotation text in panel a


def _rgb(path):
    """RGB on white (the renders carry an alpha channel)."""
    im = Image.open(path).convert('RGBA')
    return Image.alpha_composite(Image.new('RGBA', im.size, (255, 255, 255, 255)), im).convert('RGB')


def _image(ax, im, crop=None):
    if crop: im = im.crop(crop)
    ax.imshow(np.asarray(im), interpolation='lanczos'); ax.axis('off')
    return im.size


def _ann(ax, text, xy, xytext, ha='left', va='center', color=None):
    ax.annotate(text, xy=xy, xytext=xytext, textcoords='data', fontsize=FS_ANN, color=color or C_TEXT, ha=ha, va=va,
                linespacing=1.05, annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=C_TEXT, lw=0.5, shrinkA=1.5, shrinkB=0))


# ------------------------------------------------------------------------------------------------------------ panel a
# Pixels of the original renders (1672 x 1022 and 1770 x 1030). CROP: (left, top, right, bottom); the left crop is the one of
# the source slide. MASK: the targeting ligand (black chevron) at the right end of the origami is not part of this device and
# is painted white. SCALE_BAR_PX: 10 nm in the active-state render (from the scale bar of the source slide).
CROP = {'inactive': (538, 20, 1672, 800), 'active': (776, 95, 1770, 990)}
MASK = {'inactive': (1541, 262, 1672, 382), 'active': (1644, 495, 1770, 614)}
SCALE_BAR_PX, SCALE_BAR_NM = 202.5, 10
ARROW_H = 0.42                                                       # height of the arrow between the two states (in)
TITLE_H = 0.17                                                       # space for a state title above its picture (in)


def _render(path, state):
    im = _rgb(path)
    l, t, r, b = MASK[state]
    im.paste((255, 255, 255), (l, t, r, b))
    return im.crop(CROP[state])


def _aspect(state):
    l, t, r, b = CROP[state]
    return (r - l) / (b - t)


def a_height(w):
    """Height (in) of panel a (both pictures, the arrow and the second title) when the pictures are `w` inches wide."""
    return w / _aspect('inactive') + ARROW_H + TITLE_H + w / _aspect('active')


def draw_a(fig, x0, y0, w, renders):
    """Inactive state above, active state below, arrow '+ sRNA mimic' between them. Pictures `w` inches wide, the bottom of
    the lower picture at y0 (figure inches). renders: (inactive, active) image files; missing files give placeholders."""
    W, H = fig.get_size_inches()
    h_in, h_act = w / _aspect('inactive'), w / _aspect('active')
    y_act = y0; y_in = y0 + h_act + TITLE_H + ARROW_H
    a1 = fig.add_axes([x0 / W, y_in / H, w / W, h_in / H])
    a2 = fig.add_axes([x0 / W, y_act / H, w / W, h_act / H])
    paths = dict(zip(('inactive', 'active'), renders or (None, None)))
    for ax, st in ((a1, 'inactive'), (a2, 'active')):
        p = paths[st]
        if p and os.path.exists(p):
            _image(ax, _render(p, st))
            ax._drawn = True
        else:
            if p: print(f'Fig. 6a: render not found ({p}); drawing a placeholder.')
            placeholder(ax, f'Fig. 6a {st}-state render', text=f'{st.capitalize()}-state render')
            ax._drawn = False
    for ax, t in ((a1, 'Inactive: high FRET'), (a2, 'Active: low FRET')):
        ax.text(0.5, 1.03, t, transform=ax.transAxes, ha='center', va='bottom', fontsize=BASE_FS, fontweight='bold', color=C_TEXT)

    if a1._drawn:                                                     # labels in render pixels (after cropping)
        ox, oy = CROP['inactive'][:2]
        _ann(a1, 'Displaceable strand', (642 - ox, 400 - oy), (548 - ox, 742 - oy), ha='left', va='center')
        _ann(a1, 'DNA origami', (1380 - ox, 540 - oy), (1545 - ox, 742 - oy), ha='center', va='center')
        _ann(a1, 'Cy5*', (1075 - ox, 395 - oy), (1560 - ox, 590 - oy), ha='center', va='center')
    if a2._drawn:
        ox, oy = CROP['active'][:2]
        _ann(a2, 'Cy3*', (1128 - ox, 238 - oy), (936 - ox, 150 - oy), ha='center', va='center')
        _ann(a2, 'Displaced strand', (1480 - ox, 215 - oy), (1536 - ox, 120 - oy), ha='center', va='center')
        _ann(a2, 'Cy5', (1286 - ox, 662 - oy), (1460 - ox, 950 - oy), ha='left', va='center')
        xb, yb = 800 - ox, 955 - oy                                   # scale bar, bottom left
        a2.plot([xb, xb + SCALE_BAR_PX], [yb, yb], color=C_TEXT, lw=1.2, solid_capstyle='butt', clip_on=False)
        a2.text(xb + SCALE_BAR_PX / 2, yb - 22, f'{SCALE_BAR_NM} nm', ha='center', va='bottom', fontsize=FS_ANN, color=C_TEXT)

    ar = fig.add_axes([x0 / W, (y_act + h_act + TITLE_H) / H, w / W, ARROW_H / H]); ar.axis('off'); ar.set_xlim(0, 1); ar.set_ylim(0, 1)
    ar.annotate('', xy=(0.5, 0.04), xytext=(0.5, 0.96), arrowprops=dict(arrowstyle='-|>', color=C_TEXT, lw=0.8, mutation_scale=7))
    ar.text(0.55, 0.5, '+ sRNA mimic', ha='left', va='center', fontsize=FS_ANN, color=C_TEXT)
    return a1, a2


# ------------------------------------------------------------------------------------------------------------ panel c
def _mathregular(v):
    m, e = f'{v:.0e}'.split('e')
    return rf'$\mathregular{{{int(m)}\times10^{{{int(e)}}}}}$'


C_UNITS = 'log2'           # default (as in the manuscript); the figure scripts pass config.FIG6C_UNITS
C_YLIM = {'log10': (-1.6, -0.8), 'log_axis': (1e-2, 2e-1), 'linear': (0, 0.12), 'dB': (-4.5, 1.5), 'log2': (-1.5, 0.5)}
_POW = {0: '0', 2: r'$\mathregular{2\times10^{-2}}$', 4: r'$\mathregular{4\times10^{-2}}$', 6: r'$\mathregular{6\times10^{-2}}$',
        8: r'$\mathregular{8\times10^{-2}}$', 10: r'$\mathregular{1\times10^{-1}}$', 12: r'$\mathregular{1.2\times10^{-1}}$'}


def _minus(v, fmt='{:g}'):
    return fmt.format(v).replace('-', '\u2212')


def draw_c(ax, df, units=None):
    """FRET efficiency E on chip and in tube (n = 3; bars mean ± SD, points single replicates). units:
    'log10'    y = log10(E); mean ± SD of the log values; bars start at the bottom of the axis
    'log_axis' E on a logarithmic axis; geometric mean with the SD of log10(E) (asymmetric error bars)
    'linear'   E on a linear axis from 0, ticks written as m x 10^e
    'dB'       10 log10(E / E_origami) of the same setting (chip or tube); 0 dB = origami alone
    'log2'     log2(E / E_origami) of the same setting; -1 = half of origami alone"""
    units = units or C_UNITS
    if units not in C_YLIM:
        raise ValueError(f"Fig. 6c units must be one of {list(C_YLIM)}, not {units!r}")
    x = np.arange(len(CATS)); w = 0.36
    lo, hi = C_YLIM[units]
    for dx, setting, col, lab in ((-w / 2, 'Chip', COL_CHIP, 'On chip'), (w / 2, 'Tube', COL_TUBE, 'In tube')):
        d = df[df.setting == setting].copy()
        ref = d.loc[d.condition == 'Origami', 'fret'].mean()
        d['y'] = {'log10': np.log10(d['fret']), 'log_axis': np.log10(d['fret']), 'linear': d['fret'],
                  'dB': 10 * np.log10(d['fret'] / ref), 'log2': np.log2(d['fret'] / ref)}[units]
        g = d.groupby('condition')['y']
        m, sd = g.mean().reindex(CATS).values, g.std(ddof=1).reindex(CATS).values
        if units == 'log_axis':                                       # back-transform: geometric mean, asymmetric SD
            m, err = 10 ** m, [10 ** m - 10 ** (m - sd), 10 ** (m + sd) - 10 ** m]
        else:
            err = sd
        base = lo if units in ('log10', 'log_axis') else 0
        ax.bar(x + dx, m - base, width=w, bottom=base, color=col, edgecolor='none', zorder=2, label=lab)
        ax.errorbar(x + dx, m, yerr=err, fmt='none', ecolor=C_TEXT, elinewidth=0.7, capsize=2.0, capthick=0.7, zorder=4)
        if SHOW_REPLICATES:
            for i, c in enumerate(CATS):
                v = d[d.condition == c].sort_values('replicate')['y'].values
                if units == 'log_axis': v = 10 ** v
                ax.scatter(x[i] + dx + np.linspace(-0.09, 0.09, len(v)), v, s=5, facecolor='white', edgecolor=C_TEXT,
                           linewidth=0.45, zorder=5)
    if units == 'log_axis':
        ax.set_yscale('log')
        major = [1e-2, 2e-2, 5e-2, 1e-1, 2e-1]
        ax.yaxis.set_major_locator(mt.FixedLocator(major))
        ax.yaxis.set_major_formatter(mt.FixedFormatter([_mathregular(v) for v in major]))
        ax.yaxis.set_minor_locator(mt.LogLocator(base=10, subs=np.arange(2, 10)))
        ax.yaxis.set_minor_formatter(mt.NullFormatter())
    elif units == 'linear':
        ax.yaxis.set_major_locator(mt.MultipleLocator(0.02))
        ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: _POW.get(int(round(v * 100)), '')))
    else:
        ax.yaxis.set_major_locator(mt.MultipleLocator({'log10': 0.2, 'dB': 1, 'log2': 0.5}[units]))
        ax.yaxis.set_major_formatter(mt.FuncFormatter(lambda v, _: _minus(v, '{:.1f}' if units == 'log10' else '{:g}')))
        if units in ('dB', 'log2'):
            ax.axhline(0, color=C_TEXT, lw=0.6, zorder=3)
    ax.set_ylim(lo, hi)
    ax.set_ylabel({'log10': r'log$_{\mathregular{10}}$(FRET efficiency)', 'log_axis': 'FRET efficiency', 'linear': 'FRET efficiency',
                   'dB': 'FRET change (dB)', 'log2': r'log$_{\mathregular{2}}$ fold change'}[units])
    ax.set_xticks(x); ax.set_xticklabels(CAT_LABELS, linespacing=1.0); ax.tick_params(axis='x', length=0)
    ax.set_xlim(-0.6, len(CATS) - 0.4)
    ax.grid(axis='y', which='major', ls=':', lw=0.4, color='#CCCCCC', alpha=0.8, zorder=0); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc='upper left', fontsize=BASE_FS - 1, handlelength=1.0, handleheight=0.8, borderaxespad=0.1,
              ncol=2, columnspacing=1.0)


# ------------------------------------------------------------------------------------------------------------ panel b
def _photo_axes(ax, photo):
    if photo and not os.path.exists(photo):
        print(f'Fig. 6b: photograph not found ({photo}); drawing a placeholder.'); photo = None
    if photo:
        im = Image.open(photo).convert('RGB'); a = im.size[0] / im.size[1]
        if abs(a - PHOTO_ASPECT) > 0.0005:                             # centre-crop to the box so the layout never changes
            if a > PHOTO_ASPECT: d = int(round((im.size[0] - PHOTO_ASPECT * im.size[1]) / 2)); im = im.crop((d, 0, im.size[0] - d, im.size[1]))
            else: d = int(round((im.size[1] - im.size[0] / PHOTO_ASPECT) / 2)); im = im.crop((0, d, im.size[0], im.size[1] - d))
        _image(ax, im)
    else:
        placeholder(ax, 'Fig. 6b chip photograph', text='Chip photograph')


def _letters(fig, pos):
    W, H = fig.get_size_inches()
    for t, (x, y) in zip('abc', pos):
        fig.text(x / W, y / H, t, fontsize=10, fontweight='bold', color=C_TEXT, ha='left', va='top')


# ------------------------------------------------------------------------------------------------------------ figure
W_FIG = 17.1 * CM                                                    # double-column width (in)
X_A, GAP_AB, GAP_B_LAYOUT_PHOTO = 0.20, 0.30, 0.08
W_C, X_C_LABELS = 1.80, 0.62                                         # plot width; room for the y tick labels and label of c
Y_BOT = 0.36                                                         # room for the category labels of c


def figure6(df, photo=None, renders=None, c_units=None):
    """df: FRET replicates from loaders.fret(); photo: chip photograph (None = placeholder);
    renders: (inactive, active) nanodevice render files (None or missing files = placeholders); c_units: see draw_c."""
    lay = Image.open(LAYOUT).convert('RGB'); r_lay = lay.size[0] / lay.size[1]
    k_b = 1 / r_lay + 1 / PHOTO_ASPECT                                # column b height = k_b * width + gap
    k_a = 1 / _aspect('inactive') + 1 / _aspect('active')             # column a height = k_a * width + arrow + title
    # widths: X_A + w_a + GAP_AB + w_b + X_C_LABELS + W_C + 0.03 = W_FIG, and both columns equally high
    rest = W_FIG - X_A - GAP_AB - X_C_LABELS - W_C - 0.03
    w_b = (rest + (ARROW_H + TITLE_H - GAP_B_LAYOUT_PHOTO) / k_a) / (1 + k_b / k_a)
    w_a = rest - w_b
    h = k_b * w_b + GAP_B_LAYOUT_PHOTO
    top = Y_BOT + h
    H = top + TITLE_H + 0.06
    fig = plt.figure(figsize=(W_FIG, H), facecolor='white')
    def ax_in(l, b, w, hh): return fig.add_axes([l / W_FIG, b / H, w / W_FIG, hh / H])

    draw_a(fig, X_A, Y_BOT, w_a, renders)
    x_b = X_A + w_a + GAP_AB
    h_lay, h_ph = w_b / r_lay, w_b / PHOTO_ASPECT
    _image(ax_in(x_b, top - h_lay, w_b, h_lay), lay)
    _photo_axes(ax_in(x_b, Y_BOT, w_b, h_ph), photo)
    x_c = x_b + w_b + X_C_LABELS
    c = ax_in(x_c, Y_BOT, W_C, h); draw_c(c, df, c_units)
    y_let = H - 0.02
    _letters(fig, ((0.02, y_let), (x_b - 0.16, y_let), (x_c - X_C_LABELS + 0.02, y_let)))
    return fig
