"""
fig5_panels.py - Figure 5 (pump + mixer) as ONE tight, aligned composite and as INDIVIDUAL panels (a, b, d, e).

Every panel is drawn by the same function whether it goes into the composite or is saved alone.
Layout (composite): a, b and e have the same width and height; a and c share top and bottom edges; b, d and e share top and bottom
edges (d = the 3 x 4 frame block, with its time labels below like tick labels). Panel c is padded with its own grey to span the right-hand width.

Frequencies are the PROGRAMMED phase frequencies, 1000 / ESI delay (ms): 25 ms -> 40.0 Hz, 33 -> 30.3, 48 -> 20.8, 83 -> 12.0, 167 -> 6.0.
One cycle = 6 phases (cycle frequency = phase frequency / 6).  The delivered rate was lower (every phase lasted 17-25 ms longer than
programmed); flow rates therefore use the measured elapsed time of each run (pump_analysis.py).
Panel c is the device schematic (assets/panel_c_schematic.png, PANEL_C).
Panel d (video frames of the chip) is a grid of placeholders by default, each box the size of one frame in the manuscript
(width/height = FRAME_ASPECT); composite(..., frames=True) draws the frames from the videos instead (needs the mixing videos).
Panel e: with the videos, line = primary analysis and band = min-max over the 48 analysis settings; from the result tables
(no videos) only the primary line can be drawn.
SAVING: the savefig lines at the bottom are commented out - uncomment what you need.
"""
import os, sys, numpy as np, cv2
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch2_style import *
from figio import placeholder
from matplotlib.lines import Line2D
import mixing_analysis as MA
import pump_analysis as PA

PANEL_C = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets', 'panel_c_schematic.png')   # schematic as in the manuscript (with t = 7.7 s)
# Scale bar of panel d. The displayed crop of each video frame spans FRAME_WIDTH_MM. Calibrated against the layout in panel c (1 mm bar):
# inlet centre-to-centre spacing 1.35 mm -> 54.6 px/mm and left-channel-to-inlet distance 1.64 mm -> 59.3 px/mm in the 600 dpi composite,
# mean 57.0 px/mm over a 298 px frame -> 5.23 mm. Replace with the microscope calibration (objective + camera pixel size) when available.
FRAME_WIDTH_MM = 5.23
SCALE_BAR_MM = 1.0
FRAME_ASPECT = 298 / 293          # width/height of one frame of panel d in the manuscript (crop of the aligned videos)
TIMES = (0, 5, 15, 25)            # frames of panel d, seconds after sealing
XMAX_E = 30                       # panel e: first 30 s (whole record: Fig. S11)
PS = 9


def hz_label(cond): return MA.hz(cond)


# ------------------------------------------------------------------------------------------------ panels
def _xaxis(ax, ticks):
    ax.set_xticks(ticks); ax.set_xticklabels([f'{v:.1f}' for v in ticks]); ax.set_xlim(ticks[0] - 5, ticks[-1] + 5)
    ax.set_xlabel('Phase frequency (Hz)'); ax.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', alpha=0.6)


def draw_a(ax, S):
    for bar in (0.5, 1.0, 1.5):
        s = S[S.bar == bar].sort_values('delay_ms'); ok = s.flow.notna()
        ax.errorbar(s.f_prog[ok], s.flow[ok], yerr=s.flow_sd[ok], color=PA.PCOL[bar], marker='o', ms=3.4, lw=1.5, capsize=2, elinewidth=0.7)
    hs = [Line2D([], [], color=PA.PCOL[b], marker='o', ms=3.4, lw=1.5, label=f'{b:g} bar') for b in (0.5, 1.0, 1.5)]
    ax.set_ylim(0, 5); ax.set_ylabel('Flow rate (µL min$^{-1}$)'); _xaxis(ax, sorted(S.f_prog.unique()))
    ax.legend(handles=hs, title='Actuation pressure', title_fontsize=BASE_FS - 1.3, frameon=False, fontsize=BASE_FS - 1.3, loc='upper left', handlelength=1.5, labelspacing=0.3, borderaxespad=0.2)


def draw_b(ax, S):
    for bar in (0.5, 1.0, 1.5):
        s = S[S.bar == bar].sort_values('delay_ms')
        ax.errorbar(s.f_prog, s.stroke, yerr=s.stroke_sd, color=PA.PCOL[bar], marker='o', ms=3.4, lw=1.5, capsize=2, elinewidth=0.7)
    ax.set_ylim(12, 22); ax.set_ylabel('Volume per cycle (nL)'); _xaxis(ax, sorted(S.f_prog.unique()))
    hs = [Line2D([], [], color=PA.PCOL[b], marker='o', ms=3.4, lw=1.5, label=f'{b:g} bar') for b in (0.5, 1.0, 1.5)]
    ax.legend(handles=hs, title='Actuation pressure', title_fontsize=BASE_FS - 1.3, frameon=False, fontsize=BASE_FS - 1.3, loc='lower left', handlelength=1.5, labelspacing=0.3, borderaxespad=0.2)


def common_crop(res):
    return min(r['G']['xr'] + 30 for r in res), min(r['G']['h'] - 40 for r in res)


def draw_d(fig, gs, res):
    """3 rows x 4 frames on a prepared GridSpec `gs` (3 x len(TIMES))."""
    cw, ch = common_crop(res); axes = []
    for i, r in enumerate(res):
        c = MA.CONDS[r['cond']]; cap = cv2.VideoCapture(os.path.join(MA.DATA_DIR, c['folder'], c['video']))
        for j, T in enumerate(TIMES):
            ax = fig.add_subplot(gs[i, j]); axes.append(ax)
            cap.set(cv2.CAP_PROP_POS_FRAMES, int(round((r['t_seal'] + T) * r['fps']))); ok, f = cap.read()
            ax.imshow(cv2.cvtColor(f[:ch, :cw], cv2.COLOR_BGR2RGB)); ax.axis('off')
            if i == 2: ax.text(0.5, -0.04, f't = {T} s', transform=ax.transAxes, ha='center', va='top', fontsize=BASE_FS - 0.5)
            if i == 2 and j == 0:                                   # 1 mm scale bar, lower left of the first frame of the last row
                L = cw * SCALE_BAR_MM / FRAME_WIDTH_MM; x0, y0 = 0.15 * cw, 0.86 * ch
                ax.plot([x0, x0 + L], [y0, y0], color='#111111', lw=1.6, solid_capstyle='butt')
                ax.text(x0 + L / 2, y0 - 0.03 * ch, f'{SCALE_BAR_MM:g} mm', ha='center', va='bottom', fontsize=BASE_FS - 2, fontweight='bold', color='#111111')
            if j == 0: ax.text(-0.05, 0.5, hz_label(r['cond']), transform=ax.transAxes, rotation=90, ha='right', va='center', color=MA.COLORS[r['cond']], fontweight='bold', fontsize=BASE_FS)
        cap.release()
    return axes


def _scale_bar_axes(ax):
    """1 mm scale bar in axes fractions, same place as on the frames (assumes the frame spans FRAME_WIDTH_MM)."""
    L = SCALE_BAR_MM / FRAME_WIDTH_MM; x0, y0 = 0.15, 0.14
    ax.plot([x0, x0 + L], [y0, y0], color='#111111', lw=1.6, solid_capstyle='butt', transform=ax.transAxes)
    ax.text(x0 + L / 2, y0 + 0.03, f'{SCALE_BAR_MM:g} mm', transform=ax.transAxes, ha='center', va='bottom', fontsize=BASE_FS - 2,
            fontweight='bold', color='#111111')


def draw_d_placeholder(fig, gs, conds):
    """Same grid, labels and scale bar as draw_d, with an empty box per frame (insert the video frames yourself)."""
    axes = []
    for i, cond in enumerate(conds):
        for j, T in enumerate(TIMES):
            ax = fig.add_subplot(gs[i, j]); axes.append(ax)
            placeholder(ax, f'Fig. 5d frame {hz_label(cond)}, t = {T} s', text='', show_size=False, aspect=FRAME_ASPECT)
            if i == len(conds) - 1: ax.text(0.5, -0.04, f't = {T} s', transform=ax.transAxes, ha='center', va='top', fontsize=BASE_FS - 0.5)
            if i == len(conds) - 1 and j == 0: _scale_bar_axes(ax)
            if j == 0: ax.text(-0.05, 0.5, hz_label(cond), transform=ax.transAxes, rotation=90, ha='right', va='center', color=MA.COLORS[cond], fontweight='bold', fontsize=BASE_FS)
    return axes


def _frames_available(res, frames):
    return bool(frames) and all('G' in r for r in res)


def draw_e(ax, res, xmax=XMAX_E):
    MA.draw_curves(ax, res, xmax); ax.set_xticks(np.arange(0, xmax + 1, 10 if xmax <= 30 else 20)); MA.curve_legend(ax, res)


def _panel_c_rgb():
    """PANEL_C as RGB with its transparency resolved. cv2.imread drops the alpha channel, so transparent pixels came out with the
    colour stored under them (black in a PNG exported from PowerPoint); here they take the schematic's own background colour."""
    from PIL import Image
    a = np.asarray(Image.open(PANEL_C).convert('RGBA')).astype(np.float32)
    h, w = a.shape[:2]
    opaque = np.where(a[h // 2, :, 3] == 255)[0]                     # background = opaque pixel near the left edge, middle row
    bg = a[h // 2, opaque[opaque >= 3][0] if (opaque >= 3).any() else (opaque[0] if opaque.size else 3), :3]
    al = a[..., 3:4] / 255.0
    return np.round(a[..., :3] * al + bg * (1 - al)).astype(np.uint8), np.round(bg).astype(np.uint8)


def draw_c(ax, aspect=None):
    """Schematic; if `aspect` (width/height of the target box) is wider than the picture, extend the grey background left and right
    (the drawing keeps its size and stays centred)."""
    im, bg = _panel_c_rgb(); h, w = im.shape[:2]
    if aspect and aspect * h > w:
        pad = int(round((aspect * h - w) / 2)); side = np.empty((h, pad, 3), np.uint8); side[:] = bg
        im = np.concatenate([side, im, side], axis=1)
    ax.imshow(im); ax.axis('off')


# ------------------------------------------------------------------------------------------------ figures
def prepare():
    res = [MA.analyse(c) for c in MA.CONDS]; rows = [MA.summarise(r) for r in res]
    for r, row in zip(res, rows): r['t50_report'] = row['t50']
    return res, PA.summarise(PA.load())


def _ax(fig, x, y, w, h):
    W, H = fig.get_size_inches(); return fig.add_axes([x / W, y / H, w / W, h / H])


def _letter(fig, s, x, y):
    W, H = fig.get_size_inches(); fig.text(x / W, y / H, s, fontsize=PS + 1, fontweight='bold', va='bottom', ha='left', color=C_TEXT)


def composite(res, S, frames=False):
    """Layout in inches. a, b and e have the same width and height (WP x HP); a/c share top and bottom edges; b/d/e share top and bottom edges.
    frames=False: panel d = placeholders; frames=True: video frames (needs the videos)."""
    W = 7.2; xl, WP, HP = 0.52, 1.86, 1.50
    y0 = 0.50; y1 = y0 + HP; ya = y1 + 0.60; H = ya + HP + 0.14
    xr = xl + WP + 0.18; wr = W - 0.04 - xr                    # right-hand region
    fig = plt.figure(figsize=(W, H))
    axa, axb = _ax(fig, xl, ya, WP, HP), _ax(fig, xl, y0, WP, HP); axc = _ax(fig, xr, ya, wr, HP)
    draw_a(axa, S); draw_b(axb, S); draw_c(axc, aspect=wr / HP)
    real = _frames_available(res, frames)
    if frames and not real: print('Fig. 5d: the mixing videos are needed for frames=True; drawing placeholders.')
    cw, ch = common_crop(res) if real else (FRAME_ASPECT, 1.0); fw = HP / 3 * cw / ch            # frame width for a block exactly HP high
    gap = 0.012; bw = 4 * fw + 3 * gap * fw; dx0 = xr + 0.22
    gsd = gridspec.GridSpec(3, len(TIMES), figure=fig, left=dx0 / W, right=(dx0 + bw) / W, bottom=y0 / H, top=y1 / H, wspace=gap, hspace=gap * 3)
    draw_d(fig, gsd, res) if real else draw_d_placeholder(fig, gsd, [r['cond'] for r in res])
    ex0 = W - 0.04 - WP; axe = _ax(fig, ex0, y0, WP, HP); draw_e(axe, res)           # e: same size as a and b, flush right
    top = lambda y: y + HP + 0.06
    _letter(fig, 'a', 0.03, top(ya)); _letter(fig, 'b', 0.03, top(y0)); _letter(fig, 'c', xr, top(ya))
    _letter(fig, 'd', xr, top(y0)); _letter(fig, 'e', ex0 - 0.50, top(y0))
    print(f'composite {W:.2f} x {H:.2f} in; d block {bw:.2f} in wide, gap d-e axes {ex0 - dx0 - bw:.2f} in')
    return fig


def panels(res, S, frames=False):
    """Individual panels with exactly the same axes size as in the composite (a, b, e: 1.86 x 1.50 in)."""
    out = {}; WP, HP = 1.86, 1.50; L, B, R, T = 0.52, 0.50, 0.06, 0.08
    def one(draw, top_extra=0.0):
        f = plt.figure(figsize=(L + WP + R, B + HP + T + top_extra)); ax = f.add_axes([L / f.get_size_inches()[0], B / f.get_size_inches()[1], WP / f.get_size_inches()[0], HP / f.get_size_inches()[1]]); draw(ax); return f
    out['a_pump_flow'] = one(lambda ax: draw_a(ax, S)); out['b_pump_stroke'] = one(lambda ax: draw_b(ax, S)); out['e_mixing_index'] = one(lambda ax: draw_e(ax, res), 0.22)
    real = _frames_available(res, frames)
    cw, ch = common_crop(res) if real else (FRAME_ASPECT, 1.0); fw = HP / 3 * cw / ch; bw = 4 * fw * 1.036; Wd, Hd = bw + 0.30, HP + 0.30
    f = plt.figure(figsize=(Wd, Hd)); gs = gridspec.GridSpec(3, 4, figure=f, left=0.24 / Wd, right=(0.24 + bw) / Wd, bottom=0.26 / Hd, top=(0.26 + HP) / Hd, wspace=0.012, hspace=0.036)
    draw_d(f, gs, res) if real else draw_d_placeholder(f, gs, [r['cond'] for r in res]); out['d_timelapse'] = f
    return {k: out[k] for k in ('a_pump_flow', 'b_pump_stroke', 'd_timelapse', 'e_mixing_index')}


if __name__ == '__main__':
    res, S = prepare(); fig_all = composite(res, S); fig_panels = panels(res, S)
    # ---- SAVE (uncomment what you need)
    # fig_all.savefig('Fig5_pump_and_mixer.jpeg', dpi=600, bbox_inches='tight', pad_inches=0.02)
    # for name, f in fig_panels.items(): f.savefig(f'Fig5{name}.jpeg', dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.show()
