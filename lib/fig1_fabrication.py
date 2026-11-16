"""Fig. 1 - fabrication workflow for multilayer push-up valve chips from laser-SLA molds (drawing code).

make_figure() returns the figure. It draws with matplotlib defaults + STYLE inside an rc_context, so the result is the
same whether it runs alone, after other figures, or in a Jupyter kernel, and the global style is left untouched.

Panels a-c: mold design, printing, washing and post-curing (cross-sections perpendicular to each channel).
Panels e-g: cross-sections along the control channel, so the control line crosses the flow channel and extends past it
on both sides.  Colours: flow layer teal, control layer gold (shared thesis palette).
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Polygon, Circle, FancyArrowPatch


# ---------------------------------------------------------------- shared style (Arial, Nature-style)
SLATE = '#2C3E50'; TEAL = '#298c8c'; TEAL_D = '#2d635b'; TEAL_T = '#59a89c'
GOLD = '#d9c077'; GOLD_D = '#8c773a'; OVER = '#98996b'; RED = '#E74C3C'
MOLD = '#8a96a0'; PANEL = '#F3F5F6'; GLASS = '#cfd9e0'; VIOLET = '#7B4FA8'
FONT = next((f for f in ('Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans')
             if f in {x.name for x in mpl.font_manager.fontManager.ttflist}), 'DejaVu Sans')
# font.family names the face itself (not the generic 'sans-serif'): a generic family is resolved when the figure is
# drawn, i.e. after the rc_context below has closed, and could then pick up another face
STYLE = {'font.family': [FONT], 'font.sans-serif': [FONT], 'text.color': SLATE,
         'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic',
         'mathtext.bf': FONT + ':bold', 'mathtext.default': 'it',
         'savefig.dpi': 600, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'}



def tint(c, a):
    return tuple(a * np.array(mcolors.to_rgb(c)) + (1 - a) * np.ones(3))


MM = 1 / 25.4
W, H = 180, 112                    # figure size in mm; every panel axis is drawn in mm
FS_T, FS_L, FS_S = 7.2, 6.0, 5.6   # panel title, label, secondary label (pt)



def make_figure():
    with mpl.rc_context():
        mpl.rcdefaults()
        mpl.rcParams.update(STYLE)

        fig = plt.figure(figsize=(W * MM, H * MM))
        bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis('off')


        def panel(x, y, w, h, letter, title):
            ax = fig.add_axes([x / W, y / H, w / W, h / H]); ax.set_xlim(0, w); ax.set_ylim(0, h); ax.axis('off')
            ax.add_patch(Rectangle((0, 0), w, h, fc=PANEL, ec='none', zorder=0))
            ax.text(1.6, h - 1.8, letter, fontsize=9, fontweight='bold', va='top', ha='left')
            ax.text(5.6, h - 2.5, title, fontsize=FS_T, va='top', ha='left')
            return ax


        def rect(ax, x, y, w, h, fc, ec='none', lw=0.5, z=2, ls='-', alpha=1.0):
            ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z, ls=ls, alpha=alpha))


        def semi(ax, xc, y0, r, fc, ec='none', lw=0.5, z=3, ls='-', up=True, n=60):
            t = np.linspace(0, np.pi, n); s = 1 if up else -1
            ax.add_patch(Polygon(np.c_[xc + r * np.cos(t), y0 + s * r * np.sin(t)], closed=True, fc=fc, ec=ec, lw=lw, zorder=z, ls=ls))


        def arrow(ax, p0, p1, c=SLATE, lw=0.8, ms=5, z=6, style='-|>'):
            ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=c, lw=lw, zorder=z,
                                         shrinkA=0, shrinkB=0))


        def txt(ax, x, y, s, fs=FS_L, c=SLATE, ha='center', va='center', **kw):
            ax.text(x, y, s, fontsize=fs, color=c, ha=ha, va=va, zorder=8, linespacing=1.15, **kw)


        def squiggle(ax, x0, y0, y1, amp=0.55, cycles=2.5, c=RED, lw=0.8, z=6):
            yy = np.linspace(y0, y1, 60)
            ax.plot(x0 + amp * np.sin(2 * np.pi * cycles * (yy - y0) / (y1 - y0)), yy, color=c, lw=lw, zorder=z, solid_capstyle='round')


        def dots(ax, x0, x1, y, r=0.28, c=VIOLET, step=2.2):
            for xx in np.arange(x0, x1 + 1e-6, step):
                ax.add_patch(Circle((xx, y), r, fc=c, ec='none', zorder=7))


        # ---------------------------------------------------------------- layout (mm)
        PW1, PH, PY1 = 56, 48, 60                    # row 1: three panels
        PW2, PY2 = 41, 4                             # row 2: four panels
        g1 = (W - 3 * PW1) / 2; g2 = (W - 4 * PW2) / 3
        x1 = [i * (PW1 + g1) for i in range(3)]
        x2 = [i * (PW2 + g2) for i in range(4)]

        # ================================================================ a  CAD design
        ax = panel(x1[0], PY1, PW1, PH, 'a', 'CAD design')
        # flow layer: rounded channel, drawn as the raised mold feature
        rect(ax, 6, 28, 28, 2.6, tint(MOLD, 0.25), ec=MOLD, lw=0.6)
        semi(ax, 20, 30.6, 7, 'none', ec=TEAL, lw=1.0, ls=(0, (3, 2)))
        arrow(ax, (13, 26.3), (27, 26.3), lw=0.6, ms=4, style='<|-|>')
        txt(ax, 20, 24.0, '$W_f$', c=TEAL_D)
        arrow(ax, (29.8, 30.6), (29.8, 37.6), lw=0.6, ms=4, style='<|-|>')
        txt(ax, 31.4, 34.1, '$H_f$', c=TEAL_D, ha='left')
        txt(ax, 37.5, 35.6, 'Flow layer', c=TEAL_D, ha='left', fontweight='bold')
        txt(ax, 37.5, 32.3, 'rounded\nchannel', fs=FS_S, ha='left')
        # control layer: rectangular channel
        rect(ax, 6, 10, 28, 2.6, tint(MOLD, 0.25), ec=MOLD, lw=0.6)
        rect(ax, 13, 12.6, 14, 7, 'none', ec=GOLD_D, lw=1.0, ls=(0, (3, 2)))
        arrow(ax, (13, 8.2), (27, 8.2), lw=0.6, ms=4, style='<|-|>')
        txt(ax, 20, 5.6, '$W_c$', c=GOLD_D)
        arrow(ax, (29.8, 12.6), (29.8, 19.6), lw=0.6, ms=4, style='<|-|>')
        txt(ax, 31.4, 16.1, '$H_c$', c=GOLD_D, ha='left')
        txt(ax, 37.5, 17.6, 'Control layer', c=GOLD_D, ha='left', fontweight='bold')
        txt(ax, 37.5, 14.2, 'rectangular\nchannel', fs=FS_S, ha='left')

        # ================================================================ b  laser SLA printing
        ax = panel(x1[1], PY1, PW1, PH, 'b', 'Laser SLA printing')
        tx0, tx1, tyb, tyt = 4, 38, 14, 25               # resin tank
        rect(ax, tx0 + 0.5, tyb + 0.5, tx1 - tx0 - 1, 9.5, tint(VIOLET, 0.13), z=1)           # resin
        for (xa, ya, xb, yb) in [(tx0, tyt, tx0, tyb), (tx0, tyb, tx1, tyb), (tx1, tyb, tx1, tyt)]:
            ax.plot([xa, xb], [ya, yb], color=SLATE, lw=1.1, zorder=4, solid_capstyle='round')
        rect(ax, 12, 30, 18, 2, SLATE, z=5)                                                    # build platform
        rect(ax, 20, 32, 2, 5, SLATE, z=5)                                                     # arm
        arrow(ax, (26.8, 33.2), (26.8, 38.0), lw=0.8, ms=4.5)
        txt(ax, 28.4, 35.6, 'lifts after\neach layer', fs=FS_S, ha='left')
        nb, nf, lh = 3, 5, 1.7                                                                 # printed part, layer by layer
        yy = 30.0
        for k in range(nb + nf + 1):
            wl = 14 if k < nb else 6
            last = (k == nb + nf)
            rect(ax, 21 - wl / 2, yy - lh, wl, lh, tint(VIOLET, 0.55) if last else MOLD, ec='white', lw=0.35, z=3)
            yy -= lh
        # laser, mirror, beam
        rect(ax, 3, 4, 8, 5, SLATE, z=5); rect(ax, 10.6, 5.7, 0.9, 1.6, VIOLET, z=6)
        ax.plot([11.5, 21, 21], [6.5, 6.5, 14.0], color=VIOLET, lw=1.5, zorder=4, solid_capstyle='butt')
        ax.add_patch(Polygon([[19.9, 5.4], [22.1, 7.6], [21.6, 8.1], [19.4, 5.9]], closed=True, fc=SLATE, ec='none', zorder=6))
        ax.add_patch(Circle((21, 14.4), 1.3, fc=VIOLET, ec='none', alpha=0.35, zorder=5))
        ax.add_patch(Circle((21, 14.4), 0.55, fc=VIOLET, ec='none', zorder=6))
        txt(ax, 7, 1.6, '405 nm laser', fs=FS_S, va='bottom')
        # labels with leaders (kept clear of the beam)
        for (yl, s, xe) in [(31.0, 'build\nplatform', 30.2), (21.6, 'printed mold,\nlayer by layer', 25.0),
                            (17.0, 'resin tank', 38.0), (9.4, 'laser cures\neach layer', None)]:
            if xe is not None:
                ax.plot([xe, 40.3], [yl, yl], color=SLATE, lw=0.5, zorder=4)
                ax.add_patch(Circle((xe, yl), 0.35, fc=SLATE, ec='none', zorder=5))
            txt(ax, 40.8, yl, s, fs=FS_S, ha='left')
        ax.plot([22.5, 27.5, 40.3], [13.6, 9.4, 9.4], color=SLATE, lw=0.5, zorder=4)

        # ================================================================ c  wash and post-cure
        ax = panel(x1[2], PY1, PW1, PH, 'c', 'Wash and post-cure')
        # 1 IPA wash
        rect(ax, 5.5, 14.5, 18, 12.0, tint(TEAL, 0.14), z=1)
        for (xa, ya, xb, yb) in [(5, 30, 5, 14), (5, 14, 24, 14), (24, 14, 24, 30)]:
            ax.plot([xa, xb], [ya, yb], color=SLATE, lw=1.0, zorder=4, solid_capstyle='round')
        rect(ax, 8.5, 17.5, 12, 1.5, MOLD, z=3); rect(ax, 10.0, 19.0, 2.6, 2.0, MOLD, z=3); semi(ax, 16.5, 19.0, 1.6, MOLD)
        xs = np.linspace(8.8, 20.2, 80)
        for yo in (22.6, 23.6):
            ax.plot(xs, yo + 0.35 * np.sin(2 * np.pi * (xs - 8.8) / 3.8), color=TEAL_D, lw=0.6, zorder=4)
        for (bx, by, br) in [(21.6, 20.6, 0.5), (22.3, 22.6, 0.4), (21.2, 24.6, 0.45)]:
            ax.add_patch(Circle((bx, by), br, fc='white', ec=TEAL_D, lw=0.4, zorder=4))
        txt(ax, 14.5, 25.4, 'IPA', c=TEAL_D, fontweight='bold')
        txt(ax, 14.5, 11.6, 'IPA wash', fontweight='bold')
        txt(ax, 14.5, 8.8, 'removes uncured resin', fs=FS_S, c='#5c6b78')
        arrow(ax, (25.8, 21), (30.4, 21), lw=0.9, ms=5)
        # 2 UV post-cure (+ optional bake)
        rect(ax, 32, 30, 19, 3.6, SLATE, z=4); txt(ax, 41.5, 31.8, 'UV', c='white', fontweight='bold')
        for (xa, xb) in [(34.5, 34.0), (38.0, 37.6), (41.5, 41.5), (45.0, 45.4), (48.5, 49.0)]:
            arrow(ax, (xa, 29.6), (xb, 21.4), c=VIOLET, lw=0.9, ms=4)
        rect(ax, 32, 14, 19, 2, tint(SLATE, 0.8), ec=SLATE, lw=0.5, z=3)
        rect(ax, 34, 16, 15, 1.5, MOLD, z=3); rect(ax, 36, 17.5, 2.6, 2.0, MOLD, z=3); semi(ax, 44.2, 17.5, 1.6, MOLD)
        squiggle(ax, 33.2, 16.5, 21.5); squiggle(ax, 49.8, 16.5, 21.5)
        txt(ax, 41.5, 11.6, 'UV post-cure', fontweight='bold')
        txt(ax, 41.5, 8.8, '± thermal bake', fs=FS_S, c='#5c6b78')
        txt(ax, 28, 3.6, 'prevents PDMS cure inhibition', fs=FS_S, c='#5c6b78', style='italic')

        # ================================================================ d  PDMS casting and curing
        ax = panel(x2[0], PY2, PW2, PH, 'd', 'PDMS casting')
        rect(ax, 4, 30, 33, 9, tint(TEAL, 0.32), z=1)                          # thick flow layer (7:1)
        rect(ax, 4, 27.5, 33, 2.5, MOLD, z=2); semi(ax, 20.5, 30, 6, MOLD)
        txt(ax, 20.5, 24.6, 'flow layer', c=TEAL_D, fontweight='bold')
        txt(ax, 20.5, 22.0, 'cast, 7:1', fs=FS_S)
        rect(ax, 4, 15, 33, 6, tint(GOLD, 0.45), z=1)                          # thin control layer (20:1)
        rect(ax, 4, 12.5, 33, 2.5, MOLD, z=2); rect(ax, 15, 15, 11, 4.2, MOLD, z=2)
        txt(ax, 20.5, 10.0, 'control layer', c=GOLD_D, fontweight='bold')
        txt(ax, 20.5, 7.4, 'spin-coated, 20:1', fs=FS_S)
        for xs_ in (8, 12, 16):
            squiggle(ax, xs_, 1.4, 4.8, c=RED, amp=0.45, cycles=1.5)
        txt(ax, 19.5, 3.1, 'thermal cure', fs=FS_S, ha='left')

        # ================================================================ e  release, plasma and bond (section along the control channel)
        ax = panel(x2[1], PY2, PW2, PH, 'e', 'Plasma bonding')
        rect(ax, 2, 3, 37, 2.5, MOLD, z=2)                                     # control mold
        rect(ax, 2, 5.5, 37, 7, tint(GOLD, 0.45), z=1); rect(ax, 5, 5.5, 31, 3.5, MOLD, z=2)
        rect(ax, 2, 18, 37, 15, tint(TEAL, 0.32), z=1); semi(ax, 20.5, 18, 4.8, PANEL, z=3)     # released flow layer, channel facing down
        dots(ax, 3.5, 38, 13.3); dots(ax, 3.5, 38, 17.2)
        txt(ax, 3.0, 15.3, 'O$_2$ plasma', fs=FS_S, c=VIOLET, ha='left')
        arrow(ax, (35.4, 16.9), (35.4, 14.1), lw=1.0, ms=5)
        txt(ax, 20.5, 35.4, 'flow layer (released)', fs=FS_S, c=TEAL_D)
        txt(ax, 20.5, 1.2, 'control layer on mold', fs=FS_S, c=GOLD_D)

        # ================================================================ f  peel and punch
        ax = panel(x2[2], PY2, PW2, PH, 'f', 'Peel and punch ports')
        rect(ax, 2, 5, 37, 7, tint(GOLD, 0.45), z=1); rect(ax, 5, 5, 31, 3.5, PANEL, z=2)       # control groove open to the bottom
        rect(ax, 2, 12, 37, 15, tint(TEAL, 0.32), z=1); semi(ax, 20.5, 12, 4.2, PANEL, z=3)      # flow channel (closed)
        rect(ax, 5.5, 8.5, 2, 18.5, PANEL, z=3); rect(ax, 33.5, 8.5, 2, 18.5, PANEL, z=3)        # punched ports
        rect(ax, 6.1, 21.0, 0.8, 11.0, SLATE, z=5); rect(ax, 4.5, 32.0, 4.0, 2.0, SLATE, z=5)    # punch
        arrow(ax, (11.0, 33.8), (11.0, 29.6), lw=1.0, ms=5)
        txt(ax, 20.5, 22.4, 'ports', fs=FS_S, c='#5c6b78')
        ax.plot([17.6, 8.2], [22.4, 22.4], color='#5c6b78', lw=0.4, zorder=5)
        ax.plot([23.4, 32.9], [22.4, 22.4], color='#5c6b78', lw=0.4, zorder=5)
        txt(ax, 20.5, 1.8, 'mold removed', fs=FS_S, c='#5c6b78')

        # ================================================================ g  bond to glass (final device)
        ax = panel(x2[3], PY2, PW2, PH, 'g', 'Bonding to glass')
        rect(ax, 2, 3, 37, 3, GLASS, ec=MOLD, lw=0.4, z=1); txt(ax, 20.5, 4.5, 'glass', fs=FS_S)
        rect(ax, 2, 6, 37, 7, tint(GOLD, 0.45), z=1); rect(ax, 5, 6, 31, 3.5, PANEL, z=2)        # sealed control channel
        rect(ax, 16.3, 9.5, 8.4, 3.5, OVER, z=2, alpha=0.85)                                      # membrane at the valve
        rect(ax, 2, 13, 37, 15, tint(TEAL, 0.32), z=1); semi(ax, 20.5, 13, 4.2, PANEL, z=3)       # flow channel
        rect(ax, 5.5, 9.5, 2, 18.5, PANEL, z=3); rect(ax, 33.5, 9.5, 2, 18.5, PANEL, z=3)
        dots(ax, 3.5, 37.5, 6.0, r=0.25)
        txt(ax, 20.5, 7.75, 'control channel', fs=FS_S, c=GOLD_D)
        txt(ax, 20.5, 14.9, 'flow', fs=FS_S, c=TEAL_D)
        txt(ax, 15.4, 11.3, 'membrane', fs=5.0, ha='right')
        arrow(ax, (6.5, 33.8), (6.5, 29.2), lw=1.0, ms=5)
        txt(ax, 9.3, 31.6, 'pressure', fs=FS_S, ha='left')

        # ---------------------------------------------------------------- arrows between steps and between rows
        for xa, xb, yc in [(x1[0] + PW1, x1[1], PY1 + PH / 2), (x1[1] + PW1, x1[2], PY1 + PH / 2)] + \
                          [(x2[i] + PW2, x2[i + 1], PY2 + PH / 2) for i in range(3)]:
            arrow(bg, (xa + 0.7, yc), (xb - 0.7, yc), lw=1.0, ms=6)
        xc = x1[2] + PW1 / 2; xd = x2[0] + PW2 / 2; ym = (PY1 + PY2 + PH) / 2
        bg.plot([xc, xc], [PY1 - 0.3, ym], color=SLATE, lw=1.0, solid_capstyle='butt')
        bg.plot([xc, xd], [ym, ym], color=SLATE, lw=1.0, solid_capstyle='butt')
        arrow(bg, (xd, ym), (xd, PY2 + PH + 0.3), lw=1.0, ms=6)
        bg.text((xc + xd) / 2, ym + 1.0, 'printed molds', fontsize=FS_S, color='#5c6b78', ha='center', va='bottom', style='italic')

    # mathtext fonts are looked up when the figure is drawn (displayed or saved), i.e. after the rc_context has closed;
    # keep them set globally (identical to lib/ch2_style.py, so the other figures are unaffected)
    mpl.rcParams.update({k: v for k, v in STYLE.items() if k.startswith('mathtext.')})
    return fig
