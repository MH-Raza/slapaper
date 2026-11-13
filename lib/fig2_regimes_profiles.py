"""
Figure 2: valve regimes (schematic, panel a) and representative mold / replica profiles (panels b, c) from the
FINAL profilometry pipeline (profilometry_analysis.py). Panel a is the schematic of Figure_1_Rev_V2.ipynb (idealised dome,
normalised units); its legend now says 'schematic' because the real flow channels are three-point-arc segments (Section 2.1).
Panels b, c: C1 (rectangle, W_c = 500 um) and F3 (H_f = 100 um). Profiles: each chip's median registered scan line; line =
mean of the 3 chips, band = SD between chips. Bars: mean +/- SD over 3 chips, open circles = single chips. No edge clipping.
Reference for the widths: the design silhouette at 10 % of the height (C1 500 um, F3 337.4 um).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from ch2_style import *
import profilometry_analysis as PA

C_CAD = C_TEXT; C_MEMBRANE = '#B8860B'; C_GLASS = '#EAECEE'; C_PDMS_BG = '#F4F6F8'
C_FLUID_0 = '#FFFFFF'; C_FLUID_A = '#FBF0E1'; C_FLUID_M = '#F2C793'

# ---- panel a: copied from the notebook (idealised schematic, normalised units) -------------------------------------------
W_f_s      = 2.0    # normalised flow channel width
H_f_s      = 1.0    # normalised flow channel height
W_c_narrow = 0.8    # bridge regime W_c
W_c_wide   = 2.0    # normal regime W_c
H_c_s      = 0.9    # control channel depth


def get_membrane(x, state, W_c):
    """
    Compute PDMS membrane deflection profile.

    States:
      open        : flat at y=0 across full width
      bridge      : partial deflection, peak stays below flow roof
      bridge_high : max pressure, peak hard-capped below roof (leakage gap)
      deflect     : normal regime partial deflection
      close       : full conformal seal, membrane matches flow roof
    """
    roof = np.where(
        np.abs(x) <= W_f_s / 2,
        H_f_s * np.sqrt(np.maximum(1 - (x / (W_f_s / 2)) ** 2, 0)),
        0)
    y    = np.zeros_like(x)

    if state == 'open':
        return y

    mask   = np.abs(x) <= W_c / 2
    x_norm = x[mask] / (W_c / 2)

    if state == 'bridge':
        y[mask] = (H_f_s * 0.65) * (1 - x_norm ** 2) ** 2

    elif state == 'bridge_high':
        raw      = (H_f_s * 0.90) * (1 - x_norm ** 2) ** 1.5
        y[mask]  = np.minimum(raw, roof[mask] * 0.82)   # hard gap preserved

    elif state == 'deflect':
        y[mask] = (H_f_s * 0.82) * (1 - x_norm ** 2) ** 2

    elif state == 'close':
        y    = np.copy(roof)
        y    = np.where(np.abs(x) > W_c / 2, 0, y)
        return np.maximum(y - 0.03, 0)

    limit = np.maximum(roof - 0.04, 0)
    k     = 25.0
    y     = np.minimum(y, limit) - np.log(1 + np.exp(-k * np.abs(y - limit))) / k
    return np.maximum(y, 0)


def draw_schematic(ax, W_c, state, col_idx, is_top_row):
    """Draw one schematic panel for Panel a."""
    x    = np.linspace(-1.8, 1.8, 3000)
    roof = np.where(
        np.abs(x) <= W_f_s / 2,
        H_f_s * np.sqrt(np.maximum(1 - (x / (W_f_s / 2)) ** 2, 0)),
        0)

    # Backgrounds
    ax.fill_between([-1.8, 1.8], -H_c_s,  1.5,   color=C_PDMS_BG, zorder=0)
    ax.fill_between([-1.8, 1.8], -1.05,  -H_c_s, color=C_GLASS,   zorder=0)
    ax.axhline(-H_c_s, color='#BDC3C7', lw=0.5, zorder=1)

    # Fluid fill
    fc = (C_FLUID_0 if state == 'open' else
          C_FLUID_M if state in ['bridge_high', 'close'] else
          C_FLUID_A)
    ax.fill_between(x, 0, roof, where=np.abs(x) <= W_f_s / 2,
                    color=C_FLUID_0, zorder=1)
    ax.fill_between([-W_c / 2, W_c / 2], -H_c_s, 0, color=fc, zorder=1)

    # CAD outlines
    arc_x = np.linspace(-W_f_s / 2, W_f_s / 2, 500)
    arc_y = H_f_s * np.sqrt(1 - (arc_x / (W_f_s / 2)) ** 2)
    ax.plot(arc_x, arc_y,
            color=C_FLOW, ls='--', lw=0.9, dashes=(3, 2), zorder=2)
    ax.plot([-W_c/2, -W_c/2, W_c/2, W_c/2],
            [-H_c_s, 0, 0, -H_c_s],
            color=C_CTRL, ls='--', lw=0.9, dashes=(3, 2), zorder=2)

    # PDMS membrane
    mem = get_membrane(x, state, W_c)
    if state == 'open':
        ax.plot([-1.8, 1.8], [0, 0], color=C_MEMBRANE, lw=2.2,
                solid_capstyle='butt', zorder=4, clip_on=True)
    else:
        ax.plot(x, mem, color=C_MEMBRANE, lw=2.2,
                solid_capstyle='round', solid_joinstyle='round', zorder=4)

    # Conformal seal fill
    if state == 'close':
        ax.fill_between(x, 0, mem, where=(np.abs(x) <= W_c / 2),
                        color=C_FLUID_M, alpha=0.55, zorder=3)

    # Pressure annotation
    if state == 'open':
        ax.text(0, -0.60, r'$P=0$', color=C_TEXT,
                fontsize=BASE_FS, ha='center', va='center')
    else:
        p_txt = r'$P_{act}$' if state in ['bridge', 'deflect'] else r'$P_{max}$'
        y0    = -0.58 if state in ['bridge', 'deflect'] else -0.62
        dy    = 0.36  if state in ['bridge', 'deflect'] else 0.40
        hw    = 0.09  if state in ['bridge', 'deflect'] else 0.11
        hl    = 0.12  if state in ['bridge', 'deflect'] else 0.15
        ax.arrow(0, y0, 0, dy, head_width=hw, head_length=hl,
                 fc=C_TEXT, ec=C_TEXT, lw=0.7, length_includes_head=True, zorder=5)
        ax.text(0, y0 - 0.06, p_txt, color=C_TEXT,
                fontsize=BASE_FS, ha='center', va='top')

    # Dimension annotations (open state only)
    if state == 'open':
        xd = -W_f_s / 2 - 0.22
        ax.annotate('', xy=(xd, 0), xytext=(xd, H_f_s),
                    arrowprops=dict(arrowstyle='<->', color=C_FLOW, lw=0.6, mutation_scale=5))
        ax.text(xd - 0.08, H_f_s / 2, r'$H_f$',
                color=C_FLOW, fontsize=BASE_FS, fontweight='bold', ha='right', va='center')
        yd = -0.46
        ax.annotate('', xy=(-W_c / 2, yd), xytext=(W_c / 2, yd),
                    arrowprops=dict(arrowstyle='<->', color=C_TEXT, lw=0.6, mutation_scale=5))
        ax.text(0, yd + 0.06, r'$W_c$', color=C_TEXT,
                fontsize=BASE_FS, fontweight='bold', ha='center', va='bottom')

    # State-specific labels
    if state == 'bridge_high':
        mem_peak = float(mem[len(mem) // 2]) if len(mem) > 0 else 0.7
        gap_y    = (mem_peak + H_f_s) / 2
        ax.annotate('Leakage\nGap', xy=(0.08, gap_y), xytext=(1.10, 0.92),
                    arrowprops=dict(arrowstyle='->', color=C_TEXT, lw=0.7,
                                    mutation_scale=6),
                    fontsize=BASE_FS - 0.5, color=C_TEXT,
                    ha='center', va='bottom', fontweight='bold')
    elif state == 'close':
        ax.text(0, 0.35, 'Conformal\nSeal', color=C_TEXT,
                fontsize=BASE_FS - 0.5, fontweight='bold',
                ha='center', va='center', linespacing=1.3, zorder=6)

    ax.set_aspect('equal')
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.20, 1.42)
    ax.axis('off')

    if is_top_row:
        headers = [r'Open State ($P=0$)',
                   r'Actuated ($P_{act}$)',
                   r'High Pressure ($P_{max}$)']
        ax.text(0, 1.44, headers[col_idx],
                fontsize=BASE_FS + 0.5, color=C_TEXT,
                fontweight='bold', ha='center', va='bottom')


def fig2(R, F, store):
    fig = plt.figure(figsize=(7.2, 9.0))
    gs1 = gridspec.GridSpec(5, 1, height_ratios=[1.15, 0.12, 0.85, 0.12, 0.55], hspace=0.20, left=0.09, right=0.97, top=0.965, bottom=0.055)
    gs_a = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs1[0], wspace=0.04, hspace=-0.15)
    axes_a = [[fig.add_subplot(gs_a[r, c]) for c in range(3)] for r in range(2)]
    cfg = [(W_c_narrow, 'open', 0, True, 0), (W_c_narrow, 'bridge', 1, True, 0), (W_c_narrow, 'bridge_high', 2, True, 0),
           (W_c_wide, 'open', 0, False, 1), (W_c_wide, 'deflect', 1, False, 1), (W_c_wide, 'close', 2, False, 1)]
    for (Wc, state, col, top, row) in cfg: draw_schematic(axes_a[row][col], Wc, state, col, top)
    for ri, label in enumerate(['Bridge Regime\n' + r'($W_c < H_f$)', 'Normal Regime\n' + r'($W_c \geq H_f$)']):
        axes_a[ri][0].text(-2.32, 0.10, label, va='center', ha='center', fontsize=BASE_FS, color=C_TEXT, fontweight='bold', rotation=90,
                           transform=axes_a[ri][0].transData, clip_on=False)
    axes_a[0][0].text(-2.78, 1.58, 'a', fontsize=10, fontweight='bold', color=C_TEXT, va='top', ha='left', transform=axes_a[0][0].transData, clip_on=False)
    axes_a[1][1].legend(handles=[Line2D([0], [0], color=C_FLOW, lw=1.2, ls='--', dashes=(3, 2), label='Flow channel (schematic)'),
                                 Line2D([0], [0], color=C_CTRL, lw=1.2, ls='--', dashes=(3, 2), label='Control channel (schematic)'),
                                 Line2D([0], [0], color=C_MEMBRANE, lw=2.2, label='PDMS membrane'),
                                 patches.Patch(facecolor=C_GLASS, edgecolor='#AAAAAA', lw=0.5, label='Glass substrate')],
                        loc='upper center', bbox_to_anchor=(0.50, 0), ncol=4, frameon=False, fontsize=BASE_FS - 0.5, handlelength=2.0,
                        handletextpad=0.5, columnspacing=1.2)
    # ---- panel b: profiles ------------------------------------------------------------------------------------------------
    st_b = strip_axes(fig, gs1[1])
    gs_b = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs1[2], wspace=0.22)
    ax_b = [fig.add_subplot(gs_b[0]), fig.add_subplot(gs_b[1])]
    for ax, ch, ccol, ttl, ylim, yt in ((ax_b[0], 'C1', C_CTRL, r'Rectangular profile ($W_c$ = 500 µm)', (-15, 265), [0, 100, 200]),
                                        (ax_b[1], 'F3', C_FLOW, r'Circular-segment profile ($H_f$ = 100 µm)', (-10, 150), [0, 50, 100])):
        ax.plot(PA.GRID, PA.cad_profile(ch, PA.GRID), color=C_CAD, ls='--', lw=0.9, dashes=(3, 2), zorder=5)
        for mat, col in (('Mold', C_MOLD), ('PDMS', ccol)):
            P = np.vstack([np.median(store[(ch, mat, r)], axis=0) for r in (1, 2, 3)]); m, s = P.mean(0), P.std(0, ddof=1)
            ax.plot(PA.GRID, m, color=col, lw=1.0, zorder=4 if mat == 'Mold' else 3); ax.fill_between(PA.GRID, m - s, m + s, color=col, alpha=0.16, lw=0)
        ax.set_title(ttl, fontsize=BASE_FS, fontweight='bold', pad=4); ax.set_xlim(-430, 430); ax.set_ylim(*ylim); ax.set_yticks(yt)
        ax.set_xlabel('Width (µm)', labelpad=2); ax.grid(axis='y', linestyle=':', linewidth=0.4, color='#CCCCCC', alpha=0.5, zorder=0); ax.tick_params(length=2.5, pad=2)
    ax_b[0].set_ylabel('Height (µm)', labelpad=2); ax_b[0].text(-0.18, 1.12, 'b', transform=ax_b[0].transAxes, fontsize=10, fontweight='bold', va='top', ha='left')
    st_b.legend(handles=[Line2D([0], [0], color=C_CAD, ls='--', lw=0.9, dashes=(3, 2), label='Design silhouette'),
                         Line2D([0], [0], color=C_MOLD, lw=1.0, label='SLA mold'), Line2D([0], [0], color=C_CTRL, lw=1.0, label='PDMS, control'),
                         Line2D([0], [0], color=C_FLOW, lw=1.0, label='PDMS, flow')], loc='center', ncol=4, frameon=False, fontsize=BASE_FS - 0.5,
                handlelength=1.6, handletextpad=0.4, columnspacing=1.6)
    # ---- panel c: height and footprint against the design ---------------------------------------------------------------------
    st_c = strip_axes(fig, gs1[3])
    gs_c = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs1[4], wspace=0.30)
    ax_c = [fig.add_subplot(gs_c[0]), fig.add_subplot(gs_c[1])]
    chs = ['C1', 'F3']; ccols = [C_CTRL, C_FLOW]; x2 = np.array([0, 1]); bw = 0.28
    for ax, key, ylabel, cadv, ylim_, yt in ((ax_c[0], 'H', 'Height (µm)', [PA.CAD[c][0] for c in chs], (0, 260), [0, 100, 200]),
                                             (ax_c[1], 'W10', 'Footprint $W_{10}$ (µm)', [PA.CAD_W[c][10] for c in chs], (0, 820), [0, 200, 400, 600, 800])):
        for j, (ch, cc) in enumerate(zip(chs, ccols)):
            for mat, off, col in (('Mold', -bw / 2, C_MOLD), ('PDMS', bw / 2, cc)):
                v = R[(R.channel == ch) & (R.material == mat)].sort_values('chip')[key].values
                ax.bar(x2[j] + off, v.mean(), bw, yerr=v.std(ddof=1), color=col, alpha=0.88, error_kw={'lw': 0.8, 'capsize': 3.0, 'capthick': 0.8}, zorder=3)
                ax.scatter(np.full(3, x2[j] + off) + np.array([-0.06, 0, 0.06]), v, s=6, facecolor='white', edgecolor=C_TEXT, lw=0.5, zorder=6)
            ax.plot([x2[j] - bw * 1.1, x2[j] + bw * 1.1], [cadv[j]] * 2, color=C_CAD, lw=1.2, ls='--', zorder=5)
        ax.set_xticks(x2); ax.set_xticklabels([r'C1  ($W_c$ = 500 µm)', r'F3  ($H_f$ = 100 µm)'], fontsize=BASE_FS - 0.5)
        ax.set_ylabel(ylabel, labelpad=2); ax.set_ylim(ylim_); ax.set_yticks(yt); ax.set_xlim(-0.55, 1.55)
        ax.grid(axis='y', linestyle=':', linewidth=0.4, color='#CCCCCC', alpha=0.5, zorder=0); ax.tick_params(length=2.5, pad=2)
    ax_c[0].text(-0.20, 1.18, 'c', transform=ax_c[0].transAxes, fontsize=10, fontweight='bold', va='top', ha='left')
    st_c.legend(handles=[Line2D([0], [0], color=C_CAD, ls='--', lw=1.2, label='Design'), patches.Patch(facecolor=C_MOLD, alpha=0.88, label='SLA mold'),
                         patches.Patch(facecolor=C_CTRL, alpha=0.88, label='PDMS, control'), patches.Patch(facecolor=C_FLOW, alpha=0.88, label='PDMS, flow'),
                         Line2D([0], [0], marker='o', ls='', mfc='white', mec=C_TEXT, mew=0.5, ms=3.2, label='Single chip')],
                loc='center', ncol=5, frameon=False, fontsize=BASE_FS - 0.5, handlelength=1.4, handletextpad=0.4, columnspacing=1.3)
    return fig


if __name__ == '__main__':
    R, store = PA.run(); F = PA.per_feature(R)
    f2 = fig2(R, F, store)
    # f2.savefig('Fig2_regimes_and_profiles.jpeg', dpi=600, bbox_inches='tight')    # save figure (uncomment to export)
    plt.show()
