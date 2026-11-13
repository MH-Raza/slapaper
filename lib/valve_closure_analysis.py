"""
Valve-closure analysis (25 geometries x 5 pressures = 125 runs) and figures.

DATA    one tab-separated log per run, 10 Hz, latin-1:  time (s) | OB1 ch1 target/read (flow-driving pressure, nominal
        100 mbar) | OB1 ch2 target/read (control-line supply, mbar) | MFS2 target/read (flow, uL/min).
        File name 'F<i> C<j> <P> mbar.txt'.  H_f: F1-F5 = 350/200/100/50/25 um;  W_c: C1-C5 = 500/300/200/50/25 um.
PROTOCOL (experimenter + verified in the traces)  the control supply is constant for the whole 180 s run; the valve is
        switched by an external solenoid that is not logged: open 0-30 s, closed 30-60 s, open, closed, open, closed
        (3 cycles; closing edges at 30.8/90.9/151.0 s, reopening at 60/120 s in the 74 runs that close cleanly).
        Some logs contain a restarted acquisition (time resets); the last continuous acquisition is the valid one.

METHOD
  * Analysis windows are the ones chosen a priori by the experimenter to avoid switching transients and sensor noise:
    open 10-25, 70-85, 130-145 s; closed 40-55, 100-115, 160-175 s (5 s before the next switch, >= 9 s after closing).
    Alternative window sets are used ONLY as a sensitivity analysis (Supplementary Fig. S-V4).
  * Closure per cycle = 100 x (1 - Q_closed / Q_open), Q = window mean of the MFS2 reading; each cycle uses its own open
    plateau. Values are NOT clipped to 0-100 (small negative values of permanently open crossings are drift and noise).
  * Replication: 3 actuation cycles of ONE device per geometry = technical replicates (cycle-to-cycle repeatability);
    they are not independent devices. Mean +/- SD (ddof = 1) is reported; no test treats cycles as independent devices.
  * Sealing: mean closure >= 95%. The threshold pressure of a geometry lies in (highest non-sealing tested pressure,
    P_min]; it is therefore reported as an interval, with left-censoring (sealed already at 500 mbar) and right-censoring
    (not sealed at 1500 mbar) shown as such. Pressure resolution is 250 mbar.
  * Data checks on every file: logged control setpoint and read-back must match the file name, drive pressure must be
    ~100 mbar, duration 180 s; time resets are handled as above; flow-sensor readings at the logging ceiling (8.00 uL/min)
    in the open windows are flagged (closure is then a lower bound).
  * 'F3 C4 1000 mbar.txt' in the data folder was acquired at 1250 mbar (logged setpoint 1250, read-back 1262 mbar) and
    is excluded; the run supplied separately (FILE_OVERRIDES) has setpoint 1000 (read-back 1009), 100 mbar drive,
    the standard 30 s schedule and 180 s duration, and replaces it. That it comes from the same device cannot be verified.
  * As-fabricated aspect ratio AR = W_c,fab / H_f,fab from the profilometry of the PDMS replicas: W_c,fab = full width at
    half maximum (W50) of the control channel, H_f,fab = depth of the flow channel; chip means, SD propagated in
    quadrature (the valve chips themselves were not profiled; this is the process spread). The width definition is
    varied in Supplementary Fig. S-V3.
"""
import os, re, glob
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from ch2_style import *

DATA_DIR = r'./Valve data'
EXTRA_FILES = {('F3', 'C4', 1000): r'./F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt'}     # replaces the mislabelled run
EXCLUDE = {'F3 C4 1000 mbar.txt'}
PROF_PER_FILE = r'./results/profilometry_per_file.csv'
OUT_DIR = r'./results'

FLOW_H = {'F1': 350, 'F2': 200, 'F3': 100, 'F4': 50, 'F5': 25}
CTRL_W = {'C1': 500, 'C2': 300, 'C3': 200, 'C4': 50, 'C5': 25}
OPEN_WIN = [(10, 25), (70, 85), (130, 145)]
CLOSED_WIN = [(40, 55), (100, 115), (160, 175)]
ALT_WINDOWS = {'last 10 s of each phase': ([(20, 30), (80, 90), (140, 150)], [(50, 60), (110, 120), (170, 180)]),
               'conservative (15 s margin)': ([(5, 15), (65, 75), (125, 135)], [(45, 55), (105, 115), (165, 175)])}
SENSOR_CEILING = 8.00
SEAL = 95.0


def load(path):
    df = pd.read_csv(path, sep='\t', encoding='latin-1')
    df = df.loc[:, ~df.columns.str.startswith('Unnamed')]
    df.columns = ['t', 'p1_target', 'p1_read', 'p2_target', 'p2_read', 'q_target', 'q']
    resets = np.where(np.diff(df.t.values) < 0)[0]
    if len(resets):
        df = df.iloc[resets[-1] + 1:].reset_index(drop=True)
    df['t'] = df.t - df.t.iloc[0]
    return df, len(resets)


def parse(name):
    m = re.match(r'^(F\d) (C\d) (\d+) mbar\.txt$', name)
    return (m.group(1), m.group(2), int(m.group(3))) if m else None


def run_files(data_dir=DATA_DIR, extra=EXTRA_FILES):
    files = {}
    for f in sorted(glob.glob(os.path.join(data_dir, '*.txt'))):
        p = parse(os.path.basename(f))
        if p: files[p] = f
    for k, f in extra.items():
        files[k] = f
    return files


def cycle_closure(df, ow, cw):
    t, q = df.t.values, df.q.values
    out = []
    for (a, b), (c, d) in zip(ow, cw):
        qo = q[(t > a) & (t < b)].mean(); qc = q[(t > c) & (t < d)].mean()
        out.append((qo, qc, 100 * (1 - qc / qo)))
    return out


def analyse(data_dir=DATA_DIR, extra=EXTRA_FILES):
    cyc, qc = [], []
    for (F, C, P), f in sorted(run_files(data_dir, extra).items()):
        df, nres = load(f); t, q = df.t.values, df.q.values
        tgt = float(df.p2_target[df.p2_target > 0].median()); rd = float(df.p2_read[df.t > 5].median()); dr = float(df.p1_read[df.t > 5].median())
        valid = abs(tgt - P) <= 1 and abs(rd - P) <= 25 and abs(dr - 100) <= 5 and 179 < t[-1] < 181
        sat = float(np.mean([np.mean(q[(t > a) & (t < b)] >= SENSOR_CEILING - 0.01) for a, b in OPEN_WIN]))
        qc.append(dict(F=F, C=C, P=P, file=os.path.basename(f), target_logged=tgt, readback=rd, drive=dr, duration_s=float(t[-1]),
                       time_resets=nres, valid=valid, open_window_frac_at_ceiling=sat))
        for k, (o, c) in enumerate(zip(cycle_closure(df, OPEN_WIN, CLOSED_WIN), [None] * 3)):
            alt = {n: cycle_closure(df, *w)[k][2] for n, w in ALT_WINDOWS.items()}
            cyc.append(dict(F=F, C=C, Hf=FLOW_H[F], Wc=CTRL_W[C], P=P, cycle=k + 1, Q_open=o[0], Q_closed=o[1], closure=o[2],
                            valid=valid, sat_frac=sat, **{f'closure_{n}': v for n, v in alt.items()}))
    D = pd.DataFrame(cyc); Q = pd.DataFrame(qc)
    # the excluded mislabelled run is documented in the QC table only
    bad = os.path.join(data_dir, 'F3 C4 1000 mbar.txt')
    if os.path.exists(bad):
        df, nres = load(bad)
        Q = pd.concat([Q, pd.DataFrame([dict(F='F3', C='C4', P=1000, file='F3 C4 1000 mbar.txt (EXCLUDED)',
                       target_logged=float(df.p2_target[df.p2_target > 0].median()), readback=float(df.p2_read[df.t > 5].median()),
                       drive=float(df.p1_read[df.t > 5].median()), duration_s=float(df.t.iloc[-1]), time_resets=nres, valid=False,
                       open_window_frac_at_ceiling=np.nan)])], ignore_index=True)
    V = D[D.valid]
    agg = dict(closure=('closure', 'mean'), sd=('closure', lambda v: v.std(ddof=1)), cmin=('closure', 'min'), cmax=('closure', 'max'),
               Q_open=('Q_open', 'mean'), sat_frac=('sat_frac', 'first'))
    for n in ALT_WINDOWS: agg[f'closure_{n}'] = (f'closure_{n}', 'mean')
    G = V.groupby(['F', 'C', 'Hf', 'Wc', 'P']).agg(**agg).reset_index()
    return D, Q, G


def audit_run(F, C, P, data_dir=DATA_DIR, extra=EXTRA_FILES):
    """Print the numbers behind one run: file, logged pressures, window means of the flow, closure of each cycle."""
    f = run_files(data_dir, extra)[(F, C, P)]; df, nres = load(f); t, q = df.t.values, df.q.values
    print(f'{os.path.basename(f)}: control setpoint {df.p2_target[df.p2_target > 0].median():.0f} mbar, read-back {df.p2_read[df.t > 5].median():.0f} mbar, '
          f'drive {df.p1_read[df.t > 5].median():.0f} mbar, {t[-1]:.1f} s, time resets {nres}')
    cl = []
    for k, ((a, b), (c, d)) in enumerate(zip(OPEN_WIN, CLOSED_WIN)):
        qo, qc = q[(t > a) & (t < b)].mean(), q[(t > c) & (t < d)].mean(); cl.append(100 * (1 - qc / qo))
        print(f'  cycle {k+1}: Q_open[{a}-{b} s] = {qo:.3f}, Q_closed[{c}-{d} s] = {qc:.3f} uL/min -> closure {cl[-1]:.1f}%')
    print(f'  mean {np.mean(cl):.1f}% (SD {np.std(cl, ddof=1):.1f}) -> {"seals (>= 95%)" if np.mean(cl) >= SEAL else "does not seal (< 95%)"}')


def pmin_table(G, col='closure'):
    rows = []
    for (F, C), g in G.groupby(['F', 'C']):
        g = g.set_index('P').sort_index()
        sealed = [p for p in PRESSURES if p in g.index and g.loc[p, col] >= SEAL]
        pm = min(sealed) if sealed else np.nan
        below = [p for p in PRESSURES if p in g.index and (np.isnan(pm) or p < pm)]
        lower = max(below) if below else np.nan                      # highest tested non-sealing pressure (NaN: sealed at 500)
        rows.append(dict(F=F, C=C, Hf=FLOW_H[F], Wc=CTRL_W[C], P_min=pm, P_low=lower, closure_1500=g.loc[1500, col],
                         closure_500=g.loc[500, col], not_sealed_1500=bool(np.isnan(pm)),
                         borderline=bool(((g[col] - SEAL).abs() < 1.0).any())))
    return pd.DataFrame(rows)


def add_aspect_ratio(T, prof_path=PROF_PER_FILE, wdef='W50'):
    R = pd.read_csv(prof_path); R = R[R.material == 'PDMS']
    m = R.groupby('channel')[[wdef, 'H']].mean(); s = R.groupby('channel')[[wdef, 'H']].std(ddof=1)
    T = T.copy(); W, Ws = m.loc[T.C, wdef].values, s.loc[T.C, wdef].values; H, Hs = m.loc[T.F, 'H'].values, s.loc[T.F, 'H'].values
    T['Wc_fab'], T['Hf_fab'] = W, H
    T['AR_fab'] = W / H; T['AR_fab_sd'] = T.AR_fab * np.hypot(Ws / W, Hs / H); T['AR_design'] = T.Wc / T.Hf
    return T


def separation(T, wdefs=('W10', 'W20', 'W50', 'W80'), prof_path=PROF_PER_FILE):
    """Bracket in which sealing at 1500 mbar switches on, for all flow heights and for H_f <= 200 um."""
    rows = []
    for wd in wdefs:
        A = add_aspect_ratio(T, prof_path, wd)
        for lab, sub in (('all', A), ('Hf<=200', A[A.Hf <= 200])):
            s = ~sub.not_sealed_1500
            lo = sub[~s].AR_fab.max(); hi = sub[s].AR_fab.min()
            rows.append(dict(width_def=wd, subset=lab, largest_unsealed_AR=lo, smallest_sealed_AR=hi, separable=bool(lo < hi)))
    return pd.DataFrame(rows)


# =========================================================================================== FIGURE 4 (main text)
def _flow_cad_silhouette(h, w=350, n=240):
    """Three-point-arc cross-section (chord w, sagitta h) seen from above: (x, z)."""
    R = h / 2 + w ** 2 / (8 * h); cy = h - R; xmax = R if h > w / 2 else w / 2
    x = np.linspace(-xmax, xmax, n); z = np.maximum(cy + np.sqrt(np.maximum(R ** 2 - x ** 2, 0)), 0)
    return x, z, 2 * xmax


def fig4(G, T, data_dir=DATA_DIR, extra=EXTRA_FILES):
    """Same layout and theme as Figure_2Rev_Final: a matrix schematic, b normal vs bridge, c threshold vs aspect ratio,
    d closure maps; the AR/threshold panel uses as-fabricated dimensions and shows the pressure interval."""
    fig = plt.figure(figsize=(7.2, 4.8), dpi=150)
    LEFT, RIGHT, TOP, BOTTOM, VGAP = 0.07, 0.93, 0.95, 0.14, 0.015
    w_a = 0.222; row1_h = w_a * 7.2 / 4.8                       # panel a is square
    ROW1_TOP, ROW1_BOTTOM = TOP, TOP - row1_h; ROW2_TOP, ROW2_BOTTOM = ROW1_BOTTOM - VGAP, BOTTOM
    COL_A = (LEFT, LEFT + w_a); COL_B = (0.386, 0.566); COL_C = (0.650, 0.935)

    def plabel(xl, yt, t):
        fig.text(xl - 0.028, yt + 0.012, t, fontsize=10, fontweight='bold', va='bottom', ha='left', color=C_TEXT)

    # ------------------------------------------------------------------ a  matrix schematic (true CAD cross-sections)
    gs_A = fig.add_gridspec(2, 2, left=COL_A[0], right=COL_A[1], top=ROW1_TOP, bottom=ROW1_BOTTOM, width_ratios=[1, 4.5],
                            height_ratios=[1, 4.5], wspace=0.03, hspace=0.03)
    ax_m = fig.add_subplot(gs_A[1, 1]); ax_t = fig.add_subplot(gs_A[0, 1], sharex=ax_m); ax_l = fig.add_subplot(gs_A[1, 0], sharey=ax_m)
    for sp in ax_m.spines.values(): sp.set_visible(False)
    ax_m.set_xticks([]); ax_m.set_yticks([])
    GC = np.arange(1, 6) * 1000; MIN_C = 55
    ax_m.add_patch(patches.FancyBboxPatch((0, 0), 6000, 6000, boxstyle='square,pad=0', facecolor=C_PDMS_BG, edgecolor='#BBBBBB', lw=0.5, zorder=0))
    cw = list(CTRL_W.values()); fh = list(FLOW_H.values()); fw = [_flow_cad_silhouette(h)[2] for h in fh]
    for y, w in zip(GC, cw):
        ax_m.add_patch(patches.Rectangle((0, y - max(w, MIN_C) / 2), 6000, max(w, MIN_C), lw=0, facecolor=C_CTRL, alpha=0.82, zorder=1))
    for x, wf in zip(GC, fw):
        ax_m.add_patch(patches.Rectangle((x - wf / 2, 0), wf, 6000, lw=0, facecolor=C_FLOW, alpha=0.82, zorder=2))
        for y, w in zip(GC, cw):
            ax_m.add_patch(patches.Rectangle((x - wf / 2, y - max(w, MIN_C) / 2), wf, max(w, MIN_C), lw=0, facecolor=C_FLOW, alpha=0.42, zorder=3))
    BSZ = 680
    ax_m.add_patch(patches.Rectangle((4000 - BSZ / 2, 1000 - BSZ / 2), BSZ, BSZ, fill=False, edgecolor=C_NORMAL, lw=1.2, zorder=10))
    ax_m.text(4000 + BSZ / 2 + 140, 1000, 'Normal\n(b)', color=C_NORMAL, fontsize=4.8, va='center', ha='left', fontweight='bold', zorder=10)
    ax_m.add_patch(patches.Rectangle((1000 - BSZ / 2, 5000 - BSZ / 2), BSZ, BSZ, fill=False, edgecolor=C_BRIDGE, lw=1.2, zorder=10))
    ax_m.text(1000 + BSZ / 2 + 140, 5000, 'Bridge\n(b)', color=C_BRIDGE, fontsize=4.8, va='center', ha='left', fontweight='bold', zorder=10)
    ax_m.set_xlim(0, 6000); ax_m.set_ylim(6000, 0)
    ax_m.text(3000, 6400, r'$H_f$ (µm)', color=C_FLOW, fontsize=BASE_FS, fontweight='bold', ha='center', va='top')
    ax_m.text(1.12, 0.5, r'$W_c$ (µm)', color=C_CTRL, fontsize=BASE_FS, fontweight='bold', va='center', ha='center', rotation=90,
              transform=ax_m.transAxes, clip_on=False)
    sb0 = 6000 - 500 - 600; sb1 = sb0 + 500; sby = 5200; sbh = 85
    ax_m.plot([sb0, sb1], [sby, sby], color=C_TEXT, lw=1.3, solid_capstyle='butt', zorder=12)
    ax_m.plot([sb0, sb0], [sby - sbh, sby + sbh], color=C_TEXT, lw=1.0, zorder=12); ax_m.plot([sb1, sb1], [sby - sbh, sby + sbh], color=C_TEXT, lw=1.0, zorder=12)
    ax_m.text((sb0 + sb1) / 2, sby + 200, '500 µm', ha='center', va='top', fontsize=5.0, color=C_TEXT, zorder=12)
    ax_t.axis('off'); ax_l.axis('off')
    ax_t.plot([-300, 5300], [0, 0], color=C_TEXT, lw=0.8); ax_l.plot([0, 0], [-300, 5300], color=C_TEXT, lw=0.8)
    ax_t.set_ylim(0, 430); ax_l.set_xlim(320, 0)
    for xc, hv, lb in zip(GC, fh, FLOW_H.keys()):
        xs, zs, _ = _flow_cad_silhouette(hv); ax_t.fill_between(xs + xc, 0, zs, color=C_FLOW, alpha=0.65, lw=0)
        ax_t.text(xc, hv + 25, lb, color=C_FLOW, fontsize=5.5, ha='center', va='bottom', fontweight='bold')
    for yc, wv, lb in zip(GC, cw, CTRL_W.keys()):
        bw = max(wv, 30); ax_l.fill_betweenx([yc - bw / 2, yc - bw / 2, yc + bw / 2, yc + bw / 2], 0, [0, 190, 190, 0], color=C_CTRL, alpha=0.65)
        ax_l.text(220, yc, lb, color=C_CTRL, fontsize=5.5, ha='right', va='center', fontweight='bold')
    plabel(COL_A[0], ROW1_TOP, 'a')

    # ------------------------------------------------------------------ b  normal vs bridge
    ax_B = fig.add_axes([COL_B[0], ROW1_BOTTOM, COL_B[1] - COL_B[0], ROW1_TOP - ROW1_BOTTOM])
    ax_B.axhline(0, color='#AAAAAA', lw=0.4, zorder=1)
    for (F, C), col, lab, z in ((('F4', 'C1'), C_NORMAL, r'Normal ($H_f$=50, $W_c$=500 µm)', 4), (('F1', 'C5'), C_BRIDGE, r'Bridge ($H_f$=350, $W_c$=25 µm)', 3)):
        g = G[(G.F == F) & (G.C == C)].sort_values('P')
        ax_B.errorbar(g.P, g.closure, yerr=g.sd, color=col, lw=1.0, marker='o', ms=3.5, capsize=2, elinewidth=0.7, zorder=z, label=lab)
    ax_B.axhline(SEAL, color=C_FLOW, ls=':', lw=0.8, zorder=1)
    ax_B.set_xlim(400, 1600); ax_B.set_ylim(-15, 115); ax_B.set_xticks(PRESSURES)
    ax_B.set_xlabel('Actuation pressure (mbar)', labelpad=3); ax_B.set_ylabel('% Flow channel closed', labelpad=3)
    ax_B.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', zorder=0)
    ax_B.legend(loc='lower right', bbox_to_anchor=(0.98, 0.20), frameon=False, fontsize=4.8, handlelength=1.6, handletextpad=0.5, labelspacing=0.32, borderpad=0)
    plabel(COL_B[0], ROW1_TOP, 'b')

    # ------------------------------------------------------------------ c  threshold pressure vs as-fabricated aspect ratio
    # Only tested pressures appear on the y axis (500-1500 mbar). A filled marker is the lowest tested pressure at which the geometry
    # reached >= 95% closure; an open marker (drawn at 1500 mbar, the last pressure tested) is a geometry that never did.
    ax_C = fig.add_axes([COL_C[0], ROW1_BOTTOM, COL_C[1] - COL_C[0], ROW1_TOP - ROW1_BOTTOM])
    hcol = {25: PRESSURE_COLORS[500], 50: PRESSURE_COLORS[750], 100: PRESSURE_COLORS[1000], 200: PRESSURE_COLORS[1250], 350: PRESSURE_COLORS[1500]}
    hedge = {25: '#A8893A', 50: '#A8893A', 100: '#5E7F4F', 200: '#1F6F6B', 350: '#0F3F3F'}
    hmark = {350: 'o', 200: 's', 100: '^', 50: 'D', 25: 'v'}
    ax_C.axvspan(0.25, 1.0, color='#EFEFEF', zorder=0); ax_C.axvline(1.0, color='#666666', ls='--', lw=1.0, zorder=1)
    for _, r in T.iterrows():
        sealed = not r.not_sealed_1500; y = r.P_min if sealed else 1500; ec = hedge[r.Hf]
        ax_C.errorbar(r.AR_fab, y, xerr=r.AR_fab_sd, fmt='none', ecolor=ec, elinewidth=0.4, alpha=0.45, capsize=0, zorder=2)
        ax_C.scatter(r.AR_fab, y, s=22, marker=hmark[r.Hf], facecolor=hcol[r.Hf] if sealed else 'white', edgecolors=ec, lw=0.9, zorder=4)
    ax_C.text(0.50, 800, 'Bridge', color='#555555', fontsize=5.6, ha='center', va='center', fontstyle='italic')
    ax_C.text(14.0, 800, 'Normal', color='#555555', fontsize=5.6, ha='center', va='center', fontstyle='italic')
    ax_C.set_xscale('log'); ax_C.set_xlim(0.25, 42); ax_C.set_ylim(400, 1600)
    ax_C.set_xticks([0.5, 1, 2, 5, 10, 20]); ax_C.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v_, _: '{:g}'.format(v_)))
    ax_C.tick_params(axis='x', which='minor', bottom=False)
    ax_C.set_yticks(PRESSURES); ax_C.set_yticklabels([str(p) for p in PRESSURES])
    ax_C.set_xlabel(r'$W_{c,fab}\,/\,H_{f,fab}$', labelpad=2); ax_C.set_ylabel(r'$P_{min}$ (mbar)', labelpad=3)
    ax_C.grid(True, which='major', axis='both', ls=':', lw=0.4, alpha=0.40, zorder=0)
    # one compact legend in the empty upper-right of the axes
    h1 = [Line2D([], [], marker=hmark[h], ls='', mfc=hcol[h], mec=hedge[h], mew=0.9, ms=3.8, label=f'{h}') for h in (350, 200, 100, 50, 25)]
    h1.append(Line2D([], [], marker='o', ls='', mfc='white', mec='#555555', mew=0.9, ms=3.8, label='no seal'))
    ax_C.legend(handles=h1, title=r'$H_f$ (µm)', title_fontsize=5.6, loc='upper right', bbox_to_anchor=(0.985, 0.985), frameon=False, fontsize=5.4,
                handlelength=1.0, handletextpad=0.4, labelspacing=0.4, borderpad=0.1, alignment='left')
    plabel(COL_C[0], ROW1_TOP, 'c')

    # ------------------------------------------------------------------ d  closure maps
    gs_D = fig.add_gridspec(1, 5, left=LEFT, right=RIGHT, top=ROW2_TOP, bottom=ROW2_BOTTOM, wspace=0.10)
    axes_D = [fig.add_subplot(gs_D[0, i]) for i in range(5)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('w2gold', [(0, '#ffffff'), (0.20, '#f2e6b8'), (0.55, '#dfc46a'), (1.0, '#d9c077')], N=512)
    Hs, Ws = sorted(FLOW_H.values()), sorted(CTRL_W.values())
    ar = {(r.Hf, r.Wc): r.AR_fab for _, r in T.iterrows()}
    for idx, p in enumerate(PRESSURES):
        ax = axes_D[idx]; ax.set_facecolor('#e0e0e0'); sub = G[G.P == p]
        M = np.full((5, 5), np.nan)
        for i, h in enumerate(Hs):
            for j, w in enumerate(Ws):
                v = sub[(sub.Hf == h) & (sub.Wc == w)].closure
                if len(v): M[i, j] = v.iloc[0]
        ax.imshow(np.clip(M, 0, 100), cmap=cmap, vmin=0, vmax=100, origin='lower', extent=(0, 5, 0, 5), aspect='equal', zorder=1)
        for k in range(6):
            ax.axhline(k, color='#D0D0D0', lw=0.4, zorder=2); ax.axvline(k, color='#D0D0D0', lw=0.4, zorder=2)
        ax.add_patch(patches.Rectangle((0, 0), 5, 5, fill=False, edgecolor='#3A3A3A', lw=0.65, clip_on=False, zorder=5))
        for i in range(5):
            for j in range(5):
                v = M[i, j]
                if np.isnan(v):
                    ax.add_patch(patches.Rectangle((j, i), 1, 1, fc='#e0e0e0', ec='none', hatch='////', zorder=3))
                elif v >= SEAL:
                    ax.add_patch(patches.Rectangle((j, i), 1, 1, fill=False, edgecolor=C_FLOW, lw=0.9, zorder=4))
        # boundary W_c,fab = H_f,fab (cells with AR < 1 on the left/below), drawn on the cell edges
        for i, h in enumerate(Hs):
            for j, w in enumerate(Ws):
                a = ar[(h, w)]
                if j < 4 and (a < 1) != (ar[(h, Ws[j + 1])] < 1):
                    ax.plot([j + 1, j + 1], [i, i + 1], color=C_TEXT, ls='--', lw=0.9, zorder=6, alpha=0.8)
                if i < 4 and (a < 1) != (ar[(Hs[i + 1], w)] < 1):
                    ax.plot([j, j + 1], [i + 1, i + 1], color=C_TEXT, ls='--', lw=0.9, zorder=6, alpha=0.8)
        ax.set_xlim(0, 5); ax.set_ylim(0, 5)
        ax.set_title(f'{p} mbar', fontsize=BASE_FS, pad=4, fontweight='bold', color=C_TEXT)
        ax.set_xticks(np.arange(5) + 0.5); ax.set_yticks(np.arange(5) + 0.5); ax.minorticks_off(); ax.tick_params(axis='both', which='major', length=2, pad=2)
        ax.set_xticklabels([str(w) for w in Ws])
        if idx == 0: ax.set_yticklabels([str(h) for h in Hs]); ax.set_ylabel(r'$H_f$ (µm)', labelpad=4)
        else: ax.set_yticklabels([])
        ax.set_xlabel(r'$W_c$ (µm)' if idx == 2 else '', labelpad=4)
        for sp in ax.spines.values(): sp.set_visible(False)
    fig.canvas.draw()
    ref = axes_D[2].get_position(); cw_, ch_ = ref.width, ref.height
    for ax in axes_D:
        pos = ax.get_position(); ax.set_position([pos.x0 + (pos.width - cw_) / 2, ROW2_BOTTOM + ((ROW2_TOP - ROW2_BOTTOM) - ch_) / 2, cw_, ch_])
    plabel(LEFT, axes_D[0].get_position().y1, 'd')
    pl, pf = axes_D[4].get_position(), axes_D[0].get_position()
    cax = fig.add_axes([pl.x1 + 0.008, pf.y0, 0.010, pf.height])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize(0, 100)); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax); cb.set_label('Closure (%)', fontsize=BASE_FS - 1, labelpad=4, fontweight='bold')
    cb.set_ticks([0, 25, 50, 75, 100]); cb.ax.tick_params(labelsize=BASE_FS - 2, length=2, width=0.5)
    cb.ax.axhline(y=95, color=C_FLOW, lw=1.2, xmin=0.0, xmax=0.55, zorder=5); cb.outline.set_linewidth(0.5)
    fig.legend(handles=[Line2D([0], [0], color='none', marker='s', markersize=7, markeredgecolor=C_FLOW, markerfacecolor='#d9c077', markeredgewidth=0.9, label=r'$\geq\!95\%$ Closure'),
                        Line2D([0], [0], color=C_TEXT, ls='--', lw=1.0, label=r'$W_{c,fab} = H_{f,fab}$')],
               loc='lower center', bbox_to_anchor=(0.50, 0.125), ncol=2, frameon=False, fontsize=BASE_FS, handlelength=1.6, handletextpad=0.5, columnspacing=1.5)
    return fig


# =========================================================================================== SUPPLEMENTARY FIGURES
def figSV_traces(data_dir=DATA_DIR, extra=EXTRA_FILES):
    """Flow of every run (normalised to its open plateau) with the analysis windows; rows F1-F5, columns C1-C5."""
    files = run_files(data_dir, extra)
    fig, axs = plt.subplots(5, 5, figsize=(7.2, 8.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.935, bottom=0.065, wspace=0.07, hspace=0.12)
    for i, F in enumerate(FLOW_H):
        for j, C in enumerate(CTRL_W):
            ax = axs[i, j]
            for a, b in OPEN_WIN: ax.axvspan(a, b, color=C_FLOW, alpha=0.10, lw=0)
            for a, b in CLOSED_WIN: ax.axvspan(a, b, color=C_CTRL, alpha=0.18, lw=0)
            for P in PRESSURES:
                df, _ = load(files[(F, C, P)]); qo = np.mean([df.q[(df.t > a) & (df.t < b)].mean() for a, b in OPEN_WIN])
                ax.plot(df.t, df.q / qo, color=PRESSURE_COLORS[P], lw=0.45, alpha=0.95)
            ax.axhline(0, color='#AAAAAA', lw=0.3); ax.set_xlim(0, 180); ax.set_ylim(-0.25, 1.45)
            ax.text(0.97, 0.96, f'{F}×{C}\n{FLOW_H[F]}/{CTRL_W[C]} µm', transform=ax.transAxes, fontsize=4.8, va='top', ha='right')
            ax.tick_params(labelsize=BASE_FS - 2.2, length=2, pad=1)
            if j == 0: ax.set_ylabel(r'$Q\,/\,Q_{open}$', fontsize=BASE_FS - 0.8)
            if i == 4: ax.set_xlabel('Time (s)', fontsize=BASE_FS - 0.8)
    hs = [Line2D([], [], color=PRESSURE_COLORS[p], lw=1.0, label=f'{p} mbar') for p in PRESSURES]
    hs += [patches.Patch(fc=C_FLOW, alpha=0.2, label='open window'), patches.Patch(fc=C_CTRL, alpha=0.3, label='closed window')]
    fig.legend(handles=hs, loc='upper center', ncol=7, frameon=False, fontsize=BASE_FS - 1.4, bbox_to_anchor=(0.54, 0.995), handlelength=1.5, columnspacing=1.0)
    return fig


def figSV_cycles(D):
    """Closure vs pressure for all 25 geometries: the three cycles (open) and their mean (line)."""
    V = D[D.valid]
    fig, axs = plt.subplots(5, 5, figsize=(7.2, 8.0), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.935, bottom=0.065, wspace=0.07, hspace=0.12)
    for i, F in enumerate(FLOW_H):
        for j, C in enumerate(CTRL_W):
            ax = axs[i, j]; g = V[(V.F == F) & (V.C == C)]
            ax.axhline(SEAL, color=C_FLOW, lw=0.5, ls=':'); ax.axhline(0, color='#AAAAAA', lw=0.3)
            mean = g.groupby('P').closure.mean()
            ax.plot(mean.index, mean.values, color=C_TEXT, lw=0.6, zorder=2)
            for P in PRESSURES:
                gp = g[g.P == P]
                ax.scatter(np.full(len(gp), P) + np.array([-25, 0, 25])[:len(gp)], gp.closure, s=6, marker=PRESSURE_MARKERS[P], facecolor='none',
                           edgecolor=PRESSURE_COLORS[P] if P != 500 else '#A8893A', lw=0.6, zorder=3)
            ax.set_xlim(350, 1650); ax.set_ylim(-15, 108); ax.set_xticks([500, 1000, 1500])
            ax.text(0.04, 0.05, f'{F}×{C}', transform=ax.transAxes, fontsize=4.8, va='bottom', ha='left')
            ax.tick_params(labelsize=BASE_FS - 2.2, length=2, pad=1)
            if j == 0: ax.set_ylabel('Closure (%)', fontsize=BASE_FS - 0.8)
            if i == 4: ax.set_xlabel('Pressure (mbar)', fontsize=BASE_FS - 0.8)
    fig.legend(handles=[Line2D([], [], marker='o', ls='', mfc='none', mec=C_TEXT, ms=3.4, label='single cycle'), Line2D([], [], color=C_TEXT, lw=0.8, label='mean of 3 cycles'),
                        Line2D([], [], color=C_FLOW, lw=0.6, ls=':', label='95% criterion')], loc='upper center', ncol=3, frameon=False,
               fontsize=BASE_FS - 1.3, bbox_to_anchor=(0.54, 0.995))
    return fig


def figSV_ar(G, T, prof_path=PROF_PER_FILE):
    """All 125 closure values against the design and the as-fabricated aspect ratio; sensitivity of the sealing bracket to the
    width definition used for W_c,fab."""
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, left=0.075, right=0.975, top=0.94, bottom=0.085, wspace=0.28, hspace=0.45)   # right margin: panel-d label was clipped at 0.99
    A = add_aspect_ratio(T, prof_path, 'W50'); arf = {(r.Hf, r.Wc): r for _, r in A.iterrows()}
    fmt = ticker.FuncFormatter(lambda v, _: '{:g}'.format(v))
    for k, (key, ttl) in enumerate((('design', 'Design aspect ratio'), ('fab', 'As-fabricated aspect ratio (W50)'))):
        ax = fig.add_subplot(gs[0, k]); rng = np.random.default_rng(3)
        for p in PRESSURES:
            s = G[G.P == p]
            x = np.array([arf[(r.Hf, r.Wc)].AR_design if key == 'design' else arf[(r.Hf, r.Wc)].AR_fab for _, r in s.iterrows()])
            ax.scatter(x * 10 ** rng.uniform(-0.03, 0.03, len(x)), s.closure, s=14, marker=PRESSURE_MARKERS[p], color=PRESSURE_COLORS[p],
                       edgecolors='white', lw=0.3, alpha=0.9, zorder=3, label=f'{p}')
        ax.axhline(SEAL, color=C_FLOW, ls=':', lw=0.8); ax.axvline(1, color='#555555', ls='--', lw=1.0)
        ax.set_xscale('log'); ax.set_xlim(0.05, 30); ax.set_ylim(-12, 108); ax.set_xticks([0.1, 0.5, 1, 2, 5, 20])
        ax.xaxis.set_major_formatter(fmt); ax.tick_params(axis='x', which='minor', bottom=False)
        ax.set_xlabel(r'$W_c\,/\,H_f$' if key == 'design' else r'$W_{c,fab}\,/\,H_{f,fab}$'); ax.set_ylabel('Closure (%)')
        ax.set_title(ttl, fontsize=BASE_FS, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'ab'[k], x=-0.17)
        if k == 0:
            ax.legend(title='Pressure (mbar)', title_fontsize=5.3, fontsize=5.0, loc='center left', bbox_to_anchor=(0.0, 0.55), ncol=1, frameon=True,
                      framealpha=0.92, edgecolor='#CCCCCC', handletextpad=0.2, labelspacing=0.2)
    wds = ('W10', 'W20', 'W50', 'W80')
    ax = fig.add_subplot(gs[1, 0])
    for n, wd in enumerate(wds):
        B = add_aspect_ratio(T, prof_path, wd); s_ = ~B.not_sealed_1500
        ax.scatter(B.AR_fab[s_], np.full(s_.sum(), n) + 0.13, marker='o', color=C_FLOW, s=9)
        ax.scatter(B.AR_fab[~s_], np.full((~s_).sum(), n) - 0.13, marker='o', facecolor='none', edgecolor=C_BRIDGE, s=9)
    ax.set_xscale('log'); ax.set_xlim(0.3, 30); ax.set_yticks(range(4)); ax.set_yticklabels(wds); ax.invert_yaxis()
    ax.set_xticks([0.5, 1, 2, 5, 20]); ax.xaxis.set_major_formatter(fmt); ax.tick_params(axis='x', which='minor', bottom=False)
    ax.set_xlabel(r'$W_{c,fab}\,/\,H_{f,fab}$'); ax.set_ylabel(r'Width used for $W_{c,fab}$')
    ax.set_title('Sealed (teal) and not sealed (red, open) at 1500 mbar', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'c', x=-0.17)
    S = separation(T, prof_path=prof_path); ax = fig.add_subplot(gs[1, 1]); k = 0; labs = []
    for _, r in S.iterrows():
        sep = r.separable; col = C_FLOW if sep else C_BRIDGE
        ax.plot([r.largest_unsealed_AR, r.smallest_sealed_AR], [k, k], color=col, lw=2.2, solid_capstyle='butt')
        ax.scatter([r.largest_unsealed_AR], [k], marker='|', color=C_BRIDGE, s=34, zorder=3); ax.scatter([r.smallest_sealed_AR], [k], marker='|', color=C_TEXT, s=34, zorder=3)
        labs.append(f"{r.width_def}, " + ('all $H_f$' if r.subset == 'all' else r'$H_f$ ≤ 200 µm')); k += 1
    ax.set_yticks(range(k)); ax.set_yticklabels(labs, fontsize=BASE_FS - 1.6); ax.invert_yaxis(); ax.axvline(1, color='#555555', ls='--', lw=0.8)
    ax.set_xscale('log'); ax.set_xlim(0.5, 2.2); ax.set_xticks([0.5, 0.7, 1, 1.5, 2.0]); ax.xaxis.set_major_formatter(fmt); ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel('Aspect ratio'); ax.set_title('Bracket where sealing starts (teal) or overlap (red)', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'd', x=-0.36)   # tick meaning is in the caption
    return fig


def figSV_qc(G, Q, T):
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.55], left=0.07, right=0.99, top=0.93, bottom=0.04, wspace=0.40, hspace=0.62)
    ax = fig.add_subplot(gs[0, 0]); q = Q.copy(); q['bad'] = ~q.valid
    ax.scatter(q.P[~q.bad], q.readback[~q.bad], s=8, color=C_FLOW, label='used'); ax.scatter(q.P[q.bad], q.readback[q.bad], s=22, marker='x', color=C_BRIDGE, label='excluded')
    ax.plot([400, 1600], [400, 1600], color=C_TEXT, lw=0.5, ls=':'); ax.set_xlabel('File-name pressure (mbar)'); ax.set_ylabel('Logged read-back (mbar)')
    ax.legend(fontsize=5.5, loc='upper left'); ax.set_title('Setpoint check', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'a', x=-0.25)
    ax = fig.add_subplot(gs[0, 1]); s = G.sort_values('Q_open').reset_index(drop=True)
    ax.scatter(range(len(s)), s.Q_open, s=6, color=C_FLOW); ax.axhline(SENSOR_CEILING, color=C_BRIDGE, lw=0.7, ls='--')
    ax.text(2, SENSOR_CEILING + 0.15, 'logging ceiling 8.00 µL/min', fontsize=5.0, color=C_BRIDGE)
    ax.set_ylim(0, 9); ax.set_xlabel('Run (sorted)'); ax.set_ylabel(r'Open-phase flow (µL/min)'); ax.set_title('Flow-sensor ceiling', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'b', x=-0.25)
    ax = fig.add_subplot(gs[0, 2]); cols = {'last 10 s of each phase': C_CTRL, 'conservative (15 s margin)': C_FLOW}
    for n, c in cols.items():
        ax.scatter(G.closure, G[f'closure_{n}'], s=7, color=c, alpha=0.8, label=n)
    ax.plot([-10, 105], [-10, 105], color=C_TEXT, lw=0.5, ls=':'); ax.axhline(SEAL, color=C_TEXT, lw=0.4, ls=':'); ax.axvline(SEAL, color=C_TEXT, lw=0.4, ls=':')
    ax.set_xlabel('Closure, primary windows (%)'); ax.set_ylabel('Closure, alternative windows (%)'); ax.legend(fontsize=4.8, loc='lower right')
    ax.set_title('Window sensitivity (125 runs)', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'c', x=-0.25)
    ax = fig.add_subplot(gs[1, :])
    rows = []
    for n, col in (('primary', 'closure'), ('last 10 s of each phase', 'closure_last 10 s of each phase'), ('conservative (15 s margin)', 'closure_conservative (15 s margin)')):
        P = pmin_table(G, col); P['set'] = n; rows.append(P)
    PM = pd.concat(rows); piv = PM.pivot_table(index=['Hf', 'Wc'], columns='set', values='P_min', aggfunc='first')
    diff = piv[(piv.nunique(axis=1, dropna=False) > 1)]
    ax.axis('off')
    ax.set_title('Geometries whose P_min depends on the window set (mbar)', fontsize=BASE_FS - 0.8, loc='left', pad=4, fontweight='bold'); panel_label(ax, 'd', x=-0.02, y=1.0)
    cell = [[f'{int(i[0])} / {int(i[1])}'] + ['no seal' if np.isnan(r[c]) else f'{int(r[c])}' for c in ('primary', 'last 10 s of each phase', 'conservative (15 s margin)')] for i, r in diff.iterrows()]
    if cell:
        tb = ax.table(cellText=cell, colLabels=['H_f / W_c (µm)', 'primary windows', 'last 10 s of phase', 'conservative'], loc='upper left', cellLoc='center', colWidths=[0.2, 0.2, 0.2, 0.2])
        tb.auto_set_font_size(False); tb.set_fontsize(5.6); tb.scale(1, 1.0)
        for (r_, c_), cellobj in tb.get_celld().items():
            cellobj.set_linewidth(0.3); cellobj.set_edgecolor('#BBBBBB')
            if r_ == 0: cellobj.set_text_props(fontweight='bold'); cellobj.set_facecolor('#F0F0F0')
    return fig


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    D, Q, G = analyse()
    T = add_aspect_ratio(pmin_table(G))
    D.to_csv(os.path.join(OUT_DIR, 'valve_closure_per_cycle.csv'), index=False); Q.to_csv(os.path.join(OUT_DIR, 'valve_file_QC.csv'), index=False)
    G.to_csv(os.path.join(OUT_DIR, 'valve_closure_per_run.csv'), index=False); T.to_csv(os.path.join(OUT_DIR, 'valve_Pmin_aspect_ratio.csv'), index=False)
    separation(T).to_csv(os.path.join(OUT_DIR, 'aspect_ratio_separation.csv'), index=False)
    f4 = fig4(G, T); s1 = figSV_traces(); s2 = figSV_cycles(D); s3 = figSV_ar(G, T); s4 = figSV_qc(G, Q, T)
    # ---- figure saving (uncomment to export)
    # f4.savefig('Fig4_valve_closure.jpeg', dpi=600, bbox_inches='tight')
    # s1.savefig('FigS_valve_traces.jpeg', dpi=600, bbox_inches='tight')
    # s2.savefig('FigS_valve_closure_vs_pressure.jpeg', dpi=600, bbox_inches='tight')
    # s3.savefig('FigS_valve_aspect_ratio.jpeg', dpi=600, bbox_inches='tight')
    # s4.savefig('FigS_valve_QC.jpeg', dpi=600, bbox_inches='tight')
    plt.show()
