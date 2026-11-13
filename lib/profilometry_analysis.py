"""
Profilometry of laser-SLA molds and PDMS replicas (KLA Zeta-20): processing, statistics and figures.

INPUT   60 CSV exports  'C1 CM1.csv', 'C1 CM1 PDMS.csv', ... 'F5 FM3 PDMS.csv' in DATA_DIR.
        col 0 = row index, col 1 = lateral position x (um), cols 2.. = individual scan lines z(x) (um), 62-63 lines.
        CM1-3 / FM1-3 = three separately printed control / flow mold chips; the PDMS replica with the same number was cast
        from that chip (confirmed by the experimenter). Each chip carries all five features (C1-C5 or F1-F5).
DESIGN  Control C1-C5: rectangles, H = 200 um, W = 500/300/200/50/25 um.
        Flow F1-F5: 'three-point arc' (Fusion 360) through the two ends of a 350 um base chord and an apex of height
        H = 350/200/100/50/25 um, extruded -> circular segment of chord 350 um and sagitta H. R = H/2 + 350^2/(8H).
        For H > 175 um (F1, F2) this is a MAJOR arc: the circle is wider than the chord (F1: 437.5 um at 131 um height),
        i.e. F1 is a bulb with an overhang, F2 is nearly a semicircle. Seen from above (what a profiler measures) the
        design 'silhouette' of F1 is therefore 437.5 um wide, not 350 um. All offsets below are taken against this
        silhouette at the SAME threshold fraction as the measured width.

PROCESSING (per scan line -> per file -> per chip; nothing is synthesised, replaced or interpolated into the data)
  1  PDMS replicas (valleys) are inverted so that every feature is a raised profile.
  2  Levelling: least-squares line through the substrate on BOTH sides of the feature (more than 15 um outside the
     footprint = contiguous region above 5% of the peak of a 61-point running median), subtracted; two iterations.
     If one wall does not reach 5% of the peak (one replica) the footprint is taken at 15%; flagged in the QC columns.
  3  Height H = median of the levelled raw points over the central 20% of the half-height width (plateau of a rectangle,
     apex of a rounded feature). A maximum would pick up sidewall spikes of the PDMS replicas and noise.
  4  Widths W10, W20, W50, W80 = distance between the outermost crossings of 10/20/50/80% of H around the feature centre,
     on a 15-point (~10 um) running median (edge preserving), linear interpolation (pixel pitch 0.685 um).
     W10 = footprint (primary lateral measure); W50 = full width at half maximum (primary width for aspect ratios).
  5  Lines are measured one by one (the feature centre drifts up to ~20 um across a scan) and the file value is the
     median over lines. Lines deviating from the file median by more than max(5 scaled MAD, 5%) are excluded.
  6  The UNIT OF REPLICATION IS THE CHIP (n = 3). Scan lines and the five features on one chip are sub-samples of one
     print. Per feature: mean, SD, t-based 95% CI over the 3 chips. Pooled quantities use chip means (df = 2). Dependence
     on height uses a randomised-block regression  dW = a*H + chip effect  (chip as fixed block; balanced design).

OUTPUT  tables in results/, figure functions at the bottom (the lines that save figures are commented out).
"""
import os, re
import numpy as np, pandas as pd
from scipy import stats
from scipy.ndimage import median_filter
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from ch2_style import *

DATA_DIR = r'./Data'
OUT_DIR = r'./results'
CAD = {'C1': (200, 500), 'C2': (200, 300), 'C3': (200, 200), 'C4': (200, 50), 'C5': (200, 25),
       'F1': (350, 350), 'F2': (200, 350), 'F3': (100, 350), 'F4': (50, 350), 'F5': (25, 350)}   # (H, W or chord) um
CHANNELS = list(CAD)
FRACS = {10: 0.10, 20: 0.20, 50: 0.50, 80: 0.80}
GRID = np.arange(-480, 480.5, 1.0)


# =========================================================================================== CAD geometry
def cad_profile(ch, x):
    """Design silhouette z(x) seen from above (um). Rectangles: H over |x| <= W/2. Flow: circular segment through
    (-W/2, 0), (0, H), (W/2, 0); a major arc (H > W/2) extends to |x| = R and drops vertically to the substrate."""
    H, W = CAD[ch]
    z = np.zeros_like(x, dtype=float)
    if ch.startswith('C'):
        z[np.abs(x) <= W / 2] = H
    else:
        R = H / 2 + W ** 2 / (8 * H); cy = H - R
        xmax = R if H > W / 2 else W / 2
        m = np.abs(x) <= xmax
        z[m] = np.maximum(cy + np.sqrt(np.maximum(R ** 2 - x[m] ** 2, 0)), 0)
    return z


def _cad_width(ch, f):
    H, W = CAD[ch]
    if ch.startswith('C'):
        return float(W)
    x = np.arange(-300, 300.0001, 0.02); z = cad_profile(ch, x)
    m = z >= f * H
    return float(x[m].max() - x[m].min())


def _cad_area(ch):
    x = np.arange(-300, 300.0001, 0.02); return float(np.trapezoid(cad_profile(ch, x), x))


CAD_A = {ch: _cad_area(ch) for ch in CHANNELS}                                   # design cross-sectional area (um^2)
CAD_W = {ch: {k: _cad_width(ch, f) for k, f in FRACS.items()} for ch in CHANNELS}      # silhouette width at each threshold
CAD_MAXW = {ch: (CAD[ch][1] if ch.startswith('C') else _cad_width(ch, 1e-6)) for ch in CHANNELS}


# =========================================================================================== per-line processing
def contiguous(mask, i):
    if not mask[i]:
        return None
    L = i
    while L > 0 and mask[L - 1]: L -= 1
    R = i
    while R < len(mask) - 1 and mask[R + 1]: R += 1
    return L, R


def crossing(x, z, level, L, R):
    xl = x[L] if L == 0 else x[L - 1] + (level - z[L - 1]) * (x[L] - x[L - 1]) / (z[L] - z[L - 1])
    xr = x[R] if R == len(x) - 1 else x[R] + (level - z[R]) * (x[R + 1] - x[R]) / (z[R + 1] - z[R])
    return xl, xr


def level_line(x, z):
    z = z - np.median(np.sort(z)[:len(z) // 5])
    info = {'two_sided': False, 'tilt_um_per_mm': np.nan, 'sub_rms': np.nan, 'footprint_15': False}
    for it in range(2):
        zb = median_filter(z, size=61, mode='nearest')
        i = int(np.argmax(zb)); h = zb[i]
        ok = False
        for frac in (0.05, 0.15):
            reg = contiguous(zb > frac * h, i)
            if reg is None:
                return None, info
            L, R = reg
            sl, sr = x < x[L] - 15, x > x[R] + 15
            if sl.sum() >= 40 and sr.sum() >= 40:
                ok = True; info['footprint_15'] |= (frac == 0.15); break
        if ok:
            p = np.polyfit(x[sl | sr], z[sl | sr], 1); info['two_sided'] = True
            if it == 0: info['tilt_um_per_mm'] = p[0] * 1000
            z = z - np.polyval(p, x)
        elif (sl | sr).sum() >= 40:
            z = z - np.median(z[sl | sr])
        else:
            return None, info
        info['sub_rms'] = float(np.std(z[sl | sr]))
    return z, info


def measure_line(x, z):
    zb = median_filter(z, size=61, mode='nearest')
    i = int(np.argmax(zb)); hb = zb[i]
    reg = contiguous(zb > 0.5 * hb, i)
    if reg is None or hb <= 0:
        return None
    a, b = crossing(x, zb, 0.5 * hb, *reg)
    xc0, w0 = 0.5 * (a + b), b - a
    core = np.abs(x - xc0) <= 0.1 * w0
    if core.sum() < 5:
        return None
    H = float(np.median(z[core]))
    zs = median_filter(z, size=15, mode='nearest')
    j0 = int(np.argmin(np.abs(x - xc0)))
    out = {'H': H}
    reg5 = contiguous(zs >= 0.05 * H, j0)                          # cross-sectional area of the levelled raw profile
    if reg5 is not None:
        a5, b5 = crossing(x, zs, 0.05 * H, *reg5); mm = (x >= a5) & (x <= b5)
        out['A'] = float(np.trapezoid(z[mm], x[mm]))
    else:
        out['A'] = np.nan
    for k, f in FRACS.items():
        reg = contiguous(zs >= f * H, j0)
        if reg is None:
            out[f'W{k}'] = np.nan; continue
        l, r = crossing(x, zs, f * H, *reg)
        out[f'W{k}'] = r - l
        if k == 50: out['xc'] = 0.5 * (l + r)
    out.setdefault('xc', xc0)
    return out


def process_file(path, invert):
    df = pd.read_csv(path, header=None)
    x = df.iloc[:, 1].values.astype(float); Z = df.iloc[:, 2:].values.astype(float)
    rows, lev = [], []
    for k in range(Z.shape[1]):
        z = Z[:, k]
        if np.isnan(z).any():
            continue
        zl, info = level_line(x, -z if invert else z)
        if zl is None:
            continue
        m = measure_line(x, zl)
        if m is None:
            continue
        m.update(info); m['line'] = k; rows.append(m); lev.append(zl)
    L = pd.DataFrame(rows)
    keep = pd.Series(True, index=L.index)
    for c in ('H', 'W10', 'W50'):
        med = L[c].median(); mad = 1.4826 * (L[c] - med).abs().median()
        keep &= (L[c] - med).abs() <= max(5 * mad, 0.05 * abs(med))
    Lk = L[keep]
    reg_lines = np.vstack([np.interp(GRID, x - L.loc[i, 'xc'], lev[i]) for i in Lk.index])
    res = {c: float(Lk[c].median()) for c in ('H', 'W10', 'W20', 'W50', 'W80', 'A')}
    res.update({f'{c}_line_sd': float(Lk[c].std(ddof=1)) for c in ('H', 'W10', 'W50')})
    res.update({'lines_total': int(Z.shape[1]), 'lines_used': int(keep.sum()),
                'two_sided_frac': float(Lk['two_sided'].mean()), 'footprint15_frac': float(Lk['footprint_15'].mean()),
                'tilt_med_um_per_mm': float(np.nanmedian(Lk['tilt_um_per_mm'])), 'substrate_rms_um': float(np.nanmedian(Lk['sub_rms'])),
                'centre_drift_um': float(np.polyfit(Lk['line'], Lk['xc'], 1)[0] * (Lk['line'].max() - Lk['line'].min()))})
    return res, reg_lines, L, keep


def parse(fname):
    m = re.match(r'^(C\d|F\d)\s+((?:CM|FM)\d)(?:\s+(Mold|PDMS))?\.csv$', fname)
    return None if not m else (m.group(1), m.group(2), m.group(3) or 'Mold')        # CM1/FM1 mold files carry no suffix


def run(data_dir=DATA_DIR):
    recs, store = [], {}
    for f in sorted(os.listdir(data_dir)):
        p = parse(f)
        if p is None:
            continue
        ch, rep, mat = p
        res, lines, _, _ = process_file(os.path.join(data_dir, f), invert=(mat == 'PDMS'))
        res.update({'channel': ch, 'chip': int(rep[-1]), 'material': mat, 'file': f, 'H_cad': CAD[ch][0], 'W_cad': CAD[ch][1],
                    'layer': 'control' if ch[0] == 'C' else 'flow'})
        for k in FRACS:
            res[f'Wcad{k}'] = CAD_W[ch][k]; res[f'dW{k}'] = res[f'W{k}'] - CAD_W[ch][k]
        res['dH'] = res['H'] - CAD[ch][0]
        res['A_cad'] = CAD_A[ch]; res['A_ratio'] = res['A'] / CAD_A[ch]
        res['Weq'] = res['A'] / res['H']; res['Weq_cad'] = CAD_A[ch] / CAD[ch][0]; res['dWeq'] = res['Weq'] - res['Weq_cad']
        recs.append(res); store[(ch, mat, int(rep[-1]))] = lines
    return pd.DataFrame(recs), store


# =========================================================================================== statistics
def ci_t(v, level=0.95):
    v = np.asarray(v, float); v = v[~np.isnan(v)]; n = len(v)
    if n < 2:
        return v.mean(), np.nan, n, np.nan
    s = v.std(ddof=1)
    return v.mean(), s, n, stats.t.ppf(0.5 + level / 2, n - 1) * s / np.sqrt(n)


def per_feature(R):
    rows = []
    for (ch, mat), g in R.groupby(['channel', 'material']):
        d = {'channel': ch, 'material': mat, 'n_chips': len(g), 'H_cad': g.H_cad.iloc[0], 'W_cad': g.W_cad.iloc[0]}
        for c in ['H', 'W10', 'W20', 'W50', 'W80', 'A', 'A_ratio', 'Weq', 'dH', 'dW10', 'dW20', 'dW50', 'dW80', 'dWeq']:
            m, s, n, h = ci_t(g[c]); d.update({f'{c}_mean': m, f'{c}_sd': s, f'{c}_ci95': h})
        rows.append(d)
    return pd.DataFrame(rows)


def block_regression(y, h, chip):
    """y = a*h + block(chip) by least squares (chip = fixed block). Returns slope, intercept (mean of block effects),
    t-based 95% CIs, r (partial, within-chip), p and a function giving the CI half-width of the fitted mean line."""
    y, h, chip = map(np.asarray, (y, h, chip)); ids = np.unique(chip)
    X = np.column_stack([h] + [(chip == i).astype(float) for i in ids])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta; df = len(y) - X.shape[1]; s2 = (res ** 2).sum() / df
    cov = s2 * np.linalg.inv(X.T @ X); t = stats.t.ppf(0.975, df)
    a, b = beta[0], beta[1:].mean()
    var_b = cov[1:, 1:].sum() / len(ids) ** 2
    se_icpt = np.sqrt(var_b)
    # partial correlation of h and y after removing chip means
    hy = np.column_stack([h, y]); hy = hy - np.array([hy[chip == i].mean(0) for i in chip])
    r = np.corrcoef(hy[:, 0], hy[:, 1])[0, 1]
    p = 2 * stats.t.sf(abs(a / np.sqrt(cov[0, 0])), df)
    def band(x):
        x = np.atleast_1d(x); out = []
        for xv in x:
            v = np.concatenate([[xv], np.full(len(ids), 1 / len(ids))])
            out.append(t * np.sqrt(v @ cov @ v))
        return np.array(out)
    return dict(slope=a, slope_ci=t * np.sqrt(cov[0, 0]), icpt=b, icpt_ci=t * se_icpt, r=r, p=p, df=df, n=len(y), band=band)


def pooled_offsets(R, frac=10, col=None):
    """Rectangular C1-C3 footprint offset, chip-level pooling (mean of the three channels on each chip -> n = 3 chips)."""
    out = {}; col = col or f'dW{frac}'
    for mat in ('Mold', 'PDMS'):
        g = R[(R.material == mat) & R.channel.isin(['C1', 'C2', 'C3'])]
        chipm = g.groupby('chip')[col].mean()
        m, s, n, h = ci_t(chipm.values)
        sl = block_regression(g[col], g.W_cad, g.chip)
        out[mat] = dict(mean=m, sd=s, n=n, ci95=h, chip_means=chipm.values, slope_vs_width=sl['slope'], slope_p=sl['p'])
    return out


def flow_fits(R, frac=10, col=None):
    out = {}; col = col or f'dW{frac}'
    for mat in ('Mold', 'PDMS'):
        for lab, chans in (('F1-F4', ['F1', 'F2', 'F3', 'F4']), ('F1-F5', ['F1', 'F2', 'F3', 'F4', 'F5'])):
            g = R[(R.material == mat) & R.channel.isin(chans)]
            out[(mat, lab)] = block_regression(g[col], g.H_cad, g.chip)
    return out


def transfer(R):
    """Replica minus mold height, paired by chip number; the layer effect is the mean over the 5 channels of each chip
    (3 chip-level values, df = 2)."""
    M = R[R.material == 'Mold'].set_index(['channel', 'chip']).H; P = R[R.material == 'PDMS'].set_index(['channel', 'chip']).H
    D = (P - M).rename('dH').reset_index(); D['layer'] = np.where(D.channel.str.startswith('C'), 'control', 'flow')
    D['dH_rel_pct'] = 100 * D.dH.values / M.reindex(pd.MultiIndex.from_frame(D[['channel', 'chip']])).values
    feat = D.groupby('channel').dH.agg(['mean', 'std', 'count']).reset_index()
    feat['ci95'] = stats.t.ppf(0.975, feat['count'] - 1) * feat['std'] / np.sqrt(feat['count'])
    lay = {}
    for name in ('control', 'flow'):
        g = D[D.layer == name]; cm = g.groupby('chip').dH.mean()
        m, s, n, h = ci_t(cm.values); t = stats.ttest_1samp(cm.values, 0)
        mold_chip = R[(R.material == 'Mold') & (R.layer == name)].groupby('chip').H.mean().values
        pdms_chip = R[(R.material == 'PDMS') & (R.layer == name)].groupby('chip').H.mean().values
        w = stats.ttest_ind(pdms_chip, mold_chip, equal_var=False)                   # pairing-free comparison
        lay[name] = dict(mean=m, sd=s, n=n, ci95=h, p_paired=t.pvalue, chip_means=cm.values,
                         rel_pct=100 * m / mold_chip.mean(), p_unpaired=w.pvalue,
                         mean_features=g.groupby('channel').dH.mean().mean())
    return D, feat, lay


def f5_single_layer(R, frac=10, col=None):
    out = {}; col = col or f'dW{frac}'
    for mat in ('Mold', 'PDMS'):
        g = R[(R.material == mat) & (R.channel == 'F5')]
        out[mat] = ci_t(g[col].values)
    return out


def pooled_flow(R, col='dW10', chans=('F1', 'F2', 'F3', 'F4')):
    out = {}
    for mat in ('Mold', 'PDMS'):
        g = R[(R.material == mat) & R.channel.isin(chans)]
        out[mat] = ci_t(g.groupby('chip')[col].mean().values)
    return out


def flow_offset_diagnostics(R, frac=10):
    """Why a straight-line 'slope' of the footprint offset against the designed height is NOT reported as a result.
    Mold, F1-F4, offset against the design silhouette (chip as block): lack-of-fit of the straight line against H as a
    four-level factor, slopes with one feature left out, slope against the 350 um chord instead of the silhouette, and
    paired (by chip) differences between neighbouring features."""
    col = f'dW{frac}'; mold = R[R.material == 'Mold']
    g = mold[mold.channel.isin(['F1', 'F2', 'F3', 'F4'])]
    y, chip, H = g[col].values, g.chip.values, g.H_cad.values.astype(float)
    blk = np.column_stack([(chip == i).astype(float) for i in (1, 2, 3)])
    def rss(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; return (r ** 2).sum(), len(y) - np.linalg.matrix_rank(X)
    r0, d0 = rss(blk); r1, d1 = rss(np.column_stack([H, blk]))
    r2, d2 = rss(np.column_stack([(H == l).astype(float) for l in sorted(set(H))[1:]] + [blk]))
    F_lof = ((r1 - r2) / (d1 - d2)) / (r2 / d2); F_h = ((r0 - r2) / (d0 - d2)) / (r2 / d2)
    out = [dict(test='lack of fit of a straight line (F1-F4)', stat=f'F({d1-d2},{d2}) = {F_lof:.1f}', p=stats.f.sf(F_lof, d1 - d2, d2)),
           dict(test='H_f as a four-level factor (F1-F4)', stat=f'F({d0-d2},{d2}) = {F_h:.1f}', p=stats.f.sf(F_h, d0 - d2, d2))]
    b = block_regression(y, H, chip)
    out.append(dict(test='slope, F1-F4 (um/um)', stat=f"{b['slope']:.3f} +/- {b['slope_ci']:.3f} (95% CI, {b['df']} d.f.)", p=b['p']))
    for drop in ('F1', 'F2', 'F3', 'F4'):
        s_ = mold[mold.channel.isin([c for c in ('F1', 'F2', 'F3', 'F4') if c != drop])]
        bb = block_regression(s_[col], s_.H_cad, s_.chip)
        out.append(dict(test=f'slope without {drop}', stat=f"{bb['slope']:.3f} +/- {bb['slope_ci']:.3f} ({bb['df']} d.f.)", p=bb['p']))
    bc = block_regression(g.W10 - 350, g.H_cad, g.chip)
    out.append(dict(test='slope against the 350 um chord (not the silhouette)', stat=f"{bc['slope']:.3f} +/- {bc['slope_ci']:.3f} ({bc['df']} d.f.)", p=bc['p']))
    pm_ = mold.pivot_table(index='chip', columns='channel', values=col)
    for a_, b_ in (('F1', 'F2'), ('F2', 'F3'), ('F3', 'F4'), ('F4', 'F5')):
        d = pm_[a_] - pm_[b_]; t = stats.ttest_1samp(d, 0)
        out.append(dict(test=f'mean offset {a_} minus {b_} (paired over chips, um)', stat=f"{d.mean():+.1f} +/- {stats.t.ppf(.975, 2) * d.std(ddof=1) / np.sqrt(3):.1f} (2 d.f.)", p=t.pvalue))
    return pd.DataFrame(out)


def threshold_sensitivity(R):
    rows = []
    for k in FRACS:
        po = pooled_offsets(R, k); ff = flow_fits(R, k)[('Mold', 'F1-F4')]
        rows.append({'threshold_pct': k, 'rect_offset_mold': po['Mold']['mean'], 'rect_offset_mold_ci': po['Mold']['ci95'],
                     'rect_offset_pdms': po['PDMS']['mean'], 'rect_offset_pdms_ci': po['PDMS']['ci95'],
                     'F1F4_slope': ff['slope'], 'F1F4_slope_ci': ff['slope_ci'], 'F1F4_icpt': ff['icpt'], 'F1F4_icpt_ci': ff['icpt_ci'],
                     'F1F4_r': ff['r'], 'F1F4_df': ff['df']})
    po = pooled_offsets(R, col='dWeq'); ff = flow_fits(R, col='dWeq')[('Mold', 'F1-F4')]
    rows.append({'threshold_pct': 'area (W_eq)', 'rect_offset_mold': po['Mold']['mean'], 'rect_offset_mold_ci': po['Mold']['ci95'],
                 'rect_offset_pdms': po['PDMS']['mean'], 'rect_offset_pdms_ci': po['PDMS']['ci95'],
                 'F1F4_slope': ff['slope'], 'F1F4_slope_ci': ff['slope_ci'], 'F1F4_icpt': ff['icpt'], 'F1F4_icpt_ci': ff['icpt_ci'],
                 'F1F4_r': ff['r'], 'F1F4_df': ff['df']})
    return pd.DataFrame(rows)


# =========================================================================================== helper for legends
def _strip(fig, spec, title, handles=None, ncol=3, stat_lines=(), letter=None):
    ax = strip_axes(fig, spec)
    ax.text(0, 0.99, title, fontsize=BASE_FS, fontweight='bold', va='top', ha='left', color=C_TEXT)
    if letter: ax.text(-0.15, 0.99, letter, fontsize=10, fontweight='bold', va='top', ha='left')
    if handles:
        ax.legend(handles=handles, loc='center left', bbox_to_anchor=(-0.01, 0.50), ncol=ncol, frameon=False, fontsize=BASE_FS - 1.4,
                  handlelength=1.4, handletextpad=0.35, columnspacing=0.9, borderpad=0, labelspacing=0.25)
    for k, s in enumerate(stat_lines):
        ax.text(0, 0.02 + 0.215 * (len(stat_lines) - 1 - k), s, fontsize=BASE_FS - 1.4, va='bottom', ha='left', color=C_TEXT)
    return ax


# =========================================================================================== FIGURE 3 (main text)
def fig3(R, F, PO, FIT, T_D, T_F, T_L):
    """Profilometry only: footprint offset of control (a) and flow (b) features, and mold -> PDMS height transfer (c)."""
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.40, 1], width_ratios=[1, 1, 1.08], left=0.068, right=0.992, top=0.985,
                          bottom=0.145, wspace=0.34, hspace=0.07)
    rng = np.random.default_rng(1)
    mk = {'Mold': ('s', C_MOLD, C_MOLD), 'PDMS': ('o', None, None)}
    # ---------------- a
    colp = C_CTRL; dx = {'Mold': -9, 'PDMS': 9}
    ax = fig.add_subplot(gs[1, 0])
    for mat in ('Mold', 'PDMS'):
        g = R[(R.material == mat) & (R.layer == 'control')]; f = F[(F.material == mat) & F.channel.str.startswith('C')]
        col = C_MOLD if mat == 'Mold' else colp
        ax.scatter(g.W_cad + dx[mat] + rng.uniform(-3, 3, len(g)), g.dW10, s=7, marker='s' if mat == 'Mold' else 'o', facecolor='none',
                   edgecolor=col, lw=0.5, zorder=2)
        ax.errorbar(f.W_cad + dx[mat], f.dW10_mean, yerr=f.dW10_ci95, fmt='s' if mat == 'Mold' else 'o', ms=3.8, color=col,
                    mec=col if mat == 'Mold' else '#8C6D1F', capsize=1.8, lw=0.7, zorder=3)
        pm = PO[mat]
        ax.hlines(pm['mean'], 170, 535, color=col if mat == 'Mold' else '#8C6D1F', lw=0.7, ls='--', zorder=1)
    ax.axhline(0, color=C_TEXT, lw=0.5, ls=':')
    ax.set_xlim(0, 545); ax.set_ylim(-20, 300)
    ax.set_xlabel(r'Designed control width $W_c$ (µm)'); ax.set_ylabel(r'Footprint offset $\Delta W_{10}$ (µm)')
    ax.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', alpha=0.6, zorder=0)
    hd = [Line2D([], [], marker='s', ls='', mfc=C_MOLD, mec=C_MOLD, ms=3.8, label='Mold'),
          Line2D([], [], marker='o', ls='', mfc=colp, mec='#8C6D1F', ms=3.8, label='PDMS'),
          Line2D([], [], marker='o', ls='', mfc='none', mec=C_TEXT, ms=3.4, mew=0.5, label='Single chip')]
    _strip(fig, gs[0, 0], 'Control layer (rectangular)', hd, 3,
           [f"C1–C3 pooled: mold +{PO['Mold']['mean']:.0f} ± {PO['Mold']['ci95']:.0f} µm (dashed)",
            f"PDMS +{PO['PDMS']['mean']:.0f} ± {PO['PDMS']['ci95']:.0f} µm; mean ± 95% CI over 3 chips"], letter='a')
    # ---------------- b   (descriptive: per-feature mean +/- 95% CI over 3 chips and the single chips; NO fitted line:
    #                         the offsets are not a straight-line function of H_f, see Supplementary Fig. S4e,f)
    ax = fig.add_subplot(gs[1, 1]); dxb = {'Mold': -6, 'PDMS': 6}
    for mat in ('Mold', 'PDMS'):
        g = R[(R.material == mat) & (R.layer == 'flow')]; f = F[(F.material == mat) & F.channel.str.startswith('F')]
        col = C_MOLD if mat == 'Mold' else C_FLOW
        ax.scatter(g.H_cad + dxb[mat] + rng.uniform(-2.5, 2.5, len(g)), g.dW10, s=7, marker='s' if mat == 'Mold' else 'o',
                   facecolor='none', edgecolor=col, lw=0.5, zorder=2)
        ax.errorbar(f.H_cad + dxb[mat], f.dW10_mean, yerr=f.dW10_ci95, fmt='s' if mat == 'Mold' else 'o', ms=3.8, color=col,
                    capsize=1.8, lw=0.7, zorder=3)
    for ch in ('F5', 'F4', 'F3', 'F2', 'F1'):
        ax.text(CAD[ch][0], 322, ch, ha='center', va='top', fontsize=BASE_FS - 1.2, color=C_FLOW, fontweight='bold')
    ax.axhline(0, color=C_TEXT, lw=0.5, ls=':')
    ax.set_xlim(0, 385); ax.set_ylim(-60, 330)
    ax.set_xlabel(r'Designed flow height $H_f$ (µm)'); ax.set_ylabel(r'Footprint offset $\Delta W_{10}$ (µm)')
    ax.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', alpha=0.6, zorder=0)
    hd = [Line2D([], [], marker='s', ls='', mfc=C_MOLD, mec=C_MOLD, ms=3.8, label='Mold'),
          Line2D([], [], marker='o', ls='', mfc=C_FLOW, mec=C_FLOW, ms=3.8, label='PDMS'),
          Line2D([], [], marker='o', ls='', mfc='none', mec=C_TEXT, ms=3.4, mew=0.5, label='Single chip')]
    fm_ = F[(F.material == 'Mold') & F.channel.isin(['F1', 'F2', 'F3', 'F4'])].dW10_mean; f5 = f5_single_layer(R)
    _strip(fig, gs[0, 1], 'Flow layer (circular segment)', hd, 3,
           [f"mold, F1–F4 (2+ layers): +{fm_.min():.0f} to +{fm_.max():.0f} µm, not monotonic in $H_f$",
            f"single layer (F5): mold +{f5['Mold'][0]:.0f} ± {f5['Mold'][3]:.0f} µm (95% CI, 3 chips)"], letter='b')
    # ---------------- c
    ax = fig.add_subplot(gs[1, 2])
    labs = CHANNELS; xpos = np.arange(10) + np.r_[np.zeros(5), np.ones(5) * 0.8]
    for xi, ch in zip(xpos, labs):
        colm = '#8C6D1F' if ch[0] == 'C' else C_TEAL_D; colo = C_CTRL if ch[0] == 'C' else C_FLOW
        d = T_D[T_D.channel == ch].sort_values('chip').dH.values
        ax.scatter(xi + np.array([-0.17, 0, 0.17])[:len(d)], d, s=8, facecolor='none', edgecolor=colo, lw=0.6, zorder=2)
        r = T_F[T_F.channel == ch].iloc[0]
        ax.errorbar(xi, r['mean'], yerr=r['ci95'], fmt='_', ms=8, mew=1.3, color=colm, capsize=2, lw=0.8, zorder=3)
    for name, sl, col in (('control', slice(0, 5), '#8C6D1F'), ('flow', slice(5, 10), C_TEAL_D)):
        ax.hlines(T_L[name]['mean'], xpos[sl][0] - 0.45, xpos[sl][-1] + 0.45, color=col, lw=0.8, ls='--', zorder=1)
    ax.axhline(0, color=C_TEXT, lw=0.5)
    ax.set_xticks(xpos); ax.set_xticklabels(labs); ax.set_xlim(-0.7, xpos[-1] + 0.7)
    ax.set_xlabel('Feature')
    lo = min(T_D.dH.min(), (T_F['mean'] - T_F['ci95']).min()) - 5; hi = max(T_D.dH.max(), (T_F['mean'] + T_F['ci95']).max()) + 5
    ax.set_ylim(max(lo, -200), min(hi, 120))
    ax.set_ylabel(r'$H_{PDMS} - H_{mold}$ (µm)')
    ax.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', alpha=0.6, zorder=0)
    hd = [Line2D([], [], marker='o', ls='', mfc='none', mec=C_TEXT, ms=3.4, mew=0.5, label='Chip pair'),
          Line2D([], [], marker='_', ls='', color=C_TEXT, ms=7, mew=1.3, label='Mean ± 95% CI'),
          Line2D([], [], color=C_TEXT, ls='--', lw=0.8, label='Layer mean')]
    _strip(fig, gs[0, 2], 'Height transfer, mold → PDMS', hd, 3,
           [f"control {T_L['control']['mean']:+.1f} ± {T_L['control']['ci95']:.1f} µm; flow {T_L['flow']['mean']:+.1f} ± {T_L['flow']['ci95']:.1f} µm".replace('-', '−'),
            "layer means over 3 chips (paired), 95% CI"], letter='c')
    return fig


# =========================================================================================== SUPPLEMENTARY FIGURES
def figS_individual(R, store, group):
    """Individual profiles of every mold and PDMS replica of one layer: all scan lines (thin) registered on the feature
    centre, file median (bold), design silhouette (dashed)."""
    chans = [c for c in CHANNELS if c.startswith(group)]
    fig, axs = plt.subplots(5, 6, figsize=(7.2, 6.9), sharex='row', sharey='row')
    fig.subplots_adjust(left=0.105, right=0.995, top=0.915, bottom=0.07, wspace=0.09, hspace=0.34)
    for i, ch in enumerate(chans):
        H, W = CAD[ch]; half = max(R[R.channel == ch].W10.max(), CAD_MAXW[ch]) / 2 + 45
        for j, (mat, rep) in enumerate([(m, r) for m in ('Mold', 'PDMS') for r in (1, 2, 3)]):
            ax = axs[i, j]; lines = store[(ch, mat, rep)]
            col = C_MOLD if mat == 'Mold' else (C_CTRL if group == 'C' else C_FLOW)
            for ln in lines: ax.plot(GRID, ln, color=col, lw=0.3, alpha=0.30)
            ax.plot(GRID, np.median(lines, axis=0), color=C_TEXT if mat == 'Mold' else ('#8C6D1F' if group == 'C' else C_TEAL_D), lw=0.9)
            ax.plot(GRID, cad_profile(ch, GRID), color=C_TEXT, lw=0.6, ls=(0, (3, 2)))
            r_ = R[(R.channel == ch) & (R.material == mat) & (R.chip == rep)].iloc[0]
            ax.text(0.03, 0.97, f"H {r_.H:.0f}\nW10 {r_.W10:.0f}\nW50 {r_.W50:.0f}", transform=ax.transAxes, fontsize=4.6, va='top', ha='left')
            ax.set_xlim(-half, half); ax.set_ylim(-0.08 * max(H, 40), max(H, 40) * 1.5)
            ax.axhline(0, color='#AAAAAA', lw=0.3)
            ax.tick_params(labelsize=BASE_FS - 2, length=2, pad=1)
            if i == 0: ax.set_title(f"{'Mold' if mat == 'Mold' else 'PDMS'} chip {rep}", fontsize=BASE_FS - 0.5, pad=3)
            if j == 0:
                ax.set_ylabel('$z$ (µm)', fontsize=BASE_FS - 0.5)
                ax.text(-0.52, 0.5, ch, transform=ax.transAxes, fontsize=BASE_FS + 1, fontweight='bold', va='center', ha='center', color=C_FLOW if group == 'F' else '#8C6D1F')
            if i == 4: ax.set_xlabel('$x$ (µm)', fontsize=BASE_FS - 0.5)
    hd = [Line2D([], [], color=C_MOLD, lw=0.6, alpha=0.6, label='single scan lines (mold)'),
          Line2D([], [], color=C_CTRL if group == 'C' else C_FLOW, lw=0.6, alpha=0.6, label='single scan lines (PDMS, inverted)'),
          Line2D([], [], color=C_TEXT, lw=0.9, label='file median'), Line2D([], [], color=C_TEXT, lw=0.6, ls=(0, (3, 2)), label='design (CAD)')]
    fig.legend(handles=hd, loc='upper center', ncol=4, frameon=False, fontsize=BASE_FS - 1.2, bbox_to_anchor=(0.54, 0.995), handlelength=1.6)
    return fig


def figS_cad(F, store):
    """Design silhouettes (dashed) with mean registered profiles of the three molds and the three PDMS replicas."""
    fig, axs = plt.subplots(2, 5, figsize=(7.2, 3.5), sharex=False)
    fig.subplots_adjust(left=0.07, right=0.995, top=0.86, bottom=0.12, wspace=0.28, hspace=0.5)
    for k, ch in enumerate(CHANNELS):
        ax = axs.ravel()[k]; H, W = CAD[ch]
        for mat in ('Mold', 'PDMS'):
            ps = np.vstack([np.median(store[(ch, mat, r)], axis=0) for r in (1, 2, 3)])
            col = C_MOLD if mat == 'Mold' else (C_CTRL if ch[0] == 'C' else C_FLOW)
            ax.fill_between(GRID, ps.min(0), ps.max(0), color=col, alpha=0.22, lw=0)
            ax.plot(GRID, ps.mean(0), color=col if mat == 'Mold' else col, lw=1.0)
        ax.plot(GRID, cad_profile(ch, GRID), color=C_TEXT, lw=0.7, ls=(0, (3, 2)))
        ax.set_xlim(-400, 400); ax.set_ylim(-0.06 * H, 1.45 * H)
        ttl = f"{ch}: " + (f"{W} × {H} µm" if ch[0] == 'C' else rf"chord 350, $H_f$ {H} µm")
        ax.set_title(ttl, fontsize=BASE_FS - 0.8, pad=3)
        if k % 5 == 0: ax.set_ylabel('$z$ (µm)')
        if k >= 5: ax.set_xlabel('$x$ (µm)')
        ax.tick_params(labelsize=BASE_FS - 1.8)
        if ch == 'F1':
            ax.annotate('design bulb:\n437.5 µm wide at\n131 µm height', xy=(-218, 131), xytext=(-380, 480), fontsize=BASE_FS - 2.3,
                        arrowprops=dict(arrowstyle='-', lw=0.4, color=C_TEXT), ha='left', va='top')
    hd = [Line2D([], [], color=C_MOLD, lw=1.0, label='mold, mean of 3 chips (band = range)'),
          Line2D([], [], color=C_CTRL, lw=1.0, label='PDMS control, inverted'), Line2D([], [], color=C_FLOW, lw=1.0, label='PDMS flow, inverted'),
          Line2D([], [], color=C_TEXT, lw=0.7, ls=(0, (3, 2)), label='design silhouette seen from above')]
    fig.legend(handles=hd, loc='upper center', ncol=4, frameon=False, fontsize=BASE_FS - 1.4, bbox_to_anchor=(0.54, 0.995), handlelength=1.6)
    return fig


def figS_qc(R, F, TS, store, data_dir=DATA_DIR):
    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 3, left=0.07, right=0.99, top=0.93, bottom=0.09, wspace=0.42, hspace=0.45)
    # a levelling
    ax = fig.add_subplot(gs[0, 0])
    df = pd.read_csv(os.path.join(data_dir, 'F1 FM1 Mold.csv'), header=None); x = df.iloc[:, 1].values
    tl = np.array([abs(level_line(x, df.iloc[:, k].values)[1]['tilt_um_per_mm']) for k in range(2, df.shape[1])])
    z = df.iloc[:, 2 + int(np.nanargmin(np.abs(tl - np.nanpercentile(tl, 90))))].values
    zl, info = level_line(x, z); m = measure_line(x, zl)
    ax.plot(x, z - np.median(np.sort(z)[:len(z) // 5]), color=C_MOLD, lw=0.6, label='raw (offset only)')
    ax.plot(x, zl, color=C_TEXT, lw=0.6, label='levelled (two-sided line)')
    for f, ls_, lb in ((0.1, '-', 'W10'), (0.5, '--', 'W50')): ax.axhline(f * m['H'], color=C_FLOW, lw=0.5, ls=ls_)
    ax.axhline(m['H'], color=C_BRIDGE, lw=0.5)
    ax.set_xlabel(r'$x$ (µm)'); ax.set_ylabel(r'$z$ (µm)'); ax.set_ylim(-30, m['H'] * 1.45)
    ax.legend(loc='upper left', fontsize=BASE_FS - 1.8, borderaxespad=0, handlelength=1.4)
    ax.text(880, m['H'] + 8, 'H', color=C_BRIDGE, fontsize=BASE_FS - 1.5, ha='right'); ax.text(880, 0.5 * m['H'] + 8, '50%', color=C_FLOW, fontsize=BASE_FS - 1.5, ha='right')
    ax.text(880, 0.1 * m['H'] + 8, '10%', color=C_FLOW, fontsize=BASE_FS - 1.5, ha='right')
    ax.set_title(f"F1 FM1 mold, one line (tilt {info['tilt_um_per_mm']:.0f} µm/mm)", fontsize=BASE_FS - 0.5, loc='left', pad=4)
    panel_label(ax, 'a', x=-0.2)
    # b drift
    ax = fig.add_subplot(gs[0, 1])
    for f, col in (('F5 FM3 Mold.csv', C_FLOW), ('C3 CM3 PDMS.csv', C_CTRL), ('C1 CM1.csv', C_MOLD)):
        mat = parse(f)[2]; _, _, L, keep = process_file(os.path.join(data_dir, f), invert=(mat == 'PDMS'))
        ax.plot(L.line[keep], L.xc[keep] - L.xc[keep].median(), lw=0.8, color=col, label=f[:-4])
    ax.axhline(0, color=C_TEXT, lw=0.5, ls=':'); ax.set_xlabel('Scan line'); ax.set_ylabel('Centre − median (µm)')
    ax.legend(loc='upper right', fontsize=BASE_FS - 1.8, borderaxespad=0, handlelength=1.4)
    ax.set_title('Feature centre across scan lines', fontsize=BASE_FS - 0.5, loc='left', pad=4); panel_label(ax, 'b', x=-0.2)
    # c chip effect on heights
    ax = fig.add_subplot(gs[0, 2]); ctr = R[(R.layer == 'control') & (R.material == 'Mold')]
    for chip, col, mk in ((1, '#D4B96A', 'o'), (2, C_CTRL, 's'), (3, '#8C6D1F', '^')):
        g = ctr[ctr.chip == chip]; ax.scatter(g.channel.str[1].astype(int) + (chip - 2) * 0.18, g.H, s=14, marker=mk, color=col, label=f'chip {chip}', zorder=3)
    ax.axhline(200, color=C_TEXT, lw=0.5, ls='--'); ax.set_xticks(range(1, 6)); ax.set_xticklabels(['C1', 'C2', 'C3', 'C4', 'C5'])
    ax.set_ylabel('Mold height (µm)'); ax.set_ylim(150, 245); ax.legend(ncol=3, loc='upper center', fontsize=BASE_FS - 1.8, borderaxespad=0.1, handletextpad=0.2, columnspacing=0.8)
    ax.set_title('Chip-to-chip variation (design 200 µm)', fontsize=BASE_FS - 0.5, loc='left', pad=4); panel_label(ax, 'c', x=-0.2)
    # d threshold sensitivity
    TS = TS[TS.threshold_pct.apply(lambda v: not isinstance(v, str))].copy(); TS['threshold_pct'] = TS.threshold_pct.astype(int)
    ax = fig.add_subplot(gs[1, 0]); th = TS.threshold_pct.values
    ax.errorbar(th - 1, TS.rect_offset_mold, yerr=TS.rect_offset_mold_ci, fmt='s-', ms=3, color=C_MOLD, capsize=1.6, lw=0.7, label='C1–C3 mold')
    ax.errorbar(th + 1, TS.rect_offset_pdms, yerr=TS.rect_offset_pdms_ci, fmt='o-', ms=3, color='#8C6D1F', capsize=1.6, lw=0.7, label='C1–C3 PDMS')
    ax.axhline(0, color=C_TEXT, lw=0.5, ls=':'); ax.set_xticks(th); ax.set_xlabel('Width threshold (% of H)'); ax.set_ylabel(r'Pooled offset (µm)')
    ax.legend(loc='upper right', fontsize=BASE_FS - 1.8, borderaxespad=0.1); ax.set_title('Rectangles vs threshold', fontsize=BASE_FS - 0.5, loc='left', pad=4); panel_label(ax, 'd', x=-0.2)
    ax = fig.add_subplot(gs[1, 1])                                   # e: reference definition
    mold = R[R.material == 'Mold']
    for lab, col_, mk, off in (('vs design silhouette (437.5 µm for F1)', C_FLOW, 'o', -4), ('vs 350 µm chord', C_MOLD, 's', 4)):
        xs, ym, yc = [], [], []
        for ch in ('F5', 'F4', 'F3', 'F2', 'F1'):
            v = mold[mold.channel == ch].W10.values - (CAD_W[ch][10] if 'silhouette' in lab else 350)
            m_, s_, n_, h_ = ci_t(v); xs.append(CAD[ch][0] + off); ym.append(m_); yc.append(h_)
        ax.errorbar(xs, ym, yerr=yc, fmt=mk + '-', ms=3, color=col_, capsize=1.6, lw=0.7, label=lab)
    ax.set_xlabel(r'Designed flow height $H_f$ (µm)'); ax.set_ylabel(r'Mold offset $\Delta W_{10}$ (µm)'); ax.set_ylim(-10, 320)
    ax.legend(loc='lower right', fontsize=BASE_FS - 2.2, borderaxespad=0.1, handlelength=1.4)
    ax.set_title('Offset depends on the reference', fontsize=BASE_FS - 0.5, loc='left', pad=4); panel_label(ax, 'e', x=-0.2)
    ax = fig.add_subplot(gs[1, 2])                                   # f: leave-one-feature-out slopes
    labs, sl_, ci_ = ['F1–F4'], [], []
    b0 = block_regression(mold[mold.channel.isin(['F1', 'F2', 'F3', 'F4'])].dW10, mold[mold.channel.isin(['F1', 'F2', 'F3', 'F4'])].H_cad, mold[mold.channel.isin(['F1', 'F2', 'F3', 'F4'])].chip)
    sl_.append(b0['slope']); ci_.append(b0['slope_ci'])
    for drop in ('F1', 'F2', 'F3', 'F4'):
        g = mold[mold.channel.isin([c for c in ('F1', 'F2', 'F3', 'F4') if c != drop])]
        bb = block_regression(g.dW10, g.H_cad, g.chip); labs.append(f'–{drop}'); sl_.append(bb['slope']); ci_.append(bb['slope_ci'])
    ypos = np.arange(len(labs))
    ax.errorbar(sl_, ypos, xerr=ci_, fmt='o', ms=3.4, color=C_TEAL_D, capsize=1.8, lw=0.7)
    ax.axvline(0, color=C_TEXT, lw=0.5, ls=':'); ax.set_yticks(ypos); ax.set_yticklabels(labs); ax.invert_yaxis()
    ax.set_xlabel('Slope of $\\Delta W_{10}$ vs $H_f$ (µm/µm)')
    ax.set_title('Slope with one feature left out', fontsize=BASE_FS - 0.5, loc='left', pad=4); panel_label(ax, 'f', x=-0.3)
    return fig


# =========================================================================================== driver
def analyse(data_dir=DATA_DIR):
    R, store = run(data_dir)
    F = per_feature(R); PO = pooled_offsets(R); FIT = flow_fits(R); T_D, T_F, T_L = transfer(R)
    TS = threshold_sensitivity(R)
    return R, store, F, PO, FIT, T_D, T_F, T_L, TS


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    R, store, F, PO, FIT, T_D, T_F, T_L, TS = analyse()
    R.to_csv(os.path.join(OUT_DIR, 'profilometry_per_file.csv'), index=False)
    F.to_csv(os.path.join(OUT_DIR, 'profilometry_per_feature.csv'), index=False)
    T_D.to_csv(os.path.join(OUT_DIR, 'height_transfer_pairs.csv'), index=False)
    TS.to_csv(os.path.join(OUT_DIR, 'threshold_sensitivity.csv'), index=False)
    flow_offset_diagnostics(R).to_csv(os.path.join(OUT_DIR, 'flow_offset_diagnostics.csv'), index=False)
    f3 = fig3(R, F, PO, FIT, T_D, T_F, T_L)
    s1 = figS_individual(R, store, 'C'); s2 = figS_individual(R, store, 'F'); s3 = figS_cad(F, store); s4 = figS_qc(R, F, TS, store)
    # ---- figure saving (uncomment to export)
    # f3.savefig('Fig3_profilometry.jpeg', dpi=600, bbox_inches='tight')
    # s1.savefig('FigS_profiles_control.jpeg', dpi=600, bbox_inches='tight')
    # s2.savefig('FigS_profiles_flow.jpeg', dpi=600, bbox_inches='tight')
    # s3.savefig('FigS_design_vs_measured.jpeg', dpi=600, bbox_inches='tight')
    # s4.savefig('FigS_profilometry_QC.jpeg', dpi=600, bbox_inches='tight')
    plt.show()
