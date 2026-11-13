"""
mixing_analysis.py - homogenisation of dye in the closed-loop peristaltic mixer from bright-field video.

DATA   three aligned/cropped videos (ESI delay 33, 83, 167 ms; 25 fps) and, per condition, the TIFF stills of the same
       field (blank chip, water, pure dye, 1:1 premix).  ONE video per condition (n = 1).

WHY THE EARLIER RED-CHANNEL RMI WAS NOT USED (checked on these videos)
  * Red channel: pure dye (0.5 % FD&C Blue 1) and the 0.25 % premix are both at the black floor (median 0) in the TIFFs; 13-15 % of
    loop pixels sit at the floor in the first frame and the dye-mass balance of the sealed loop is violated by 69-87 %.
  * Std of the raw red over the bounding box of the loop is dominated by PDMS background (the channel is ~23 % of the frame), and its
    value changes by up to +/-40 % of the printed t90 with the box size.
  * The TIFFs cannot be used as pixel-wise references for the videos (video/TIFF background ratio 1.09-1.13, MPEG-4 re-encode,
    one set of stills for all three videos).

METHOD (video only; no TIFF intensity is used - the TIFFs only supply the channel mask)
  1. Channel mask: red-channel difference water-minus-dye TIFF, thresholded; inlet/outlet stubs and the arm carrying the pump valves
     are removed; 2 x 2 binning; only pixels fully inside the channel are kept.
  2. Signal: green channel (the only channel in its linear range at this dye concentration).  Absorbance relative to the FINAL state,
        L(x,t) = -log10( I(x,t) / I_final(x) ),    I_final = median of the last REF_S seconds,
     after 0.5 s temporal averaging and a K x K normalised box smoothing restricted to the channel.  Static structure (walls, texture,
     illumination, valves) cancels pixel by pixel because both frames come from the same video.
  3. Index (PRIMARY: concentration-based, after Johnson et al. 2002 and Bossink et al. 2022)
        M(t) = 1 - sigma_L(t) / sigma_ref ,   sigma_L = std of L over the selected channel pixels,
     sigma_ref = std of the ideal unmixed state with the same mean (half the fluid at 0, half at twice the mean). In units of L this equals the
     absorbance of the mean concentration, taken here from the initially clear pixels (their absorbance after mixing); no calibration series was
     recorded.  M(0) is therefore ~0.1-0.2, not 0: the real initial field is less heterogeneous than an ideal step.  The variant with the
     MEASURED initial std, M = 1 - sigma_L(t)/sigma_L(0), is kept as a secondary output (reference='measured').
  4. Pump onset = first time the channel changes persistently (frame-difference energy), not a hand-set t_start.  The loop is still being
     filled for ~0.5-1 s after onset (total dye in the loop rises); t = 0 of closed-loop mixing is the end of that step (t_seal),
     and sigma(0) is the std in the first 0.6 s after t_seal.
  5. The curve is only shown up to (end - REF_S): later frames are part of the reference and trivially approach M = 1.
  6. Times are reported as t50 / t75 / t90 (first time M reaches 0.5 / 0.75 / 0.9 and stays within 0.03 of it) AND as M at fixed times, WITHIN THE
     FIRST T_WINDOW = 30 s after sealing (the final-state reference still comes from the end of the whole record; charts show 30 s only).
     The primary analysis excludes the bottom arm (a trapped bubble sits on it in all three videos).  Sensitivity grid: K in {5,7,9,11},
     REF_S in {2,3,5} s, arm exclusion {none, bottom}, edge erosion {0, 3} px.
  LIMITS  the final state must be uniform for M to mean "mixed": true for 33 and 83 ms (residual 2-4 % of the start), NOT for 167 ms
  (residual 21-30 %), whose t50 is therefore indicative only.  n = 1: times have no replicate uncertainty, only the sensitivity range.
"""
import os, sys, glob, csv
import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.signal import welch, detrend
from scipy.ndimage import uniform_filter1d
from ch2_style import *

DATA_DIR = './Mixing data'
OUT_DIR = './results_mixing'
CONDS = {'33 ms':  dict(folder='33 ms - Aligned_Cropped_Data',  video='33ms_cropped_Mixing_video.mp4'),
         '83 ms':  dict(folder='83 ms - Aligned_Cropped_Data',  video='83ms_cropped_Mixing_video.mp4'),
         '167 ms': dict(folder='167 ms - Aligned_Cropped_Data', video='167ms_cropped_Mixing_video.mp4')}
COLORS = {'33 ms': '#0F4C4C', '83 ms': '#2A8A8A', '167 ms': '#79BDB4'}      # teal ramp (flow-layer colour), dark = highest frequency
LSTYLE = {'33 ms': (0, ()), '83 ms': (0, (5, 1.5)), '167 ms': (0, (3, 1, 1, 1))}
PHASE_HZ = {'33 ms': 30.3, '83 ms': 12.0, '167 ms': 6.0}    # PROGRAMMED phase frequency = 1000 / ESI delay (ms); one cycle = 6 phases
PHASE_S = {'33 ms': 0.0563, '83 ms': 0.1007, '167 ms': 0.2001}  # seconds per phase DELIVERED in the pump runs (elapsed time / (6 x cycles)); used only in Fig. S9b
K_PRIMARY, REF_PRIMARY = 7, 3.0


def hz(cond):
    return f'{PHASE_HZ[cond]:.1f} Hz'

GRID_K, GRID_REF = (5, 7, 9, 11), (2.0, 3.0, 5.0)
AVG_HALF_S = 0.25            # temporal averaging half-window (s)
SAMPLE_HZ = 5                # curve sampling rate (s^-1)
LEVELS = (0.5, 0.75, 0.9)
SUSTAIN = 0.03
T_WINDOW = 30.0           # reporting and chart window after sealing (s)
FIXED_T = (5, 10, 20, 30)


# ------------------------------------------------------------------------------------------------ geometry
def _find(folder, key):
    f = [p for p in glob.glob(os.path.join(DATA_DIR, folder, '*')) if key.lower() in os.path.basename(p).lower() and p.lower().endswith('.tiff')]
    return sorted(f, key=len)[0]


def geometry(cond):
    c = CONDS[cond]
    A = cv2.imread(_find(c['folder'], 'Stream_A')).astype(np.float32); Bw = cv2.imread(_find(c['folder'], 'Stream_B')).astype(np.float32)
    blank = cv2.imread(_find(c['folder'], 'Blank_chip')); cap = cv2.VideoCapture(os.path.join(DATA_DIR, c['folder'], c['video']))
    h, w = int(cap.get(4)), int(cap.get(3)); cap.release()
    A, Bw, blank = A[:h, :w], Bw[:h, :w], blank[:h, :w]
    d = cv2.GaussianBlur(Bw[:, :, 2] - A[:, :, 2], (0, 0), 2.0); mask = d > 0.5 * np.percentile(d, 98)
    mask = cv2.morphologyEx(mask.astype(np.uint8), cv2.MORPH_OPEN, np.ones((5, 5), np.uint8)) > 0
    rows, cols = mask.sum(1), mask.sum(0)
    vc, hr = np.where(cols > 0.45 * rows.max())[0], np.where(rows > 0.45 * cols.max())[0]
    xl, xr, yt, yb = int(vc[0]), int(vc[-1]), int(hr[0]), int(hr[-1])
    g = cv2.GaussianBlur(cv2.cvtColor(blank, cv2.COLOR_BGR2GRAY).astype(np.float32), (0, 0), 1.5)
    sy = np.abs(cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=3))[:, xr + 8:w].mean(1)
    pk = sorted(np.where((sy > 0.45 * sy.max()) & (np.r_[sy[1:], 0] <= sy) & (np.r_[0, sy[:-1]] <= sy))[0].tolist())
    vy0, vy1 = pk[0], pk[-1]
    valves = [(pk[1] + 3, pk[2] - 3), (pk[5] + 3, pk[6] - 3), (pk[9] + 3, pk[10] - 3)]            # interiors of the three valve rectangles
    ring = mask.copy(); ring[:yt] = False; ring[yb + 1:] = False; ring[:, :xl - 3] = False
    vz = np.zeros_like(ring); vz[vy0 - 12:vy1 + 13, int(0.5 * w):] = True; vz &= ring; ring_noV = ring & ~vz                # arm along the three pump valves
    hb, wb = h // 2, w // 2
    ring_b = cv2.resize(ring.astype(np.float32), (wb, hb), interpolation=cv2.INTER_AREA) > 0.99
    vz_b = cv2.resize(vz.astype(np.float32), (wb, hb), interpolation=cv2.INTER_AREA) > 0.5
    ys, xs = np.where(ring_b)
    aw = 30                                                                                 # arm width (binned px)
    yt_b, yb_b, xl_b, xr_b = yt // 2, yb // 2, xl // 2, xr // 2
    arm = np.full(len(ys), -1); arm[ys < yt_b + aw] = 0; arm[ys > yb_b - aw] = 1
    mid = (arm < 0); arm[mid & (xs < xl_b + aw)] = 2; arm[(arm < 0) & (xs > xr_b - aw)] = 3                # 0 top, 1 bottom, 2 left, 3 right (above/below valves)
    arm[vz_b[ys, xs]] = 4                                                                                  # 4 = right arm along the pump valves (kept for the dye inventory only)
    edge = cv2.distanceTransform(ring_b.astype(np.uint8), cv2.DIST_L2, 3)[ys, xs]
    return dict(h=h, w=w, hb=hb, wb=wb, ys=ys, xs=xs, arm=arm, edge=edge, valves=valves, xr=xr)


# ------------------------------------------------------------------------------------------------ video
def read_video(cond, G):
    c = CONDS[cond]; cap = cv2.VideoCapture(os.path.join(DATA_DIR, c['folder'], c['video']))
    fps = cap.get(cv2.CAP_PROP_FPS); n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    px = np.zeros((n, len(G['ys'])), np.uint8); vm = np.zeros((n, 3), np.float32); frames = {}
    vx0, vx1 = G['xr'] + 15, G['w'] - 6; i = 0
    keep = set()                                                                            # frames kept for the picture panel are returned by request
    while True:
        ok, f = cap.read()
        if not ok:
            break
        gch = f[:, :, 1].astype(np.float32); gb = cv2.resize(gch, (G['wb'], G['hb']), interpolation=cv2.INTER_AREA)
        px[i] = np.clip(np.round(gb[G['ys'], G['xs']]), 0, 255).astype(np.uint8)
        for j, (a, b) in enumerate(G['valves']):
            vm[i, j] = f[a:b, vx0:vx1].mean()
        i += 1
    cap.release()
    return px[:i], vm[:i], fps


def smooth_norm(v, G, k):
    img = np.zeros((G['hb'], G['wb']), np.float32); m = np.zeros_like(img); img[G['ys'], G['xs']] = v; m[G['ys'], G['xs']] = 1
    return (cv2.blur(img, (k, k)) / np.clip(cv2.blur(m, (k, k)), 1e-6, None))[G['ys'], G['xs']]


def detect_onset(px, fps, G):
    X = px[:int(min(len(px), 25 * fps))].astype(np.float32); ref = X[:int(0.3 * fps)].mean(0)
    d = np.array([np.sqrt(np.mean(smooth_norm(x - ref, G, 7) ** 2)) for x in X]); t = np.arange(len(X)) / fps
    b = d[:int(0.3 * fps)]; thr = max(b.mean() + 6 * max(b.std(), 1e-3), 1.5 * b.max())
    for j in np.where(d > thr)[0]:
        if (d[j:j + int(0.4 * fps)] > thr).all():
            return float(t[j]), t, d, thr
    return None, t, d, thr


def detect_seal(px, fps, onset, G):
    """End of the fast filling step: the mean absorbance of the WHOLE loop (relative to the final state) rises while dye is pumped in
    from the inlet and is constant once the loop is sealed.  t_seal = first time after the fastest rise at which the rise rate has fallen
    below 10 % of its peak.  This is t = 0 of closed-loop mixing."""
    n = len(px); I_f = np.median(px[-int(3 * fps):].astype(np.float32), axis=0)
    inv = (-np.log10(np.clip(px.astype(np.float32), 1, None) / np.clip(I_f, 1, None))).mean(1)
    t = np.arange(n) / fps; m = uniform_filter1d(inv, size=max(1, int(0.4 * fps))); sl = np.gradient(m, t)
    i0 = int(onset * fps); ip = i0 + int(np.argmax(sl[i0:i0 + int(8 * fps)])); j = ip
    while j < n - 1 and sl[j] > 0.10 * sl[ip]:
        j += 1
    return float(t[j]), t, m


def sigma_curve(px, fps, onset, G, k, ref_s):
    """sigma_L over every channel pixel (and the sample times) for one smoothing/reference setting; subsets are chosen later."""
    n = len(px); I_f = np.median(px[-int(ref_s * fps):].astype(np.float32), axis=0)
    idx = np.arange(0, n, max(1, int(fps / SAMPLE_HZ))); Ls = []
    for i in idx:
        a, b = max(0, i - int(AVG_HALF_S * fps)), min(n, i + int(AVG_HALF_S * fps) + 1)
        r = smooth_norm(px[a:b].astype(np.float32).mean(0) / np.clip(I_f, 1, None), G, k)
        Ls.append((-np.log10(np.clip(r, 1e-3, None))).astype(np.float32))
    return idx / fps, np.array(Ls)


def index_from(t, Ls, t0, sel, end_t, ref_s, reference='step'):
    """M = 1 - sigma_L(t)/sigma_ref.  reference='step': ideal unmixed state with the same mean (scale from the initially clear pixels; PRIMARY);
    reference='measured': sigma_L over the first 0.6 s after t0 (variant).  Curve = t0 ... end - ref_s."""
    sig = Ls[:, sel].std(1); m0 = (t >= t0) & (t <= t0 + 0.6)
    if reference == 'measured':
        s0 = sig[m0].mean()
    else:
        L0 = Ls[m0].mean(0)[sel]; q = np.percentile(L0, 2); s0 = -np.median(L0[L0 < 0.5 * q])
    M = 1 - sig / s0; keep = (t >= t0) & (t <= end_t - ref_s)
    return t[keep] - t0, M[keep], s0


def t_level(t, M, level, after=0.0):
    sel = np.where(t >= after)[0]; tt, mm = t[sel], M[sel]
    for j in np.where(mm >= level)[0]:
        if mm[j:].min() >= level - SUSTAIN:
            return float(tt[0]) if j == 0 else float(tt[j - 1] + (level - mm[j - 1]) / (mm[j] - mm[j - 1] + 1e-12) * (tt[j] - tt[j - 1]))
    return None


def analyse(cond):
    G = geometry(cond); px, vm, fps = read_video(cond, G); onset, ton, don, thr = detect_onset(px, fps, G); end_t = len(px) / fps
    t_seal, t_inv, inv = detect_seal(px, fps, onset, G)
    res = dict(cond=cond, G=G, fps=fps, onset=onset, t_seal=t_seal, inventory=(t_inv, inv), end_t=end_t, vm=vm, onset_trace=(ton, don, thr), curves={}, curves_alt={}, noise={})
    for k in GRID_K:
        for ref in GRID_REF:
            t, Ls = sigma_curve(px, fps, onset, G, k, ref)
            for arm_ex in (False, True):
                for ero in (0, 3):
                    sel = G['arm'] != 4                                                     # pump-valve arm excluded from sigma (valves flicker)
                    if arm_ex: sel &= G['arm'] != 1
                    if ero: sel &= G['edge'] >= ero
                    tt, M, s0 = index_from(t, Ls, t_seal, sel, end_t, ref)
                    res['curves'][(k, ref, arm_ex, ero)] = (tt, M, s0)
                    res['curves_alt'][(k, ref, arm_ex, ero)] = index_from(t, Ls, t_seal, sel, end_t, ref, reference='measured')
            if (k, ref) == (K_PRIMARY, REF_PRIMARY):
                res['L_primary'] = (t, Ls)
    # noise floor of sigma_L: two disjoint 1-s windows inside the last 3 s
    n = len(px); A = px[-int(3 * fps):-int(2 * fps)].astype(np.float32).mean(0); B_ = px[-int(1 * fps):].astype(np.float32).mean(0)
    dL = -np.log10(np.clip(smooth_norm(A, G, K_PRIMARY) / np.clip(smooth_norm(B_, G, K_PRIMARY), 1, None), 1e-3, None))
    res['noise_L'] = float(dL[(G['arm'] != 1) & (G['arm'] != 4)].std() / np.sqrt(2))
    res['px_first'] = None
    return res


def _win(c):
    t, M = c[0], c[1]; k = t <= T_WINDOW; return t[k], M[k]


def summarise(res):
    t, M, s0 = res['curves'][(K_PRIMARY, REF_PRIMARY, True, 0)]
    out = dict(cond=res['cond'], fps=res['fps'], duration_s=res['end_t'], onset_s=res['onset'], t_seal_s=res['t_seal'], after_seal_s=res['end_t'] - res['t_seal'],
               sigma_ref=s0, noise_floor_over_sigma_ref=res['noise_L'] / s0, M_start=float(np.interp(0.0, t, M)), M_last=float(M[-5:].mean()), usable_until_s=float(t[-1]))
    allc = list(res['curves'].values())
    for T in FIXED_T:
        out[f'M@{T}s'] = float(np.interp(T, t, M)) if T <= t[-1] else None                 # None = beyond the end of the record minus the reference window
        vals = [float(np.interp(T, c[0], c[1])) for c in allc if T <= c[0][-1]]
        out[f'M@{T}s_min'] = min(vals) if vals else None; out[f'M@{T}s_max'] = max(vals) if vals else None
    for lv in LEVELS:
        v = [t_level(*_win(c), lv) for c in allc]; p = t_level(*_win((t, M)), lv); ok = [x for x in v if x is not None]
        rep_ = len(ok) / len(v) >= 0.9 and bool(ok) and p is not None and max(ok) <= 1.5 * p                    # reported only if reached in >= 90 % of the settings
        out[f't{int(lv*100)}'] = p if rep_ else None; out[f't{int(lv*100)}_min'] = min(ok) if (ok and rep_) else None; out[f't{int(lv*100)}_max'] = max(ok) if (ok and rep_) else None
        out[f't{int(lv*100)}_n_reached'] = f'{len(ok)}/{len(v)}'
        out[f'measured_sigma0_t{int(lv*100)}'] = t_level(*_win(res['curves_alt'][(K_PRIMARY, REF_PRIMARY, True, 0)]), lv)
    return out


# ------------------------------------------------------------------------------------------------ figures
def draw_curves(ax, results, xmax=T_WINDOW):
    """Mixing index of the three runs: line = primary analysis, band = min-max over the 48 analysis settings (NOT a confidence interval)."""
    for r in results:
        c = r['cond']; col = COLORS[c]; t, M, s0 = r['curves'][(K_PRIMARY, REF_PRIMARY, True, 0)]
        grid = np.linspace(0, t[-1], 600); lo = np.full(600, np.inf); hi = np.full(600, -np.inf)
        for tt, mm, _ in r['curves'].values():
            y = np.interp(grid, tt, mm, right=np.nan); lo = np.fmin(lo, y); hi = np.fmax(hi, y)
        ok = np.isfinite(lo) & np.isfinite(hi)
        ax.fill_between(grid[ok], lo[ok], hi[ok], color=col, alpha=0.16, lw=0); ax.plot(t, M, color=col, lw=1.5, ls=LSTYLE[c])
    ax.axhline(0.5, color=C_TEXT, lw=0.6, ls=':'); ax.axhline(0.9, color=C_TEXT, lw=0.6, ls=':')
    ax.set_xlim(0, xmax); ax.set_ylim(-0.1, 1.02); ax.set_xticks(np.arange(0, xmax + 1, 10)); ax.set_xlabel('Time from loop sealing (s)'); ax.set_ylabel('Mixing index $M$')
    ax.grid(True, ls=':', lw=0.35, color='#CCCCCC', alpha=0.7)


def curve_legend(ax, results, fs=None):
    """One-row legend ABOVE the axes (never touches data)."""
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=COLORS[r['cond']], lw=1.5, ls=LSTYLE[r['cond']], label=hz(r['cond'])) for r in results]
    ax.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, fontsize=fs or BASE_FS - 1.4, handlelength=2.0, handletextpad=0.4, columnspacing=1.0, borderaxespad=0.2)


def fig_diagnostics(results):
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8)); fig.subplots_adjust(left=0.08, right=0.99, top=0.93, bottom=0.10, wspace=0.50, hspace=0.62)
    ax = axs[0, 0]
    for r in results:
        t, d, thr = r['onset_trace']; ax.semilogy(t, d, color=COLORS[r['cond']], lw=0.8, label=hz(r['cond'])); ax.plot([r['onset']], [thr], 'v', color=COLORS[r['cond']], ms=4)
    ax.set_xlim(0, 8); ax.set_xlabel('Time from start of video (s)'); ax.set_ylabel('RMS change of channel (counts)'); ax.legend(frameon=False, fontsize=BASE_FS - 1.5); ax.set_title('Pump onset (triangles)', fontsize=BASE_FS, loc='left')
    ax = axs[0, 1]
    for r in results:
        vm, fps = r['vm'], r['fps']; w0 = int((r['onset'] + 4) * fps); hp = vm[w0:] - uniform_filter1d(vm[w0:], int(1.5 * fps), axis=0)
        f, P = welch(detrend(hp[:, 0]), fs=fps, nperseg=min(len(hp), int(20 * fps)), noverlap=int(10 * fps)); ax.semilogy(f, P / P.max(), color=COLORS[r['cond']], lw=0.8, label=hz(r['cond']))
    for r in results:
        ax.axvline(1 / (6 * PHASE_S[r['cond']]), color=COLORS[r['cond']], lw=0.8, ls='--')                       # delivered (pump gravimetry)
        ax.axvline(PHASE_HZ[r['cond']] / 6, color=COLORS[r['cond']], lw=0.8, ls=':')                             # programmed
    ax.set_xlim(0, 12.5); ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel('Valve-1 flicker (norm.)'); ax.set_title('Valve flicker spectrum', fontsize=BASE_FS, loc='left')
    ax = axs[0, 2]
    for r in results:
        t, M, s0 = r['curves'][(K_PRIMARY, REF_PRIMARY, True, 0)]; ax.plot(t, 1 - M, color=COLORS[r['cond']], lw=0.9)
    ax.set_yscale('log'); ax.set_ylim(0.02, 1.3); ax.set_xlim(0, T_WINDOW); ax.set_xlabel('Time from loop sealing (s)'); ax.set_ylabel('$\\sigma_L / \\sigma_{ref}$'); ax.axhline(0.1, color='k', ls=':', lw=0.6)
    ax.set_title('Residual heterogeneity', fontsize=BASE_FS, loc='left')
    def strip(ax, fn, ylabel, title):
        for i, r in enumerate(results):
            v = [fn(c) for c in r['curves'].values()]; v = [x for x in v if x is not None and np.isfinite(x)]
            ax.scatter(np.full(len(v), i) + np.random.default_rng(0).uniform(-0.15, 0.15, len(v)), v, s=6, color=COLORS[r['cond']], alpha=0.7)
        ax.set_xticks(range(3)); ax.set_xticklabels([f"{hz(r['cond'])}\n{sum(1 for c in r['curves'].values() if (fn(c) is not None and np.isfinite(fn(c))))}/{len(r['curves'])}" for r in results], fontsize=BASE_FS - 1)
        ax.set_ylabel(ylabel); ax.set_title(title, fontsize=BASE_FS, loc='left'); ax.set_xlim(-0.5, 2.5)
    strip(axs[1, 0], lambda c: t_level(*_win(c), 0.5), '$t_{50}$ (s)', '$t_{50}$ over the settings')
    strip(axs[1, 1], lambda c: float(np.interp(20.0, c[0], c[1])) if c[0][-1] >= 20.0 else None, '$M$ at 20 s', '$M$ at 20 s over the settings')
    strip(axs[1, 2], lambda c: float(np.interp(0.0, c[0], c[1])), '$M$ at sealing', '$M$ at sealing over the settings')
    for ax, l in zip(axs.ravel(), 'abcdef'): ax.text(-0.22, 1.12, l, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
    return fig


def fig_maps(results, times=(5, 10, 20, 25)):
    fig, axs = plt.subplots(3, len(times), figsize=(6.4, 4.8)); fig.subplots_adjust(left=0.06, right=0.9, top=0.95, bottom=0.04, wspace=0.03, hspace=0.06)
    for r, row in zip(results, axs):
        G = r['G']; t, Ls = r['L_primary']
        for ax, T in zip(row, times):
            ax.axis('off')
            if r['t_seal'] + T > t[-1] - REF_PRIMARY:
                continue
            i = int(np.argmin(np.abs(t - (r['t_seal'] + T)))); img = np.full((G['hb'], G['wb']), np.nan, np.float32); img[G['ys'], G['xs']] = Ls[i]
            im = ax.imshow(img, cmap='RdBu_r', vmin=-0.08, vmax=0.08)
            if r is results[0]: ax.set_title(f't = {T} s', fontsize=BASE_FS)
        row[0].text(-0.05, 0.5, hz(r['cond']), transform=row[0].transAxes, rotation=90, ha='right', va='center', color=COLORS[r['cond']], fontweight='bold', fontsize=BASE_FS + 1)
    cax = fig.add_axes([0.92, 0.25, 0.015, 0.5]); fig.colorbar(im, cax=cax).set_label('$L$ (log$_{10}$ ratio to the final state)', fontsize=BASE_FS - 1)
    return fig


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)
    results = [analyse(c) for c in CONDS]; rows = [summarise(r) for r in results]
    for r, row in zip(results, rows): r['t50_report'] = row['t50']
    with open(os.path.join(OUT_DIR, 'mixing_summary.csv'), 'w', newline='') as fh:
        keys = list(rows[0].keys()); w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
    for r in results:
        t, M, s0 = r['curves'][(K_PRIMARY, REF_PRIMARY, True, 0)]; np.savetxt(os.path.join(OUT_DIR, f"mixing_curve_{r['cond'].replace(' ', '')}.csv"), np.c_[t, M], delimiter=',', header='t_from_onset_s,M', comments='')
    for r in rows: print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    f_s1 = fig_diagnostics(results); f_s2 = fig_maps(results)
    # ---- figure saving (uncomment to export)
    # f_s1.savefig('FigS9_mixing_diagnostics.jpeg', dpi=600, bbox_inches='tight')
    # f_s2.savefig('FigS10_mixing_maps.jpeg', dpi=600, bbox_inches='tight')
    plt.show()
