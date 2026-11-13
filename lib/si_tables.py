"""Supplementary Tables S1-S8, formatted exactly as printed in the SI (strings, Unicode minus, same rounding).
Each function takes the analysis outputs from loaders.py and returns a DataFrame of strings."""
import re
import numpy as np, pandas as pd
from scipy import stats
import profilometry_analysis as PRO
import valve_closure_analysis as VAL

MINUS = '−'


def f(v, nd=1):
    """Number -> string with nd decimals and a Unicode minus; '' for NaN."""
    if v is None or (isinstance(v, float) and np.isnan(v)): return ''
    s = f'{v:.{nd}f}'
    return MINUS + s[1:] if s.startswith('-') else s     # keeps '−0.0' for small negative values, as in the SI


def pm(m, s, nd=1): return f'{f(m, nd)} ± {f(s, nd)}'


def ci95(v):
    v = np.asarray(v, float); n = len(v)
    return stats.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n)


ORDER = PRO.CHANNELS  # C1..C5, F1..F5


# ---------------------------------------------------------------------------------------------------------- S1
def table_S1(R):
    """Result for each of the 60 profilometry scans."""
    R = R.copy()
    R['_m'] = (R.material == 'PDMS').astype(int); R['_c'] = R.channel.map(ORDER.index)
    rows = []
    for _, r in R.sort_values(['_c', '_m', 'chip']).iterrows():
        mold = re.match(r'\S+\s+([CF]M\d)', r.file).group(1)                   # 'CM1' / 'FM2' (file names vary: 'C1 CM1.csv', 'C1 CM2 Mold.csv')
        rem = []
        if r.footprint15_frac > 0.5: rem.append('footprint at 15% (uneven walls)')   # most scan lines needed the 15% footprint
        if abs(r.tilt_med_um_per_mm) > 5: rem.append('tilt > 5 µm/mm')
        rows.append({'Scan': f'{r.channel} {mold} {"PDMS" if r.material == "PDMS" else "Mold"}',
                     'Lines used / total': f'{int(r.lines_used)} / {int(r.lines_total)}',
                     'Median tilt (µm/mm)': f(r.tilt_med_um_per_mm), 'H (µm)': f(r.H), 'W10 (µm)': f(r.W10),
                     'W20 (µm)': f(r.W20), 'W50 (µm)': f(r.W50), 'W80 (µm)': f(r.W80), 'Remark': '; '.join(rem)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S2
def table_S2(F):
    """Per-feature statistics over three chips (µm)."""
    F = F.copy(); F['_c'] = F.channel.map(ORDER.index); F['_m'] = (F.material == 'PDMS').astype(int)
    rows = []
    for _, r in F.sort_values(['_c', '_m']).iterrows():
        rows.append({'Feature': r.channel, 'Material': r.material,
                     'H mean ± SD (95% CI)': f'{pm(r.H_mean, r.H_sd)} ({f(r.H_ci95)})',
                     'W10 mean ± SD': pm(r.W10_mean, r.W10_sd),
                     'ΔW10 mean ± SD (95% CI)': f'{pm(r.dW10_mean, r.dW10_sd)} ({f(r.dW10_ci95)})',
                     'W50 mean ± SD': pm(r.W50_mean, r.W50_sd), 'ΔW50 mean ± SD': pm(r.dW50_mean, r.dW50_sd)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S3
def table_S3(T_D):
    """Height change from mold to PDMS replica for each chip pair (µm)."""
    rows = []
    for ch in ORDER:
        s = T_D[T_D.channel == ch].sort_values('chip')
        rows.append({'Feature': ch, **{f'Chip {int(c)}': f(v) for c, v in zip(s.chip, s.dH)},
                     'Mean ± 95% CI (2 d.f.)': pm(s.dH.mean(), ci95(s.dH))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S4
def table_S4(TS):
    """Pooled footprint offset of the rectangles C1-C3 for different width definitions (mean ± 95% CI over three chips)."""
    rows = []
    for _, r in TS.iterrows():
        lab = f'{int(float(r.threshold_pct))}% of H' if str(r.threshold_pct).replace('.', '').isdigit() else 'area-based Weq'
        rows.append({'Width threshold': lab, 'C1–C3 mold offset (µm)': pm(r.rect_offset_mold, r.rect_offset_mold_ci, 0),
                     'C1–C3 PDMS offset (µm)': pm(r.rect_offset_pdms, r.rect_offset_pdms_ci, 0)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S5
def table_S5():
    """Design silhouette widths used as reference for the offsets (µm); computed from the CAD geometry."""
    rows = []
    for ch in ORDER:
        H, W = PRO.CAD[ch]
        if ch.startswith('C'): wmax = W
        else:
            R_ = H / 2 + W ** 2 / (8 * H); wmax = 2 * R_ if H > W / 2 else W     # a major arc is widest at its centre height
        rows.append({'Feature': ch, 'Design W at 10% (µm)': f(PRO._cad_width(ch, 0.10)), 'at 20%': f(PRO._cad_width(ch, 0.20)),
                     'at 50%': f(PRO._cad_width(ch, 0.50)), 'at 80%': f(PRO._cad_width(ch, 0.80)), 'Maximum width': f(wmax)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S6
def table_S6(G):
    """Closure of all 125 valve runs (primary windows)."""
    rows = []
    for _, r in G.sort_values(['F', 'C', 'P']).iterrows():
        rows.append({'Run': f'{r.F}×{r.C} {int(r.P)} mbar',
                     'Closure mean ± SD (min–max) (%)': f'{pm(r.closure, r.sd)} ({f(r.cmin)} to {f(r.cmax)})',
                     'Open flow (µL/min)': f(r.Q_open, 2)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S7
def table_S7(G, prof_path):
    """Minimum sealing pressure and as-fabricated aspect ratio for four definitions of the width."""
    T = VAL.pmin_table(G)
    A = {w: VAL.add_aspect_ratio(T, prof_path, w).set_index(['F', 'C']) for w in ('W10', 'W20', 'W50', 'W80')}
    rows = []
    for _, r in T.sort_values(['F', 'C']).iterrows():
        k = (r.F, r.C)
        row = {'Valve': f'{r.F}×{r.C}', 'Pmin (mbar)': 'no seal at ≤ 1500' if np.isnan(r.P_min) else f'{int(r.P_min)}'}
        for w in ('W10', 'W20', 'W50', 'W80'):
            row[f'ARfab ({w})'] = pm(A[w].loc[k, 'AR_fab'], A[w].loc[k, 'AR_fab_sd'], 2)
        row['AR design'] = f(A['W50'].loc[k, 'AR_design'], 2)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------- S8
def table_S8(df):
    """Pump gravimetry: every weighing."""
    rows = []
    for _, r in df.sort_values(['bar', 'delay_ms', 'rep'], ascending=[True, True, True]).iterrows():
        rows.append({'Pressure (bar)': f'{r.bar:g}', 'ESI delay (ms)': f'{int(r.delay_ms)}', 'Programmed phase frequency (Hz)': f(r.f_prog_Hz),
                     'Replicate': f'{int(r.rep)}', 'Cycles': f'{int(r.cycles)}', 'Mass difference (mg)': f'{r.dm_mg:g}',
                     'Elapsed time (s)': f(r.t_s, 2), 'Stroke volume (nL)': f(r.stroke_nL), 'Flow rate (µL/min)': f(r.flow_uL_min, 2),
                     'Delivered phase frequency (Hz)': f(r.f_phase_Hz)})
    return pd.DataFrame(rows)


TITLES = {
    'S1': 'Table S1. Result for each of the 60 profilometry scans.',
    'S2': 'Table S2. Per-feature statistics over three chips (µm).',
    'S3': 'Table S3. Height change from mold to PDMS replica for each chip pair (µm).',
    'S4': 'Table S4. Pooled footprint offset of the resolved rectangles C1–C3 for different width definitions (mean ± 95% CI over three chips).',
    'S5': 'Table S5. Design silhouette widths used as reference for the offsets (µm).',
    'S6': 'Table S6. Closure of all 125 valve runs (primary windows).',
    'S7': 'Table S7. Minimum sealing pressure and as-fabricated aspect ratio for four definitions of the width.',
    'S8': 'Table S8. Pump gravimetry.',
}


def to_docx(tables, path):
    """All tables in one Word file (Times New Roman 9 pt, header row bold), ready to paste into the SI."""
    from docx import Document
    from docx.shared import Pt
    doc = Document(); st = doc.styles['Normal']; st.font.name = 'Times New Roman'; st.font.size = Pt(9)
    for key, df in tables.items():
        p = doc.add_paragraph(); r = p.add_run(TITLES[key]); r.bold = True
        t = doc.add_table(rows=1, cols=len(df.columns)); t.style = 'Table Grid'
        for c, name in zip(t.rows[0].cells, df.columns):
            c.text = ''; c.paragraphs[0].add_run(str(name)).bold = True
        for _, row in df.iterrows():
            cells = t.add_row().cells
            for c, v in zip(cells, row): c.text = str(v)
        doc.add_paragraph()
    cp = doc.core_properties                                  # replace python-docx's template properties (author 'python-docx', 2013)
    import datetime
    cp.author = cp.last_modified_by = 'Muhammad Hassan Raza'; cp.title = 'Supplementary Tables S1-S8'; cp.comments = ''
    cp.created = cp.modified = datetime.datetime.now(); cp.revision = 1
    doc.save(path)
