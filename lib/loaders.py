"""Data access for every figure and table script: run the analysis on the raw data when it is present, otherwise read the
result tables (results_reference/). Every loader also writes its result tables to results/ so the table scripts use the
latest numbers."""
import os, sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import config as C
import profilometry_analysis as PRO
import valve_closure_analysis as VAL
import pump_analysis as PUMP


def _have(path):
    return os.path.exists(path) and (not os.path.isdir(path) or len(os.listdir(path)) > 0)


def _raw(name):
    """True if the raw dataset is complete, False if it is absent; an incomplete dataset raises an error that lists what
    is missing (the analysis would otherwise run on part of the data without saying so)."""
    import datacheck
    st = datacheck.DATASETS[name][1]()
    if st['ok']:
        return True
    if st['found'] == 0 and 'not found' in st['note']:
        return False
    raise FileNotFoundError(f'{datacheck.DATASETS[name][0]}: incomplete ({st["found"]}/{st["expected"]}). {st["note"]} '
                            f'Complete the folder, or move it away to use the result tables.')


def raw_status(show=True):
    """Which raw datasets are complete (paths from config.py). Returns {'profilometry', 'valves', 'pump', 'mixing', 'fret': bool}.
    For the detailed check (packages, missing files, what each figure will use) run datacheck.check()."""
    import datacheck
    st = {k: fn() for k, (_, fn) in datacheck.DATASETS.items()}
    if show:
        for k, (desc, _) in datacheck.DATASETS.items():
            print(f'{"found  " if st[k]["ok"] else "missing"}  {desc}' + (f'  ({st[k]["note"]})' if st[k]['note'] else ''))
    return {k: v['ok'] for k, v in st.items()}


def _csv(name):
    """Latest result table: results/ (written by a previous raw-data run) before results_reference/."""
    for d in (C.RESULTS, C.RESULTS_REFERENCE):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


def _need(path, what):
    raise FileNotFoundError(f'{what} needs the raw data at {path} (see config.py).')


# ----------------------------------------------------------------------------------------------------------- profilometry
def profilometry(need_raw=False):
    """Returns dict R, store, F, PO, FIT, T_D, T_F, T_L, TS. store (scan lines) is None when read from the result tables."""
    if _raw('profilometry'):
        R, store = PRO.run(C.PROFILOMETRY_DIR); src = 'raw data'
        R.to_csv(os.path.join(C.RESULTS, 'profilometry_per_file.csv'), index=False)
    elif need_raw:
        _need(C.PROFILOMETRY_DIR, 'This figure')
    else:
        R, store = pd.read_csv(_csv('profilometry_per_file.csv')), None; src = 'result table'
    F = PRO.per_feature(R); PO = PRO.pooled_offsets(R); FIT = PRO.flow_fits(R); T_D, T_F, T_L = PRO.transfer(R)
    TS = PRO.threshold_sensitivity(R)
    F.to_csv(os.path.join(C.RESULTS, 'profilometry_per_feature.csv'), index=False)
    T_D.to_csv(os.path.join(C.RESULTS, 'height_transfer_pairs.csv'), index=False)
    TS.to_csv(os.path.join(C.RESULTS, 'threshold_sensitivity.csv'), index=False)
    print(f'profilometry: {len(R)} scans from the {src}')
    return dict(R=R, store=store, F=F, PO=PO, FIT=FIT, T_D=T_D, T_F=T_F, T_L=T_L, TS=TS)


# ----------------------------------------------------------------------------------------------------------- valves
def valves(need_raw=False):
    """Returns dict D (per cycle), Q (file checks), G (per run), T (P_min and as-fabricated aspect ratio), prof_path."""
    if _raw('valves'):
        D, Q, G = VAL.analyse(C.VALVE_DIR, C.VALVE_EXTRA); src = 'raw data'
        for df, n in ((D, 'valve_closure_per_cycle.csv'), (Q, 'valve_file_QC.csv'), (G, 'valve_closure_per_run.csv')):
            df.to_csv(os.path.join(C.RESULTS, n), index=False)
    elif need_raw:
        _need(C.VALVE_DIR, 'This figure')
    else:
        D = pd.read_csv(_csv('valve_closure_per_cycle.csv')); Q = pd.read_csv(_csv('valve_file_QC.csv'))
        G = pd.read_csv(_csv('valve_closure_per_run.csv')); src = 'result table'
    prof_path = _csv('profilometry_per_file.csv')
    T = VAL.add_aspect_ratio(VAL.pmin_table(G), prof_path)
    T.to_csv(os.path.join(C.RESULTS, 'valve_Pmin_aspect_ratio.csv'), index=False)
    VAL.separation(T, prof_path=prof_path).to_csv(os.path.join(C.RESULTS, 'aspect_ratio_separation.csv'), index=False)
    print(f'valves: {len(G)} runs from the {src}')
    return dict(D=D, Q=Q, G=G, T=T, prof_path=prof_path)


# ----------------------------------------------------------------------------------------------------------- pump
def pump():
    """Returns dict df (every weighing) and S (per condition)."""
    if _raw('pump'):
        df = PUMP.load(C.PUMP_XLSX); src = 'raw data'
    else:
        df = pd.read_csv(_csv('pump_per_replicate.csv')); src = 'result table'
    S = PUMP.summarise(df)
    df.to_csv(os.path.join(C.RESULTS, 'pump_per_replicate.csv'), index=False); S.to_csv(os.path.join(C.RESULTS, 'pump_summary.csv'), index=False)
    print(f'pump: {len(df)} weighings from the {src}')
    return dict(df=df, S=S)


# ----------------------------------------------------------------------------------------------------------- FRET (Fig. 6)
FRET_CONDITIONS = ['Origami + sRNA mimic', 'Origami', 'Origami + random sequence']
FRET_CELLS = {'Tube': ('E', 'F', 'G', (8, 9, 10)), 'Chip': ('U', 'V', 'W', (10, 11, 12))}   # three replicates per condition


def fret():
    """Every FRET-efficiency replicate (long table: setting, condition, replicate, fret). Read from the workbook cells,
    not from its charts: the Excel charts plotted the first replicate as the bar height (with the SD of all three as the
    error bar), the figure uses the mean."""
    if _raw('fret'):
        import openpyxl
        ws = openpyxl.load_workbook(C.FRET_XLSX, data_only=True)['FRET eff']
        rows = []
        for setting, (*cols, reps) in FRET_CELLS.items():
            for col, cond in zip(cols, FRET_CONDITIONS):
                for i, r in enumerate(reps, 1):
                    rows.append({'setting': setting, 'condition': cond, 'replicate': i, 'fret': float(ws[f'{col}{r}'].value)})
        df = pd.DataFrame(rows); src = 'raw data'
    else:
        df = pd.read_csv(_csv('fret_replicates.csv')); src = 'result table'
    df.to_csv(os.path.join(C.RESULTS, 'fret_replicates.csv'), index=False)
    print(f'FRET: {len(df)} values from the {src}')
    return df


# ----------------------------------------------------------------------------------------------------------- mixing
def mixing(need_raw=False):
    """Mixing runs for Fig. 5e. With the videos: the full video analysis (line + band over the 48 analysis settings; ~1-2 min).
    Without them: the primary mixing curves from the result tables (line only). Returns one dict per condition."""
    import numpy as np
    import mixing_analysis as MIX
    key = (MIX.K_PRIMARY, MIX.REF_PRIMARY, True, 0)
    if _raw('mixing'):
        MIX.DATA_DIR = C.MIXING_DIR
        res = [MIX.analyse(c) for c in MIX.CONDS]; rows = [MIX.summarise(r) for r in res]
        for r, row in zip(res, rows): r['t50_report'] = row['t50']
        pd.DataFrame(rows).to_csv(os.path.join(C.RESULTS, 'mixing_summary.csv'), index=False)
        for r in res:
            t, M, _ = r['curves'][key]
            pd.DataFrame({'t_from_onset_s': t, 'M': M}).to_csv(os.path.join(C.RESULTS, f"mixing_curve_{r['cond'].replace(' ', '')}.csv"), index=False)
        print('mixing: video analysis of the raw data')
        return res
    if need_raw:
        _need(C.MIXING_DIR, 'The mixing video analysis')
    res = []
    for cond in MIX.CONDS:
        d = pd.read_csv(_csv(f"mixing_curve_{cond.replace(' ', '')}.csv"))
        # column 0 is time from loop sealing (s); the header name 't_from_onset_s' is kept from the original export
        res.append(dict(cond=cond, curves={key: (d.iloc[:, 0].values, d.iloc[:, 1].values, np.nan)}, from_table=True))
    print('mixing: curves from the result tables (the min-max band of Fig. 5e needs the videos)')
    return res
