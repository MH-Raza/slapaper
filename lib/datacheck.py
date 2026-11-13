"""Check the Python packages and the data before running the figures and tables.

check() prints, for every dataset, what was found against what the code expects, then for every figure and table whether it
will be built from the raw data, from the result tables (results_reference/), or not at all. Nothing is modified.
"""
import os, re, glob, importlib
import config as C

PACKAGES = [('numpy', 'numpy'), ('pandas', 'pandas'), ('scipy', 'scipy'), ('matplotlib', 'matplotlib'), ('PIL', 'pillow'),
            ('cv2', 'opencv-python'), ('openpyxl', 'openpyxl'), ('statsmodels', 'statsmodels'), ('docx', 'python-docx')]

CHANNELS = ['C1', 'C2', 'C3', 'C4', 'C5', 'F1', 'F2', 'F3', 'F4', 'F5']
PRESSURES = [500, 750, 1000, 1250, 1500]
MIX_CONDS = {'33 ms': ('33 ms - Aligned_Cropped_Data', '33ms_cropped_Mixing_video.mp4'),
             '83 ms': ('83 ms - Aligned_Cropped_Data', '83ms_cropped_Mixing_video.mp4'),
             '167 ms': ('167 ms - Aligned_Cropped_Data', '167ms_cropped_Mixing_video.mp4')}
PUMP_SHEETS = ['0.5 mm spacing 0.5 bar', '1 bar', '1.5 bar']


def packages(show=True):
    out = {}
    for mod, pip in PACKAGES:
        try:
            m = importlib.import_module(mod); out[pip] = getattr(m, '__version__', 'installed')
        except ImportError:
            out[pip] = None
    if show:
        print('Python packages')
        for k, v in out.items():
            print(f'  {"ok     " if v else "MISSING"}  {k:<15s} {v or "pip install " + k}')
    return out


def _profilometry():
    d = C.PROFILOMETRY_DIR
    if not os.path.isdir(d):
        return dict(ok=False, found=0, expected=60, note=f'folder not found: {d}')
    pat = re.compile(r'^(C\d|F\d)\s+((?:CM|FM)(\d))(?:\s+(Mold|PDMS))?\.csv$')
    have = set()
    for f in os.listdir(d):
        m = pat.match(f)
        if m: have.add((m.group(1), int(m.group(3)), m.group(4) or 'Mold'))
    want = {(c, k, mat) for c in CHANNELS for k in (1, 2, 3) for mat in ('Mold', 'PDMS')}
    miss = sorted(want - have)
    note = '' if not miss else 'missing: ' + ', '.join(f'{c} {"CM" if c[0] == "C" else "FM"}{k} {m}' for c, k, m in miss[:8]) + \
                               (' ...' if len(miss) > 8 else '')
    return dict(ok=not miss, found=len(have & want), expected=60, note=note)


def _valves():
    d = C.VALVE_DIR
    if not os.path.isdir(d):
        return dict(ok=False, found=0, expected=125, note=f'folder not found: {d}')
    names = set(os.path.basename(p) for p in glob.glob(os.path.join(d, '*.txt')))
    want = {f'F{i} C{j} {p} mbar.txt' for i in range(1, 6) for j in range(1, 6) for p in PRESSURES} - {'F3 C4 1000 mbar.txt'}
    extra = list(C.VALVE_EXTRA.values())[0]
    miss = sorted(want - names) + ([] if os.path.exists(extra) else [os.path.basename(extra) + ' (replacement run, in data/)'])
    note = '' if not miss else 'missing: ' + ', '.join(miss[:6]) + (' ...' if len(miss) > 6 else '')
    return dict(ok=not miss, found=len(want & names) + os.path.exists(extra), expected=125, note=note)


def _pump():
    p = C.PUMP_XLSX
    if not os.path.exists(p):
        return dict(ok=False, found=0, expected=1, note=f'file not found: {p}')
    import openpyxl
    titles = [ws.title for ws in openpyxl.load_workbook(p, read_only=True).worksheets[:3]]
    ok = titles == PUMP_SHEETS
    return dict(ok=ok, found=1, expected=1, note='' if ok else f'first three sheets must be {PUMP_SHEETS}, found {titles}')


def _mixing():
    d = C.MIXING_DIR
    if not os.path.isdir(d):
        return dict(ok=False, found=0, expected=3, note=f'folder not found: {d}')
    good, notes = 0, []
    for cond, (folder, video) in MIX_CONDS.items():
        f = os.path.join(d, folder)
        tif = [os.path.basename(p).lower() for p in glob.glob(os.path.join(f, '*')) if p.lower().endswith('.tiff')]
        need = [video] if not os.path.exists(os.path.join(f, video)) else []
        need += [k + '*.tiff' for k in ('stream_a', 'stream_b', 'blank_chip') if not any(k in t for t in tif)]
        if need: notes.append(f'{folder}: missing {", ".join(need)}')
        else: good += 1
    return dict(ok=good == 3, found=good, expected=3, note='; '.join(notes))


def _fret():
    p = C.FRET_XLSX
    if not os.path.exists(p):
        return dict(ok=False, found=0, expected=1, note=f'file not found: {p}')
    import openpyxl
    ok = 'FRET eff' in openpyxl.load_workbook(p, read_only=True).sheetnames
    return dict(ok=ok, found=1, expected=1, note='' if ok else "sheet 'FRET eff' not found")


DATASETS = {'profilometry': ('60 profilometry CSV files', _profilometry), 'valves': ('125 valve logs', _valves),
            'pump': ('pump gravimetry workbook', _pump), 'mixing': ('3 mixing videos + stills', _mixing),
            'fret': ('FRET workbook', _fret)}

# output -> list of (dataset, 'raw' = required / 'table' = result tables used without it / 'part' = partial without it)
NEEDS = {
    'Fig. 1': [], 'Fig. 2': [('profilometry', 'raw')], 'Fig. 3': [('profilometry', 'table')],
    'Fig. 4': [('valves', 'table'), ('profilometry', 'table')], 'Fig. 5': [('pump', 'table'), ('mixing', 'part')],
    'Fig. 6': [('fret', 'table')],
    'Fig. S1': [('profilometry', 'raw')], 'Fig. S2': [('profilometry', 'raw')], 'Fig. S3': [('profilometry', 'raw')],
    'Fig. S4': [('valves', 'raw')], 'Fig. S5': [('valves', 'table')], 'Fig. S6': [('valves', 'table'), ('profilometry', 'table')],
    'Tables S1-S5': [('profilometry', 'table')], 'Tables S6-S7': [('valves', 'table'), ('profilometry', 'table')],
    'Table S8': [('pump', 'table')],
}
PART_NOTE = {'mixing': 'without the mixing videos: panel d = placeholders, panel e without the min-max band'}


def check(show=True):
    packages(show)
    st = {k: fn() for k, (_, fn) in DATASETS.items()}
    if show:
        print('\nData (paths in config.py)')
        for k, (desc, _) in DATASETS.items():
            s = st[k]
            print(f'  {"ok     " if s["ok"] else ("PARTIAL" if s["found"] else "MISSING")}  {desc:<28s} {s["found"]}/{s["expected"]}'
                  + (f'   {s["note"]}' if s['note'] else ''))
        print('\nWhat will be built')
        for out, needs in NEEDS.items():
            if not needs:
                msg = 'ok (no data needed)'
            elif any(m == 'raw' and not st[d]['ok'] for d, m in needs):
                msg = 'NEEDS RAW DATA: ' + ', '.join(DATASETS[d][0] for d, m in needs if m == 'raw' and not st[d]['ok'])
            else:
                src = ['raw data' if st[d]['ok'] else 'result tables' for d, m in needs if m != 'part']
                msg = 'ok from ' + ' + '.join(dict.fromkeys(src))
                part = [PART_NOTE[d] for d, m in needs if m == 'part' and not st[d]['ok']]
                if part: msg += '; ' + '; '.join(part)
            print(f'  {out:<13s} {msg}')
        ph = getattr(C, 'FIG6_PHOTO', None)
        rd = [r for r in (getattr(C, 'FIG6_RENDERS', None) or ()) if r and os.path.exists(r)]
        print('\nImages: Fig. 6a nanodevice renders ' + (f'found ({len(rd)} of 2)' if rd else
              'not found: placeholders are drawn (see README, FIG6_RENDERS in config.py)') + '; Fig. 6b chip photograph ' +
              (f'found ({os.path.basename(ph)})' if ph and os.path.exists(ph) else 'not found: a placeholder is drawn (FIG6_PHOTO in config.py)') +
              '; Fig. 5d frames are taken from the mixing videos.')
    return st
