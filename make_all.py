"""Build every table (CSV + one Word file) and every figure that the available data allow. Figures are only displayed
unless the save lines in the individual figure scripts are uncommented; set SAVE_FIGURES = True here to export them all."""
import os, sys, runpy, traceback
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import config as C
import loaders as L
import si_tables as ST
from figio import save_tiff

SAVE_FIGURES = False

# ---------------------------------------------------------------- tables
pro, val, pmp = L.profilometry(), L.valves(), L.pump()
tables = {'S1': ST.table_S1(pro['R']), 'S2': ST.table_S2(pro['F']), 'S3': ST.table_S3(pro['T_D']), 'S4': ST.table_S4(pro['TS']),
          'S5': ST.table_S5(), 'S6': ST.table_S6(val['G']), 'S7': ST.table_S7(val['G'], val['prof_path']), 'S8': ST.table_S8(pmp['df'])}
for k, t in tables.items():
    t.to_csv(os.path.join(C.TABLES, f'Table{k}.csv'), index=False, encoding='utf-8-sig')
ST.to_docx(tables, os.path.join(C.TABLES, 'SI_Tables_S1-S8.docx'))
print('tables written to', C.TABLES)

# ---------------------------------------------------------------- figures
here = os.path.dirname(os.path.abspath(__file__))
scripts = ['Fig1_fabrication.py', 'Fig2_regimes_profiles.py', 'Fig3_profilometry.py', 'Fig4_valve_closure.py', 'Fig5_pump_mixer.py', 'Fig6_origami.py',
           'FigS1_profiles_control.py', 'FigS2_profiles_flow.py', 'FigS3_design_vs_measured.py', 'FigS4_valve_traces.py',
           'FigS5_closure_vs_pressure.py', 'FigS6_aspect_ratio.py']
plt.show = lambda *a, **k: None
for s in scripts:
    try:
        ns = runpy.run_path(os.path.join(here, s), run_name='__main__')
        out = s.split('_')[0] + '.tif'
        if SAVE_FIGURES: save_tiff(ns['fig'], os.path.join(C.FIGURES, out))
        print('OK     ', s)
    except FileNotFoundError as e:
        print('SKIPPED', s, '-', e)
    except Exception:
        print('FAILED ', s); traceback.print_exc()
    plt.close('all')
