"""Table S7 - minimum sealing pressure and as-fabricated aspect ratio (four width definitions).
Data: valve + profilometry (raw or result tables).
Writes tables/TableS7.csv and prints the table. Run from anywhere: python TableS7_Pmin_aspect_ratio.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

v = L.valves()
tab = ST.table_S7(v['G'], v['prof_path'])
tab.to_csv(os.path.join(C.TABLES, 'TableS7.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S7']); print(tab.to_string(index=False))
