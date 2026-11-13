"""Table S1 - result for each of the 60 profilometry scans.
Data: profilometry (raw scans, or results/profilometry_per_file.csv).
Writes tables/TableS1.csv and prints the table. Run from anywhere: python TableS1_profilometry_scans.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

d = L.profilometry()
tab = ST.table_S1(d['R'])
tab.to_csv(os.path.join(C.TABLES, 'TableS1.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S1']); print(tab.to_string(index=False))
