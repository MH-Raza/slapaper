"""Table S4 - pooled C1-C3 footprint offset for different width definitions.
Data: profilometry (raw or result table).
Writes tables/TableS4.csv and prints the table. Run from anywhere: python TableS4_width_definitions.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

d = L.profilometry()
tab = ST.table_S4(d['TS'])
tab.to_csv(os.path.join(C.TABLES, 'TableS4.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S4']); print(tab.to_string(index=False))
