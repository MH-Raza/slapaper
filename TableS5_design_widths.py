"""Table S5 - design silhouette widths used as reference for the offsets.
Data: CAD geometry only (no data).
Writes tables/TableS5.csv and prints the table. Run from anywhere: python TableS5_design_widths.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

tab = ST.table_S5()
tab.to_csv(os.path.join(C.TABLES, 'TableS5.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S5']); print(tab.to_string(index=False))
