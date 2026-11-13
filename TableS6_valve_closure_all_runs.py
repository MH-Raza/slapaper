"""Table S6 - closure of all 125 valve runs.
Data: valve logs (or results/valve_closure_per_run.csv).
Writes tables/TableS6.csv and prints the table. Run from anywhere: python TableS6_valve_closure_all_runs.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

v = L.valves()
tab = ST.table_S6(v['G'])
tab.to_csv(os.path.join(C.TABLES, 'TableS6.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S6']); print(tab.to_string(index=False))
