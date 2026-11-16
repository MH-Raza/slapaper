"""Table S8 - pump gravimetry, every weighing.
Data: pump xlsx (or results/pump_per_replicate.csv).
Writes tables/TableS8.csv and prints the table. Run from anywhere: python TableS8_pump_gravimetry.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

p = L.pump()
tab = ST.table_S8(p['df'])
tab.to_csv(os.path.join(C.TABLES, 'TableS8.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S8']); print(tab.to_string(index=False))
