"""Table S3 - height change from mold to PDMS replica for each chip pair.
Data: profilometry (raw or result table).
Writes tables/TableS3.csv and prints the table. Run from anywhere: python TableS3_height_transfer.py"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
import pandas as pd
import config as C
import loaders as L
import si_tables as ST

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20); pd.set_option('display.max_rows', 200)

d = L.profilometry()
tab = ST.table_S3(d['T_D'])
tab.to_csv(os.path.join(C.TABLES, 'TableS3.csv'), index=False, encoding='utf-8-sig')   # utf-8-sig: opens correctly in Excel
print(ST.TITLES['S3']); print(tab.to_string(index=False))
