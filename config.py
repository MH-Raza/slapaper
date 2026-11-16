"""Paths and switches for the code. Edit this file only.

Raw data (put them under data/ or change the paths):
    data/Data/                         60 profilometry CSV exports ('C1 CM1.csv', 'C1 CM1 PDMS.csv', ... 'F5 FM3 PDMS.csv')
    data/Valve data/                   125 valve logs ('F1 C1 500 mbar.txt', ...)
    data/F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt   replacement run for F3 x C4 at 1000 mbar
    data/pump_gravimetry.xlsx           pump gravimetry
    data/Mixing data/                  three aligned mixing videos + stills (33, 83, 167 ms folders)

If a raw-data folder is missing, the scripts that can work from the result tables use RESULTS_REFERENCE instead
(the tables produced from the raw data for the manuscript). Figures that draw raw scans, logs or video frames
(Fig. 2, Fig. 5, Figs S1-S4) need the raw data.
"""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, 'data')

PROFILOMETRY_DIR = os.path.join(DATA, 'Data')
VALVE_DIR = os.path.join(DATA, 'Valve data')
VALVE_EXTRA = {('F3', 'C4', 1000): os.path.join(DATA, 'F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt')}
PUMP_XLSX = os.path.join(DATA, 'pump_gravimetry.xlsx')
MIXING_DIR = os.path.join(DATA, 'Mixing data')
FRET_XLSX = os.path.join(DATA, 'FRET_Chip_and_tube.xlsx')    # sheet 'FRET eff': tube E8:G10, chip U10:W12

# ---- images
FIG6_PHOTO = os.path.join(ROOT, 'assets', 'fig6b_chip_photo.jpg')   # origami chip photograph for Fig. 6b; None = placeholder
FIG6_RENDERS = (os.path.join(ROOT, 'assets', 'fig6a_nanodevice_inactive.png'),   # Fig. 6a nanodevice renders (inactive, active);
                os.path.join(ROOT, 'assets', 'fig6a_nanodevice_active.png'))     # missing files = placeholders (see README)
FIG6C_UNITS = 'log2'       # Fig. 6c: log2 fold change vs origami alone (manuscript); also 'log10', 'log_axis', 'linear', 'dB'
FIG5_VIDEO_FRAMES = True   # Fig. 5d frames taken from the mixing videos (data/Mixing data/); placeholders only if the videos are missing

RESULTS = os.path.join(ROOT, 'results')                       # written by the scripts
RESULTS_REFERENCE = os.path.join(ROOT, 'results_reference')   # tables used for the manuscript (fallback)
FIGURES = os.path.join(ROOT, 'figures')                       # target of the (commented) savefig lines
TABLES = os.path.join(ROOT, 'tables')                         # CSV + Word output of the table scripts

for d in (RESULTS, FIGURES, TABLES):
    os.makedirs(d, exist_ok=True)
