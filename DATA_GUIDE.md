# Data guide

The raw data are in `data/`, as shown below. To use another location, change the paths in `config.py`. Names must match: the code finds files by name. `notebooks/00_Check_data.ipynb` checks all of this and lists anything missing.

```
data/
├── Data/                                        60 profilometry CSV files
├── Valve data/                                  124 valve logs (125 runs with the replacement file below)
├── F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt   (included)
├── pump_gravimetry.xlsx                          pump gravimetry
├── Mixing data/                                 three folders, one per mixing run
│   ├── 33 ms - Aligned_Cropped_Data/
│   ├── 83 ms - Aligned_Cropped_Data/
│   └── 167 ms - Aligned_Cropped_Data/
└── FRET_Chip_and_tube.xlsx                      (included)
```

## 1. Profilometry: `data/Data/` (60 files)

**Used by:** Figs 2, 3, S1, S2, S3, S6; Tables S1–S4 and S7. Figs 2 and S1–S3 cannot be drawn without these files.

**Names:** `<feature> <mold><n>[ Mold| PDMS].csv`.

- Features are C1–C5 (control molds CM1–CM3) and F1–F5 (flow molds FM1–FM3).
- Each feature has three mold chips and their three PDMS replicas: 10 features × 3 chips × 2 materials = 60 files.
- Examples: `C1 CM1.csv` (mold; the suffix may be left out), `C1 CM2 Mold.csv`, `C1 CM1 PDMS.csv`, `F5 FM3 PDMS.csv`.
- Extra spaces are accepted.

**Content:** KLA Zeta-20 CSV export with no header row.

- Column 1: row index.
- Column 2: lateral position x (µm).
- Columns 3 onwards: one scan line each, z (µm); usually 62–63 lines.
- PDMS replicas are valleys; the code inverts them.

## 2. Valve closure: `data/Valve data/` (124 files) + the replacement file

**Used by:** Figs 4, S4, S5, S6; Tables S6 and S7. Fig. S4 cannot be drawn without these files.

**Names:** `F<i> C<j> <P> mbar.txt`, with i, j = 1–5 and P = 500, 750, 1000, 1250, 1500 (for example `F1 C1 500 mbar.txt`).

- `F3 C4 1000 mbar.txt` was recorded at 1250 mbar, so it is excluded. If it is in the folder it is only listed in the QC table.
- The replacement run `F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt` is already in `data/`.

**Content:** tab-separated Fluigent log, latin-1, about 10 Hz, 180 s. There are seven columns:

1. Time [s]
2. OB1 ch1 target
3. OB1 ch1 read (drive pressure, 100 mbar)
4. OB1 ch2 target
5. OB1 ch2 read (control pressure, mbar)
6. MFS2 target
7. MFS2 read (flow, µL/min)

**Checks the code makes on every file:**

- the setpoint must match the pressure in the file name;
- the drive pressure must be 100 mbar;
- the run must last 180 s.

Restarted acquisitions are handled: the last continuous record is used.

## 3. Pump: `data/pump_gravimetry.xlsx`

**Used by:** Fig. 5a, b; Table S8.

**Sheets:** the first three sheets must be named `0.5 mm spacing 0.5 bar`, `1 bar` and `1.5 bar`; a notes sheet may follow.

**Each sheet:**

- Row 1: the ESI delays (25, 33, 48, 83, 167 ms) across the columns.
- Rows labelled in column A:
  - `Weight before (1)` … `(3)` and `Weight after (1)` … `(3)`, in mg;
  - three `Actual time` rows, as `mm:ss:cs`;
  - `Cycles`.

## 4. Mixing: `data/Mixing data/` (optional)

**Used by:** Fig. 5d (the 12 frames are taken from the videos) and Fig. 5e (the min–max band over the 48 analysis settings).

Without these folders, panel d is left as placeholders and panel e is drawn from `results_reference/mixing_curve_*.csv` (primary curve only).

| Folder | Video | Stills (.tiff) |
|---|---|---|
| `33 ms - Aligned_Cropped_Data` | `33ms_cropped_Mixing_video.mp4` | names containing `Stream_A`, `Stream_B`, `Blank_chip` |
| `83 ms - Aligned_Cropped_Data` | `83ms_cropped_Mixing_video.mp4` | same |
| `167 ms - Aligned_Cropped_Data` | `167ms_cropped_Mixing_video.mp4` | same |

The full video analysis takes 1–2 min.

## 5. FRET: `data/FRET_Chip_and_tube.xlsx` (included)

**Used by:** Fig. 6c. On the `FRET eff` sheet:

- tube replicates are cells E8:G10;
- chip replicates are cells U10:W12;
- columns are origami + sRNA mimic, origami, origami + random sequence.

## Chip images

| Panel | Source | Box |
|---|---|---|
| Fig. 6a, nanodevice renders (inactive, active) | `assets/fig6a_nanodevice_inactive.png`, `assets/fig6a_nanodevice_active.png` (`FIG6_RENDERS` in `config.py`; not in the public repository) | cropped by the code; placeholders of the same size if missing |
| Fig. 6b, origami chip photograph | `assets/fig6b_chip_photo.jpg` (included; `FIG6_PHOTO` in `config.py`) | width/height 1.38; the photo is centre-cropped to it |
| Fig. 5d, 12 video frames | the mixing videos (`FIG5_VIDEO_FRAMES = True`) | 298 × 293 px each at 600 dpi |
