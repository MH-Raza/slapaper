# Stereolithography-enabled multilayer soft lithographic microfluidic valve networks

Analysis and figure code for:

> M. H. Raza, M. Ebrahimian, M. Olowe, S. Aravamudhan, S. Desai, S. Oliveira and R. Zadegan, Stereolithography-enabled multilayer soft lithographic microfluidic valve networks, submitted to *Lab on a Chip*, 2026.

The code reproduces every main-text figure (Figs 1–6), every Supplementary Figure (Figs S1–S6) and every Supplementary Table (Tables S1–S8) from the raw measurements:

- optical profilometry of laser-SLA molds and PDMS replicas;
- valve closure for 25 geometries × 5 actuation pressures (125 runs);
- peristaltic-pump gravimetry;
- mixing videos of the recirculating mixer;
- FRET of a DNA origami nanodevice on chip and in tube.

Each output has a Jupyter notebook (`notebooks/`) and an identical Python script (repository root).

## Installation

Python 3.9 or newer.

```bash
git clone https://github.com/SBMI-LAB/slapaper.git
cd slapaper
pip install -r requirements.txt
```

or, with conda:

```bash
conda env create -f environment.yml
conda activate sla-soft-litho
```

Tested with Python 3.13, NumPy 2.5, pandas 3.0, SciPy 1.18, Matplotlib 3.11, OpenCV 5.0 and statsmodels 0.15 (Linux).

## Data

**Included in this repository:**

| Path | Content |
|---|---|
| `data/FRET_Chip_and_tube.xlsx` | FRET efficiencies, three replicates per condition, on chip and in tube |
| `data/F3_C4_-_1000_millibars_valve_-_100_millibars_flow.txt` | replacement valve run F3 × C4 at 1000 mbar |
| `results_reference/` | result tables computed from the raw data (per scan, per feature, per valve run, per pump replicate, mixing curves) |
| `assets/` | pump/mixer schematic (Fig. 5c), chip layout and chip photograph (Fig. 6b) |

**Raw data** (profilometry scans, valve logs, pump gravimetry and mixing videos) are in `data/`, organised as described in [DATA_GUIDE.md](DATA_GUIDE.md).

**Without the raw data:**

- These outputs are built from `results_reference/`: Figs 1, 3, 4, 5 and 6, Figs S5 and S6, and Tables S1–S8.
- These outputs need the raw data: Figs 2 and S1–S4, the video frames of Fig. 5d, and the min–max band of Fig. 5e.

## Usage

**Jupyter:** open `notebooks/` and run, in order:

1. `00_Check_data.ipynb` lists which data were found and what each output will use. It changes nothing.
2. `00_All_figures_and_tables.ipynb` builds all tables and figures.
3. Alternatively, run `Fig1.ipynb` … `FigS6.ipynb` and `TableS1.ipynb` … `TableS8.ipynb` for one output each. The first cell of each notebook lists the data it needs.

**Command line:**

```bash
python Fig4_valve_closure.py      # one figure
python make_all.py                # all tables and figures
```

**Outputs:**

- Tables are written to `tables/`: one CSV per table and `SI_Tables_S1-S8.docx`.
- Figures are displayed only. To export one, uncomment the `save_tiff(...)` line near the end of its script or notebook. The export is a 600 dpi RGB TIFF in `figures/`. `make_all.py` exports all figures when `SAVE_FIGURES = True`.

## Outputs and data

| Output | Script / notebook | Data |
|---|---|---|
| Fig. 1 fabrication workflow | `Fig1_fabrication.py` | none |
| Fig. 2 valve regimes and profiles | `Fig2_regimes_profiles.py` | profilometry (raw) |
| Fig. 3 dimensional fidelity | `Fig3_profilometry.py` | profilometry (raw or tables) |
| Fig. 4 valve closure | `Fig4_valve_closure.py` | valve logs + profilometry (raw or tables) |
| Fig. 5 pump and recirculating mixer | `Fig5_pump_mixer.py` | pump gravimetry (raw or tables); mixing videos for panel d |
| Fig. 6 DNA nanodevice on chip | `Fig6_origami.py` | FRET workbook (included); nanodevice renders for panel a (not included, see below) |
| Figs S1, S2 control and flow scans | `FigS1_profiles_control.py`, `FigS2_profiles_flow.py` | profilometry (raw) |
| Fig. S3 design vs measured | `FigS3_design_vs_measured.py` | profilometry (raw) |
| Fig. S4 valve traces | `FigS4_valve_traces.py` | valve logs (raw) |
| Fig. S5 closure vs pressure | `FigS5_closure_vs_pressure.py` | valve logs (raw or tables) |
| Fig. S6 aspect ratio | `FigS6_aspect_ratio.py` | valve logs + profilometry (raw or tables) |
| Tables S1–S4 | `TableS1_…` – `TableS4_…` | profilometry (raw or tables) |
| Table S5 | `TableS5_design_widths.py` | none (CAD dimensions) |
| Tables S6, S7 | `TableS6_…`, `TableS7_…` | valve logs (+ profilometry for S7) |
| Table S8 | `TableS8_pump_gravimetry.py` | pump gravimetry (raw or tables) |

If a raw-data folder is present but incomplete, the code stops and lists the missing files. It never analyses a partial data set without saying so.

## Repository layout

```
config.py            paths and figure switches (the only file you normally edit)
make_all.py          builds every table and figure
Fig*.py, Table*.py   one script per output
notebooks/           the same as Jupyter notebooks, plus 00_Check_data and 00_All_figures_and_tables
lib/                 analysis and plotting modules
  profilometry_analysis.py   scan processing, edge detection, width/height statistics
  valve_closure_analysis.py  closure per cycle and per run, minimum sealing pressure, aspect-ratio separation
  pump_analysis.py           stroke volume and flow rate from gravimetry
  mixing_analysis.py         mixing index from the videos
  loaders.py, datacheck.py   data loading (raw data or result tables) and checks
  fig*_panels.py, ch2_style.py, figio.py   figures, shared style, TIFF export
  si_tables.py               Supplementary Tables (CSV + Word)
assets/              images placed in Figs 5 and 6
data/                FRET workbook and replacement valve run; put the raw data here
results_reference/   result tables of the raw-data analysis
```

`results/`, `figures/` and `tables/` are created when the code runs.

## Settings

All settings are in `config.py`:

- `FIG6_RENDERS`: the two nanodevice renders of Fig. 6a (inactive, active). Missing files give placeholder boxes of the same size.
- `FIG6_PHOTO`: the chip photograph placed in Fig. 6b. `None` draws a placeholder box.
- `FIG6C_UNITS`: the y-axis of Fig. 6c. The manuscript uses `'log2'`, the log₂ fold change relative to origami alone. The other options are `'log10'`, `'log_axis'`, `'linear'` and `'dB'`, all described in `lib/fig6_panels.draw_c`.
- `FIG5_VIDEO_FRAMES`: whether the frames of Fig. 5d are taken from the mixing videos.

## Notes

- **Fig. 6a renders:** the 3D renders of the DNA origami nanodevice are adapted from Zadegan & Hughes (2017), © American Chemical Society, and are therefore not distributed here. The code draws placeholder boxes of the same size in their place, so the layout and every other panel are unchanged. To reproduce the published panel, put the two renders in `assets/` as `fig6a_nanodevice_inactive.png` and `fig6a_nanodevice_active.png`. The code crops them, removes the targeting ligand, and draws the labels and the 10 nm scale bar.

- **Fonts:** figures use Arial when it is installed. Otherwise a metric-compatible sans-serif is used, and a `findfont: cursive` message is harmless.
- **Mixing analysis:** `lib/mixing_analysis.py` takes 1–2 min per video.

## Citation

If you use this code, please cite the paper above. Citation metadata are in [CITATION.cff](CITATION.cff).

## License

- **Code:** MIT License (see [LICENSE](LICENSE)).
- **Data and images** in `data/`, `results_reference/` and `assets/`: CC BY 4.0. The Fig. 6a renders are not covered and are not included.

## Contact

The corresponding authors are Samuel Oliveira and Reza Zadegan, Department of Nanoengineering, North Carolina A&T State University. Please report issues through the GitHub issue tracker.
