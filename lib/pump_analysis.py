"""
pump_analysis.py - gravimetric characterisation of the three-valve peristaltic pump (0.5 mm valve spacing).

DATA   pump_gravimetry.xlsx, sheets '0.5 mm spacing 0.5 bar' / '1 bar' / '1.5 bar': per ESI delay (25, 33, 48, 83, 167 ms) three replicate
       weighings (mg, before/after N cycles) and the measured elapsed time of each run ("mm:ss:cs").  Notes sheet: weights in mg.
WHAT CHANGED vs Pump.ipynb
  * Frequency (x axis) = what was PROGRAMMED:  f_phase = 1000 / delay_ms  (40.0, 30.3, 20.8, 12.0 Hz);  one cycle = 6 phases, so
    f_cycle = f_phase / 6.  It is the setting, not a measurement.
  * Time: the notebook divided by the NOMINAL run time (cycles x 6 x delay).  The recorded elapsed times are 1.2-2.0x longer (every phase
    lasts 17-25 ms longer than the programmed delay), so nominal-time flow rates were too high.  Flow rate is a measurement and uses the
    MEASURED time of the same run:   Q = V / t ,  V_stroke = V / N .  The delivered phase rate 6 N / t is also listed (f_meas_Hz).
  * SD: ddof = 1 (Pump.ipynb used ddof = 0).
  * Replicates without a recorded time contribute to the stroke volume but not to Q or f_phase (n is shown); n = 1 has no error bar.
  * 1.5 bar / 83 ms uses the spreadsheet (197, 197, 201 mg after; 151, 151, 150 before); Pump.ipynb had 198 and a missing third replicate.
  * Stroke volume V/N does not depend on the clock.  Q is V/t and f is 6N/t, so Q vs f is a straight line with slope 6 V_stroke / 60 by
    construction; its linearity is NOT independent evidence beyond the constancy of V_stroke, which is tested directly (two-way ANOVA).
"""
import re, datetime, os
import numpy as np, pandas as pd
import openpyxl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from ch2_style import *

XLSX = './pump_gravimetry.xlsx'
DENSITY = 0.9976          # g/mL, water ~20 C (0.5 % FD&C Yellow 5 solution assumed equal to water)
PHASES = 6
PRESS = {'0.5 mm spacing 0.5 bar': 0.5, '1 bar': 1.0, '1.5 bar': 1.5}
PCOL = {0.5: '#D9BC62', 1.0: '#B08D2E', 1.5: '#5E4A14'}       # gold ramp = actuation (control-layer colour), light -> dark with pressure


def parse_time(v):
    if v is None:
        return None
    if isinstance(v, datetime.timedelta):                       # Excel read '24:37:10' as h:m:s
        tot = v.total_seconds(); mm, ss, cs = int(tot // 3600), int((tot % 3600) // 60), tot % 60
        return mm * 60 + ss + cs / 100.0
    nums = re.findall(r'\d+', str(v))
    return int(nums[0]) * 60 + int(nums[1]) + int(nums[2]) / 100.0 if len(nums) == 3 else None


def load(xlsx=XLSX):
    wb = openpyxl.load_workbook(xlsx, data_only=True); rows = []
    for ws in wb.worksheets[:3]:
        bar = PRESS[ws.title]; data = [[c.value for c in r] for r in ws.iter_rows()]
        cols = {int(v): j for j, v in enumerate(data[0]) if isinstance(v, (int, float))}
        times = [r for r in data if r[0] == 'Actual time']; cyc = [r for r in data if r[0] == 'Cycles'][0]
        for d, j in cols.items():
            for k in (1, 2, 3):
                b = [r for r in data if r[0] == f'Weight before ({k})'][0][j]; a = [r for r in data if r[0] == f'Weight after ({k})'][0][j]
                if a is None or b is None:
                    continue
                t = parse_time(times[k - 1][j]); N = cyc[j]; vol_ul = (a - b) / DENSITY            # mg -> uL
                rows.append(dict(bar=bar, delay_ms=d, rep=k, cycles=N, dm_mg=a - b, vol_uL=vol_ul, t_s=t, stroke_nL=vol_ul / N * 1000,
                                 flow_uL_min=vol_ul / t * 60 if t else np.nan, f_prog_Hz=1000.0 / d, f_phase_Hz=PHASES * N / t if t else np.nan,
                                 phase_ms=t / (PHASES * N) * 1000 if t else np.nan))
    return pd.DataFrame(rows)


def summarise(df):
    g = df.groupby(['bar', 'delay_ms'])
    S = g.agg(n_vol=('stroke_nL', 'count'), stroke=('stroke_nL', 'mean'), stroke_sd=('stroke_nL', lambda v: v.std(ddof=1) if len(v) > 1 else np.nan),
              n_t=('flow_uL_min', lambda v: v.notna().sum()), flow=('flow_uL_min', 'mean'), flow_sd=('flow_uL_min', lambda v: v.std(ddof=1) if v.notna().sum() > 1 else np.nan),
              f_prog=('f_prog_Hz', 'first'), f_phase=('f_phase_Hz', 'mean'), phase_ms=('phase_ms', 'mean')).reset_index()
    S['nominal_f_Hz'] = 1000 / S.delay_ms
    return S


def anova(df):
    import statsmodels.formula.api as smf, statsmodels.api as sm
    m = smf.ols('stroke_nL ~ C(bar) + C(delay_ms)', data=df).fit(); return sm.stats.anova_lm(m, typ=2), m


def fig_pump(S):
    """Two-panel pump figure; x = programmed phase frequency (1000 / delay)."""
    fig, axs = plt.subplots(1, 2, figsize=(5.4, 2.3)); fig.subplots_adjust(left=0.10, right=0.99, top=0.95, bottom=0.2, wspace=0.34)
    for bar in (0.5, 1.0, 1.5):
        s = S[S.bar == bar].sort_values('delay_ms'); c = PCOL[bar]; ok = s.flow.notna()
        axs[0].errorbar(s.f_prog[ok], s.flow[ok], yerr=s.flow_sd[ok], color=c, marker='o', ms=3.5, lw=1.4, capsize=2, elinewidth=0.7, label=f'{bar:g} bar')
        axs[1].errorbar(s.f_prog, s.stroke, yerr=s.stroke_sd, color=c, marker='o', ms=3.5, lw=1.4, capsize=2, elinewidth=0.7)
    ticks = sorted(S.f_prog.unique())
    for ax in axs:
        ax.set_xticks(ticks); ax.set_xticklabels([f'{v:.1f}' for v in ticks]); ax.set_xlim(ticks[0] - 5, ticks[-1] + 5)
        ax.set_xlabel('Phase frequency (Hz)'); ax.grid(axis='y', ls=':', lw=0.4, color='#CCCCCC', alpha=0.6)
    axs[0].set_ylim(0, 5); axs[1].set_ylim(12, 22); axs[0].set_ylabel('Flow rate (µL min$^{-1}$)'); axs[1].set_ylabel('Volume per cycle (nL)')
    axs[0].legend(title='Actuation pressure', title_fontsize=BASE_FS - 1.5, frameon=False, fontsize=BASE_FS - 1.5, loc='upper left', handlelength=1.6)
    return fig


if __name__ == '__main__':
    df = load(); S = summarise(df)
    os.makedirs('./results_pump', exist_ok=True); df.to_csv('./results_pump/pump_per_replicate.csv', index=False); S.to_csv('./results_pump/pump_summary.csv', index=False)
    pd.set_option('display.width', 220); print(S.round(2).to_string(index=False))
    tab, m = anova(df); print(tab.round(4)); print('pooled stroke volume: %.2f +/- %.2f nL (SD of %d replicates)' % (df.stroke_nL.mean(), df.stroke_nL.std(ddof=1), len(df)))
    f = fig_pump(S)
    # f.savefig('Fig5ab_pump.jpeg', dpi=600, bbox_inches='tight')    # save figure (uncomment to export)
    plt.show()
