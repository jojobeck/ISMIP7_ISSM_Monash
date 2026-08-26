"""Plot the outputs of greene_icemask_retreat_rate.py:

  1. figures/figure_greene_retreat_map.png
     The retreat_rate_m_per_yr map from greene_icefront_retreat_1997_2021.nc
     -- Euclidean distance from each pixel that retreated between 1997 and
     2021 to the current (2021) ice edge, divided by the FULL elapsed
     period (a long-term average, diluted by any short bursts -- see (2)).

  2. figures/figure_greene_retreat_interval_maps.png
     Three panels from greene_icefront_retreat_intervals.nc: the per-pixel
     local max / median / mean retreat rate across the 23 consecutive
     Greene snapshot intervals, each pixel normalised by the actual
     interval it retreated in (not the full 1997-2021 span) -- this is
     where genuine short-burst hotspots show up that the endpoint map (1)
     smooths away. Each panel has its own colour scale (max is typically
     much larger than median/mean; a shared scale would wash out the
     other two).

  3. figures/figure_greene_retreat_distribution.png
     A proper per-interval boxplot built from the true per-pixel
     retreat-point distribution in greene_icemask_retreat_rate.csv
     (point_rate_p10/p25/median/p75/p90 columns -- real data, not a
     proxy), one box per interval, box width = interval duration and
     position = interval midpoint date (the 23 intervals are not evenly
     spaced). Overlaid: the domain-wide perimeter-proxy rate (diamond
     marker) for comparison against the true per-pixel median, and a
     reference line at the pooled all-intervals median (from the
     interval-maps NetCDF attrs).

Run from the init/ directory: python plot_greene_icemask_retreat.py
"""
import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import Normalize, ListedColormap

try:
    import cmocean
    SEQ_CMAP = cmocean.cm.amp
except ImportError:
    SEQ_CMAP = plt.get_cmap('Reds')
SEQ_CMAP = SEQ_CMAP.copy()
SEQ_CMAP.set_bad(alpha=0)

NC_ENDPOINT  = 'greene_icefront_retreat_1997_2021.nc'
NC_INTERVALS = 'greene_icefront_retreat_intervals.nc'
CSV_DIST     = 'greene_icemask_retreat_rate.csv'
figure_dir = 'figures'
os.makedirs(figure_dir, exist_ok=True)

ICE_1997_COLOR = '#cfe3ee'   # light blue: context, not the data itself
LAND_COLOR     = '#f7f7f5'
BOX_COLOR      = '#2a6f7f'   # per-pixel distribution (the real data)
PROXY_COLOR    = '#b23a3a'   # domain-wide perimeter proxy (for comparison)
REF_COLOR      = '#444444'


def plot_endpoint_map():
    ds = xr.open_dataset(NC_ENDPOINT)
    x, y = ds['x'].values, ds['y'].values
    extent = [float(x.min()), float(x.max()), float(y.min()), float(y.max())]

    ice1997 = ds['icemask_1997'].values
    ice2021 = ds['icemask_2021'].values
    rate = ds['retreat_rate_m_per_yr'].values

    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)

    bg = np.where(ice1997 == 1, 1, 0)
    ax.imshow(bg, origin='lower', extent=extent,
              cmap=ListedColormap([LAND_COLOR, ICE_1997_COLOR]),
              vmin=0, vmax=1, interpolation='nearest', zorder=0)
    ax.contour(x, y, ice2021, levels=[0.5], colors='#3a3a3a', linewidths=0.8, zorder=1)

    vmax = float(np.nanpercentile(rate, 99))
    norm = Normalize(vmin=0, vmax=vmax)
    im = ax.imshow(rate, origin='lower', extent=extent, cmap=SEQ_CMAP, norm=norm,
                   interpolation='nearest', zorder=2)

    cbar = fig.colorbar(im, ax=ax, shrink=0.75, extend='max')
    cbar.set_label('Retreat rate, long-term average (m yr$^{-1}$)')

    ax.set_title('Greene ice-front retreat rate: endpoint (1997-10-01 → 2021-03-15)\n'
                 'diluted long-term average -- see interval-maps figure for burst rates\n'
                 'background: 1997 ice extent; outline: 2021 ice front')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

    out = os.path.join(figure_dir, 'figure_greene_retreat_map.png')
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print('Saved:', out)


def plot_interval_maps():
    ds = xr.open_dataset(NC_INTERVALS)
    x, y = ds['x'].values, ds['y'].values
    extent = [float(x.min()), float(x.max()), float(y.min()), float(y.max())]

    panels = [
        ('local_max_rate_m_per_yr', 'Local max (across 23 intervals)'),
        ('local_median_rate_m_per_yr', 'Local median (across 23 intervals)'),
        ('local_mean_rate_m_per_yr', 'Local mean (across 23 intervals)'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 7), constrained_layout=True)
    for ax, (var, label) in zip(axes, panels):
        field = ds[var].values
        vmax = float(np.nanpercentile(field, 99))
        norm = Normalize(vmin=0, vmax=vmax)
        im = ax.imshow(field, origin='lower', extent=extent, cmap=SEQ_CMAP, norm=norm,
                       interpolation='nearest')
        cbar = fig.colorbar(im, ax=ax, shrink=0.7, extend='max', orientation='horizontal', pad=0.03)
        cbar.set_label('m yr$^{-1}$', fontsize=9)
        ax.set_title(label, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

    fig.suptitle('Greene ice-front retreat rate: per-pixel across the 23 consecutive intervals\n'
                'each pixel normalised by the actual interval it retreated in, not the full record',
                fontsize=12)

    out = os.path.join(figure_dir, 'figure_greene_retreat_interval_maps.png')
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print('Saved:', out)


def plot_distribution():
    df = pd.read_csv(CSV_DIST, parse_dates=['start_date', 'end_date'])
    ds_int = xr.open_dataset(NC_INTERVALS)
    pooled_median = ds_int.attrs['pooled_median_m_per_yr']

    mid_date = df['start_date'] + (df['end_date'] - df['start_date']) / 2
    width_days = (df['end_date'] - df['start_date']).dt.days.values
    positions = mdates.date2num(mid_date)

    stats = []
    for _, row in df.iterrows():
        stats.append({
            'label': '',
            'med': row['point_rate_median_m_per_yr'],
            'q1': row['point_rate_p25_m_per_yr'],
            'q3': row['point_rate_p75_m_per_yr'],
            'whislo': row['point_rate_p10_m_per_yr'],
            'whishi': row['point_rate_p90_m_per_yr'],
            'mean': row['point_rate_mean_m_per_yr'],
            'fliers': [],
        })

    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)

    ax.axhline(pooled_median, color=REF_COLOR, lw=1, ls='--', zorder=1)
    ax.text(positions[-1], pooled_median, ' pooled median (all intervals) %.1f m/yr' % pooled_median,
           color=REF_COLOR, va='bottom', ha='right', fontsize=9)

    bp = ax.bxp(stats, positions=positions, widths=width_days * 0.8,
               patch_artist=True, showmeans=True, manage_ticks=False, zorder=2)
    for box in bp['boxes']:
        box.set_facecolor(BOX_COLOR)
        box.set_alpha(0.65)
        box.set_edgecolor(BOX_COLOR)
    for element in ('whiskers', 'caps'):
        for line in bp[element]:
            line.set_color(BOX_COLOR)
    for line in bp['medians']:
        line.set_color('white')
        line.set_linewidth(1.5)

    ax.scatter(positions, df['linear_retreat_rate_m_per_yr'], marker='D', s=22,
              color=PROXY_COLOR, zorder=3, label='domain-wide perimeter-proxy rate')

    # selective direct labels: point count for the two largest-n intervals only
    for _, row in df.nlargest(2, 'n_retreat_points').iterrows():
        pos = mdates.date2num(row['start_date'] + (row['end_date'] - row['start_date']) / 2)
        ax.annotate('n=%d' % row['n_retreat_points'], xy=(pos, row['point_rate_p90_m_per_yr']),
                   xytext=(0, 6), textcoords='offset points', ha='center', fontsize=8, color=REF_COLOR)

    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.set_ylabel('Retreat rate (m yr$^{-1}$)')
    ax.set_title('Per-interval retreat-rate distribution (23 Greene snapshot intervals, 1997–2021)\n'
                'box = true per-pixel retreat-point distribution (p10–p90, median); '
                'width = actual interval duration')
    ax.legend(loc='upper left', frameon=False, fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', color='#dddddd', lw=0.6, zorder=0)
    ax.set_axisbelow(True)

    out = os.path.join(figure_dir, 'figure_greene_retreat_distribution.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print('Saved:', out)


if __name__ == '__main__':
    plot_endpoint_map()
    plot_interval_maps()
    plot_distribution()
