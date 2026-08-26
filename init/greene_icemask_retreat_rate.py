"""Ice-front retreat rate of the Greene observed ice mask (icemask_greene)
from the ISMIP7 obs MIPkit -- this is the mask used to build the historical
calving-front spclevelset (see init/hist_run_tune_CESM_WACCM.m, step
'Greene_spclevelset...').

Source: raw_data/ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc
  icemask_greene(greene_mask_time, y, x), byte, flag_values 0/1
  (0=no_ice, 1=ice), 500 m polar-stereographic grid (x/y dims of length
  12161). greene_mask_time carries CF time attributes ("days since
  1900-1-1", standard_name="time"), so xarray auto-decodes it to
  datetime64 on open -- do NOT re-interpret ds['greene_mask_time'].values
  as raw day-counts. The 24 snapshots are NOT evenly spaced (gaps range
  from ~0.5 yr to several years). NOTE: despite the variable's "reliable
  ~mid-1985 to late-2021" documentation, the actual snapshots in this file
  only span 1997-10-01 to 2021-03-15 (23.46 yr) -- there is no 1985 data
  in this file.

SPARSE point design (not a dense per-interval grid): retreat only ever
touches a small fraction of the 12161x12161 domain (thousands to tens of
thousands of pixels per interval, not 147 million), so every quantity
below is built from the actual list of (interval, x, y, rate) retreat
points, never from a dense (n_intervals, ny, nx) array. An earlier version
of this script *did* build that dense stack (~13.6 GB at float32) and got
OOM-killed on a 64 GB job -- nanmedian/nanpercentile on an array that size
need a full internal working copy to sort, roughly doubling the footprint
on top of everything else already alive. The sparse design below is both
cheaper and gives you the actual raw distribution to work with directly,
not a pre-collapsed summary.

FRONT FILTER: a retreated pixel can end up in the point population for two
very different reasons -- (a) genuine front migration (a pixel that was
already on the ice/no-ice boundary retreated), or (b) an interior collapse
event (a chunk detaches all at once; deep-interior points get attributed a
"distance to remaining ice" of tens of km, which is the SIZE of the
collapsed area, not a migration speed -- dividing by a short interval
produces wildly inflated apparent rates, e.g. tens of thousands of m/yr).
front_mask() flags pixels that were already on the boundary BEFORE they
retreated (tested at the start of the interval); every output below comes
in both an unfiltered ("all points") and a front-filtered ("is_front"
only) version. Use the front-filtered numbers for anything meant to
represent actual front migration speed (e.g. an ISSM migration_max cap) --
the unfiltered ones are still useful for total-area/extent questions, but
their max is not a speed.

Outputs:

1. greene_icemask_retreat_points.csv -- the raw data: one row per
   individual retreat point across all 23 intervals (columns: start_date,
   end_date, x, y, rate_m_per_yr, is_front). Everything else in this
   script is a summary computed FROM this; nothing here is a proxy.

2. greene_icemask_retreat_rate.csv -- one row per interval: domain-wide
   gross retreat/advance/net area (km^2) and a perimeter-normalised linear
   rate (m/yr, a cheap single-number proxy that assumes a uniform front),
   plus the true per-interval point-count and rate percentiles
   (mean/median/p10/p25/p75/p90/max/min) computed directly from (1) for
   that interval, both for all points (point_rate_*) and front-filtered
   only (front_rate_*, n_front_points) -- the temporal-distribution view
   (does retreat accelerate, is one interval an outlier).

3. greene_icefront_retreat_1997_2021.nc -- endpoint map: for every pixel
   that retreated between the first (1997) and last (2021) snapshot, its
   distance to the current (2021) ice edge divided by the FULL 23.46 yr.
   A long-term average, diluted if a location's retreat actually happened
   in a short burst -- see (4) for the corrected version. Not front-
   filtered (it's already an average, not a max, so the interior-collapse
   distortion is much smaller here than in the pooled max below).

4. greene_icefront_retreat_intervals.nc -- per-pixel maps (local
   max / median / mean retreat rate, m/yr), both unfiltered and
   "_front"-suffixed (front-filtered), built by binning the points from
   (1) onto the native 500 m grid (scipy.stats.binned_statistic_2d, not a
   dense stack) -- each pixel is normalised by whichever interval it
   actually retreated in, fixing the dilution in (3). Antarctica-wide
   pooled mean/median/max/min, both unfiltered and front-filtered, are
   saved as dataset attributes, and the interval that produced each global
   max/min is printed so an extreme value can be checked against how
   short/long that interval was before trusting it.

Run from the init/ directory: python greene_icemask_retreat_rate.py
"""
import warnings

import numpy as np
import xarray as xr
import pandas as pd
from scipy.ndimage import distance_transform_edt
from scipy.stats import binned_statistic_2d

NC_PATH = '../raw_data/ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc'
OUT_CSV_POINTS   = 'greene_icemask_retreat_points.csv'
OUT_CSV_SUMMARY  = 'greene_icemask_retreat_rate.csv'
OUT_NC_ENDPOINT  = 'greene_icefront_retreat_1997_2021.nc'
OUT_NC_INTERVALS = 'greene_icefront_retreat_intervals.nc'


def perimeter_m(ice, dx):
    """Raster ice-front perimeter: count ice/no-ice edges (4-connectivity),
    each edge has length dx (assumes square cells, dx==dy)."""
    edges_x = np.count_nonzero(ice[:, 1:] != ice[:, :-1])
    edges_y = np.count_nonzero(ice[1:, :] != ice[:-1, :])
    return (edges_x + edges_y) * dx


def front_mask(ice):
    """Boolean mask, True where `ice` is ice-covered AND at least one
    4-connected neighbour is not ice -- i.e. pixels already sitting on the
    ice/no-ice boundary. Domain edges are treated as bordering not-ice
    (correct here: the AIS grid edges border open ocean).

    Used to filter retreat points to genuine FRONT retreat -- excludes
    both (a) deep-interior points from a sudden interior collapse, and
    (b) the newly-formed edge of that same collapse, since neither was
    adjacent to not-ice before the event. Only pixels that were already
    part of the pre-existing coastline qualify, however far or fast the
    front then moved from that starting point."""
    not_ice = ~ice
    adjacent_not_ice = np.zeros_like(ice, dtype=bool)
    adjacent_not_ice[:-1, :] |= not_ice[1:, :]
    adjacent_not_ice[1:, :]  |= not_ice[:-1, :]
    adjacent_not_ice[:, :-1] |= not_ice[:, 1:]
    adjacent_not_ice[:, 1:]  |= not_ice[:, :-1]
    return ice & adjacent_not_ice


def main():
    ds = xr.open_dataset(NC_PATH)
    mask = ds['icemask_greene']  # (greene_mask_time, y, x)

    x = ds['x'].values
    y = ds['y'].values
    dx = float(abs(x[1] - x[0]))
    dy = float(abs(y[1] - y[0]))
    assert dx == dy, 'non-square cells -- distance transform below assumes dx==dy'
    cell_km2 = (dx * dy) / 1e6

    dates = pd.to_datetime(ds['greene_mask_time'].values)
    nT = len(dates)
    print('Greene ice-mask snapshots: %d, %s -> %s' % (nT, dates[0].date(), dates[-1].date()))
    print('Grid: dx=%.0f m, dy=%.0f m, cell area=%.4f km^2\n' % (dx, dy, cell_km2))

    ice_first = None
    prev_ice = None
    prev_date = None
    prev_area = None
    rows = []
    point_chunks = []

    for i in range(nT):
        sl = mask.isel(greene_mask_time=i).values  # loads just this (y, x) slice
        ice = (sl == 1)
        area_km2 = float(ice.sum()) * cell_km2

        if i == 0:
            ice_first = ice.copy()
        if i == nT - 1:
            ice_last = ice.copy()
            last_date = dates[i]

        if prev_ice is not None:
            dt_years = (dates[i] - prev_date).days / 365.25

            retreated = prev_ice & ~ice
            retreat_km2 = float(retreated.sum()) * cell_km2
            advance_km2 = float(np.logical_and(~prev_ice, ice).sum()) * cell_km2
            net_km2 = area_km2 - prev_area
            perim_m = perimeter_m(prev_ice, dx)
            linear_retreat_m_per_yr = (retreat_km2 * 1e6) / perim_m / dt_years

            # Sparse extraction -- only the pixels that actually retreated.
            # is_front marks the subset that was already on the ice/no-ice
            # boundary at the START of the interval (see front_mask
            # docstring) -- genuine front retreat, excluding both the
            # interior and the newly-formed edge of a sudden interior
            # collapse.
            front_prev = front_mask(prev_ice)
            yi, xi = np.where(retreated)
            n_pts = int(yi.size)
            if n_pts > 0:
                dist_to_front_i = distance_transform_edt(~ice) * dx
                rate_vals = (dist_to_front_i[yi, xi] / dt_years).astype(np.float32)
                is_front = front_prev[yi, xi]
                point_chunks.append(pd.DataFrame({
                    'start_date': prev_date.date(), 'end_date': dates[i].date(),
                    'x': x[xi].astype(np.float32), 'y': y[yi].astype(np.float32),
                    'rate_m_per_yr': rate_vals,
                    'is_front': is_front,
                }))
                pt_mean = float(rate_vals.mean())
                pt_median = float(np.median(rate_vals))
                pt_p10, pt_p25, pt_p75, pt_p90 = (float(v) for v in
                    np.percentile(rate_vals, [10, 25, 75, 90]))
                pt_max = float(rate_vals.max())
                pt_min = float(rate_vals.min())

                front_vals = rate_vals[is_front]
                n_front = int(front_vals.size)
                if n_front > 0:
                    fr_mean = float(front_vals.mean())
                    fr_median = float(np.median(front_vals))
                    fr_p10, fr_p25, fr_p75, fr_p90 = (float(v) for v in
                        np.percentile(front_vals, [10, 25, 75, 90]))
                    fr_max = float(front_vals.max())
                    fr_min = float(front_vals.min())
                else:
                    fr_mean = fr_median = fr_p10 = fr_p25 = fr_p75 = fr_p90 = fr_max = fr_min = np.nan
            else:
                pt_mean = pt_median = pt_p10 = pt_p25 = pt_p75 = pt_p90 = pt_max = pt_min = np.nan
                n_front = 0
                fr_mean = fr_median = fr_p10 = fr_p25 = fr_p75 = fr_p90 = fr_max = fr_min = np.nan

            rows.append({
                'start_date': prev_date.date(), 'end_date': dates[i].date(),
                'dt_years': dt_years,
                'retreat_km2': retreat_km2, 'advance_km2': advance_km2,
                'net_km2': net_km2,
                'retreat_rate_km2_per_yr': retreat_km2 / dt_years,
                'advance_rate_km2_per_yr': advance_km2 / dt_years,
                'net_rate_km2_per_yr': net_km2 / dt_years,
                'front_perimeter_km': perim_m / 1e3,
                'linear_retreat_rate_m_per_yr': linear_retreat_m_per_yr,
                'n_retreat_points': n_pts,
                'point_rate_mean_m_per_yr': pt_mean,
                'point_rate_median_m_per_yr': pt_median,
                'point_rate_p10_m_per_yr': pt_p10,
                'point_rate_p25_m_per_yr': pt_p25,
                'point_rate_p75_m_per_yr': pt_p75,
                'point_rate_p90_m_per_yr': pt_p90,
                'point_rate_max_m_per_yr': pt_max,
                'point_rate_min_m_per_yr': pt_min,
                # front-filtered (is_front only) -- excludes interior-collapse artifacts
                'n_front_points': n_front,
                'front_rate_mean_m_per_yr': fr_mean,
                'front_rate_median_m_per_yr': fr_median,
                'front_rate_p10_m_per_yr': fr_p10,
                'front_rate_p25_m_per_yr': fr_p25,
                'front_rate_p75_m_per_yr': fr_p75,
                'front_rate_p90_m_per_yr': fr_p90,
                'front_rate_max_m_per_yr': fr_max,
                'front_rate_min_m_per_yr': fr_min,
            })
            print('%s -> %s  dt=%5.2fy  retreat=%9.1f km2 (%6.1f km2/yr, %6.1f m/yr linear)  '
                  'n_pts=%7d (front=%6d)  all: mean=%6.1f max=%9.1f  front: mean=%6.1f max=%7.1f m/yr' % (
                      prev_date.date(), dates[i].date(), dt_years,
                      retreat_km2, retreat_km2 / dt_years, linear_retreat_m_per_yr,
                      n_pts, n_front, pt_mean, pt_max, fr_mean, fr_max))

        prev_ice, prev_date, prev_area = ice, dates[i], area_km2

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV_SUMMARY, index=False)
    print('\nSaved: %s (%d intervals)\n' % (OUT_CSV_SUMMARY, len(df)))

    points_df = pd.concat(point_chunks, ignore_index=True)
    points_df.to_csv(OUT_CSV_POINTS, index=False)
    print('Saved: %s (%d individual retreat points, all intervals)\n' %
          (OUT_CSV_POINTS, len(points_df)))

    print('--- Per-interval domain-wide linear retreat rate (m/yr, perimeter-normalised) ---')
    print('Mean  : %.1f' % df['linear_retreat_rate_m_per_yr'].mean())
    print('Median: %.1f' % df['linear_retreat_rate_m_per_yr'].median())
    imax = df['linear_retreat_rate_m_per_yr'].idxmax()
    imin = df['linear_retreat_rate_m_per_yr'].idxmin()
    print('Max   : %.1f  (%s -> %s)' %
          (df.loc[imax, 'linear_retreat_rate_m_per_yr'], df.loc[imax, 'start_date'], df.loc[imax, 'end_date']))
    print('Min   : %.1f  (%s -> %s)\n' %
          (df.loc[imin, 'linear_retreat_rate_m_per_yr'], df.loc[imin, 'start_date'], df.loc[imin, 'end_date']))

    # ---------------------------------------------------------------------
    # 3. Endpoint map: 1997 -> 2021, long-term average (diluted by any
    #    short bursts). Only two masks involved -- always memory-cheap.
    # ---------------------------------------------------------------------
    total_years = (last_date - dates[0]).days / 365.25
    retreated_total = ice_first & ~ice_last
    dist_to_2021_front = distance_transform_edt(~ice_last).astype(np.float32) * dx
    distance_retreated_m = np.where(retreated_total, dist_to_2021_front, np.nan).astype(np.float32)
    endpoint_rate_m_per_yr = distance_retreated_m / total_years

    print('--- Endpoint (1997->2021) map, long-term average -- diluted, see (4) for true max/min ---')
    print('Retreated area : %.1f km2' % (retreated_total.sum() * cell_km2))
    print('Mean rate      : %.2f m/yr\n' % np.nanmean(endpoint_rate_m_per_yr))

    ds_endpoint = xr.Dataset(
        {
            'retreat_rate_m_per_yr': (('y', 'x'), endpoint_rate_m_per_yr),
            'distance_retreated_m': (('y', 'x'), distance_retreated_m),
            'icemask_1997': (('y', 'x'), ice_first.astype('int8')),
            'icemask_2021': (('y', 'x'), ice_last.astype('int8')),
        },
        coords={'x': x, 'y': y},
        attrs={
            'title': 'Greene ice-front retreat rate, first-to-last available snapshot (long-term average)',
            'start_date': str(dates[0].date()), 'end_date': str(last_date.date()),
            'elapsed_years': total_years,
            'method': ('Distance from each pixel that retreated between the first and last '
                       'snapshot to the nearest still-ice (last-snapshot) pixel, divided by the '
                       'FULL elapsed period. This dilutes any short retreat bursts -- use '
                       'greene_icefront_retreat_intervals.nc for true local max rates.'),
            'source': 'raw_data/ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc, icemask_greene',
        },
    )
    ds_endpoint.to_netcdf(OUT_NC_ENDPOINT)
    print('Saved: %s\n' % OUT_NC_ENDPOINT)

    # ---------------------------------------------------------------------
    # 4. Multi-interval per-pixel maps, built from the SPARSE points_df via
    #    binned_statistic_2d -- NOT from a dense (n_intervals, ny, nx)
    #    stack (that's what OOM-killed the earlier version of this script).
    #    Pooled Antarctica-wide mean/median/max/min come straight off the
    #    raw points, not off any grid.
    # ---------------------------------------------------------------------
    rate_all = points_df['rate_m_per_yr'].values
    print('--- Pooled retreat rate, all %d intervals, all %d points (m/yr) ---' % (nT - 1, rate_all.size))
    print('Mean  : %.2f' % rate_all.mean())
    print('Median: %.2f' % np.median(rate_all))
    print('Max   : %.2f' % rate_all.max())
    print('Min   : %.2f' % rate_all.min())

    i_max = int(points_df['rate_m_per_yr'].idxmax())
    i_min = int(points_df['rate_m_per_yr'].idxmin())
    print('Global max occurred in interval %s -> %s -- check this isn\'t an anomalously '
          'short/noisy interval before trusting it' %
          (points_df.loc[i_max, 'start_date'], points_df.loc[i_max, 'end_date']))
    print('Global min occurred in interval %s -> %s\n' %
          (points_df.loc[i_min, 'start_date'], points_df.loc[i_min, 'end_date']))

    # ---------------------------------------------------------------------
    # Front-filtered pooled stats: same population, restricted to is_front
    # points only (see front_mask docstring) -- excludes interior-collapse
    # artifacts, so this is the number to use for something like an ISSM
    # migration_max cap, NOT the unfiltered pooled stats above.
    # ---------------------------------------------------------------------
    front_df = points_df[points_df['is_front']]
    rate_front = front_df['rate_m_per_yr'].values
    print('--- Pooled retreat rate, FRONT-FILTERED (on the pre-existing ice edge only), '
          '%d points (of %d total, %.1f%%) (m/yr) ---' %
          (rate_front.size, rate_all.size, 100.0 * rate_front.size / rate_all.size))
    print('Mean  : %.2f' % rate_front.mean())
    print('Median: %.2f' % np.median(rate_front))
    print('Max   : %.2f' % rate_front.max())
    print('Min   : %.2f' % rate_front.min())

    fi_max = int(front_df['rate_m_per_yr'].idxmax())
    fi_min = int(front_df['rate_m_per_yr'].idxmin())
    print('Front-filtered global max occurred in interval %s -> %s' %
          (points_df.loc[fi_max, 'start_date'], points_df.loc[fi_max, 'end_date']))
    print('Front-filtered global min occurred in interval %s -> %s\n' %
          (points_df.loc[fi_min, 'start_date'], points_df.loc[fi_min, 'end_date']))

    # bin edges from cell centres (uniform 500 m grid, asserted dx==dy above)
    x_edges = np.append(x - dx / 2, x[-1] + dx / 2)
    y_edges = np.append(y - dy / 2, y[-1] + dy / 2)

    # NOTE argument order: binned_statistic_2d(sample_x, sample_y, values,
    # bins=[bins_x, bins_y]) returns statistic.shape == (len(bins_x)-1,
    # len(bins_y)-1). We want a (ny, nx) result, so the physical y-values
    # go in as the function's "x" argument together with y_edges first.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)  # empty bins -> NaN, expected
        local_max, _, _, _ = binned_statistic_2d(
            points_df['y'], points_df['x'], points_df['rate_m_per_yr'],
            statistic='max', bins=[y_edges, x_edges])
        local_median, _, _, _ = binned_statistic_2d(
            points_df['y'], points_df['x'], points_df['rate_m_per_yr'],
            statistic='median', bins=[y_edges, x_edges])
        local_mean, _, _, _ = binned_statistic_2d(
            points_df['y'], points_df['x'], points_df['rate_m_per_yr'],
            statistic='mean', bins=[y_edges, x_edges])

        local_max_front, _, _, _ = binned_statistic_2d(
            front_df['y'], front_df['x'], front_df['rate_m_per_yr'],
            statistic='max', bins=[y_edges, x_edges])
        local_median_front, _, _, _ = binned_statistic_2d(
            front_df['y'], front_df['x'], front_df['rate_m_per_yr'],
            statistic='median', bins=[y_edges, x_edges])
        local_mean_front, _, _, _ = binned_statistic_2d(
            front_df['y'], front_df['x'], front_df['rate_m_per_yr'],
            statistic='mean', bins=[y_edges, x_edges])

    ds_intervals = xr.Dataset(
        {
            'local_max_rate_m_per_yr': (('y', 'x'), local_max.astype('float32')),
            'local_median_rate_m_per_yr': (('y', 'x'), local_median.astype('float32')),
            'local_mean_rate_m_per_yr': (('y', 'x'), local_mean.astype('float32')),
            'local_max_rate_front_m_per_yr': (('y', 'x'), local_max_front.astype('float32')),
            'local_median_rate_front_m_per_yr': (('y', 'x'), local_median_front.astype('float32')),
            'local_mean_rate_front_m_per_yr': (('y', 'x'), local_mean_front.astype('float32')),
        },
        coords={'x': x, 'y': y},
        attrs={
            'title': 'Greene ice-front retreat rate, per-pixel across the 23 consecutive '
                     'snapshot intervals (each pixel normalised by the interval it actually '
                     'retreated in, not the full record -- fixes the dilution in the '
                     'endpoint (1997->2021) map). Built from greene_icemask_retreat_points.csv '
                     'via binned_statistic_2d on the native 500 m grid. "_front" variables are '
                     'restricted to points already on the ice/no-ice boundary before they '
                     'retreated (front_mask) -- use these, not the unfiltered ones, for anything '
                     'meant to represent front migration speed (e.g. an ISSM migration_max cap) '
                     '-- the unfiltered ones include interior-collapse artifacts that inflate '
                     'the max by orders of magnitude (see pooled_max_m_per_yr vs '
                     'pooled_max_front_m_per_yr).',
            'n_intervals': nT - 1,
            'pooled_mean_m_per_yr': float(rate_all.mean()),
            'pooled_median_m_per_yr': float(np.median(rate_all)),
            'pooled_max_m_per_yr': float(rate_all.max()),
            'pooled_min_m_per_yr': float(rate_all.min()),
            'pooled_mean_front_m_per_yr': float(rate_front.mean()),
            'pooled_median_front_m_per_yr': float(np.median(rate_front)),
            'pooled_max_front_m_per_yr': float(rate_front.max()),
            'pooled_min_front_m_per_yr': float(rate_front.min()),
            'source': 'raw_data/ISMIP7/AIS/obs/mipkit/AntarcticaObsISMIP7-v1.2.nc, icemask_greene',
        },
    )
    ds_intervals.to_netcdf(OUT_NC_INTERVALS)
    print('Saved: %s' % OUT_NC_INTERVALS)


if __name__ == '__main__':
    main()
