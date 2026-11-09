"""Checks for the pure functions in tas_fires.py (no network needed)."""
import datetime as dt
import math

import numpy as np
import pytest

import tas_fires as tf


def test_reflectance_scaling_and_nodata():
    r = tf.to_reflectance(np.array([[-999, 0, 2500, 10000]], dtype="int16"))
    assert np.isnan(r[0, 0])
    assert r[0, 1] == 0.0
    assert r[0, 2] == pytest.approx(0.25)
    assert r[0, 3] == pytest.approx(1.0)


def test_nbr_and_ndvi():
    c = {"nir": np.array([0.4, 0.2, np.nan]), "swir2": np.array([0.1, 0.3, 0.1]),
         "red": np.array([0.05, 0.1, 0.1])}
    n = tf.nbr(c)
    assert n[0] == pytest.approx(0.6)
    assert n[1] == pytest.approx(-0.2)
    assert np.isnan(n[2])
    assert tf.ndvi(c)[0] == pytest.approx(0.35 / 0.45)


def test_norm_diff_zero_denominator_is_nan():
    assert np.isnan(tf.norm_diff(np.array([0.0]), np.array([0.0]))[0])


def test_severity_classes():
    d = np.array([-0.3, 0.05, 0.10, 0.2, 0.3, 0.5, 0.7, np.nan])
    assert tf.classify_dnbr(d).tolist() == [0, 0, 1, 1, 2, 3, 4, -1]


def test_windows():
    ign, end = dt.date(2019, 1, 15), dt.date(2019, 3, 24)
    assert tf.pre_window(ign) == (dt.date(2018, 9, 17), dt.date(2019, 1, 14))
    assert tf.immediate_window(end)[0] == dt.date(2019, 3, 25)
    assert tf.summer_window(2020) == (dt.date(2019, 12, 1), dt.date(2020, 3, 31))
    assert tf.first_summer_after(end) == 2020          # Dec 2019 - Mar 2020
    assert tf.first_summer_after(dt.date(2025, 6, 24)) == 2026
    assert tf.first_summer_after(dt.date(2019, 12, 15)) == 2021   # out in Dec: skip that summer
    assert tf.latest_summer(dt.date(2026, 11, 8)) == 2026
    assert tf.latest_summer(dt.date(2026, 2, 1)) == 2025
    assert tf.pre_summer_year(dt.date(2019, 1, 15)) == 2018     # Dec 2017 - Mar 2018
    assert tf.pre_summer_year(dt.date(2019, 12, 30)) == 2019    # Dec 2018 - Mar 2019
    assert tf.pre_summer_year(dt.date(2025, 12, 4)) == 2025


def test_grid_snaps_and_coarsens():
    g = tf.grid_for_bounds((1000.0, -5000.0, 2010.0, -4000.0), buffer_m=0, res=20)
    assert g.x0 % 40 == 0 and g.y0 % 40 == 0
    assert g.width % 2 == 0 and g.height % 2 == 0
    l, b, r, t = g.bounds
    assert l <= 1000 and r >= 2010 and b <= -5000 and t >= -4000
    c = g.coarsen(2)
    assert c.res == 40 and c.width == g.width // 2 and c.bounds == g.bounds


def test_true_north_in_tasmania():
    # Tasmania is ~14 degrees east of the EPSG:3577 central meridian (132 E): grid north and true
    # north differ by several degrees there, and not at all on the central meridian.
    from pyproj import Transformer
    to_xy = Transformer.from_crs("EPSG:4326", tf.CRS, always_xy=True)
    a = tf.true_north_angle(*to_xy.transform(146.5, -42.0))
    assert 3 < abs(a) < 10
    assert abs(tf.true_north_angle(*to_xy.transform(132.0, -25.0))) < 0.1


def test_lonlat_to_albers():
    # Hobart, roughly: EPSG:3577 x about 1.3e6 m, y about -4.8e6 m
    from shapely.geometry import Point
    p = tf.to_albers(Point(147.33, -42.88))
    assert 1.1e6 < p.x < 1.5e6
    assert -5.0e6 < p.y < -4.6e6


def test_group_by_date_least_cloudy_first():
    items = [{"properties": {"datetime": "2019-01-01T00:10:00Z", "eo:cloud_cover": 50}},
             {"properties": {"datetime": "2019-01-01T00:10:05Z", "eo:cloud_cover": 30}},
             {"properties": {"datetime": "2019-01-11T00:10:00Z", "eo:cloud_cover": 5}}]
    g = tf.group_by_date(items)
    assert [x[0] for x in g] == ["2019-01-11", "2019-01-01"]
    assert len(g[1][1]) == 2 and g[1][2] == pytest.approx(40)
