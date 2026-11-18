"""
Unit tests for the pure (no-I/O) helper functions in polaris.helpers
and polaris.config. These don't require any FITS data or SExtractor
installation, so they're safe to run in CI.
"""

import numpy as np

from polaris.helpers import make_radius_map, annulus_median, peak_offset, get_patch, is_valid_cutout


def test_make_radius_map_center_is_zero():
    r_map = make_radius_map(5, 5)
    assert r_map[2, 2] == 0
    assert r_map.shape == (5, 5)


def test_make_radius_map_symmetric():
    r_map = make_radius_map(7, 7)
    assert np.isclose(r_map[0, 3], r_map[6, 3])
    assert np.isclose(r_map[3, 0], r_map[3, 6])


def test_annulus_median_basic():
    data = np.ones((21, 21)) * 5.0
    r_map = make_radius_map(21, 21)
    med = annulus_median(data, r_map, 2, 5)
    assert med == 5.0


def test_annulus_median_empty_annulus_returns_nan():
    data = np.ones((5, 5))
    r_map = make_radius_map(5, 5)
    med = annulus_median(data, r_map, 100, 200)
    assert np.isnan(med)


def test_peak_offset_centered_source():
    data = np.zeros((9, 9))
    data[4, 4] = 100  # exact center for a 9x9 array
    assert peak_offset(data) == 0


def test_peak_offset_off_center_source():
    data = np.zeros((9, 9))
    data[0, 4] = 100  # 4 pixels above center
    assert peak_offset(data) == 4


def test_get_patch_matches_first_containing_window():
    patches = [
        (0.0, 3.25, 16.0, 18.0, 200, "Core"),
        (8.0, 25.0, 11.0, 13.25, 1000, "Outer_1"),
    ]
    size, folder = get_patch(1.0, 17.0, patches)
    assert (size, folder) == (200, "Core")


def test_get_patch_no_match_returns_none():
    patches = [(0.0, 3.25, 16.0, 18.0, 200, "Core")]
    size, folder = get_patch(50.0, 5.0, patches)
    assert size is None
    assert folder is None


def test_is_valid_cutout_all_finite():
    data = np.ones((10, 10))
    assert is_valid_cutout(data, threshold=0.2)


def test_is_valid_cutout_too_many_nans():
    data = np.ones((10, 10))
    data[:5, :] = np.nan  # 50% NaN
    assert not is_valid_cutout(data, threshold=0.2)
