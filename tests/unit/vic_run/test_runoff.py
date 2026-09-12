import pytest

from vic import ffi
from vic import lib as vic_lib


BOTTOM_LAYER = 2


def make_soil():
    """Return a bottom soil layer with simple, traceable hydraulic values."""
    soil = ffi.new("soil_con_struct *")
    soil.depth[BOTTOM_LAYER] = 1.0
    soil.max_moist[BOTTOM_LAYER] = 400.0
    soil.resid_moist[BOTTOM_LAYER] = 40.0
    soil.bubble[BOTTOM_LAYER] = 20.0
    soil.expt[BOTTOM_LAYER] = 5.0
    return soil


def head_exchange(soil, liquid, groundwater_head, *, ice=0.0,
                  exchange_length=2.0, conductivity=8.0,
                  max_drain_fraction=1.0, runoff_steps_per_day=4):
    return vic_lib.calc_groundwater_exchange_head(
        soil, BOTTOM_LAYER, liquid, ice, 1.0, groundwater_head,
        exchange_length, conductivity, max_drain_fraction,
        runoff_steps_per_day)


def test_head_exchange_sign_and_equilibrium():
    soil = make_soil()

    # At 50% liquid saturation, pressure head is -0.4 m and soil head
    # at the 1 m reference depth is -1.4 m.
    assert head_exchange(soil, 200.0, -1.4) == pytest.approx(0.0)
    assert head_exchange(soil, 200.0, -2.4) == pytest.approx(1.0)
    assert head_exchange(soil, 200.0, -0.4) == pytest.approx(-1.0)


def test_head_exchange_scales_with_runoff_substeps():
    soil = make_soil()

    four_steps = head_exchange(
        soil, 200.0, -2.4, runoff_steps_per_day=4)
    eight_steps = head_exchange(
        soil, 200.0, -2.4, runoff_steps_per_day=8)

    assert eight_steps == pytest.approx(four_steps / 2.0)


def test_downward_exchange_respects_residual_storage():
    soil = make_soil()
    soil.resid_moist[BOTTOM_LAYER] = 90.0

    exchange = head_exchange(
        soil, 100.0, -101.8, exchange_length=1.0, conductivity=100.0,
        max_drain_fraction=0.5, runoff_steps_per_day=1)

    # resid_moist is already expressed in mm, leaving 10 mm drainable.
    assert exchange == pytest.approx(10.0)


def test_downward_exchange_only_removes_liquid_water():
    soil = make_soil()
    soil.resid_moist[BOTTOM_LAYER] = 120.0

    exchange = head_exchange(
        soil, 100.0, -101.8, ice=50.0, exchange_length=1.0,
        conductivity=100.0, max_drain_fraction=1.0,
        runoff_steps_per_day=1)

    assert exchange == pytest.approx(30.0)


def test_head_and_depth_formulations_are_equivalent():
    soil = make_soil()

    depth_exchange = vic_lib.calc_groundwater_exchange(
        soil, BOTTOM_LAYER, 200.0, 0.0, 2.0, 1.0, 8.0, 1.0, 4)
    head_exchange_value = head_exchange(
        soil, 200.0, -2.0, exchange_length=1.0)

    assert depth_exchange == pytest.approx(head_exchange_value)


def test_head_exchange_rejects_invalid_exchange_length():
    soil = make_soil()

    assert head_exchange(
        soil, 200.0, -2.4, exchange_length=0.0) == pytest.approx(0.0)
