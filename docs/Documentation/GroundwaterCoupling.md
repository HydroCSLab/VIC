# Groundwater coupling interface

VIC keeps its standard ARNO baseflow formulation unless the environment
variable `VIC_GW_FORMULATION` is set. The supported values are `head` and
`scheidegger`; any other nonempty value is rejected when runoff is computed. This makes
the coupling behavior opt in and preserves existing VIC runs when the
coupling variables are absent.

The exchange returned by the coupling routines is in millimeters for the
current runoff substep. Positive exchange leaves the soil column and enters
the groundwater model. Negative exchange enters the soil column from
groundwater. The value is applied to the bottom soil layer and is accumulated
in `OUT_GW_EXCHANGE`, whose aggregation type is a sum. The signed exchange is
also included in the VIC water balance.

Coupling requires `LAKES=FALSE`. VIC rejects lake-enabled coupling before the
first model step because exchange accounting during lake/wetland area changes
is not supported. Invalid or nonfinite exchange inputs and calculations abort
the run; a valid zero conductivity or equilibrium gradient still gives zero flow.

In `scheidegger` mode, VIC places the groundwater table `VIC_GW_TEST_GAP_M`
below the soil base and computes exchange from the bottom-layer pressure head.
In `head` mode, an external head requires `VIC_GW_EXCHANGE_LENGTH_M` and
exactly one of `VIC_GW_HEAD_FILE` or `VIC_GW_HEAD_OFFSET_M`. Both spatial and
scalar offsets are relative to the fixed soil base. `VIC_GW_REFERENCE_DEPTH`
selects the soil hydraulic reference (`base`, the default, or bottom-layer
`node`); it does not change the external head's datum.

For historical head/depth equivalence tests, if no spatial head is loaded and
both scalar head and exchange length are absent, head mode uses the same
groundwater depth as `scheidegger` mode. The exchange length is then the
distance from the selected soil reference to that groundwater table. Supplying
only the scalar head or only the length is an error. The VIC-MF6 framework
always supplies a spatial head file and an explicit exchange length.

`VIC_GW_HEAD_FILE` is a headerless text file containing one record per active
cell:

```text
latitude longitude head_offset_m
```

Blank lines, lines beginning with `#`, and trailing inline comments are
allowed. The image driver matches records by the model cell coordinates,
ignores records outside the local MPI domain, rejects duplicates and invalid
numbers, and requires every local active cell to have a matching record. A
spatial file and scalar head offset must not be supplied together. Spatial
head loading is supported by the image driver; the classic driver can use
the scalar head or the historical fallback.

The following optional controls are validated when coupling is enabled:

| Variable | Meaning | Default |
| --- | --- | --- |
| `VIC_GW_TEST_GAP_M` | positive groundwater-table depth below the soil base for depth mode and the head-mode fallback | `1` |
| `VIC_GW_TEST_KA_SCALE` | scale the aquifer conductivity used by the coupling test setup | `1` |
| `VIC_GW_MAX_DRAIN_FRACTION` | maximum fraction of liquid soil water that can leave in one substep | `0.01` |
| `VIC_GW_REFERENCE_DEPTH` | soil hydraulic reference (`base` or `node`); external head datum stays fixed | `base` |

The limiter uses the soil column's residual moisture directly in millimeters
and permits only liquid water to leave. Ice remains part of the available
storage calculation, but cannot be exported as liquid groundwater exchange.

The existing calendar guard still requires a complete final day by default.
For a deliberately bounded subdaily coupling run, set
`VIC_ALLOW_PARTIAL_DAY=1`.

## Verification

Build the image driver from the repository root:

```bash
make -C vic/drivers/image -j2
```

Run the focused source tests with the Python driver extension on the import
path:

```bash
PYTHONPATH="$PWD/vic/drivers/python" \
  python3 -m pytest -q \
  tests/unit/vic_run/test_runoff.py \
  tests/unit/shared/test_read_groundwater_head_file.py \
  tests/unit/shared/test_make_dmy.py
```

The coupled Stehekin acceptance workflow is maintained in the VIC-MF6
repository. It supplies a spatial head file, checks both flow directions,
and verifies the VIC and MODFLOW water budgets.
