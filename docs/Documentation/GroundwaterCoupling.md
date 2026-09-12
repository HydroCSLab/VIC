# Groundwater coupling interface

VIC keeps its standard ARNO baseflow formulation unless the environment
variable `VIC_GW_FORMULATION` is set. The supported values are `head` and
`scheidegger`; any other value is rejected during initialization. This makes
the coupling behavior opt in and preserves existing VIC runs when the
coupling variables are absent.

The exchange returned by the coupling routines is in millimeters for the
current runoff substep. Positive exchange leaves the soil column and enters
the groundwater model. Negative exchange enters the soil column from
groundwater. The value is applied to the bottom soil layer and is accumulated
in `OUT_GW_EXCHANGE`, whose aggregation type is a sum. The signed exchange is
also included in the VIC water balance.

In `scheidegger` mode, VIC computes exchange from the bottom-layer pressure
head and the supplied aquifer reference depth. In `head` mode, VIC computes
the same Darcy-style flux from an external groundwater head. The head mode
requires `VIC_GW_EXCHANGE_LENGTH_M` and either `VIC_GW_HEAD_FILE` or the
scalar `VIC_GW_HEAD_OFFSET_M`. A scalar offset is interpreted relative to the
soil-base interface. The optional `VIC_GW_REFERENCE_DEPTH` value is `base`
(the default) or `node`.

`VIC_GW_HEAD_FILE` is a headerless text file containing one record per active
cell:

```text
latitude longitude head_offset_m
```

Blank lines, lines beginning with `#`, and trailing inline comments are
allowed. The image driver matches records by the model cell coordinates,
ignores records outside the local MPI domain, rejects duplicates and invalid
numbers, and requires every local active cell to have a matching record. A
spatial file takes precedence over the scalar head offset.

The following optional controls are validated as finite numeric values:

| Variable | Meaning | Default |
| --- | --- | --- |
| `VIC_GW_TEST_KA_SCALE` | scale the aquifer conductivity used by the coupling test setup | `1` |
| `VIC_GW_MAX_DRAIN_FRACTION` | maximum fraction of liquid soil water that can leave in one substep | `0.01` |
| `VIC_GW_REFERENCE_DEPTH` | reference for the external head (`base` or `node`) | `base` |

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
