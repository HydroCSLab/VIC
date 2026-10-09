import subprocess
import sys

import pytest

from vic import ffi
from vic import lib as vic_lib


def test_read_groundwater_head_file_matches_local_cells(tmp_path):
    head_file = tmp_path / "groundwater_heads.txt"
    head_file.write_text(
        "# latitude longitude head_offset_m\n"
        "47.00 -119.00 99.0\n"
        "48.5000000001 -120.5000000001 -0.75 # local cell 2\n"
        "\n"
        "48.25 -120.75 0.25\n"
    )

    cells = ffi.new("soil_con_struct[2]")
    cells[0].lat = 48.50
    cells[0].lng = -120.50
    cells[1].lat = 48.25
    cells[1].lng = -120.75
    nloaded = vic_lib.read_groundwater_head_file(
        str(head_file).encode(), cells, 2)

    assert nloaded == 2
    assert all(cell.gw_head_offset_valid for cell in cells)
    assert [cell.gw_head_offset_m for cell in cells] == [-0.75, 0.25]
    assert vic_lib.read_groundwater_head_file(
        str(head_file).encode(), ffi.NULL, 0) == 0


@pytest.mark.parametrize("rows, second_latitude, message", [
    ("48 -120 1\n48 -120 2\n", 49, "Duplicate groundwater head"),
    ("48 -120 1\n", 49, "No groundwater head found"),
    ("48 -120 1\n", 48, "matched more than one"),
])
def test_head_file_requires_unique_complete_mapping(
        tmp_path, rows, second_latitude, message):
    head_file = tmp_path / "heads.txt"
    head_file.write_text(rows)
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; from vic import ffi, lib; "
         "cells = ffi.new('soil_con_struct[2]'); "
         "cells[0].lat = 48; cells[1].lat = float(sys.argv[2]); "
         "cells[0].lng = cells[1].lng = -120; "
         "lib.read_groundwater_head_file(sys.argv[1].encode(), cells, 2)",
         str(head_file), str(second_latitude)], capture_output=True, text=True)
    assert result.returncode == 1
    assert message in result.stderr
