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

    first_cell = ffi.new("soil_con_struct *")
    first_cell.lat = 48.25
    first_cell.lng = -120.75
    nloaded = vic_lib.read_groundwater_head_file(
        str(head_file).encode(), first_cell, 1)

    assert nloaded == 1
    assert first_cell.gw_head_offset_valid
    assert first_cell.gw_head_offset_m == 0.25

    second_cell = ffi.new("soil_con_struct *")
    second_cell.lat = 48.50
    second_cell.lng = -120.50
    nloaded = vic_lib.read_groundwater_head_file(
        str(head_file).encode(), second_cell, 1)

    assert nloaded == 1
    assert second_cell.gw_head_offset_valid
    assert second_cell.gw_head_offset_m == -0.75
