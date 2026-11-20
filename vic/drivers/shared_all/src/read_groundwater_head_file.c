/******************************************************************************
 * @section DESCRIPTION
 *
 * Read spatial groundwater hydraulic-head offsets for VIC grid cells.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

typedef struct {
    double lat;
    double lng;
    size_t index;
} groundwater_head_cell_index;

static int
compare_groundwater_head_cells(const void *first, const void *second)
{
    const groundwater_head_cell_index *a;
    const groundwater_head_cell_index *b;

    a = (const groundwater_head_cell_index *) first;
    b = (const groundwater_head_cell_index *) second;
    if (a->lat < b->lat) {
        return -1;
    }
    if (a->lat > b->lat) {
        return 1;
    }
    if (a->lng < b->lng) {
        return -1;
    }
    if (a->lng > b->lng) {
        return 1;
    }
    return 0;
}

static size_t
lower_bound_groundwater_latitude(const groundwater_head_cell_index *cells,
                                 size_t ncells, double latitude,
                                 double coordinate_tolerance)
{
    size_t lower;
    size_t upper;

    lower = 0;
    upper = ncells;
    while (lower < upper) {
        size_t middle;

        middle = lower + (upper - lower) / 2;
        if (latitude - cells[middle].lat > coordinate_tolerance) {
            lower = middle + 1;
        }
        else {
            upper = middle;
        }
    }
    return lower;
}

/******************************************************************************
 * @brief Read a whitespace-delimited latitude, longitude, head-offset file.
 *
 * Records for cells owned by other MPI ranks are ignored. Every local active
 * cell must have exactly one matching record. The returned count is therefore
 * equal to ncells when the function succeeds.
 *****************************************************************************/
size_t
read_groundwater_head_file(const char      *path,
                           soil_con_struct *soil_con,
                           size_t           ncells)
{
    const double coordinate_tolerance = 1.e-8;
    FILE        *stream;
    char         line[MAXSTRING];
    size_t       line_number;
    size_t       nloaded;
    size_t       i;
    groundwater_head_cell_index *cell_index;

    if (path == NULL || path[0] == '\0') {
        log_err("VIC_GW_HEAD_FILE must name a readable file.");
    }
    if (soil_con == NULL && ncells > 0) {
        log_err("Cannot load groundwater heads into a null soil array.");
    }

    cell_index = malloc(ncells * sizeof(*cell_index));
    if (cell_index == NULL && ncells > 0) {
        log_err("Unable to allocate groundwater-head coordinate index for "
                "%zu VIC cells.", ncells);
    }
    for (i = 0; i < ncells; i++) {
        if (!isfinite(soil_con[i].lat) || !isfinite(soil_con[i].lng)) {
            log_err("Non-finite coordinates for local VIC cell %zu.", i);
        }
        cell_index[i].lat = soil_con[i].lat;
        cell_index[i].lng = soil_con[i].lng;
        cell_index[i].index = i;
    }
    if (ncells > 1) {
        qsort(cell_index, ncells, sizeof(*cell_index),
              compare_groundwater_head_cells);
    }

    stream = fopen(path, "r");
    if (stream == NULL) {
        log_err("Unable to open VIC_GW_HEAD_FILE: %s", path);
    }

    line_number = 0;
    nloaded = 0;
    while (fgets(line, sizeof(line), stream) != NULL) {
        char   *text;
        char    trailing;
        double  latitude;
        double  longitude;
        double  head_offset_m;
        int     nfields;
        size_t  matches;

        line_number++;
        if (strchr(line, '\n') == NULL && !feof(stream)) {
            log_err("VIC_GW_HEAD_FILE line exceeds %d characters at %s:%zu.",
                    MAXSTRING - 1, path, line_number);
        }

        text = line;
        while (*text == ' ' || *text == '\t') {
            text++;
        }
        if (*text == '\0' || *text == '\n' || *text == '#') {
            continue;
        }

        nfields = sscanf(text, "%lf %lf %lf %c", &latitude, &longitude,
                         &head_offset_m, &trailing);
        if (nfields < 3 || (nfields == 4 && trailing != '#')) {
            log_err("Expected 'latitude longitude head_offset_m' at %s:%zu.",
                    path, line_number);
        }
        if (!isfinite(latitude) || !isfinite(longitude) ||
            !isfinite(head_offset_m)) {
            log_err("Non-finite groundwater-head value at %s:%zu.", path,
                    line_number);
        }

        matches = 0;
        /* Compare differences directly to preserve the original tolerance. */
        i = lower_bound_groundwater_latitude(
            cell_index, ncells, latitude, coordinate_tolerance);
        for (; i < ncells &&
             cell_index[i].lat - latitude <= coordinate_tolerance; i++) {
            size_t local_index;

            local_index = cell_index[i].index;
            if (fabs(cell_index[i].lat - latitude) <= coordinate_tolerance &&
                fabs(cell_index[i].lng - longitude) <= coordinate_tolerance) {
                if (soil_con[local_index].gw_head_offset_valid) {
                    log_err("Duplicate groundwater head for VIC cell "
                            "lat=%.12f lon=%.12f at %s:%zu.",
                            soil_con[local_index].lat,
                            soil_con[local_index].lng, path,
                            line_number);
                }
                soil_con[local_index].gw_head_offset_m = head_offset_m;
                soil_con[local_index].gw_head_offset_valid = true;
                matches++;
                nloaded++;
            }
        }

        if (matches > 1) {
            log_err("Groundwater-head record %s:%zu matched more than one "
                    "local VIC cell.", path, line_number);
        }
    }

    if (ferror(stream)) {
        log_err("Error while reading VIC_GW_HEAD_FILE: %s", path);
    }
    fclose(stream);
    free(cell_index);

    for (i = 0; i < ncells; i++) {
        if (!soil_con[i].gw_head_offset_valid) {
            log_err("No groundwater head found for local VIC cell "
                    "lat=%.12f lon=%.12f in %s.",
                    soil_con[i].lat, soil_con[i].lng, path);
        }
    }

    return nloaded;
}
