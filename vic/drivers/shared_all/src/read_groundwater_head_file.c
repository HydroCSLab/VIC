/******************************************************************************
 * @section DESCRIPTION
 *
 * Read spatial groundwater hydraulic-head offsets for VIC grid cells.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

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

    if (path == NULL || path[0] == '\0') {
        log_err("VIC_GW_HEAD_FILE must name a readable file.");
    }
    if (soil_con == NULL && ncells > 0) {
        log_err("Cannot load groundwater heads into a null soil array.");
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
        for (i = 0; i < ncells; i++) {
            if (fabs(soil_con[i].lat - latitude) <= coordinate_tolerance &&
                fabs(soil_con[i].lng - longitude) <= coordinate_tolerance) {
                if (soil_con[i].gw_head_offset_valid) {
                    log_err("Duplicate groundwater head for VIC cell "
                            "lat=%.12f lon=%.12f at %s:%zu.",
                            soil_con[i].lat, soil_con[i].lng, path,
                            line_number);
                }
                soil_con[i].gw_head_offset_m = head_offset_m;
                soil_con[i].gw_head_offset_valid = true;
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

    for (i = 0; i < ncells; i++) {
        if (!soil_con[i].gw_head_offset_valid) {
            log_err("No groundwater head found for local VIC cell "
                    "lat=%.12f lon=%.12f in %s.",
                    soil_con[i].lat, soil_con[i].lng, path);
        }
    }

    return nloaded;
}
