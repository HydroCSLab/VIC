/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate signed water exchange across the lower VIC soil boundary.
 *
 * Positive exchange removes water from VIC and recharges groundwater.
 * Negative exchange adds groundwater to the VIC soil column. Return values
 * are millimetres per runoff substep.
 *****************************************************************************/

#include <vic_run.h>

static bool
valid_exchange_inputs(const soil_con_struct *soil_con,
                      size_t                 lindex,
                      double                 liquid_moist_mm,
                      double                 ice_moist_mm,
                      double                 conductivity_mm_day,
                      double                 max_drain_fraction,
                      size_t                 runoff_steps_per_day)
{
    return soil_con != NULL &&
           lindex < MAX_LAYERS &&
           runoff_steps_per_day > 0 &&
           isfinite(liquid_moist_mm) &&
           isfinite(ice_moist_mm) &&
           isfinite(conductivity_mm_day) &&
           isfinite(max_drain_fraction) &&
           soil_con->depth[lindex] > 0.0 &&
           soil_con->max_moist[lindex] > 0.0 &&
           soil_con->bubble[lindex] > 0.0 &&
           soil_con->expt[lindex] > 3.0 &&
           conductivity_mm_day >= 0.0 &&
           max_drain_fraction > 0.0 &&
           max_drain_fraction <= 1.0;
}

static double
bottom_layer_pressure_head(const soil_con_struct *soil_con,
                           size_t                 lindex,
                           double                 liquid_moist_mm)
{
    const double min_relative_saturation = 1.e-8;
    double       relative_saturation;
    double       campbell_exponent;

    relative_saturation =
        liquid_moist_mm / soil_con->max_moist[lindex];
    relative_saturation =
        fmax(min_relative_saturation, fmin(1.0, relative_saturation));

    /*
     * VIC stores expt = 3 + 2 / lambda, while the Campbell retention
     * exponent is 1 / lambda. Bubbling pressure is a positive value in
     * centimetres; unsaturated pressure head is negative.
     */
    campbell_exponent = 0.5 * (soil_con->expt[lindex] - 3.0);
    return -0.01 * soil_con->bubble[lindex] *
           pow(relative_saturation, -campbell_exponent);
}

static double
limit_downward_exchange(const soil_con_struct *soil_con,
                        size_t                 lindex,
                        double                 liquid_moist_mm,
                        double                 ice_moist_mm,
                        double                 max_drain_fraction,
                        double                 exchange_mm)
{
    double available_above_residual_mm;
    double available_liquid_mm;
    double maximum_exchange_mm;

    if (exchange_mm <= 0.0) {
        return exchange_mm;
    }

    /*
     * resid_moist and liquid_moist_mm are both stored in millimetres.
     * Frozen water contributes to total storage above residual moisture,
     * but only liquid water can cross the interface.
     */
    available_above_residual_mm =
        fmax(0.0, liquid_moist_mm + ice_moist_mm -
             soil_con->resid_moist[lindex]);
    available_liquid_mm = fmax(0.0, liquid_moist_mm);
    maximum_exchange_mm =
        fmin(available_above_residual_mm, available_liquid_mm);
    maximum_exchange_mm =
        fmin(maximum_exchange_mm,
             max_drain_fraction * available_liquid_mm);

    return fmin(exchange_mm, maximum_exchange_mm);
}

/******************************************************************************
 * @brief Calculate exchange from a prescribed groundwater-table depth.
 *
 * Implements the below-soil-column formulation adapted from Scheidegger et
 * al. (2021), Eq. 3, and Niu et al. (2007). Depths are positive downward
 * from land surface.
 *****************************************************************************/
double
calc_groundwater_exchange(const soil_con_struct *soil_con,
                          size_t                 lindex,
                          double                 liquid_moist_mm,
                          double                 ice_moist_mm,
                          double                 groundwater_depth_m,
                          double                 hydraulic_reference_depth_m,
                          double                 aquifer_ksat_mm_day,
                          double                 max_drain_fraction,
                          size_t                 runoff_steps_per_day)
{
    const double eps = 1.e-12;
    double       pressure_head_m;
    double       exchange_mm_day;
    double       exchange_mm;

    if (!valid_exchange_inputs(soil_con, lindex, liquid_moist_mm,
                               ice_moist_mm, aquifer_ksat_mm_day,
                               max_drain_fraction,
                               runoff_steps_per_day) ||
        !isfinite(groundwater_depth_m) ||
        !isfinite(hydraulic_reference_depth_m) ||
        groundwater_depth_m <= hydraulic_reference_depth_m + eps) {
        return 0.0;
    }

    pressure_head_m =
        bottom_layer_pressure_head(soil_con, lindex, liquid_moist_mm);
    exchange_mm_day =
        -aquifer_ksat_mm_day *
        (-groundwater_depth_m -
         (pressure_head_m - hydraulic_reference_depth_m)) /
        (groundwater_depth_m - hydraulic_reference_depth_m);
    exchange_mm = exchange_mm_day / (double) runoff_steps_per_day;

    if (!isfinite(exchange_mm)) {
        return 0.0;
    }

    return limit_downward_exchange(
        soil_con, lindex, liquid_moist_mm, ice_moist_mm,
        max_drain_fraction, exchange_mm);
}

/******************************************************************************
 * @brief Calculate exchange from an externally supplied groundwater head.
 *
 * Hydraulic heads use a land-surface-relative datum:
 *
 *     H_soil = pressure_head - hydraulic_reference_depth
 *     q = K_interface (H_soil - H_groundwater) / exchange_length
 *****************************************************************************/
double
calc_groundwater_exchange_head(const soil_con_struct *soil_con,
                               size_t                 lindex,
                               double                 liquid_moist_mm,
                               double                 ice_moist_mm,
                               double                 hydraulic_reference_depth_m,
                               double                 groundwater_head_m,
                               double                 exchange_length_m,
                               double                 interface_ksat_mm_day,
                               double                 max_drain_fraction,
                               size_t                 runoff_steps_per_day)
{
    const double eps = 1.e-12;
    double       pressure_head_m;
    double       soil_head_m;
    double       exchange_mm_day;
    double       exchange_mm;

    if (!valid_exchange_inputs(soil_con, lindex, liquid_moist_mm,
                               ice_moist_mm, interface_ksat_mm_day,
                               max_drain_fraction,
                               runoff_steps_per_day) ||
        !isfinite(hydraulic_reference_depth_m) ||
        !isfinite(groundwater_head_m) ||
        !isfinite(exchange_length_m) ||
        exchange_length_m <= eps) {
        return 0.0;
    }

    pressure_head_m =
        bottom_layer_pressure_head(soil_con, lindex, liquid_moist_mm);
    soil_head_m = pressure_head_m - hydraulic_reference_depth_m;
    exchange_mm_day =
        interface_ksat_mm_day *
        (soil_head_m - groundwater_head_m) / exchange_length_m;
    exchange_mm = exchange_mm_day / (double) runoff_steps_per_day;

    if (!isfinite(exchange_mm)) {
        return 0.0;
    }

    return limit_downward_exchange(
        soil_con, lindex, liquid_moist_mm, ice_moist_mm,
        max_drain_fraction, exchange_mm);
}
