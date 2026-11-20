/******************************************************************************
* @section DESCRIPTION
*
* Calculate infiltration and runoff from the surface, gravity driven drainage
* between all soil layers, and generates baseflow from the bottom layer.
******************************************************************************/

#include <errno.h>

#include <vic_run.h>

static double
get_environment_double(const char *name,
                       double      default_value)
{
    const char *text;
    char       *end;
    double      value;

    text = getenv(name);
    if (text == NULL || text[0] == '\0') {
        return default_value;
    }

    errno = 0;
    value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value)) {
        log_err("%s must be a finite numeric value; received '%s'.", name,
                text);
    }

    return value;
}

/******************************************************************************
* @brief    Calculate infiltration and runoff from the surface, gravity driven
*           drainage between all soil layers, and generates baseflow from the
*           bottom layer.
******************************************************************************/
int
runoff(cell_data_struct  *cell,
       energy_bal_struct *energy,
       soil_con_struct   *soil_con,
       double             ppt,
       double            *frost_fract,
       int                Nnodes)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    size_t                     lindex;
    size_t                     time_step;
    int                        last_index;
    int                        tmplayer;
    int                        fidx;
    int                        ErrorFlag;
    double                     A, frac;
    double                     tmp_runoff;
    double                     inflow;
    double                     resid_moist[MAX_LAYERS]; // residual moisture (mm)
    double                     org_moist[MAX_LAYERS]; // total soil moisture (liquid and frozen) at beginning of this function (mm)
    double                     avail_liq[MAX_LAYERS][MAX_FROST_AREAS]; // liquid soil moisture available for evap/drainage (mm)
    double                     liq[MAX_LAYERS]; // current liquid soil moisture (mm)
    double                     ice[MAX_LAYERS]; // current frozen soil moisture (mm)
    double                     moist[MAX_LAYERS]; // current total soil moisture (liquid and frozen) (mm)
    double                     max_moist[MAX_LAYERS]; // maximum storable moisture (liquid and frozen) (mm)
    double                     Ksat[MAX_LAYERS];
    double                     Q12[MAX_LAYERS - 1];
    double                     Dsmax;
    double                     tmp_inflow;
    double                     tmp_moist;
    double                     tmp_moist_for_runoff[MAX_LAYERS];
    double                     tmp_liq;
    double                     dt_inflow;
    double                     dt_runoff;
    double                     runoff[MAX_FROST_AREAS];
    double                     tmp_dt_runoff[MAX_FROST_AREAS];
    double                     baseflow[MAX_FROST_AREAS];
    double                     gw_exchange[MAX_FROST_AREAS];
    double                     dt_baseflow;
    double                     dt_gw_exchange;
    const char                *gw_env;
    const char                *gw_formulation;
    const char                *gw_head_env;
    const char                *gw_length_env;
    bool                       gw_enabled;
    bool                       gw_head_formulation;
    double                     gw_test_gap_m = 1.0;
    double                     gw_test_ka_scale = 1.0;
    double                     gw_max_drain_fraction = 0.01;
    double                     soil_column_depth_m = 0.0;
    double                     gw_test_wtd_m = 0.0;
    double                     gw_reference_depth_m = 0.0;
    double                     gw_external_head_m = 0.0;
    double                     gw_exchange_length_m = 0.0;
    double                     rel_moist;
    double                     evap[MAX_LAYERS][MAX_FROST_AREAS];
    double                     sum_liq;
    double                     evap_fraction;
    double                     evap_sum;
    layer_data_struct         *layer;
    layer_data_struct          tmp_layer;
    unsigned short             runoff_steps_per_dt;

    /** Set Temporary Variables **/
    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        resid_moist[lindex] = soil_con->resid_moist[lindex];
        max_moist[lindex] = soil_con->max_moist[lindex];
        Ksat[lindex] = soil_con->Ksat[lindex] /
                       global_param.runoff_steps_per_day;
    }

    /** Allocate and Set Values for Soil Sublayers **/
    layer = cell->layer;

    gw_formulation = getenv("VIC_GW_FORMULATION");
    gw_enabled = gw_formulation != NULL && gw_formulation[0] != '\0';
    gw_head_formulation = false;

    if (gw_enabled) {
        if (strcmp(gw_formulation, "head") == 0) {
            gw_head_formulation = true;
        }
        else if (strcmp(gw_formulation, "scheidegger") != 0) {
            log_err("VIC_GW_FORMULATION must be 'scheidegger' or 'head'.");
        }

        gw_test_gap_m = get_environment_double("VIC_GW_TEST_GAP_M", 1.0);
        gw_test_ka_scale =
            get_environment_double("VIC_GW_TEST_KA_SCALE", 1.0);
        gw_max_drain_fraction =
            get_environment_double("VIC_GW_MAX_DRAIN_FRACTION", 0.01);

        if (gw_test_gap_m <= 0.0) {
            log_err("VIC_GW_TEST_GAP_M must be greater than zero.");
        }
        if (gw_test_ka_scale <= 0.0) {
            log_err("VIC_GW_TEST_KA_SCALE must be greater than zero.");
        }
        if (gw_max_drain_fraction <= 0.0 ||
            gw_max_drain_fraction > 1.0) {
            log_err("VIC_GW_MAX_DRAIN_FRACTION must be in (0, 1].");
        }

        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            soil_column_depth_m += soil_con->depth[lindex];
        }
        gw_test_wtd_m = soil_column_depth_m + gw_test_gap_m;
        gw_reference_depth_m = soil_column_depth_m;

        gw_env = getenv("VIC_GW_REFERENCE_DEPTH");
        if (gw_env != NULL && strcmp(gw_env, "base") != 0) {
            if (strcmp(gw_env, "node") == 0) {
                gw_reference_depth_m =
                    soil_column_depth_m -
                    0.5 * soil_con->depth[options.Nlayer - 1];
            }
            else {
                log_err("VIC_GW_REFERENCE_DEPTH must be 'base' or 'node'.");
            }
        }

        if (gw_head_formulation) {
            gw_head_env = getenv("VIC_GW_HEAD_OFFSET_M");
            gw_length_env = getenv("VIC_GW_EXCHANGE_LENGTH_M");

            if (soil_con->gw_head_offset_valid) {
                if (gw_head_env != NULL) {
                    log_err("Do not specify both VIC_GW_HEAD_FILE and "
                            "VIC_GW_HEAD_OFFSET_M.");
                }
                if (gw_length_env == NULL) {
                    log_err("VIC_GW_EXCHANGE_LENGTH_M is required with "
                            "VIC_GW_HEAD_FILE.");
                }
                gw_external_head_m =
                    -soil_column_depth_m + soil_con->gw_head_offset_m;
                gw_exchange_length_m = get_environment_double(
                    "VIC_GW_EXCHANGE_LENGTH_M", 0.0);
            }
            else if (gw_head_env == NULL && gw_length_env == NULL) {
                gw_external_head_m = -gw_test_wtd_m;
                gw_exchange_length_m =
                    gw_test_wtd_m - gw_reference_depth_m;
            }
            else if (gw_head_env == NULL || gw_length_env == NULL) {
                log_err("VIC_GW_HEAD_OFFSET_M and VIC_GW_EXCHANGE_LENGTH_M "
                        "must be specified together.");
            }
            else {
                gw_external_head_m = -soil_column_depth_m +
                                     get_environment_double(
                                         "VIC_GW_HEAD_OFFSET_M", 0.0);
                gw_exchange_length_m = get_environment_double(
                    "VIC_GW_EXCHANGE_LENGTH_M", 0.0);
            }

            if (gw_exchange_length_m <= 0.0) {
                log_err("VIC_GW_EXCHANGE_LENGTH_M must be greater than zero.");
            }
        }
    }


    cell->runoff = 0;
    cell->baseflow = 0;
    cell->gw_exchange = 0;
    cell->asat = 0;

    runoff_steps_per_dt = global_param.runoff_steps_per_day /
                          global_param.model_steps_per_day;

    for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
        baseflow[fidx] = 0;
        gw_exchange[fidx] = 0;
    }

    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        evap[lindex][0] = layer[lindex].evap / (double) runoff_steps_per_dt;
        org_moist[lindex] = layer[lindex].moist;
        layer[lindex].moist = 0;
        if (evap[lindex][0] > 0) { // if there is positive evaporation
            sum_liq = 0;
            // compute available soil moisture for each frost sub area.
            for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
                avail_liq[lindex][fidx] =
                    (org_moist[lindex] - layer[lindex].ice[fidx] -
                     resid_moist[lindex]);
                if (avail_liq[lindex][fidx] < 0) {
                    avail_liq[lindex][fidx] = 0;
                }
                sum_liq += avail_liq[lindex][fidx] *
                           frost_fract[fidx];
            }
            // compute fraction of available soil moisture that is evaporated
            if (sum_liq > 0) {
                evap_fraction = evap[lindex][0] / sum_liq;
            }
            else {
                evap_fraction = 1.0;
            }
            // distribute evaporation between frost sub areas by percentage
            evap_sum = evap[lindex][0];
            for (fidx = (int)options.Nfrost - 1; fidx >= 0; fidx--) {
                evap[lindex][fidx] = avail_liq[lindex][fidx] * evap_fraction;
                avail_liq[lindex][fidx] -= evap[lindex][fidx];
                evap_sum -= evap[lindex][fidx] * frost_fract[fidx];
            }
        }
        else {
            for (fidx = (int)options.Nfrost - 1; fidx > 0; fidx--) {
                evap[lindex][fidx] = evap[lindex][0];
            }
        }
    }

    for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
        /** ppt = amount of liquid water coming to the surface **/
        inflow = ppt;

        /**************************************************
           Initialize Variables
        **************************************************/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            /** Set Layer Liquid Moisture Content **/
            liq[lindex] = org_moist[lindex] - layer[lindex].ice[fidx];

            /** Set Layer Frozen Moisture Content **/
            ice[lindex] = layer[lindex].ice[fidx];
        }

        /******************************************************
           Runoff Based on Soil Moisture Level of Upper Layers
        ******************************************************/

        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_moist_for_runoff[lindex] = (liq[lindex] + ice[lindex]);
        }
        compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, inflow, &A,
                                &(runoff[fidx]));

        // save dt_runoff based on initial runoff estimate,
        // since we will modify total runoff below for the case of completely saturated soil
        tmp_dt_runoff[fidx] = runoff[fidx] /
                              (double) runoff_steps_per_dt;

        /**************************************************
           Compute Flow Between Soil Layers ()
        **************************************************/

        dt_inflow = inflow / (double) runoff_steps_per_dt;

        Dsmax = soil_con->Dsmax / global_param.runoff_steps_per_day;

        for (time_step = 0; time_step < runoff_steps_per_dt; time_step++) {
            inflow = dt_inflow;

            /*************************************
               Compute Drainage between Sublayers
            *************************************/

            for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
                /** Brooks & Corey relation for hydraulic conductivity **/

                if ((tmp_liq = liq[lindex] - evap[lindex][fidx]) <
                    resid_moist[lindex]) {
                    tmp_liq = resid_moist[lindex];
                }

                if (tmp_liq > resid_moist[lindex]) {
                    Q12[lindex] = calc_Q12(Ksat[lindex], tmp_liq,
                                           resid_moist[lindex],
                                           max_moist[lindex],
                                           soil_con->expt[lindex]);
                }
                else {
                    Q12[lindex] = 0.;
                }
            }

            /**************************************************
               Solve for Current Soil Layer Moisture, and
               Check Versus Maximum and Minimum Moisture Contents.
            **************************************************/

            last_index = 0;
            for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
                if (lindex == 0) {
                    dt_runoff = tmp_dt_runoff[fidx];
                }
                else {
                    dt_runoff = 0;
                }

                /* transport moisture for all sublayers **/

                tmp_inflow = 0.;

                /** Update soil layer moisture content **/
                liq[lindex] = liq[lindex] +
                              (inflow - dt_runoff) -
                              (Q12[lindex] + evap[lindex][fidx]);

                /** Verify that soil layer moisture is less than maximum **/
                if ((liq[lindex] + ice[lindex]) > max_moist[lindex]) {
                    tmp_inflow = (liq[lindex] + ice[lindex]) -
                                 max_moist[lindex];
                    liq[lindex] = max_moist[lindex] - ice[lindex];

                    if (lindex == 0) {
                        Q12[lindex] += tmp_inflow;
                        tmp_inflow = 0;
                    }
                    else {
                        tmplayer = lindex;
                        while (tmp_inflow > 0) {
                            tmplayer--;
                            if (tmplayer < 0) {
                                /** If top layer saturated, add to runoff **/
                                runoff[fidx] += tmp_inflow;
                                tmp_inflow = 0;
                            }
                            else {
                                /** else add excess soil moisture to next higher layer **/
                                liq[tmplayer] += tmp_inflow;
                                if ((liq[tmplayer] + ice[tmplayer]) >
                                    max_moist[tmplayer]) {
                                    tmp_inflow =
                                        ((liq[tmplayer] +
                                          ice[tmplayer]) - max_moist[tmplayer]);
                                    liq[tmplayer] = max_moist[tmplayer] -
                                                    ice[tmplayer];
                                }
                                else {
                                    tmp_inflow = 0;
                                }
                            }
                        }
                    } /** end trapped excess moisture **/
                } /** end check if excess moisture in top layer **/

                /** verify that current layer moisture is greater than minimum **/
                if (liq[lindex] < 0) {
                    /** liquid cannot fall below 0 **/
                    Q12[lindex] += liq[lindex];
                    liq[lindex] = 0;
                }
                if ((liq[lindex] + ice[lindex]) < resid_moist[lindex]) {
                    /** moisture cannot fall below minimum **/
                    Q12[lindex] +=
                        (liq[lindex] + ice[lindex]) - resid_moist[lindex];
                    liq[lindex] = resid_moist[lindex] - ice[lindex];
                }

                inflow = (Q12[lindex] + tmp_inflow);
                Q12[lindex] += tmp_inflow;

                last_index++;
            } /* end loop through soil layers */

            /**************************************************
               Compute Baseflow
            **************************************************/

            lindex = options.Nlayer - 1;
            dt_baseflow = 0.0;
            dt_gw_exchange = 0.0;

            /*
             * Coupled groundwater replaces VIC's conceptual baseflow
             * boundary. Uncoupled runs retain the standard ARNO closure.
             */
            if (!gw_enabled) {
                rel_moist =
                    (liq[lindex] - resid_moist[lindex]) /
                    (max_moist[lindex] - resid_moist[lindex]);

                frac = Dsmax * soil_con->Ds / soil_con->Ws;
                dt_baseflow = frac * rel_moist;
                if (rel_moist > soil_con->Ws) {
                    frac = (rel_moist - soil_con->Ws) /
                           (1 - soil_con->Ws);
                    dt_baseflow +=
                        Dsmax * (1 - soil_con->Ds / soil_con->Ws) *
                        pow(frac, soil_con->c);
                }
                if (dt_baseflow < 0) {
                    dt_baseflow = 0;
                }
            }

            liq[lindex] +=
                Q12[lindex - 1] -
                (evap[lindex][fidx] + dt_baseflow);

            tmp_moist = 0;

            if ((liq[lindex] + ice[lindex]) < resid_moist[lindex]) {
                dt_baseflow +=
                    (liq[lindex] + ice[lindex]) - resid_moist[lindex];
                liq[lindex] = resid_moist[lindex] - ice[lindex];
            }

            if (gw_enabled && lindex == options.Nlayer - 1) {
                if (gw_head_formulation) {
                    dt_gw_exchange =
                        calc_groundwater_exchange_head(
                            soil_con, lindex, liq[lindex], ice[lindex],
                            gw_reference_depth_m, gw_external_head_m,
                            gw_exchange_length_m,
                            soil_con->Ksat[lindex] * gw_test_ka_scale,
                            gw_max_drain_fraction,
                            global_param.runoff_steps_per_day);
                }
                else {
                    dt_gw_exchange =
                        calc_groundwater_exchange(
                            soil_con, lindex, liq[lindex], ice[lindex],
                            gw_test_wtd_m, gw_reference_depth_m,
                            soil_con->Ksat[lindex] * gw_test_ka_scale,
                            gw_max_drain_fraction,
                            global_param.runoff_steps_per_day);
                }

                /* Positive exchange removes water from the VIC soil column. */
                liq[lindex] -= dt_gw_exchange;
                gw_exchange[fidx] += dt_gw_exchange;
            }

            if ((liq[lindex] + ice[lindex]) > max_moist[lindex]) {
                /* soil moisture above maximum */
                tmp_moist = ((liq[lindex] + ice[lindex]) - max_moist[lindex]);
                liq[lindex] = max_moist[lindex] - ice[lindex];
                tmplayer = lindex;
                while (tmp_moist > 0) {
                    tmplayer--;
                    if (tmplayer < 0) {
                        /** If top layer saturated, add to runoff **/
                        runoff[fidx] += tmp_moist;
                        tmp_moist = 0;
                    }
                    else {
                        /** else if sublayer exists, add excess soil moisture **/
                        liq[tmplayer] += tmp_moist;
                        if ((liq[tmplayer] + ice[tmplayer]) >
                            max_moist[tmplayer]) {
                            tmp_moist =
                                ((liq[tmplayer] +
                                  ice[tmplayer]) - max_moist[tmplayer]);
                            liq[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
                        }
                        else {
                            tmp_moist = 0;
                        }
                    }
                }
            }

            baseflow[fidx] += dt_baseflow;
        } /* end of sub-dt time step loop */

        /** If negative baseflow, reduce evap accordingly **/
        if (baseflow[fidx] < 0) {
            layer[lindex].evap += baseflow[fidx];
            baseflow[fidx] = 0;
        }

        /** Recompute Asat based on final moisture level of upper layers **/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_moist_for_runoff[lindex] = (liq[lindex] + ice[lindex]);
        }
        compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, 0, &A,
                                &tmp_runoff);

        /** Store tile-wide values **/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            layer[lindex].moist +=
                ((liq[lindex] + ice[lindex]) * frost_fract[fidx]);
        }
        cell->asat += A * frost_fract[fidx];
        cell->runoff += runoff[fidx] * frost_fract[fidx];
        cell->baseflow += baseflow[fidx] * frost_fract[fidx];
        cell->gw_exchange +=
            gw_exchange[fidx] * frost_fract[fidx];
    }

    /** Compute water table depth **/
    wrap_compute_zwt(soil_con, cell);

    /** Recompute Thermal Parameters Based on New Moisture Distribution **/
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_layer = cell->layer[lindex];
            moist[lindex] = tmp_layer.moist;
        }

        ErrorFlag = distribute_node_moisture_properties(energy->moist,
                                                        energy->ice,
                                                        energy->kappa_node,
                                                        energy->Cs_node,
                                                        soil_con->Zsum_node,
                                                        energy->T,
                                                        soil_con->porosity_node,
                                                        soil_con->expt_node,
                                                        soil_con->bubble_node,
                                                        moist, soil_con->depth,
                                                        soil_con->soil_dens_min,
                                                        soil_con->bulk_dens_min,
                                                        soil_con->quartz,
                                                        soil_con->soil_density,
                                                        soil_con->bulk_density,
                                                        soil_con->organic,
                                                        Nnodes,
                                                        options.Nlayer,
                                                        soil_con->FS_ACTIVE);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }
    return (0);
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
compute_runoff_and_asat(soil_con_struct *soil_con,
                        double          *moist,
                        double           inflow,
                        double          *A,
                        double          *runoff)
{
    extern option_struct options;
    double               top_moist; // total moisture (liquid and frozen) in topmost soil layers (mm)
    double               top_max_moist; // maximum storable moisture (liquid and frozen) in topmost soil layers (mm)
    size_t               lindex;
    double               ex;
    double               max_infil;
    double               i_0;
    double               basis;

    top_moist = 0.;
    top_max_moist = 0.;
    for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
        top_moist += moist[lindex];
        top_max_moist += soil_con->max_moist[lindex];
    }
    if (top_moist > top_max_moist) {
        top_moist = top_max_moist;
    }

    /** A as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
    ex = soil_con->b_infilt / (1.0 + soil_con->b_infilt);
    *A = 1.0 - pow((1.0 - top_moist / top_max_moist), ex);

    max_infil = (1.0 + soil_con->b_infilt) * top_max_moist;
    i_0 = max_infil * (1.0 - pow((1.0 - *A), (1.0 / soil_con->b_infilt)));

    /** equation (3a) Wood et al. **/

    if (inflow == 0.0) {
        *runoff = 0.0;
    }
    else if (max_infil == 0.0) {
        *runoff = inflow;
    }
    else if ((i_0 + inflow) > max_infil) {
        *runoff = inflow - top_max_moist + top_moist;
    }
    /** equation (3b) Wood et al. (wrong in paper) **/
    else {
        basis = 1.0 - (i_0 + inflow) / max_infil;
        *runoff = (inflow - top_max_moist + top_moist +
                   top_max_moist *
                   pow(basis, 1.0 * (1.0 + soil_con->b_infilt)));
    }
    if (*runoff < 0.) {
        *runoff = 0.;
    }
}

/******************************************************************************
* @brief    Calculate drainage between two layers
******************************************************************************/
double
calc_Q12(double Ksat,
         double init_moist,
         double resid_moist,
         double max_moist,
         double expt)
{
    double Q12;

    Q12 = init_moist - pow(pow(init_moist - resid_moist, 1.0 - expt) -
                           Ksat /
                           pow(max_moist - resid_moist, expt) * (1.0 - expt),
                           1.0 / (1.0 - expt)) - resid_moist;

    return Q12;
}
