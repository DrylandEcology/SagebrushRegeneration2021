#!/usr/bin/env Rscript

################################################################################
#    Understanding future big sagebrush regeneration potential
#
#    Daniel R Schlaepfer (https://orcid.org/0000-0001-9973-2065)
#
################################################################################
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################



#--- Load packages, custom functions, and description of simulation experiment
source("SagebrushRegeneration_2_DescribeExperimentAnalysis.R")

pkgs <- c("reshape2", "ClustOfVar", "lme4", "MuMIn", "rSW2utils")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
# Research question:
# 3.	Can we understand model outcomes in terms of environmental site
#     characteristics like climate, soil moisture, and soils?
#   b.	Identify a parsimonious statistical relationship between climate and
#       soil moisture variables and model outcomes for the historical period
#       and site-specific soils?



#--- Specify variables
defs_q3d_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q3d_models, "analysis_group"]
models <- vars_table[defs_q3d_models, "name_original_rSFSW2"]
model_labels <- vars_table[defs_q3d_models, "label_manuscript"]

coords <- c("X_WGS84", "Y_WGS84")
defs_q3d_coords <- which(vars_table[, "name_original_rSFSW2"] %in% coords)

preds <- c(
  "SWinput_Soil_topLayers_Sand_fraction",
  "SWinput_Soil_bottomLayers_Sand_fraction",
  "SWinput_Soil_topLayers_Clay_fraction",
  "SWinput_Soil_bottomLayers_Clay_fraction",
  "MAT_C_mean",
  "MAP_mm_mean",
  "SnowOfPPT_fraction_mean",
  "Seasonality_monthlyTandPPT_PearsonCor_mean",
  "TeeriEtAl1976_NSadj_FreezeFreeGrowingPeriod_days_mean",
  "DegreeDays_Base0C_dailyTmean_Cdays_mean",
  "UNAridityIndex_Normals_none_mean",
  "MAT_C_doyRange1to250_mean",
  "MAP_mm_doyRange1to250_mean",
  "SnowOfPPT_fraction_doyRange1to250_mean",
  "Snowcover_NSadj_Peak_mmSWE_mean",
  "periodicVWCmatricMean_FirstLayer_doyRange71to100_mean",
  "TminBelow0CwithoutSpringSnow_days_mean",
  "TminBelowNeg9CwithoutSnow_days_mean",
  "TmaxAbovePos34C_days_mean",
  "RelRecharge_topLayers_DailyMax_Fraction_mean",
  "RelRecharge_bottomLayers_DailyMax_Fraction_mean",
  "RelRecharge_topLayers_DailyMax_doy_mean",
  "RelRecharge_bottomLayers_DailyMax_doy_mean",
  "ThermalSnowfreeDryPeriods_SWPcrit3000kPa_topLayers_DrySpellsAllLayers_maxDuration_days_mean",
  "ThermalSnowfreeDryPeriods_SWPcrit3000kPa_bottomLayers_DrySpellsAllLayers_maxDuration_days_mean"
)
defs_q3d_preds <- which(vars_table[, "name_original_rSFSW2"] %in% preds)
preds_labels <- vars_table[defs_q3d_preds, "label_manuscript"]

defs_q3d <- c(defs_q3d_models, defs_q3d_coords, defs_q3d_preds)




#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

fname_dats1 <- file.path(dir_res_data, "tmp_dbOut", "dats_q3d_processed.rds")

if (file.exists(fname_dats1)) {
  dats_q3d <- readRDS(file = fname_dats1)

} else {

  ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q3d.rds")

  if (file.exists(ftmp)) {
    dats_q3d <- readRDS(ftmp)

  } else {
    dats_q3d <- list(
      vals = load_rSFSW2_BSR_data_from_netCDFs(
        path = dir_prj,
        locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
        vars_table = vars_table,
        ids_vars = defs_q3d,
        sim_scenarios = sim_scens_id,
        sim_subprojects = default_subprj
      )
    )

    dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
    saveRDS(dats_q3d, file = ftmp)
  }


  #--- Calculate differences between future and historical time periods
  dats_q3d[["delta_Current_abs"]] <-
   rSW2analysis::calc_change_from_reference(
      data = dats_q3d[["vals"]],
      ref_condition = sc_hist,
      sc_hist = sc_hist,
      method = "absolute"
    )[[1]]


  #--- Calculate ensembles for dTime x RCPs
  dats_q3d[["ens"]] <- rSW2analysis::calc_cellwise_ensemble(
    data = dats_q3d[["delta_Current_abs"]],
    funs = cfun,
    na.rm = TRUE,
    verbose = TRUE
  )

  dir.create(dirname(fname_dats1), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q3d, file = fname_dats1)
}



#------------------------------------------------------------------------------#
#--- Analyse data

do_select_variables <- FALSE #TRUE #
do_dcor_estimates <- FALSE
do_fit_relationships <- TRUE
do_resid_plots <- FALSE
do_plot_relationships <- TRUE



dir_res2_q3d <- file.path(
  dir_res2,
  "Output3_EnvRelationships",
  "FutureDelta_GeneralEnvRelationships"
)

dir.create(dir_res2_q3d, recursive = TRUE, showWarnings = FALSE)

legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)

preds_all <- c(preds, paste0("delta_", preds))
var_labs_all <- c(preds_labels, paste("delta", preds_labels))
models_delta <- paste0("delta_", models)


#--- Identify best predictors and prepare data
nclusters_max_apriori <- 4L

msdata <- list()

dir_tmp <- file.path(dir_res_data, "tmp_Part3d_futdiffenvrel")
dir.create(dir_tmp, recursive = TRUE, showWarnings = FALSE)


for (ks in seq_along(models)) {
  msdata[[ks]] <- list()

  msdata[[ks]][["tag_subset"]] <- names(def_subsets)[used_subsets[ks]]


  #------ Prepare data (part 1): temporary data to transform/select variables
  msdata[[ks]][["ids_data"]] <- which(complete.cases(
    dats_q3d[["vals"]][, models[ks], sc_hist, default_subprj]
  ))


  data_hist <- reshape2::melt(
    dats_q3d[["vals"]][msdata[[ks]][["ids_data"]], preds, sc_hist, default_subprj]
  )
  colnames(data_hist)[1:2] <- c("Gridcell", "Variable")
  data_hist[, "Variable"] <- as.character(data_hist[, "Variable"])
  data_hist[, "Scenario"] <- NA

  # repeat historical data as basis for each ensemble unit
  data_hist2 <- as.data.frame(matrix(
    NA,
    nrow = length(ensemble_units) * nrow(data_hist),
    ncol = ncol(data_hist),
    dimnames = list(NULL, colnames(data_hist))
  ))
  ids <- seq_len(nrow(data_hist))

  for (ie in seq_along(ensemble_units)) {
    id_eus <- ensemble_units[ie]

    tmp <- data_hist
    tmp[, "Scenario"] <- id_eus
    data_hist2[ids + (ie - 1) * nrow(data_hist), ] <- tmp
  }

  data_ens <- reshape2::melt(dats_q3d[["ens"]][cfun, msdata[[ks]][["ids_data"]], c(models[ks], preds), , default_subprj])
  colnames(data_ens)[1:3] <- c("Gridcell", "Variable", "Scenario")
  data_ens[, "Variable"] <- paste0(
    "delta_",
    as.character(data_ens[, "Variable"])
  )
  data_ens[, "Scenario"] <- as.character(data_ens[, "Scenario"])

  tmp_data <- reshape2::dcast(
    rbind(data_hist2, data_ens),
    Gridcell + Scenario ~ Variable
  )

  # Remove variables without variation
  iuse <- apply(tmp_data[, preds_all], 2, sd, na.rm = TRUE) > 0
  tmp <- colnames(tmp_data)
  tmp_data <- tmp_data[, tmp[!(tmp %in% names(iuse[!iuse]))]]
  preds_all_used <- preds_all[iuse]
  var_labs_used <- var_labs_all[iuse]


  #--- Check for need of variable transformations
  if (FALSE) {
    plot_correlation_checks(
      data = tmp_data,
      predictors = preds_all_used,
      responses = models_delta[ks],
      path = dir_res2_q3d,
      ftag = paste0(
        "FutureDiffEnvRelationships_untransformed_",
        msdata[[ks]][["tag_subset"]]
      )
    )
  }

  #--- Transform variables (same as for `GeneralEnvRelationships`)
  transform_predictors <- list(
    var_preds_toreplace = c(
      "MAP_mm_mean",
      "TminBelowNeg9CwithoutSnow_days_mean",
      "TminBelow0CwithoutSpringSnow_days_mean",
      "MAP_mm_doyRange1to250_mean",
      "RelRecharge_topLayers_DailyMax_doy_mean",
      "RelRecharge_bottomLayers_DailyMax_doy_mean",
      "Snowcover_NSadj_Peak_mmSWE_mean",
      "SnowOfPPT_fraction_mean",
      "UNAridityIndex_Normals_none_mean",
      "TmaxAbovePos34C_days_mean",
      "RelRecharge_topLayers_DailyMax_Fraction_mean",
      "RelRecharge_bottomLayers_DailyMax_Fraction_mean",
      "ThermalSnowfreeDryPeriods_SWPcrit3000kPa_topLayers_DrySpellsAllLayers_maxDuration_days_mean",
      "ThermalSnowfreeDryPeriods_SWPcrit3000kPa_bottomLayers_DrySpellsAllLayers_maxDuration_days_mean"
    ),

    var_preds_transformed = c(
      "log_MAP_mm_mean",
      "log_inyrs_plus_TminBelowNeg9CwithoutSnow_days_mean",
      "log_inyrs_plus_TminBelow0CwithoutSpringSnow_days_mean",
      "log_MAP_mm_doyRange1to250_mean",
      "RelRecharge_topLayers_DailyMax_WaterYeardoy_mean",
      "RelRecharge_bottomLayers_DailyMax_WaterYeardoy_mean",
      "log_inyrs_plus_Snowcover_NSadj_Peak_mmSWE_mean",
      "sqrt_inyrs_plus_SnowOfPPT_fraction_mean",
      "log_UNAridityIndex_Normals_none_mean",
      "log_inyrs_plus_TmaxAbovePos34C_days_mean",
      "exp_RelRecharge_topLayers_DailyMax_Fraction_mean",
      "exp_RelRecharge_bottomLayers_DailyMax_Fraction_mean",
      "sqrt_ThermalSnowfreeDryPeriods_SWPcrit3000kPa_topLayers_DrySpellsAllLayers_maxDuration_days_mean",
      "sqrt_ThermalSnowfreeDryPeriods_SWPcrit3000kPa_bottomLayers_DrySpellsAllLayers_maxDuration_days_mean"
    ),

    var_labs_transformed = c(
      "log(MAP [mm])",
      "log(1/nyrs + Hard frost exposure [days])",
      "log(1/nyrs + Spring frost exposure [days])",
      "log(MAP [doy 1-250; mm])",
      "Shallow recharge timing (water-year doy)",
      "Deep recharge timing (water-year doy)",
      "log(1/nyrs + Peak snowpack [mm SWE])",
      "sqrt(1/nyrs + Snow/PPT)",
      "log(PPT/PET)",
      "log(1/nyrs + Heat exposure [days])",
      "exp(Shallow recharge [%])",
      "exp(Deep recharge [%])",
      "sqrt(Shallow growing season drought [days])",
      "sqrt(Deep growing season drought [days])"
    ),

    fun_transform = c(
      function(x) log(x),
      function(x) log(1 / nyrs + x),
      function(x) log(1 / nyrs + x),
      function(x) log(x),
      function(x) (x - 273) %% 365,
      function(x) (x - 273) %% 365,
      function(x) log(1 / nyrs + x),
      function(x) sqrt(1 / nyrs + x),
      function(x) log(x),
      function(x) log(1 / nyrs + x),
      function(x) exp(x),
      function(x) exp(x),
      function(x) sqrt(x),
      function(x) sqrt(x)
    )
  )

  for (it in seq_along(transform_predictors[["var_preds_toreplace"]])) {
    ids_notfinite1 <- !is.finite(
      tmp_data[, transform_predictors[["var_preds_toreplace"]][it]]
    )

    tmp_data[, transform_predictors[["var_preds_transformed"]][it]] <-
      transform_predictors[["fun_transform"]][[it]](
        tmp_data[, transform_predictors[["var_preds_toreplace"]][it]]
      )

    ids_notfinite2 <- !is.finite(
      tmp_data[, transform_predictors[["var_preds_transformed"]][it]]
    )

    # Stop if new non-finite values added
    stopifnot(sum(ids_notfinite2 & !ids_notfinite1) == 0)
  }

  if (TRUE) {
    # Remove original untransformed variables
    ids <- preds_all_used %in% transform_predictors[["var_preds_toreplace"]]
    preds_all_used <- c(preds_all_used[!ids], transform_predictors[["var_preds_transformed"]])
    var_labs_used <- c(var_labs_used[!-ids], transform_predictors[["var_labs_transformed"]])
  } else {
    preds_all_used <- c(preds_all_used, transform_predictors[["var_preds_transformed"]])
    var_labs_used <- c(var_labs_used, transform_predictors[["var_labs_transformed"]])
  }

  tmp_data <- tmp_data[, c(models_delta[ks], preds_all_used)]

  if (FALSE) {
    plot_correlation_checks(
      data = tmp_data,
      predictors = preds_all_used,
      responses = models_delta[ks],
      path = dir_res2_q3d,
      ftag = paste0(
        "FutureDiffEnvRelationships_transformed_",
        msdata[[ks]][["tag_subset"]]
      )
    )
  }


  #------ Variable selection: best representative variable from each cluster
  # on the merged ensemble dataset (ideally, done for each scenario-combination,
  # but too costly)
  if (do_select_variables) {

    #--- Determine bivariate unbiased distance correlation
    sel_preds_info <- calc_model_dcor2d(
      data = tmp_data,
      predictors = preds_all_used,
      responses = models_delta[ks]
    )


    #--- Ascendant hierarchical clustering of variables
    sel_hclust <- calc_asc_hierarch_clust(
      data = tmp_data,
      predictors = preds_all_used,
      do_stability = TRUE,
      path = dir_res2_q3d,
      ftag = paste0(
        "FutureDiffEnvRelationships_merged-ensemble_",
        msdata[[ks]][["tag_subset"]]
      )
    )


    #--- Cut tree into selected number of clusters
    sel_clusters <- ClustOfVar::cutreevar(
      sel_hclust[["hclust"]],
      k = if (isTRUE(is.na(sel_hclust[["nclusters"]]))) {
        2 * nclusters_max_apriori
      } else {
        sel_hclust[["nclusters"]]
      }
    )

    if (FALSE) {
      summary(sel_clusters)
    }

    #--- Locate best variable for each cluster
    sel_preds_info <- locate_preds_per_cluster(
      hclust_sel = sel_clusters,
      predictor_info = sel_preds_info
    )

    #--- Determine majority consensus among models
    tmp <- consensus_best_preds(
      predictor_info = sel_preds_info,
      data = tmp_data,
      nclusters = min(
        nclusters_max_apriori,
        sel_hclust[["nclusters"]],
        na.rm = TRUE
      ),
      limit_cor = 0.5
    )

    msdata[[ks]][["preds_selected"]] <- tmp[["kpreds_selected"]]


    write.csv(
      reshape2::dcast(
        reshape2::melt(tmp[["predictor_info"]]),
        Var1 ~ Var2 + Var3
      ),
      row.names = FALSE,
      file = file.path(
        dir_res2_q3d,
        paste0(
          "Table_FutureDiffEnvRelationships_selected_merged-ensemble_",
          msdata[[ks]][["tag_subset"]],
          "_nclusters", nclusters_max_apriori,
          ".csv"
        )
      )
    )


    #--- Check for residual correlations
    plot_correlation_checks(
      data = tmp_data,
      predictors = msdata[[ks]][["preds_selected"]],
      responses = models_delta[ks],
      path = dir_res2_q3d,
      ftag = paste0(
        "FutureDiffEnvRelationships_selected_merged-ensemble_",
        msdata[[ks]][["tag_subset"]]
      )
    )

  } else {
    msdata[[ks]][["preds_selected"]] <- if (ks == 1) {
      c(
        "DegreeDays_Base0C_dailyTmean_Cdays_mean",
        "delta_SnowOfPPT_fraction_doyRange1to250_mean",
        "delta_MAP_mm_mean",
        "delta_Snowcover_NSadj_Peak_mmSWE_mean"
      )
    } else if (ks == 2) {
      c(
        "delta_DegreeDays_Base0C_dailyTmean_Cdays_mean",
        "delta_MAP_mm_mean",
        "delta_RelRecharge_topLayers_DailyMax_doy_mean",
        "delta_RelRecharge_bottomLayers_DailyMax_Fraction_mean"
      )
    }
  }


  #--- Prepare data (part 3) for estimating relationships
  # on each scenario-combination as random factor
  ids_delta <- grep("delta_", msdata[[ks]][["preds_selected"]])
  preds_selected_delta <- msdata[[ks]][["preds_selected"]][ids_delta]

  preds_selected_historical <- msdata[[ks]][["preds_selected"]][-ids_delta]
  tmp <- preds_selected_historical
  ids_preds_selected_historical_transformed <- tmp %in% transform_predictors[["var_preds_transformed"]]
  if (any(ids_preds_selected_historical_transformed)) {
    ids2 <- match(
      preds_selected_historical[ids_preds_selected_historical_transformed],
      transform_predictors[["var_preds_transformed"]],
      nomatch = 0
    )
    tmp[ids_preds_selected_historical_transformed] <- transform_predictors[["var_preds_toreplace"]][ids2]
  }
  preds_selected_historical_untransformed <- tmp


  ids <- match(msdata[[ks]][["preds_selected"]], preds_all_used, nomatch = 0)
  msdata[[ks]][["labs_preds_selected"]] <- var_labs_used[ids]


  fname_dats2 <- file.path(dir_tmp, paste0("Part3d_tmp2_model", ks, ".rds"))

  if (file.exists(fname_dats2)) {
    msdata[[ks]][["tmp_data"]] <- readRDS(file = fname_dats2)

  } else {
    #--- Historical data
    if (length(preds_selected_historical_untransformed) > 0) {
      data_hist <- reshape2::melt(dats_q3d[["vals"]][, preds_selected_historical_untransformed, sc_hist, default_subprj, drop = FALSE])[, -(3:4)]
      colnames(data_hist)[1:2] <- c("Gridcell", "Variable")
      data_hist[, "Variable"] <- as.character(data_hist[, "Variable"])
      data_hist[, "Scenario"] <- NA

      id_use <- data_hist[, "Gridcell"] %in% msdata[[ks]][["ids_data"]]
      data_hist <- data_hist[id_use, ]

      # Transform
      if (any(ids_preds_selected_historical_transformed)) {
        ids2 <- match(
          preds_selected_historical[ids_preds_selected_historical_transformed],
          transform_predictors[["var_preds_transformed"]],
          nomatch = 0)

        for (it in seq_along(ids2)) {
          ids_var <- data_hist[, "Variable"] %in% transform_predictors[["var_preds_toreplace"]][ids2[it]]
          ids_notfinite1 <- !is.finite(data_hist[ids_var, "value"])

          xtmp <- transform_predictors[["fun_transform"]][[ids2[it]]](
              data_hist[ids_var, "value"]
            )

          # Stop if new non-finite values added
          ids_notfinite2 <- !is.finite(xtmp)
          stopifnot(sum(ids_notfinite2 & !ids_notfinite1) == 0)

          # Update values and variable name
          data_hist[ids_var, "value"] <- xtmp
          data_hist[ids_var, "Variable"] <-
            transform_predictors[["var_preds_transformed"]][ids2[it]]
        }
      }

      # repeat historical data as basis for each scenario-combination
      data_hist2 <- as.data.frame(matrix(
        NA,
        nrow = length(sim_scens_id[-1]) * nrow(data_hist),
        ncol = ncol(data_hist),
        dimnames = list(NULL, colnames(data_hist))
      ))
      ids <- seq_len(nrow(data_hist))

      for (isc in 1 + seq_along(sim_scens_id[-1])) {
        id_scen <- sim_scens_id[isc]

        tmp <- data_hist
        tmp[, "Scenario"] <- id_scen
        data_hist2[ids + (isc - 2) * nrow(data_hist), ] <- tmp
      }

    } else {
      data_hist2 <- NULL
    }

    #--- Scenario data
    if (length(preds_selected_delta) > 0) {
      data_scen <- reshape2::melt(dats_q3d[["delta_Current_abs"]][,  c(models[ks], sub("delta_", "", preds_selected_delta)), -1, default_subprj])
      colnames(data_scen)[1:3] <- c("Gridcell", "Variable", "Scenario")
      data_scen[, "Variable"] <- paste0("delta_", as.character(data_scen[, "Variable"]))
      data_scen[, "Scenario"] <- as.character(data_scen[, "Scenario"])

      id_use <- data_scen[, "Gridcell"] %in% msdata[[ks]][["ids_data"]]
      data_scen <- data_scen[id_use, ]

    } else {
      data_scen <- NULL
    }

    msdata[[ks]][["tmp_data"]] <- reshape2::dcast(
      rbind(data_hist2, data_scen),
      Gridcell + Scenario ~ Variable
    )

    saveRDS(msdata[[ks]][["tmp_data"]], file = fname_dats2)
  }
}


#--- (Partial) distance correlations
if (do_dcor_estimates) {
  res_dcor <- get_dcor_tests(
    data = tmp_data,
    models = models,
    model_names = model_names,
    preds_all = preds_all,
    preds_selected = preds_selected,
    R = 99,  #199 exhaust 100 GB memory
    verbose = TRUE
  )

  write.csv(res_dcor,
    row.names = FALSE,
    file = file.path(dir_res2_q3d,
      paste0("Table_DistanceCorrelations_", sc_hist, "_", msdata[[ks]][["tag_subset"]], ".csv"))
  )
}


#--- Model fitting LMM for each RCPs x dTime
if (do_fit_relationships) {
  res_fit <- list()
  tmp <- vector("list", length(ensemble_units))
  names(tmp) <- ensemble_units

  for (iy in model_names) {
    res_fit[[iy]] <- tmp
  }

  #--- Loop over RCPs x dTime
  for (ie in seq_along(ensemble_units)) {
    id_eus <- ensemble_units[ie]

    #--- Loop over models
    for (iy in seq_along(models)) {
      #--- Prepare data (part 2)
      ids_use <- grep(id_eus, msdata[[iy]][["tmp_data"]][, "Scenario"])

      tmp <- data.frame(
        y = msdata[[iy]][["tmp_data"]][ids_use, models_delta[iy]],
        GCM = rSW2analysis:::find_reqMs(
          msdata[[iy]][["tmp_data"]][ids_use, "Scenario"]
        ),
        scale(msdata[[iy]][["tmp_data"]][ids_use, msdata[[iy]][["preds_selected"]]])
      )

      ids_use2 <- complete.cases(tmp)

      #--- Prepare data (part 3)
      datax_scaled <- scale(msdata[[iy]][["tmp_data"]][ids_use[ids_use2], msdata[[iy]][["preds_selected"]]])
      datax_scale <- attr(datax_scaled, "scaled:scale")
      datax_center <- attr(datax_scaled, "scaled:center")

      m_data <- data.frame(
        y = msdata[[iy]][["tmp_data"]][ids_use[ids_use2], models_delta[iy]],
        GCM = rSW2analysis:::find_reqMs(
          msdata[[iy]][["tmp_data"]][ids_use[ids_use2], "Scenario"]
        ),
        datax_scaled
      )

      msdata[[iy]][["data_ID"]] <- msdata[[iy]][["tmp_data"]][ids_use[ids_use2], "Gridcell"]
      msdata[[ks]][["ids_data2"]] <- unique(msdata[[iy]][["data_ID"]])



      #--- Fit models
      # Model with quadratic terms and 2-way interaction plus GCMs as random intercept effects
      y_formula1 <- stats::as.formula(paste(
        "y ~",
        "(", paste0(msdata[[iy]][["preds_selected"]], collapse = "+"), ") ^ 2",
        "+",
        paste0("I(", msdata[[iy]][["preds_selected"]], "^2)", collapse = "+"),
        "+ (1 | GCM)"
      ))

      m1 <- lme4::lmer(y_formula1, data = m_data)

      # McFadden's R squared
      nullmod <- lme4::lmer(
        y ~ 1 + (1 | GCM),
        data = m_data
      )

      pseudoR2_McFadden <- as.numeric(1 - stats::logLik(m1) / stats::logLik(nullmod))

      # Nakagawa's R-squared
      pseudoR2_Nakagawa <- MuMIn::r.squaredGLMM(m1)

      # Jaeger's R-squared
      pseudoR2_Jaeger <- array(dim = c(1, 1), dimnames = list(NULL, "Rsq")) # r2glmm::r2beta(m1, partial = FALSE)

      # Root-square mean deviations (in-sample)
      rmsd <- rSW2utils::rmse(m_data[, "y"], fitted(m1), na.rm = TRUE)

      cv_rmsd <- rmsd / mean(m_data[, "y"], na.rm = TRUE)

      # Check models
      if (FALSE) {
        # Walk through code of `check_models` interactively
        m1_checks <- check_models(
          model = m1,
          preds = msdata[[iy]][["preds_selected"]],
          data = m_data,
          datay = m_data[, "y"],
          datax_scaled = datax_scaled
        )
      }

      # Save model for later
      res_fit[[model_names[iy]]][[id_eus]] <- list(
        m = m1,
        datay = m_data[, "y"],
        datax_scaled = datax_scaled,
        datax_scale = datax_scale,
        datax_center = datax_center,
        coef_standard = summary(m1)[["coefficients"]],
        pseudoR2_McFadden = pseudoR2_McFadden,
        pseudoR2_Nakagawa_marginal = pseudoR2_Nakagawa[1, "R2m"],
        pseudoR2_Nakagawa_conditional = pseudoR2_Nakagawa[1, "R2c"],
        pseudoR2_Jaeger_sgv_marginal = pseudoR2_Jaeger[1, "Rsq"],
        rmsd = rmsd,
        cv_rmsd = cv_rmsd
      )
    }


    #--- Print/plot model evaluation (part 1)
    if (do_resid_plots) {
      get_resid_plots(
        m1 = res_fit[[1]][[id_eus]][["m"]],
        m2 = res_fit[[2]][[id_eus]][["m"]],
        m_names = model_names,
        data_sp = dats_q3d[["vals"]][, coords, sc_hist, default_subprj],
        data_ID = list(msdata[[1]][["data_ID"]], msdata[[2]][["data_ID"]]),
        subsets = list(msdata[[1]][["ids_data2"]], msdata[[2]][["ids_data2"]]),
        path = dir_res2_q3d,
        ftag = paste0("FutureDiffEnvRelationships_", id_eus)
      )
    }
  }


  #--- Print/plot model evaluation (part 2)
  get_model_background(
    res = res_fit,
    vdim1 = model_names,
    vdim2 = ensemble_units,
    path = dir_res2_q3d,
    ftag = "FutureDiffEnvRelationships"
  )
}

if (do_plot_relationships) {
  #--- Loop over RCPs x dTime
  for (ie in seq_along(ensemble_units)) {
    id_eus <- ensemble_units[ie]

    tmp_fit <- lapply(res_fit, function(x) x[[id_eus]])

    plot_model_relationships(
      res_fit = tmp_fit,
      predictors = list(
        msdata[[1]][["preds_selected"]],
        msdata[[2]][["preds_selected"]]
      ),
      responses = model_names,
      xlim_probs = c(0.005, 0.995),
      ylim = c(-1, 1),
      xlabs = list(
        msdata[[1]][["labs_preds_selected"]],
        msdata[[2]][["labs_preds_selected"]]
      ),
      ylabs = model_labels,
      include_0_on_xaxis = FALSE,
      panels_by_row = FALSE,
      path = dir_res2_q3d,
      ftag = paste0(
        "FutureDiffEnvRelationships_merged-ensemble_",
        id_eus, "_",
        vtag
      )
    )
  }
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q3d)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
