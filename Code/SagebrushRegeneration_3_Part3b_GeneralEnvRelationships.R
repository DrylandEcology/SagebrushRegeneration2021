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

pkgs <- c("ClustOfVar", "reshape2", "MuMIn", "r2glmm", "rSW2utils", "spaMM")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
# Research question:
# 3.	Can we understand model outcomes in terms of environmental site
#     characteristics like climate, soil moisture, and soils?
#   b.	Identify a parsimonious statistical relationship between climate and
#       soil moisture variables and model outcomes for the historical period
#       and site-specific soils?


#--- Specify variables
defs_q3b_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q3b_models, "analysis_group"]
models <- vars_table[defs_q3b_models, "name_original_rSFSW2"]
model_labels <- vars_table[defs_q3b_models, "label_manuscript"]

coords <- c("X_WGS84", "Y_WGS84")
defs_q3b_coords <- which(vars_table[, "name_original_rSFSW2"] %in% coords)

preds <- c(
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
defs_q3b_preds <- which(vars_table[, "name_original_rSFSW2"] %in% preds)
preds_labels <- vars_table[defs_q3b_preds, "label_manuscript"]

defs_q3b <- c(defs_q3b_models, defs_q3b_coords, defs_q3b_preds)


#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))


ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q3b.rds")

if (file.exists(ftmp)) {
  dats_q3b <- readRDS(ftmp)

} else {
  dats_q3b <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q3b,
      sim_scenarios = sc_hist,
      sim_subprojects = default_subprj
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q3b, file = ftmp)
}



#------------------------------------------------------------------------------#
#--- Analyse data

do_select_variables <- FALSE #TRUE #
do_dcor_estimates <- FALSE
do_dcor2_estimates <- TRUE
do_resid_plots <- FALSE
do_fit_relationships <- TRUE
do_plot_relationships <- TRUE



dir_res2_q3b <- file.path(
  dir_res2,
  "Output3_EnvRelationships",
  "Current_GeneralEnvRelationships"
)

dir.create(dir_res2_q3b, recursive = TRUE, showWarnings = FALSE)

legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)


for (nclusters_max_apriori in c(2L, 4L)[1]) {
  # nclusters_max_apriori: 2 and 4 give the same results because
  # hclust determines 2 clusters

  #--- Identify best preditors and prepare data
  msdata <- list()

  for (ks in seq_along(models)) {
    msdata[[ks]] <- list()

    msdata[[ks]][["tag_subset"]] <- names(def_subsets)[used_subsets[ks]]

    #--- Prepare data (part 1)
    msdata[[ks]][["ids_data"]] <- which(complete.cases(
      dats_q3b[["vals"]][, models[ks], sc_hist, default_subprj]
    ))

    msdata[[ks]][["tmp_data"]] <- as.data.frame(
      dats_q3b[["vals"]][msdata[[ks]][["ids_data"]], c(models[ks], preds), sc_hist, default_subprj]
    )


    msdata[[ks]][["data_sp"]] <-
      dats_q3b[["vals"]][msdata[[ks]][["ids_data"]], coords, sc_hist, default_subprj]



    #--- Check for need of variable transformations (per geographic subset)
    if (FALSE) {
      plot_correlation_checks(
        data = msdata[[ks]][["tmp_data"]],
        predictors = preds,
        responses = models[ks],
        path = dir_res2_q3b,
        ftag = paste0(
          "GeneralEnvRelationships_untransformed_",
          sc_hist, "_",
          msdata[[ks]][["tag_subset"]]
        )
      )
    }

    #--- Transform variables
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
        "log(/nyrs +Heat exposure [days])",
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
        msdata[[ks]][["tmp_data"]][, transform_predictors[["var_preds_toreplace"]][it]]
      )

      msdata[[ks]][["tmp_data"]][, transform_predictors[["var_preds_transformed"]][it]] <-
        transform_predictors[["fun_transform"]][[it]](
          msdata[[ks]][["tmp_data"]][, transform_predictors[["var_preds_toreplace"]][it]]
        )

      ids_notfinite2 <- !is.finite(
        msdata[[ks]][["tmp_data"]][, transform_predictors[["var_preds_transformed"]][it]]
      )

      # Stop if new non-finite values added
      stopifnot(sum(ids_notfinite2 & !ids_notfinite1) == 0)
    }

    if (TRUE) {
      # Remove original untransformed variables
      ids <- preds %in% transform_predictors[["var_preds_toreplace"]]
      preds_all <- c(preds[!ids], transform_predictors[["var_preds_transformed"]])
      var_labs_all <- c(preds_labels[!-ids], transform_predictors[["var_labs_transformed"]])
    } else {
      preds_all <- c(preds, transform_predictors[["var_preds_transformed"]])
      var_labs_all <- c(preds_labels, transform_predictors[["var_labs_transformed"]])
    }

    msdata[[ks]][["tmp_data"]] <- msdata[[ks]][["tmp_data"]][, c(models[ks], preds_all)]


    if (FALSE) {
      plot_correlation_checks(
        data = msdata[[ks]][["tmp_data"]],
        predictors = preds_all,
        responses = models[ks],
        path = dir_res2_q3b,
        ftag = paste0(
          "GeneralEnvRelationships_transformed_",
          sc_hist, "_", msdata[[ks]][["tag_subset"]]
        )
      )
    }


    #--- Select variables for easy interpretation and reduced multicollinearity
    if (do_select_variables) {
      if (FALSE) {
        # PCA
        isnotna <- complete.cases(msdata[[ks]][["datax"]])
        pca_preds <- prcomp(
          x = msdata[[ks]][["datax"]][isnotna, ],
          center = TRUE,
          scale. = TRUE
        )

        summary(pca_preds)
        print(pca_preds[["rotation"]][, 1:4])
        plot(pca_preds)
        biplot(pca_preds)
      }

      #------ Variable selection: best representative variable from each cluster

      #--- Determine bivariate unbiased distance correlation
      sel_preds_info <- calc_model_dcor2d(
        data = msdata[[ks]][["tmp_data"]],
        predictors = preds_all,
        responses = models[ks]
      )

      # Retain variables with at least x bivariate unbiased distance correlations to outcome
      preds_all_gtdcorU <- dimnames(sel_preds_info)[[1]][sel_preds_info[, , "dcorU"] > 0.0]

      #--- Ascendant hierarchical clustering of variables
      sel_hclust <- calc_asc_hierarch_clust(
        data = msdata[[ks]][["tmp_data"]],
        predictors = preds_all_gtdcorU,
        do_stability = TRUE,
        path = dir_res2_q3b,
        ftag = paste0(
          "GeneralEnvRelationships_",
          sc_hist, "_", msdata[[ks]][["tag_subset"]]
        )
      )


      #--- Cut tree into selected number of clusters
      sel_clusters <- ClustOfVar::cutreevar(
        obj = sel_hclust[["hclust"]],
        k = sel_hclust[["nclusters"]]
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
        data = msdata[[ks]][["tmp_data"]],
        nclusters = min(
          nclusters_max_apriori,
          sel_hclust[["nclusters"]],
          na.rm = TRUE
        ),
        limit_cor = 0.5
      )

      preds_selected <- tmp[["kpreds_selected"]]


      write.csv(
        reshape2::dcast(
          reshape2::melt(tmp[["predictor_info"]]),
          Var1 ~ Var2 + Var3
        ),
        row.names = FALSE,
        file = file.path(
          dir_res2_q3b,
          paste0(
            "Table_ClusterInformation_",
            sc_hist, "_", msdata[[ks]][["tag_subset"]],
            "_nclusters", nclusters_max_apriori,
            ".csv"
          )
        )
      )


      #--- Check for residual correlations
      plot_correlation_checks(
        data = msdata[[ks]][["tmp_data"]],
        predictors = preds_selected,
        responses = models[ks],
        path = dir_res2_q3b,
        ftag = paste0(
          "GeneralEnvRelationships_selected_",
          sc_hist, "_", msdata[[ks]][["tag_subset"]],
          "_nclusters", nclusters_max_apriori
        )
      )


    } else {
      if (nclusters_max_apriori == 2) {
        preds_selected <- if (ks == 1) {
          c(
            "MAT_C_mean",
            "log_inyrs_plus_TminBelow0CwithoutSpringSnow_days_mean"
          )
        } else if (ks == 2) {
          c(
            "MAT_C_doyRange1to250_mean",
            "periodicVWCmatricMean_FirstLayer_doyRange71to100_mean"
          )
        }

      } else if (nclusters_max_apriori == 4) {
        preds_selected <- if (ks == 1) {
          c(
            "MAT_C_mean",
            "log_inyrs_plus_TminBelow0CwithoutSpringSnow_days_mean"
          )
        } else if (ks == 2) {
          c(
            "MAT_C_doyRange1to250_mean",
            "periodicVWCmatricMean_FirstLayer_doyRange71to100_mean",
            "sqrt_ThermalSnowfreeDryPeriods_SWPcrit3000kPa_topLayers_DrySpellsAllLayers_maxDuration_days_mean",
            "log_inyrs_plus_TminBelow0CwithoutSpringSnow_days_mean"
          )
        }
      }
    }


    #--- Prepare data (part 2)
    msdata[[ks]][["tmp_data"]] <- msdata[[ks]][["tmp_data"]][, c(models[ks], preds_selected)]

    msdata[[ks]][["preds_selected"]] <- preds_selected
    ids <- match(preds_selected, preds_all, nomatch = 0)
    msdata[[ks]][["labs_preds_selected"]] <- var_labs_all[ids]
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

    write.csv(
      res_dcor,
      row.names = FALSE,
      file = file.path(
        dir_res2_q3b,
        paste0(
          "Table_DistanceCorrelations_",
          sc_hist, "_", ,
          "_nclusters", nclusters_max_apriori,
          ".csv"
        )
      )
    )
  }


  #--- Model fitting GLM
  if (do_fit_relationships) {
    res_fit <- list()

    #--- Loop over models
    for (iy in seq_along(models)) {

      #--- Prepare data (part 2)
      tmpd <- msdata[[iy]][["tmp_data"]][, c(models[iy], msdata[[iy]][["preds_selected"]])]
      ids_use2 <- complete.cases(tmpd)

      #--- Prepare data (part 3)
      tmpy <- msdata[[iy]][["tmp_data"]][ids_use2, models[iy]]
      tmpy2 <- data.frame(
        npos = as.integer(round(tmpy * nyrs)),
        nneg = as.integer(round((1 - tmpy) * nyrs))
      )

      datax_scaled <- scale(msdata[[iy]][["tmp_data"]][ids_use2, msdata[[iy]][["preds_selected"]]])
      datax_scale <- attr(datax_scaled, "scaled:scale")
      datax_center <- attr(datax_scaled, "scaled:center")

      m_data <- data.frame(y = tmpy, tmpy2, datax_scaled)

      msdata[[iy]][["ids_data2"]] <- which(ids_use2)



      #--- Fit models
      # Model with quadratic terms and 2-way interaction
      y_formula1 <- stats::as.formula(paste(
        "cbind(npos, nneg) ~",
        "(", paste0(msdata[[iy]][["preds_selected"]], collapse = "+"), ") ^ 2",
        "+",
        paste0("I(", msdata[[iy]][["preds_selected"]], "^2)", collapse = "+")
      ))

      m1 <- stats::glm(
        y_formula1,
        family = stats::binomial(),
        data = m_data
      )

      # McFadden's R squared
      nullmod <- stats::glm(
        cbind(npos, nneg) ~ 1,
        family = stats::binomial(),
        data = m_data
      )

      pseudoR2_McFadden <- as.numeric(1 - stats::logLik(m1) / stats::logLik(nullmod))

      # Nakagawa's R-squared
      pseudoR2_Nakagawa <- MuMIn::r.squaredGLMM(m1)

      # Jaeger's R-squared
      pseudoR2_Jaeger <- try(r2glmm::r2beta(m1, partial = FALSE))
      if (inherits(pseudoR2_Jaeger, "try-error")) {
        pseudoR2_Jaeger <- array(dim = c(1, 1), dimnames = list(NULL, "Rsq"))
      }

      # Root-square mean deviations (in-sample)
      rmsd <- rSW2utils::rmse(tmpy, fitted(m1), na.rm = TRUE)

      cv_rmsd <- rmsd / mean(tmpy, na.rm = TRUE)


      # Model with added spatial correlation structure
      if (FALSE) {
        library("spaMM")
        tmp <- as.character(y_formula1)
        y_formula_sp <- as.formula(paste(
          tmp[2], tmp[1], tmp[3], "+ Matern(1 | X_WGS84 + Y_WGS84)"
        ))

        m_sp <- spaMM::fitme(
          y_formula3,
          family = binomial(),
          data = data.frame(
            m_data,
            msdata[[iy]][["data_sp"]]
          )
        )

        # Compare model_names
        anova(m1, m_sp, test = "Chisq")
        AIC(m1, m_sp)
      }


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
      res_fit[[model_names[iy]]] <- list(
        m = m1,
        datay = m_data[, "y"],
        datax_scaled = datax_scaled,
        datax_scale = datax_scale,
        datax_center = datax_center,
        coef_standard = summary(m1)[["coefficients"]],
        pseudoR2_McFadden = pseudoR2_McFadden,
        pseudoR2_Nakagawa_marginal = pseudoR2_Nakagawa["delta", "R2m"],
        pseudoR2_Nakagawa_conditional = pseudoR2_Nakagawa["delta", "R2c"],
        pseudoR2_Jaeger_sgv_marginal = pseudoR2_Jaeger[1, "Rsq"],
        rmsd = rmsd,
        cv_rmsd = cv_rmsd
      )
    }


    #--- Print/plot model evaluation
    if (do_resid_plots) {
      get_resid_plots(
        m1 = res_fit[[1]][["m"]],
        m2 = res_fit[[2]][["m"]],
        m_names = model_names,
        data_sp = dats_q3b[["vals"]][, coords, sc_hist, default_subprj],
        subsets = list(msdata[[1]][["ids_data2"]], msdata[[2]][["ids_data2"]]),
        path = dir_res2_q3b,
        ftag = paste0(
          "GeneralEnvRelationships",
          "_nclusters", nclusters_max_apriori
        )
      )
    }

    get_model_background(
      res = res_fit,
      vdim1 = model_names,
      path = dir_res2_q3b,
      ftag = paste0(
        "GeneralEnvRelationships",
        "_nclusters", nclusters_max_apriori
      )
    )
  }


  if (do_plot_relationships) {
    plot_model_relationships(
      res_fit = res_fit,
      predictors = list(
        msdata[[1]][["preds_selected"]],
        msdata[[2]][["preds_selected"]]
      ),
      responses = model_names,
      ylim = c(0, 1),
      xlabs = list(
        msdata[[1]][["labs_preds_selected"]],
        msdata[[2]][["labs_preds_selected"]]
      ),
      ylabs = model_labels,
      include_0_on_xaxis = FALSE,
      panels_by_row = nclusters_max_apriori == 2,
      path = dir_res2_q3b,
      ftag = paste0(
        "GeneralEnvRelationships_",
        sc_hist,
        "_nclusters", nclusters_max_apriori
      )
    )
  }
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q3b)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
