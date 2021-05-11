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

pkgs <- c("gbm", "mgcv", "MuMIn", "r2glmm", "rSW2utils", "spaMM")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
# Research question:
# 3.	Can we understand model outcomes in terms of environmental site
#     characteristics like climate, soil moisture, and soils?
#   a.	How do model outcomes (for the historical period and site-specific
#       soils) vary along the two predictor variables used by the
#       Shriver et al. model?

#--- Specify variables
defs_q3a_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q3a_models, "analysis_group"]
models <- vars_table[defs_q3a_models, "name_original_rSFSW2"]
model_labels <- vars_table[defs_q3a_models, "label_manuscript"]

coords <- c("X_WGS84", "Y_WGS84")
defs_q3a_coords <- which(vars_table[, "name_original_rSFSW2"] %in% coords)

preds_selected <- c(
  "MAT_C_doyRange1to250_mean",
  "periodicVWCmatricMean_FirstLayer_doyRange71to100_mean"
)
defs_q3a_preds <- which(
  vars_table[, "name_original_rSFSW2"] %in% preds_selected
)
preds_labels <- vars_table[defs_q3a_preds, "label_manuscript"]

defs_q3a <- c(defs_q3a_models, defs_q3a_coords, defs_q3a_preds)


#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q3a.rds")

if (file.exists(ftmp)) {
  dats_q3a <- readRDS(ftmp)

} else {
  dats_q3a <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q3a,
      sim_scenarios = sc_hist,
      sim_subprojects = default_subprj
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q3a, file = ftmp)
}




#------------------------------------------------------------------------------#
#--- Analyse data

do_dcor_estimates <- FALSE
do_fit_relationships <- TRUE
do_resid_plots <- FALSE
do_plot_relationships <- TRUE
do_plot_data <- TRUE




dir_res2_q3a <- file.path(
  dir_res2,
  "Output3_EnvRelationships",
  "Current_EnvRelationships_basedon_Shriver2018"
)

dir.create(dir_res2_q3a, recursive = TRUE, showWarnings = FALSE)

legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)


# Data
tmp_data <- dats_q3a[["vals"]][, c(coords, models, preds_selected), sc_hist, default_subprj]


#--- (Partial) distance correlations
if (do_dcor_estimates) {
  id_geo_joint <- names(def_subsets)[max(used_subsets)]

  res_dcor <- get_dcor_tests(
    data = tmp_data[def_subsets[[id_geo_joint]], ],
    models = models,
    model_names = model_names,
    preds_all = NULL,
    preds_selected = preds_selected,
    R = 49, # R = 99 will get killed
    path = dir_res2_q3a,
    ftag = paste0(sc_hist, "_", id_geo_joint),
    verbose = TRUE
  )
}



#--- Fit to predictors of Shriver2018-model
do_gam <- FALSE
do_glm <- !do_gam
do_m1_withInteractions <- FALSE

if (do_fit_relationships) {
  #--- Prepare data

  # for (do_glm in c(TRUE, FALSE)) {
  #   do_gam <- !do_glm
  #
  #   for (do_m1_withInteractions in c(TRUE, FALSE)) {

  #--- Model fitting (GLM or GAM)
  res_fit <- list()

  #------ Loop over models
  for (iy in seq_along(models)) {

    tmpy <- tmp_data[, models[iy]]
    ids_data <- complete.cases(tmpy)
    tmpy <- tmpy[ids_data]

    tmpy2 <- data.frame(
      npos = as.integer(round(tmpy * nyrs)),
      nneg = as.integer(round((1 - tmpy) * nyrs))
    )

    datax_unscaled <- tmp_data[ids_data, preds_selected]
    datax_scaled <- scale(datax_unscaled)
    datax_scale <- attr(datax_scaled, "scaled:scale")
    datax_center <- attr(datax_scaled, "scaled:center")

    data_sp <- tmp_data[ids_data, coords]
    m_data <- data.frame(tmpy2, datax_scaled)


    if (FALSE && iy == 2) {
      # recover Shriver2018
      predy <- p_Shriver2018(
        VWC_spring = datax_unscaled[, "periodicVWCmatricMean_FirstLayer_doyRange71to100_mean"],
        Temp_mean = datax_unscaled[, "MAT_C_doyRange1to250_mean"]
      )

      stopifnot(identical(predy, tmpy)) ## check that we use correct values

      m1rec <- stats::glm(
        formula = as.formula(paste(
          "cbind(npos, nneg) ~",
          paste0(preds_selected, collapse = "+")
        )),
        family = stats::binomial(),
        data = data.frame(tmpy2, datax_unscaled)
      )

      summary(m1rec)[["coefficients"]][, c("Estimate", "Std. Error")]
      confint(m1rec)
      # Mean and 95% confidence intervals
      # Intercept = 3.2970 (3.1889-3.4054) versus 3.306
      # VWC_spring = 2.5058 (2.3063-2.7052) versus 2.499
      # Temp_mean = -0.2882 (-0.2962--0.2803) versus = -0.289
    }


    if (FALSE) {
      library("gbm")

      data_gbm <- data.frame(y = tmpy, datax_unscaled)

      m1 <- gbm::gbm(
        as.formula(paste("y ~ ", paste(preds_selected, collapse = "+"))),
        data = data_gbm
      )
    }

    #--- Fit model
    if (FALSE) {
      mtag <- digest::digest(list(datax_scaled, tmpy2, data_sp), algo = "spookyhash")
    }

    if (do_gam) {
      # GAM with tensor shrinkage-spline predictors
      library("mgcv")

      if (do_m1_withInteractions) {
        # interactions are highly concurve
        k1 <- c(100, 10)
        k2 <- 20
        if (FALSE) {
          y_formula1 <- as.formula(paste(
            "cbind(npos, nneg) ~",
            paste0("s(", preds_selected, ", bs = 'cs', k =", k1, ")", collapse = "+"),
            "+",
            paste0("ti(", preds_selected[1], ", ", preds_selected[2], ", bs = 'cs', k =", k2, ")")
          ))
        }
        y_formula1 <- as.formula(paste(
          "cbind(npos, nneg) ~",
          paste0("te(", preds_selected[1], ", ", preds_selected[2], ", bs = 'cs', k =", k2, ")")
        ))

      } else {
        k1 <- if (iy == 1) c(150, 100) else c(50, 10)
        y_formula1 <- as.formula(paste(
          "cbind(npos, nneg) ~",
          paste0("s(", preds_selected, ", bs = 'cs', k =", k1, ")", collapse = "+")
        ))
      }

      m1 <- mgcv::bam(
        y_formula1,
        family = binomial(),
        data = m_data,
        discrete = TRUE
      )

      # Check model
      if (FALSE) {
        # https://stats.stackexchange.com/questions/179061/gam-smoother-vs-parametric-term-concurvity-difference
        # https://stats.stackexchange.com/questions/38093/how-to-deal-with-high-correlation-among-predictors-in-multiple-regression?rq=1
        # https://stats.stackexchange.com/questions/299368/gams-with-many-slightly-correlated-predictors
        check_m1_concurv <- mgcv::concurvity(m1)
        check_m1 <- mgcv::gam.check(m1)
      }

      # McFadden's R squared
      nullmod <- mgcv::gam(
        cbind(npos, nneg) ~ 1,
        family = stats::binomial(),
        data = m_data
      )

      m1_R2 <- as.numeric(1 - stats::logLik(m1) / stats::logLik(nullmod))
    }

    if (do_glm) {
      if (do_m1_withInteractions) {
        # Model with added quadratic terms and interaction
        y_formula1 <- stats::as.formula(paste(
          "cbind(npos, nneg) ~",
          paste0(preds_selected, collapse = "*"),
          "+",
          paste0("I(", preds_selected, "^2)", collapse = "*")
        ))
      } else {
        # Model from Shriver et al. 2018 GCB
        y_formula1 <- as.formula(paste(
          "cbind(npos, nneg) ~",
          paste0(preds_selected, collapse = "+")
        ))
      }

      m1 <- stats::glm(
        y_formula1,
        family = stats::binomial(),
        data = m_data
      )

      m1un <- stats::glm(
        y_formula1,
        family = stats::binomial(),
        data = data.frame(tmpy2, datax_unscaled)
      )

      # McFadden's R squared
      nullmod <- stats::glm(
        cbind(npos, nneg) ~ 1,
        family = stats::binomial(),
        data = m_data
      )

      pseudoR2_McFadden <- as.numeric(
        1 - stats::logLik(m1) / stats::logLik(nullmod)
      )


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
          data = data.frame(m_data, data_sp)
        )

        # Compare model_names
        anova(m1, m_sp, test = "Chisq")
        AIC(m1, m_sp)
      }
    }

    # Check models
    if (FALSE) {
      # Walk through code of `check_models` interactively
      m1_checks <- check_models(
        model = m1,
        preds = preds_selected,
        data = m_data,
        datay = tmp_data[ids_data, models[iy]],
        datax_scaled = datax_scaled
      )
    }

    # Save model for later
    res_fit[[model_names[iy]]] <- list(
      m = m1,
      datay = tmp_data[ids_data, models[iy]],
      datax_scaled = datax_scaled,
      datax_scale = datax_scale,
      datax_center = datax_center,
      ids_data = ids_data,
      coef_standard = summary(m1)[[if (do_glm) "coefficients" else "s.table"]],
      coef_unstandard = if (do_glm) summary(m1un)[["coefficients"]],
      pseudoR2_McFadden = pseudoR2_McFadden,
      pseudoR2_Nakagawa_marginal = pseudoR2_Nakagawa["delta", "R2m"],
      pseudoR2_Nakagawa_conditional = pseudoR2_Nakagawa["delta", "R2c"],
      pseudoR2_Jaeger_sgv_marginal = pseudoR2_Jaeger[1, "Rsq"],
      rmsd = rmsd,
      cv_rmsd = cv_rmsd
    )
  }


  #------ Print/plot model evaluation
  if (do_resid_plots) {
    get_resid_plots(
      m1 = res_fit[[1]][["m"]],
      m2 = res_fit[[2]][["m"]],
      m_names = model_names,
      data_sp = tmp_data[, coords],
      subsets = list(res_fit[[1]][["ids_data"]], res_fit[[2]][["ids_data"]]),
      path = dir_res2_q3a,
      ftag = "Shriver2018-predictors",
      device = device_type
    )
  }

  get_model_background(
    res = res_fit,
    vdim1 = model_names,
    path = dir_res2_q3a,
    ftag = "Shriver2018-predictors"
  )

  #   } # end of loop over `do_m1_withInteractions`
  # } # end of loop over `do_glm`
}


if (do_plot_relationships) {
  plot_model_relationships(
    res_fit = res_fit,
    predictors = preds_selected,
    responses = model_names,
    ylim = c(0, 1),
    xlabs = preds_labels,
    ylabs = model_labels,
    include_0_on_xaxis = TRUE,
    panels_by_row = TRUE,
    # Add Shriver2018 model
    last_plot = expression({
      if (response == "Shriver2018") {
        vm <- obj_visreg[["meta"]][["x"]]
        xr1 <- obj_visreg[["fit"]][[vm]]
        xr1 <- xr1[seq_len(length(xr1) / 3)]

        xm2 <- funtrans(iv2)(vbreaks[2])

        ymod <- if (grepl("VWC", vm)) {
          p_Shriver2018(VWC_spring = xr1, Temp_mean = xm2)
        } else if (grepl("MAT", vm)) {
          p_Shriver2018(VWC_spring = xm2, Temp_mean = xr1)
        }

        lines(xr1, ymod, col = "black")
      }
    }),
    path = dir_res2_q3a,
    ftag = paste0(
      "Shriver2018-predictors_",
      sc_hist, "_",
      if (do_glm) "_GLM" else "_GAM",
      if (do_m1_withInteractions) "_withInteractions"
    ),
    device = device_type
  )
}


#--- Data curves for predictors of Shriver2018-model
if (do_plot_data) {
  n_panels <- c(length(model_names), length(preds_selected))
  template <- array(list(), dim = n_panels)

  xlim_mtx <- ylim_mtx <- type_mtx <- template
  for (k2 in seq_len(n_panels[2])) {
    xlim_mtx[, k2] <- list(range(tmp_data[, preds_selected[k2]], na.rm = TRUE))
  }
  ylim_mtx[] <- list(c(0, 1))
  type_mtx[] <- "smoothScatter.vals"

  axlabs_mtx <- data_mtx <- template
  for (k1 in seq_len(n_panels[1])) for (k2 in seq_len(n_panels[2])) {
    axlabs_mtx[k1, k2][[1]] <- c(preds_labels[k2], model_labels[k1])
    data_mtx[k1, k2][[1]] <- cbind(
     x = dats_q3a[["vals"]][, preds_selected[k2], sc_hist, default_subprj],
     y = dats_q3a[["vals"]][, models[k1], sc_hist, default_subprj]
    )
  }

  rSW2analysis::plot_matrix_of_panels(
    n_panels = n_panels,
    type_matrix = type_mtx,
    add_1to1 = FALSE,
    data_matrix = data_mtx,
    subset = NULL,
    xlim_matrix = xlim_mtx,
    ylim_matrix = ylim_mtx,
    label_axis_matrix = axlabs_mtx,
    use_labels = "none",
    fexp_axis = 1,
    path = dir_res2_q3a,
    device = device_type,
    ftag = paste0("SmoothScatter_", sc_hist, "-", default_subprj)
  )
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q3a)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
