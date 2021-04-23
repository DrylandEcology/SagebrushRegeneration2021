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

pkgs <- c("reshape2", "lme4", "MuMIn", "rSW2utils", "spaMM")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
# Research question:
# 3.	Can we understand model outcomes in terms of environmental site
#     characteristics like climate, soil moisture, and soils?
#   b.	Identify a parsimonious statistical relationship between climate and
#       soil moisture variables and model outcomes for the historical period
#       and site-specific soils?


#--- Specify variables
defs_q3c_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q3c_models, "analysis_group"]
models <- vars_table[defs_q3c_models, "name_original_rSFSW2"]
model_labels <- vars_table[defs_q3c_models, "label_manuscript"]

coords <- c("X_WGS84", "Y_WGS84")
defs_q3c_coords <- which(vars_table[, "name_original_rSFSW2"] %in% coords)

preds <- c(
  "SWinput_Soil_topLayers_Sand_fraction",
  "SWinput_Soil_topLayers_Clay_fraction"
)
defs_q3c_preds <- which(vars_table[, "name_original_rSFSW2"] %in% preds)
preds_labels <- vars_table[defs_q3c_preds, "label_manuscript"]

defs_q3c <- c(defs_q3c_models, defs_q3c_coords, defs_q3c_preds)




#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))


ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q3c.rds")

if (file.exists(ftmp)) {
  dats_q3c <- readRDS(ftmp)

} else {
  dats_q3c <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q3c,
      sim_scenarios = sc_hist,
      sim_subprojects = tag_subprj[-1]
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q3c, file = ftmp)
}





#--- Calculate differences among fixed soil types
tmp <- attr(terms(as.formula(paste(
  "y ~ (",
  paste0(tag_subprj[-1], collapse = "+"), ")^2"))),
  "term.labels"
)
tmp <- strsplit(tmp, split = ":")
ids <- lengths(tmp) == 2
fxsoils_diffgrid <- tmp[ids]

# Data container
tmp_dim <- dim(dats_q3c[["vals"]])
tmp_names <- dimnames(dats_q3c[["vals"]])

dats_q3c[["delta_fixedSoils"]] <- array(
  NA,
  dim = c(tmp_dim[1:3], length(fxsoils_diffgrid)),
  dimnames = c(
    tmp_names[1:3],
    list(sapply(fxsoils_diffgrid, function(x) paste0(x[2], "-", x[1])))
  )
)


# Calculate differences
for (k in seq_along(fxsoils_diffgrid)) {
  fxsoil_ref <- fxsoils_diffgrid[[k]][1]
  fxsoil_target <- fxsoils_diffgrid[[k]][2]

  dats_q3c[["delta_fixedSoils"]][, , , k] <-
    dats_q3c[["vals"]][, , , fxsoils_diffgrid[[k]][2], drop = FALSE] -
    dats_q3c[["vals"]][, , , fxsoils_diffgrid[[k]][1], drop = FALSE]
}

models_diffs <- paste0("diff_", models)
preds_selected <- paste0("diff_", preds)
tmp <- paste("delta", preds_labels)
labs_preds_selected <- sub("Shallow ", "", tmp)



#------------------------------------------------------------------------------#
#--- Analyse data

do_dcor_estimates <- FALSE
do_resid_plots <- FALSE
do_fit_relationships <- TRUE
do_plot_relationships <- TRUE
do_plot_data <- TRUE


#--------- Fit model deltas to sand/clay differences among fixed soil types


dir_res2_q3c <- file.path(
  dir_res2,
  "Output3_EnvRelationships",
  "Current_SoilDifferenceRelationships"
)

dir.create(dir_res2_q3c, recursive = TRUE, showWarnings = FALSE)

legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)


#--- Prepare data (part 1)
msdata <- list()

for (ks in seq_along(models)) {
  msdata[[ks]] <- list()
  msdata[[ks]][["tag_subset"]] <- names(def_subsets)[used_subsets[ks]]


  msdata[[ks]][["ids_data"]] <- !apply(
    dats_q3c[["delta_fixedSoils"]][, c(models[ks], preds), , ],
    1,
    anyNA
  )

  tmp <- reshape2::melt(
    dats_q3c[["delta_fixedSoils"]][msdata[[ks]][["ids_data"]], c(models[ks], preds), , ]
  )
  tmp[, "Var2"] <- paste0("diff_", tmp[, "Var2"])


  tmp_soiltype <- do.call(rbind, strsplit(as.character(tmp[, "Var3"]), split = "-"))
  tmp[, "Var4"] <- tmp_soiltype[, 2] # reference fixed soil
  tmp[, "Var5"] <- tmp_soiltype[, 1] # target fixed soil for diff = target - ref

  tmp_data <- reshape2::dcast(tmp, Var1 + Var4 + Var5 ~ Var2)
  colnames(tmp_data)[1:3] <- c("Gridcell", "SoilPrj_ref", "SoilPrj_target")
  tmp_data[, "SoilDiff"] <- apply(
    tmp_data[, c("SoilPrj_ref", "SoilPrj_target")],
    1,
    paste, collapse = "-"
  )

  msdata[[ks]][["tmp_data"]] <- tmp_data

  msdata[[ks]][["datax_unscaled"]] <- msdata[[ks]][["tmp_data"]][, preds_selected]
  msdata[[ks]][["datax_scaled"]] <- scale(msdata[[ks]][["datax_unscaled"]])
  msdata[[ks]][["datax_scale"]] <- attr(msdata[[ks]][["datax_scaled"]], "scaled:scale")
  msdata[[ks]][["datax_center"]] <- attr(msdata[[ks]][["datax_scaled"]], "scaled:center")

  msdata[[ks]][["data_sp"]] <- dats_q3c[["vals"]][msdata[[ks]][["ids_data"]], coords, sc_hist, tag_subprj[2]]
}


#--- (Partial) distance correlations
if (do_dcor_estimates) {
  res_dcor <- get_dcor_tests(
    data = tmp_data,
    models = models_diffs,
    model_names = model_names,
    preds_all = NULL,
    preds_selected = preds_selected,
    R = 49, # R = 99 will get killed
    path = dir_res2_q3c,
    ftag = paste0(sc_hist, "_", tag_subset),
    verbose = TRUE
  )
}


#--- Model fitting LMM
if (do_fit_relationships) {
  library("lme4")

  res_fit <- list()

  #--- Loop over models
  for (iy in seq_along(models_diffs)) {

    m_data <- data.frame(
      y = msdata[[iy]][["tmp_data"]][, models_diffs[iy]],
      fSoilDiff = as.factor(msdata[[iy]][["tmp_data"]][, "SoilDiff"]),
      msdata[[iy]][["datax_scaled"]]
    )

    #--- Fit models

    # Model with quadratic terms and 2-way interaction
    y_formula1 <- stats::as.formula(paste(
      "y ~",
      "(", paste0(preds_selected, collapse = "+"), ") ^ 2",
      "+",
      paste0("I(", preds_selected, "^2)", collapse = "+"),
      "+ (1 | fSoilDiff)"
    ))

    m1 <- lme4::lmer(y_formula1, data = m_data)

    # McFadden's R squared
    nullmod <- lme4::lmer(
      y ~ 1 + (1 | fSoilDiff),
      data = m_data
    )

    pseudoR2_McFadden <- as.numeric(1 - stats::logLik(m1) / stats::logLik(nullmod))

    # Nakagawa's R-squared
    pseudoR2_Nakagawa <- MuMIn::r.squaredGLMM(m1)

    # Jaeger's R-squared
    pseudoR2_Jaeger <- array(dim = c(1, 1), dimnames = list(NULL, "Rsq"))  #r2glmm::r2beta(m1, partial = FALSE)

    # Root-square mean deviations (in-sample)
    rmsd <- rSW2utils::rmse(m_data[, "y"], fitted(m1), na.rm = TRUE)

    cv_rmsd <- rmsd / mean(m_data[, "y"], na.rm = TRUE)


    # Model with added spatial correlation structure
    if (FALSE) {
      library("spaMM")
      tmp <- as.character(y_formula1)
      y_formula_sp <- as.formula(paste(
        tmp[2], tmp[1], tmp[3], "+ Matern(1 | X_WGS84 + Y_WGS84)"
      ))

      m_sp <- spaMM::fitme(
        y_formula3,
        family = gaussian(),
        data = data.frame(m_data, msdata[[iy]][["data_sp"]])
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
        preds = preds_selected,
        data = m_data,
        datay = m_data[, "y"],
        datax_scaled = msdata[[iy]][["datax_scaled"]]
      )
    }

    # Save model for later
    res_fit[[model_names[iy]]] <- list(
      m = m1,
      datay = m_data[, "y"],
      datax_scaled = msdata[[iy]][["datax_scaled"]],
      datax_scale = msdata[[iy]][["datax_scale"]],
      datax_center = msdata[[iy]][["datax_center"]],
      coef_standard = summary(m1)[["coefficients"]],
      pseudoR2_McFadden = pseudoR2_McFadden,
      pseudoR2_Nakagawa_marginal = pseudoR2_Nakagawa[1, "R2m"],
      pseudoR2_Nakagawa_conditional = pseudoR2_Nakagawa[1, "R2c"],
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
      data_sp = dats_q3c[["vals"]][, coords, sc_hist, tag_subprj[2]],
      data_ID = list(
        msdata[[1]][["tmp_data"]][, "Gridcell"],
        msdata[[2]][["tmp_data"]][, "Gridcell"]
      ),
      subsets = list(msdata[[1]][["ids_data"]], msdata[[2]][["ids_data"]]),
      path = dir_res2_q3c,
      ftag = "FixedSoilDifferences"
    )
  }

  get_model_background(
    res = res_fit,
    vdim1 = model_names,
    path = dir_res2_q3c,
    ftag = "FixedSoilDifferences"
  )
}


if (do_plot_relationships) {
  plot_model_relationships(
    res_fit = res_fit,
    predictors = preds_selected,
    responses = model_names,
    ylim = range(
      msdata[[1]][["tmp_data"]][, models_diffs[1]],
      msdata[[2]][["tmp_data"]][, models_diffs[2]],
      na.rm = TRUE
    ),
    xlabs = labs_preds_selected,
    ylabs = model_labels,
    include_0_on_xaxis = FALSE,
    panels_by_row = TRUE,
    path = dir_res2_q3c,
    ftag = paste0("FixedSoilDifferences_", sc_hist)
  )
}


#------- Plot maps
if (do_plot_data) {
  # Prepare plot data
  n_panels <- c(length(fxsoils_diffgrid), length(models))
  template <- array(list(), dim = n_panels)

  zlim_mtx <- type_mtx <- label_mtx <- data_mtx <- template

  type_mtx[] <- "map.delta"

  zlim_tmp <- 0
  ids_geo <- def_subsets[["BigSage"]]

  # Prepare data
  for (k1 in seq_along(fxsoils_diffgrid)) for (k2 in seq_along(models)) {
    fxsoil_ref <- fxsoils_diffgrid[[k1]][1]
    fxsoil_target <- fxsoils_diffgrid[[k1]][2]

    sim_delta <- dats_q3c[["delta_fixedSoils"]][, models[k2], sc_hist, k1]

    zlim_tmp <- range(zlim_tmp, sim_delta[ids_geo], na.rm = TRUE)
    data_mtx[k1, k2][[1]] <- sim_delta
    label_mtx[k1, k2] <- paste0(model_names[k2], ": ",
      fxsoil_target, "-", fxsoil_ref)
  }

  zlim_mtx[] <- list(c(-1, 1) * max(abs(zlim_tmp), na.rm = TRUE))


  rSW2analysis::plot_matrix_of_panels(
    n_panels = n_panels,
    data_matrix = data_mtx,
    meta = SFSW2_prj_meta,
    subset = ids_geo,
    zlim_matrix = zlim_mtx,
    label_title_matrix = label_mtx,
    use_labels = "panel_identifier",
    type_matrix = type_mtx,
    map_legend_pos = legend_pos,
    map_extent = extent_subsets[[tag_subset]],
    path = dir_res2_q3c,
    ftag = paste(tag_subset, sc_hist, "diffsAmongFxSoils", sep = "_"),
    pborders = pborders
  )
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q3c)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
