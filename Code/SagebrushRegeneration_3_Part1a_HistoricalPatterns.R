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


pkgs <- c(
  "MuMIn", "r2glmm", "rSW2utils", "energy", "mgcv", "perturb",
  "car", "rms", "DHARMa", "visreg"
)
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#
#--- Specify variables
defs_q1a <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)

vars <- vars_table[defs_q1a, "name_original_rSFSW2"]
var_pGISSM <- vars[1]
var_pShriver2018 <- vars[2]


#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q1a.rds")

if (file.exists(ftmp)) {
  dats_q1a <- readRDS(ftmp)

} else {
  dats_q1a <- load_rSFSW2_BSR_data_from_netCDFs(
    path = dir_prj,
    locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    vars_table = vars_table,
    ids_vars = defs_q1a,
    sim_scenarios = sc_hist,
    sim_subprojects = default_subprj
  )

  dats_q1a <- drop(dats_q1a)

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q1a, file = ftmp)
}


#------------------------------------------------------------------------------#
#--- Analyse data

#-----------
#' Check that \code{\link{p_Shriver2018}} works correctly
#' Compare against Fig. 3ab of Shriver et al. 2018 GCB
if (FALSE) {
  dir_temp <- file.path(dir_res2, "Output2_ModelComparison")
  dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

  png(height = 5, width = 10, units = "in", res = 150,
    file = file.path(dir_temp, "CheckUse_Shriver2018.png"))

  par_prev <- par(mfrow = c(1, 2), mar = c(2, 3, 1, 3), mgp = c(1, 0, 0),
    tcl = 0.3, cex = 1)

  for (var in c("VWC_spring", "Temp_mean")) {
    if (var == "VWC_spring") {
      xs <- seq(0, 0.5, length.out = 100)
      ys <- p_Shriver2018(xs, 10)
    } else {
      xs <- seq(4, 15, length.out = 100)
      ys <- p_Shriver2018(0.15, xs)
    }

    plot(xs, ys, type = "l", ylim = c(0, 1),
      xlab = var, ylab = "p(Shriver2018)")

  }

  par(par_prev)
  dev.off()
}



#--- Make comparison figures
do_correlations <- TRUE
do_plot_comparison <- TRUE
do_regress_models <- TRUE
do_ecoregion_table <- TRUE


dir_res2_q1a <- file.path(
  dir_res2,
  "Output1_GeographicPatterns",
  "HistoricalPatterns"
)




if (do_correlations) {

  dir_tmp <- file.path(dir_res_data, "tmp_hist_cortests")
  dir.create(dir_tmp, recursive = TRUE, showWarnings = FALSE)

  fname_cortests <- file.path(dir_tmp, "bdcors.rds")

  if (file.exists(fname_cortests)) {
    tests_bdcor <- readRDS(fname_cortests)

  } else {
    res_lincor <- tests_bdcor <- list()

    # Linear correlation
    if (FALSE) {
      res_lincor[["BigSage"]] <- cor.test(
        dats_q1a[def_subsets[["BigSage"]], 1],
        dats_q1a[def_subsets[["BigSage"]], 2],
        method = "pearson"
      )
      # t = -38.575, df = 23200, p-value < 2.2e-16
      # 95 percent confidence interval: -0.2575630 -0.2333792
      # sample estimates: -0.2455093
    }

    res_lincor[["GreatBasin"]] <- cor.test(
      dats_q1a[def_subsets[["GreatBasin"]], 1],
      dats_q1a[def_subsets[["GreatBasin"]], 2],
      method = "pearson"
    )
    # t = -15.981, df = 5035, p-value < 2.2e-16
    # 95 percent confidence interval: -0.2458384 -0.1932682
    # sample estimates: -0.2197128

    # Brownian distance correlation:
    if (FALSE) {
      tests_bdcor[["BigSage"]] <- energy::dcor.test(
        x = dats_q1a[def_subsets[["BigSage"]], 1],
        y = dats_q1a[def_subsets[["BigSage"]], 2],
        R = 199
      )
      #	dCor independence test (permutation test)
      #
      #data:  index 1, replicates 199
      #dCor = 0.40281, p-value = 0.005
      #sample estimates:
      #      dCov       dCor    dVar(X)    dVar(Y)
      #0.05278047 0.40281206 0.12952762 0.13254983
    }

    tests_bdcor[["GreatBasin"]] <- energy::dcor.test(
      x = dats_q1a[def_subsets[["GreatBasin"]], 1],
      y = dats_q1a[def_subsets[["GreatBasin"]], 2],
      R = 199
    )
    #	dCor independence test (permutation test)
    #
    #data:  index 1, replicates 199
    #dCor = 0.29586, p-value = 0.005
    #sample estimates:
    #      dCov       dCor    dVar(X)    dVar(Y)
    #0.02976162 0.29585840 0.10220587 0.09900783

    saveRDS(tests_bdcor, file = fname_cortests)
  }
}

if (do_plot_comparison) {
  # Historical prediction map and scatterplot between GISSM and Shriver2018
  # for ks-subset and GreatBasin-subset
  n_panels <- c(2, 2)
  template <- array(list(), dim = n_panels)

  legend_pos <- switch(
    names(def_subsets)[min(used_subsets)],
    BigSage = "right",
    "left"
  )

  # plot info
  addfun_mtx <- data_mtx <- zlim_mtx <- vlim_mtx <- type_mtx <- title_mtx <-
    axlabs_mtx <- template

  # map info
  zlim_mtx[, 1] <- list(c(0, 1))
  type_mtx[, 1] <- "map.vals"
  title_mtx[1, 1] <- list(vars_table[defs_q1a[1], "label_manuscript"])
  title_mtx[2, 1] <- list(vars_table[defs_q1a[2], "label_manuscript"])

  data_mtx[1, 1][[1]] <- dats_q1a[, 1]
  data_mtx[2, 1][[1]] <- dats_q1a[, 2]


  # Outlines of ecoregions and study extent
  tmpe1 <- list(expression({
    sp::plot(pecoregs2, border = "gray", lwd = 1, add = TRUE)
  }))

  tmpe2 <- list(expression({
    sp::plot(pecoregs2, border = "gray", lwd = 1, add = TRUE)
    sp::plot(spoly_mask_Shriver2018, border = "orange", lwd = 2, add = TRUE)
  }))

  if (diff(used_subsets) == 0) {
    addfun_mtx[1, 1] <- tmpe1
    addfun_mtx[2, 1] <- tmpe2

  } else {
    addfun_mtx[1, 1] <- tmpe2
    addfun_mtx[2, 1] <- tmpe1
  }

  # scatterplot info
  vlim_mtx[, 2] <- list(c(0, 1))
  type_mtx[, 2] <- "smoothScatter.vals"
  axlabs_mtx[, 2] <- list(c(
    vars_table[defs_q1a[2], "label_manuscript"],
    vars_table[defs_q1a[1], "label_manuscript"]
  ))


  data_mtx[2, 2][[1]] <- dats_q1a[, 2:1]

  if (diff(used_subsets) == 0) {
    data_mtx[1, 2][[1]] <- dats_q1a[, 2:1]
    data_mtx[2, 2][[1]][!def_subsets[["GreatBasin"]],] <- NA
  }



  rSW2analysis::plot_matrix_of_panels(
    n_panels = n_panels,
    type_matrix = type_mtx,
    add_1to1 = TRUE,
    data_matrix = data_mtx,
    meta = SFSW2_prj_meta,
    subset = NULL,
    asp = 1,
    xlim_matrix = vlim_mtx,
    ylim_matrix = vlim_mtx,
    zlim_matrix = zlim_mtx,
    addfun_matrix = addfun_mtx,
    label_title_matrix = title_mtx,
    label_axis_matrix = axlabs_mtx,
    use_labels = "legend_title",
    map_legend_pos = legend_pos,
    map_extent = extent_subsets[[min(used_subsets)]],
    fexp_axis = 1,
    path = file.path(dir_res2_q1a),
    ftag = paste0(
      "Maps_GISSM_vs_Shriver2018_",
      sc_hist, "-", default_subprj, "_", vtag
    ),
    pborders = pborders
  )
}

# Code intended for manual, interactive use:
if (do_regress_models && interactive() && sys.nframe() == 0) {

  # Spline regression:
  library("mgcv")

  tmp_data2 <- data.frame(
    y = dats_q1a[, var_pGISSM],
    npos = as.integer(round(dats_q1a[, var_pGISSM] * nyrs)),
    nneg = as.integer(round((1 - dats_q1a[, var_pGISSM]) * nyrs)),
    x = scale(dats_q1a[, var_pShriver2018])
  )

  res_m <- mgcv::gam(
    cbind(npos, nneg) ~ s(x, bs = "cr", k = 100),
    family = binomial(),
    data = tmp_data2
  )
  res_m <- glm(
    cbind(npos, nneg) ~ x + I(x ^ 2),
    family = binomial(),
    data = tmp_data2
  )


  # Check model
  check_m_concurv <- mgcv::concurvity(res_m)
  check_m <- mgcv::gam.check(res_m)

  check_m_condicolli <- perturb::colldiag(res_m)
  check_m_vif1 <- car::vif(res_m)
  check_m_vif2 <- rms::vif(res_m)

  visreg::visreg(res_m, "x",
    scale = "response", ylim = c(0, 1),
    overlay = TRUE)

  plot(res_m)
  R2beta <- r2glmm::r2beta(res_m, partial = FALSE)
  R2GLMM_tm <- MuMIn::r.squaredGLMM(res_m)
  R2LR_tm <- MuMIn::r.squaredLR(res_m)

  sim_res <- DHARMa::simulateResiduals(res_m)
  plot(sim_res)

}

if (do_ecoregion_table) {
  f1 <- function(x) c(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )

  f2 <- function(x) {
    if (sum(complete.cases(x)) > 0) {
      tmp <- t.test(
        x[, 1], x[, 2],
        paired = TRUE
      )
      tmp[c("statistic", "parameter", "p.value")]
    } else {
      c(statistic = NA, parameter = NA, p.value = NA)
    }
  }

  # aggregate for cec-iii ecoregions
  ttmp1a <- aggregate(
    x = dats_q1a,
    by = list(rep("All", length(id_ecoregs2))),
    FUN = f1
  )


  ttmp1b <- aggregate(
    x = dats_q1a,
    by = list(id_ecoregs2),
    FUN = f1
  )

  ttmp2a <- f2(dats_q1a)

  ttmp2b <- by(
    data = dats_q1a,
    INDICES = id_ecoregs2,
    FUN = f2
  )

  ttmp2c <- do.call(rbind, ttmp2b)
  tmp2b <- matrix(
    unlist(ttmp2c),
    ncol = 3,
    dimnames = list(NULL, colnames(ttmp2c))
  )

  stopifnot(as.character(ttmp1b$Group.1) == rownames(ttmp2c))

  write.csv(
    rbind(
      cbind(ttmp1a, ttmp2a),
      cbind(ttmp1b, tmp2b)
    ),
    file = file.path(dir_res2_q1a,
      paste0(
        "Table_GISSM_vs_Shriver2018_",
        sc_hist, "-", default_subprj,
        "_EcoregionsIII_", vtag, ".csv"
      )
    ),
    row.names = FALSE
  )


  if (FALSE) {
    r1 <- rSW2st::create_raster_from_variables(
      data = dats_q1a[, 1],
      site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
    )
    plot(r1)
    plot(pecoregs2, add = TRUE)

    r2 <- rSW2st::create_raster_from_variables(
      data = dats_q1a[, 2],
      site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
    )

    plot(r2)
    plot(pecoregs2, add = TRUE)
  }
}




#------ Map of CEC-III ------
fname_map_ceciii <- file.path(
  dir_res2,
  "Output1_GeographicPatterns",
  "Maps_CEC_Ecoregions_LevelIII"
)

if (
  length(list.files(dirname(fname_map_ceciii), basename(fname_map_ceciii))) > 0
) {
  stopifnot(exists("id_ecoregs2"), exists("pecoregs2"))

  n_panels <- c(1, 1)
  template <- array(list(), dim = n_panels)

  # data
  tmp_data <- rep(1, length(id_ecoregs2)) #as.integer(factor(id_ecoregs2))
  tmp_data[!def_subsets[["BigSage"]]] <- NA

  # plot info
  addfun_mtx <- data_mtx <- zlim_mtx <- vlim_mtx <- type_mtx <- title_mtx <-
    axlabs_mtx <- legend_mtx <- template

  # map info
  zlim_mtx[, 1] <- list(c(-1, 3))
  type_mtx[, 1] <- "map.delta"

  data_mtx[1, 1][[1]] <- tmp_data
  addfun_mtx[1, 1] <- list(expression({
    sp::plot(pecoregs2, border = "darkgray", lwd = 1, add = TRUE)
    raster::text(
      x = raster::aggregate(pecoregs2, by = "NA_L3CODE"),
      labels = "NA_L3CODE",
      cex = 0.65
    )
  }))
  legend_mtx[,] <- FALSE

  rSW2analysis::plot_matrix_of_panels(
    n_panels = n_panels,
    type_matrix = type_mtx,
    data_matrix = data_mtx,
    meta = SFSW2_prj_meta,
    subset = NULL,
    asp = 1,
    zlim_matrix = zlim_mtx,
    addfun_matrix = addfun_mtx,
    label_title_matrix = title_mtx,
    label_axis_matrix = axlabs_mtx,
    use_labels = "legend_title",
    legend_matrix = legend_mtx,
    map_legend_pos = "right",
    map_extent = extent_subsets[["BigSage"]],
    fexp_axis = 1,
    path = dirname(fname_map_ceciii),
    ftag = basename(fname_map_ceciii),
    pborders = pborders
  )
}



#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q1a)




#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
