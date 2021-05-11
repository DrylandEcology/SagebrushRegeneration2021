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



#------------------------------------------------------------------------------#
# Research question:
#
# 2d. How much do models extrapolate beyond range of fitted field locations

#--- Specify variables
defs_q2d_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q2d_models, "analysis_group"]
models <- vars_table[defs_q2d_models, "name_original_rSFSW2"]
model_labels <- vars_table[defs_q2d_models, "label_manuscript"]

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
  "DegreeDays_Base0C_dailyTmean_Cdays_mean"
)

defs_q2d <- which(vars_table[, "name_original_rSFSW2"] %in% preds)



#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q2d.rds")

if (file.exists(ftmp)) {
  dats_q2d <- readRDS(ftmp)

} else {
  dats_q2d <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q2d,
      sim_scenarios = sim_scens_id,
      sim_subprojects = default_subprj
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q2d, file = ftmp)
}



#------------------------------------------------------------------------------#
#--- Load coordinates of field sites used for model development
dir_prj_raw <- file.path(dir_prj, "0_Simulation_Raw")

fname_xy_GISSM <- file.path(
  dir_prj_raw,
  "FieldSites_Schlaepfer+_2014_ECOMOD.csv"
)

if (file.exists(fname_xy_GISSM)) {
  xy_GISSM <- read.csv(file = fname_xy_GISSM)

} else {
  stop("Acquire field site coordinates from Schlaepfer et al. 2014 ECOMOD.")
}

fname_xy_Shriver2018 <- file.path(
  dir_prj_raw,
  "FieldSites_Shriver+_2018_GCB.csv"
)

if (file.exists(fname_xy_Shriver2018)) {
  xy_Shriver2018 <- read.csv(file = fname_xy_Shriver2018)

} else {
  stop("Acquire field site coordinates from Shriver et al. 2018 GCB.")
}


ids_models <- list(
  GISSM = find_siteIDs(
    grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]],
    locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    targets = xy_GISSM[, c("X_WGS84", "Y_WGS84")]
  ),
  Shriver2018 = find_siteIDs(
    grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]],
    locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    targets = xy_Shriver2018[, c("X_WGS84", "Y_WGS84")]
  )
)
stopifnot(names(ids_models) == model_names)



#------------------------------------------------------------------------------#
#------ Analyse data
extraps <- c("NT1", "NT2")

#--- Calculate degree of extrapolation: NT1 and NT2
fname_extrap <- file.path(dir_res_data, "tmp_q3a_extrap.rds")

if (file.exists(fname_extrap)) {
  dats_q2d[["extrap"]] <- readRDS(fname_extrap)

} else {
  dats_q2d[["extrap"]] <- array(
    NA,
    dim = c(
      SFSW2_prj_meta[["sim_size"]][["runsN_sites"]],
      length(extraps),
      length(sim_scens_id),
      length(def_subsets),
      length(ids_models)
    ),
    dimnames = list(
      NULL,
      extraps,
      sim_scens_id,
      names(def_subsets),
      names(ids_models)
    )
  )


  for (ks in seq_along(def_subsets)) {
    # Geographic subset
    ids_geo <- def_subsets[[ks]]

    # Historical reference data
    refdat <- dats_q2d[["vals"]][, preds, sc_hist, default_subprj]

    for (k1 in seq_along(sim_scens_id)) {
      # Scenario data onto which models are projected
      prodat <- dats_q2d[["vals"]][ids_geo, preds, sim_scens_id[k1], default_subprj]

      for (k3 in seq_along(ids_models)) {
        tmp_refdat <- refdat[ids_geo & ids_models[[k3]], ]

        # Calculate NT1: univariate degree of extrapolation
        tmp <- try(rSW2analysis::calc_NT1(
          refdat = tmp_refdat, prodat = prodat),
          silent = TRUE)

        if (!inherits(tmp, "try-error")) {
          dats_q2d[["extrap"]][ids_geo, "NT1", k1, ks, k3] <- tmp
        }

        # Calculate NT2: multivariate degree of extrapolation
        tmp <- try(
          rSW2analysis::calc_NT2(
            refdat = tmp_refdat,
            prodat = prodat
          ),
          silent = TRUE
        )

        if (!inherits(tmp, "try-error")) {
          dats_q2d[["extrap"]][ids_geo, "NT2", k1, ks, k3] <- tmp
        }
      }
    }
  }

  saveRDS(dats_q2d[["extrap"]], file = fname_extrap)
}



#--- Calculate cell-wise ensembles across GCMs x RCPs x future-times
ens_fun <-  c("min", cfun, "max")
datEns_q2d <- list()
datEns_q2d[["vals"]] <- rSW2analysis::calc_cellwise_ensemble(
  data = dats_q2d[["extrap"]],
  funs = ens_fun,
  na.rm = TRUE,
  verbose = TRUE
)





#--- Make plots: loop over different geographic subsets
do_map_extrapolation <- TRUE
do_table_extrapolation <- TRUE

legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)


dir_res2_q2d <- file.path(
  dir_res2,
  "Output2_Uncertainty"
)

if (do_map_extrapolation) {
  # Extrapolation figure
  #   - separate figures
  #      * NT1, NT2
  #      * historical, future
  #   - rows: 2 pairs by regen-model: each time periods
  #   - columns: 2 RCPs

  for (ie in seq_along(extraps)) {

    tmp_hist <- cbind(
      dats_q2d[["extrap"]][, extraps[ie], sc_hist, used_subsets[1], model_names[1]],
      dats_q2d[["extrap"]][, extraps[ie], sc_hist, used_subsets[2], model_names[2]]
    )

    tmp_fut <- list(
      datEns_q2d[["vals"]][cfun, , extraps[ie], , used_subsets[1], model_names[1]],
      datEns_q2d[["vals"]][cfun, , extraps[ie], , used_subsets[2], model_names[2]]
    )

    tmp <- c(tmp_hist, unlist(tmp_fut))
    zlim_tmp <- switch(
      EXPR = extraps[ie],
      NT1 = c(quantile(tmp, probs = 0, na.rm = TRUE), 0),
      NT2 = c(0, quantile(tmp, probs = 1, na.rm = TRUE))
    )


    #--- historical
    n_panels <- c(length(ids_models), 1)
    template <- array(list(), dim = n_panels)

    type_mtx <- label_mtx <- legend_mtx <- zlim_mtx <- template
    data_mtx <- addfun_mtx <- template

    type_mtx[] <- "map.vals"
    zlim_mtx[] <- list(zlim_tmp)
    legend_mtx[] <- FALSE
    legend_mtx[1, 1] <- TRUE
    label_mtx[, 1] <- paste0(
      model_labels, ": ",
      rSW2analysis:::cur_to_hist(sc_hist),
      " (", req_CalendarPeriods[1], ")"
    )

    tmp <- tmp_hist

    if (extraps[ie] == "NT2") {
      tmp[tmp <= 1] <- 0
    }

    data_mtx[1, ][[1]] <- tmp[, 1]
    data_mtx[2, ][[1]] <- tmp[, 2]

    if (FALSE) {
      addfun_mtx[2, ] <- list(expression(
        sp::plot(spoly_mask_Shriver2018, border = "orange", lwd = 2, add = TRUE)
      ))
    }


    rSW2analysis::plot_matrix_of_panels(
      n_panels = n_panels,
      data_matrix = data_mtx,
      meta = SFSW2_prj_meta,
      subset = NULL,
      zlim_matrix = zlim_mtx,
      label_title_matrix = label_mtx,
      use_labels = "panel_identifier",
      type_matrix = type_mtx,
      col_rev = extraps[ie] == "NT1",
      legend_matrix = legend_mtx,
      map_legend_pos = legend_pos,
      map_extent = extent_subsets[[min(used_subsets)]],
      addfun_matrix = addfun_mtx,
      path = dir_res2_q2d,
      device = device_type,
      ftag = paste(
        "Extrapolation",
        default_subprj,
        "Historical", extraps[ie], vtag,
        sep = "_"
      ),
      pborders = pborders
    )



    #--- future
    n_panels <- c(2 * length(req_dTime[-1]), length(reqCSs))
    template <- array(list(), dim = n_panels)

    data_mtx <- addfun_mtx <- template
    type_mtx <- label_mtx <- legend_mtx <- zlim_mtx <- template

    type_mtx[] <- "map.vals"
    zlim_mtx[] <- list(zlim_tmp)
    legend_mtx[] <- FALSE
    legend_mtx[1, 1] <- TRUE

    for (k1 in seq_len(n_panels[1])) for (k2 in seq_len(n_panels[2])) {
      iv <- 1 + (k1 - 1) %/% 2
      it <- 1 + (k1 - 1) %% 2
      id_eus <- grep(reqCSs[k2], ensemble_units)


      label_mtx[k1, k2] <- paste0(
        model_labels[iv], ": ",
        reqCSs[k2],
        " (", req_CalendarPeriods[-1][it], ")"
      )

      id_eus <- grep(reqCSs[k2], ensemble_units)
      tmp <- tmp_fut[[iv]][, id_eus[it]]

      if (extraps[ie] == "NT2") {
        tmp[tmp <= 1] <- 0
      }

      data_mtx[k1, k2][[1]] <- tmp

      if (FALSE && iv == 2) {
        addfun_mtx[k1, k2] <- list(expression(
          sp::plot(spoly_mask_Shriver2018, border = "orange", lwd = 2, add = TRUE)
        ))
      }
    }


    rSW2analysis::plot_matrix_of_panels(
      n_panels = n_panels,
      data_matrix = data_mtx,
      meta = SFSW2_prj_meta,
      subset = NULL,
      zlim_matrix = zlim_mtx,
      label_title_matrix = label_mtx,
      use_labels = "panel_identifier",
      type_matrix = type_mtx,
      col_rev = extraps[ie] == "NT1",
      legend_matrix = legend_mtx,
      map_legend_pos = legend_pos,
      map_extent = extent_subsets[[min(used_subsets)]],
      addfun_matrix = addfun_mtx,
      path = dir_res2_q2d,
      device = device_type,
      ftag = paste(
        "Extrapolation",
        default_subprj,
        "Future", cfun, extraps[ie], vtag,
        sep = "_"
      ),
      pborders = pborders
    )
  }
}


# Overview table:
if (do_table_extrapolation) {
  Ns_geo <- sapply(def_subsets[used_subsets], sum)
  N_ens <- length(ensemble_units)
  N_funs <- length(ens_fun)

  resp_type <- expand.grid(
    Extrapolation = extraps,
    Model = model_labels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  table_overview <- t(matrix(
    NA,
    nrow = 1 + N_funs * N_ens,
    ncol = nrow(resp_type),
    dimnames = list(
      c(
        sc_hist,
        paste0(
          rep(ensemble_units, each = N_funs), "_",
          rep(ens_fun, times = N_ens)
        )
      ),
      paste0("NoExtrapolation_", apply(resp_type, 1, paste, collapse = "."))
    )
  ))

  for (k in seq_len(ncol(table_overview))) {
    if (k > 1) {
      iens <- 1 + (k - 2) %/% N_funs
      ifun <- 1 + (k - 2) %% N_funs
      coltag <- paste0(ensemble_units[iens], "_", ens_fun[ifun])
    } else {
      iens <- ifun <- NA
      coltag <- sc_hist
    }

    for (iv in seq_along(models)) {
      if (k == 1) {
        tmp <- dats_q2d[["extrap"]][def_subsets[[used_subsets[iv]]], extraps, sc_hist, used_subsets[iv], model_names[iv]]
      } else {
        tmp <- datEns_q2d[["vals"]][ens_fun[ifun], def_subsets[[used_subsets[iv]]], extraps, ensemble_units[iens], used_subsets[iv], model_names[iv]]
      }


      table_overview[paste0("NoExtrapolation_NT1.", model_labels[iv]), coltag] <-
        sum(abs(tmp[, "NT1"]) < sqrt(.Machine$double.eps), na.rm = TRUE) / Ns_geo[iv]

      table_overview[paste0("NoExtrapolation_NT2.", model_labels[iv]), coltag] <-
        sum(tmp[, "NT2"] <= 1, na.rm = TRUE) / Ns_geo[iv]
    }
  }

  dir.create(dir_res2_q2d, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    table_overview,
    file = file.path(
      dir_res2_q2d,
      paste0(
        paste(
          "Table_Extrapolation",
          default_subprj,
          paste(ens_fun, collapse = "-"),
          sep = "_"
        ),
        ".csv"
      )
    ),
    row.names = TRUE
  )
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q2d, datEns_q2d)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
