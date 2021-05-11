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
#--- Specify variables
defs_q1b <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)

vars <- vars_table[defs_q1b, "name_original_rSFSW2"]
var_pGISSM <- vars[1]
var_pShriver2018 <- vars[2]

#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q1b.rds")

if (file.exists(ftmp)) {
  dats_q1b <- readRDS(ftmp)

} else {

  dats_q1b <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q1b,
      sim_scenarios = sim_scens_id,
      sim_subprojects = default_subprj
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q1b, file = ftmp)
}


#--- Calculate absolute differences between
# * future and historical time periods
# * end-of-century and mid-century time periods
for (method in c("absolute")) {
  dats_q1b <- c(
    dats_q1b,
    rSW2analysis::calc_change_from_reference(
      data = dats_q1b[["vals"]],
      ref_condition = c(sc_hist, req_dTime[2]),
      sc_hist = sc_hist,
      method = method
    )
  )
}


#--- Calculate cell-wise median ensembles across GCMs x RCPs x future-times
datEns_q1b <- list()
for (kn in names(dats_q1b)) {
  datEns_q1b[[kn]] <- rSW2analysis::calc_cellwise_ensemble(
    data = dats_q1b[[kn]],
    funs = cfun,
    na.rm = TRUE,
    verbose = TRUE
  )
}

datEns_q1bOld <- datEns_q1b

# This script analyses across-GCM values and not individual GCM-driven runs
rm(dats_q1b)


#------------------------------------------------------------------------------#
#--- Analyse data

do_fig_overview <- TRUE
do_fig_overview_small <- TRUE
do_rate_comparison <- TRUE
do_add_mask_Shriver2018 <- FALSE
do_ecoregion_table <- TRUE


dir_res2_q1b <- file.path(
  dir_res2,
  "Output1_GeographicPatterns",
  "FuturePatterns"
)



legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)


if (do_fig_overview) {
  # Overview figure:
  #  - rows: 2 pairs by regen-model: each future time periods
  #  - columns: 2 pairs by RCP: each 1st, values; 2nd, deltas

  n_panels <- c(2 * length(req_dTime[-1]), 2 * length(reqCSs))
  template <- array(list(), dim = n_panels)

  zlim_mtx <- type_mtx <- label_mtx <- addfun_mtx <- template
  zlim_mtx[, c(1, 3)] <- list(c(0, 1))
  zlim_mtx[, c(2, 4)] <- list(c(-1, 1))

  type_mtx[, c(1, 3)] <- "map.vals"
  type_mtx[, c(2, 4)] <- "map.delta"

  for (k1 in seq_len(n_panels[1])) {
    label_mtx[k1, c(1, 3)] <- paste0(
      vars_table[defs_q1b[1 + (k1 - 1) %/% 2], "label_manuscript"], ": ",
      reqCSs, " (", req_CalendarPeriods[-1][1 + (k1 - 1) %% 2], ")"
    )
  }

  legend_mtx <- array(FALSE, dim = n_panels)
  legend_mtx[1, 1:2] <- TRUE

  if (do_add_mask_Shriver2018) {
    addfun_mtx[3:4, c(1, 3)] <- list(expression(
      sp::plot(spoly_mask_Shriver2018, border = "orange", lwd = 2, add = TRUE)
    ))
  }


  data_mtx <- template

  for (ircp in seq_along(reqCSs)) {
    k2 <- (ircp - 1) * 2

    id_eus <- grep(reqCSs[ircp], ensemble_units)

    for (ie in seq_along(id_eus)) {
      for (iv in seq_along(vars)) {
        k1 <- (iv - 1) * length(vars) + ie

        # projections
        data_mtx[k1, k2 + 1][[1]] <- datEns_q1b[["vals"]][cfun, , vars[iv], id_eus[ie], default_subprj]

        # deltas
        data_mtx[k1, k2 + 2][[1]] = datEns_q1b[["delta_Current_abs"]][cfun, , vars[iv], id_eus[ie], default_subprj]
      }
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
    addfun_matrix = addfun_mtx,
    legend_matrix = legend_mtx,
    map_legend_pos = legend_pos,
    map_extent = extent_subsets[[min(used_subsets)]],
    path = dir_res2_q1b,
    device = device_type,
    ftag = paste(
      "Maps_GISSM_vs_Shriver2018",
      default_subprj, cfun, vtag, sep = "_"
    ),
    pborders = pborders
  )
}



if (do_fig_overview_small) {
  # Small overview figure by RCP and by future time period:
  #  - 2 rows: regen-model
  #  - 2 columns: each 1st, values; 2nd, deltas

  for (it in seq_along(req_dTime[-1])) {
    for (ircp in seq_along(reqCSs)) {

      n_panels <- c(2, 2)
      template <- array(list(), dim = n_panels)

      zlim_mtx <- type_mtx <- label_mtx <- addfun_mtx <- template
      zlim_mtx[, 1] <- list(c(0, 1))
      zlim_mtx[, 2] <- list(c(-1, 1))

      type_mtx[, 1] <- "map.vals"
      type_mtx[, 2] <- "map.delta"

      for (k1 in seq_along(vars)) {
        label_mtx[k1, 1] <- paste0(
          vars_table[defs_q1b[k1], "label_manuscript"], ": ",
          reqCSs[ircp], " (", req_CalendarPeriods[-1][it], ")"
        )
      }

      legend_mtx <- array(FALSE, dim = n_panels)
      legend_mtx[1, 1:2] <- TRUE

      if (do_add_mask_Shriver2018) {
        addfun_mtx[3:4, c(1, 3)] <- list(expression(
          sp::plot(spoly_mask_Shriver2018, border = "orange", lwd = 2, add = TRUE)
        ))
      }


      data_mtx <- template

      id_eus <- grep(
        pattern = paste0(req_dTime[1 + it], ".", reqCSs[ircp]),
        x = ensemble_units,
        value = TRUE
      )

      for (iv in seq_along(vars)) {
        # projections
        data_mtx[iv, 1][[1]] <- datEns_q1b[["vals"]][cfun, , vars[iv], id_eus, default_subprj]

        # deltas
        data_mtx[iv, 2][[1]] <- datEns_q1b[["delta_Current_abs"]][cfun, , vars[iv], id_eus, default_subprj]
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
        addfun_matrix = addfun_mtx,
        legend_matrix = legend_mtx,
        map_legend_pos = legend_pos,
        map_extent = extent_subsets[[min(used_subsets)]],
        path = dir_res2_q1b,
        device = device_type,
        ftag = paste(
          "Maps_GISSM_vs_Shriver2018",
          default_subprj, cfun, id_eus, vtag, sep = "_"
        ),
        pborders = pborders
      )

    }
  }
}



# Compare rates of change among time periods:
if (do_rate_comparison) {
  id_d40yrs <- grep(req_dTime[2], ensemble_units)
  id_d90yrs <- grep(req_dTime[3], ensemble_units)

  dt10 <- 1 / 10 * c(40, 40, 90, 90)


  n_panels <- c(length(req_dTime[-1]), length(vars))
  template <- array(list(), dim = n_panels)

  type_mtx <- axlabs_mtx <- template

  type_mtx[] <- "smoothScatter.vals"
  axlabs_mtx[] <- list(c(
    "Rate of change (1st half of century) [1 / 10 yrs]",
    "Rate of change (2nd half of century) [1 / 10 yrs]"
  ))

  data_mtx <- vlim_mtx <- label_mtx <- template
  dxr_lims <- 0
  for (iv in seq_along(vars)) for (ik in seq_along(reqCSs)) {
    label_mtx[ik, iv] <- paste0(
      vars_table[defs_q1b[iv], "label_manuscript"], ": ",
      reqCSs[ik]
    )

    dx <- cbind(
      datEns_q1b[["delta_Current_abs"]][cfun, , vars[iv], id_d40yrs, default_subprj],
      datEns_q1b[["delta_d40yrs_abs"]][cfun, , vars[iv], id_d90yrs, default_subprj]
    )
    dxr <- sweep(dx, MARGIN = 2, STATS = dt10, FUN = "/")
    dxr_lims <- range(
      dxr_lims,
      quantile(dxr, probs = c(0.005, 0.995), na.rm = TRUE)
    )

    data_mtx[ik, iv][[1]] <- dxr[, grep(reqCSs[ik], colnames(dxr))]
  }

  vlim_mtx[] <- list(dxr_lims, dxr_lims)

  rSW2analysis::plot_matrix_of_panels(
    n_panels = n_panels,
    type_matrix = type_mtx,
    add_1to1 = TRUE,
    data_matrix = data_mtx,
    subset = NULL,
    asp = 1,
    xlim_matrix = vlim_mtx,
    ylim_matrix = vlim_mtx,
    label_axis_matrix = axlabs_mtx,
    label_title_matrix = label_mtx,
    use_labels = "panel_identifier",
    fexp_axis = 1,
    path = dir_res2_q1b,
    device = device_type,
    ftag = paste(
      "RateOfChange_GISSM_vs_Shriver2018",
      default_subprj, cfun, vtag, sep = "_"
    )
  )
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

  # data
  for (it in seq_along(req_dTime[-1])) {
    for (ircp in seq_along(reqCSs)) {
      id_eus <- grep(
        pattern = paste0(req_dTime[1 + it], ".", reqCSs[ircp]),
        x = ensemble_units,
        value = TRUE
      )

      tmp_data <- datEns_q1b[["delta_Current_abs"]][cfun, , vars, id_eus, default_subprj]

      # aggregate for cec-iii ecoregions
      ttmp1a <- aggregate(
        x = tmp_data,
        by = list(rep("All", length(id_ecoregs2))),
        FUN = f1
      )


      ttmp1b <- aggregate(
        x = tmp_data,
        by = list(id_ecoregs2),
        FUN = f1
      )

      ttmp2a <- f2(tmp_data)

      ttmp2b <- by(
        data = tmp_data,
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
        file = file.path(dir_res2_q1b,
          paste0(
            "Table_GISSM_vs_Shriver2018_",
            paste(default_subprj, cfun, id_eus, sep = "_"),
            "_EcoregionsIII_", vtag, ".csv"
          )
         ),
        row.names = FALSE
      )
    }
  }

  if (FALSE) {
    r1 <- rSW2st::create_raster_from_variables(
      data = tmp_data[, 1],
      site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
    )
    plot(r1)
    plot(pecoregs2, add = TRUE)

    r2 <- rSW2st::create_raster_from_variables(
      data = tmp_data[, 2],
      site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
    )

    plot(r2)
    plot(pecoregs2, add = TRUE)
  }
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(datEns_q1b)


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
