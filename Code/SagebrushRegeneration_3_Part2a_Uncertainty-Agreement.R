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
# 3b. How much do GCM-runs agree in direction of responses

#--- Specify variables
defs_q2a <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)

vars_model <- vars_table[defs_q2a, "analysis_group"]
vars <- vars_table[defs_q2a, "name_original_rSFSW2"]
var_pGISSM <- vars[1]
var_pShriver2018 <- vars[2]



#------------------------------------------------------------------------------#
#--- Load aggregated simulation output from SOILWAT2 runs
print(paste(Sys.time(), "--", "Load simulation data"))

ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q2a.rds")

if (file.exists(ftmp)) {
  dats_q2a <- readRDS(ftmp)

} else {
  dats_q2a <- list(
    vals = load_rSFSW2_BSR_data_from_netCDFs(
      path = dir_prj,
      locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
      vars_table = vars_table,
      ids_vars = defs_q2a,
      sim_scenarios = sim_scens_id,
      sim_subprojects = default_subprj
    )
  )

  dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dats_q2a, file = ftmp)
}


#--- Calculate absolute and relative differences between
# * future and historical time periods
dats_q2a <- c(
  dats_q2a,
  rSW2analysis::calc_change_from_reference(
    data = dats_q2a[["vals"]],
    ref_condition = sc_hist,
    sc_hist = sc_hist,
    method = "absolute"
  )
)


#--- Calculate cell-wise ensembles across GCMs x RCPs x future-times
datEns_q2a <- list()
for (kn in names(dats_q2a)) {
  datEns_q2a[[kn]] <- rSW2analysis::calc_cellwise_ensemble(
    data = dats_q2a[[kn]],
    funs = c(cfun, "agreement"),
    na.rm = TRUE,
    verbose = TRUE
  )
}





#------------------------------------------------------------------------------#
#------ Analyse data
do_map_agreement <- TRUE
do_map_agreement_small <- TRUE
do_table_overview <- TRUE
do_ecoregion_table <- TRUE


legend_pos <- switch(
  names(def_subsets)[min(used_subsets)],
  BigSage = "right",
  "left"
)

dir_res2_q2a <- file.path(
  dir_res2,
  "Output2_Uncertainty"
)

if (do_map_agreement) {
  # Agreement figure
  #  - rows: 2 pairs by regen-model: each future time periods
  #  - columns: 2 RCPs

  alpha_agree_limit <- 90

  n_panels <- c(2 * length(req_dTime[-1]), length(reqCSs))
  template <- array(list(), dim = n_panels)

  type_mtx <- legend_mtx <- zlim_mtx <- template

  type_mtx[] <- "map.agreement_sign"
  zlim_mtx[] <- list(100 * c(-1, 1))
  legend_mtx[] <- FALSE
  legend_mtx[1, 1] <- TRUE

  data_mtx <- label_mtx <- addfun_mtx <- addenv_mtx <- template

  for (k1 in seq_len(n_panels[1])) for (k2 in seq_len(n_panels[2])) {
    iv <- 1 + (k1 - 1) %/% 2
    it <- 1 + (k1 - 1) %% 2

    id_eus <- grep(reqCSs[k2], ensemble_units)
    n_gcms <- sum(grepl(reqCSs[k2],
      unlist(SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]])
    ))

    label_mtx[k1, k2] <- paste0(
      vars_table[defs_q2a[iv], "label_manuscript"], ": ",
      reqCSs[k2],
      " (", req_CalendarPeriods[-1][it], ")"
    )

    data_agree <-
      100 / n_gcms * datEns_q2a[["delta_Current_abs"]]["agreement", , vars[iv], id_eus[it], default_subprj]

    data_mtx[k1, k2][[1]] <- rSW2analysis::get_data_for_agreement_sign(
      data_direction = datEns_q2a[["delta_Current_abs"]][cfun, , vars[iv], id_eus[it], default_subprj],
      data_agree = data_agree
    )

    addenv_mtx[k1, k2][[1]] <- list(
      shp_agree = rSW2st::isoline_from_raster(
        grid = rSW2st::create_raster_from_variables(
          data = data_agree,
          site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
          grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
        ),
        alpha = alpha_agree_limit
      )
    )

    addfun_mtx[k1, k2] <- list(expression(
        sp::plot(shp_agree, border = "black", lwd = 1, add = TRUE)
      )
    )
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
    legend_matrix = legend_mtx,
    map_legend_pos = legend_pos,
    map_extent = extent_subsets[[min(used_subsets)]],
    addfun_matrix = addfun_mtx,
    addenv_matrix = addenv_mtx,
    path = dir_res2_q2a,
    ftag = paste(
      "Agreement_GISSM_vs_Shriver2018",
      default_subprj, cfun,
      paste0(alpha_agree_limit, "perc-agreement"),
      vtag,
      sep = "_"
    ),
    pborders = pborders
  )
}


if (do_map_agreement_small) {
  # Agreement figure for each RCP and time period
  #  - rows: 2 pairs by regen-model
  #  - 1 column
  alpha_agree_limit <- 90

  for (it in seq_along(req_dTime[-1])) {
    for (ircp in seq_along(reqCSs)) {

      n_panels <- c(2, 1)
      template <- array(list(), dim = n_panels)

      type_mtx <- legend_mtx <- zlim_mtx <- template

      type_mtx[] <- "map.agreement_sign"
      zlim_mtx[] <- list(100 * c(-1, 1))
      legend_mtx[] <- FALSE
      legend_mtx[1, 1] <- TRUE

      data_mtx <- label_mtx <- addfun_mtx <- addenv_mtx <- template

      id_eus <- grep(
        pattern = paste0(req_dTime[1 + it], ".", reqCSs[ircp]),
        x = ensemble_units,
        value = TRUE
      )

      n_gcms <- sum(grepl(
        pattern = reqCSs[ircp],
        x = unlist(SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]])
      ))

      for (iv in seq_along(vars)) {
        label_mtx[iv, 1] <- paste0(
          vars_table[defs_q2a[iv], "label_manuscript"], ": ",
          reqCSs[ircp], " (", req_CalendarPeriods[-1][it], ")"
        )

        data_agree <-
          100 / n_gcms * datEns_q2a[["delta_Current_abs"]]["agreement", , vars[iv], id_eus, default_subprj]

        data_mtx[iv, 1][[1]] <- rSW2analysis::get_data_for_agreement_sign(
          data_direction = datEns_q2a[["delta_Current_abs"]][cfun, , vars[iv], id_eus, default_subprj],
          data_agree = data_agree
        )

        addenv_mtx[iv, 1][[1]] <- list(
          shp_agree = rSW2st::isoline_from_raster(
            grid = rSW2st::create_raster_from_variables(
              data = data_agree,
              site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
              grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
            ),
            alpha = alpha_agree_limit
          )
        )

        addfun_mtx[iv, 1] <- list(expression(
          sp::plot(shp_agree, border = "black", lwd = 1, add = TRUE)
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
        legend_matrix = legend_mtx,
        map_legend_pos = legend_pos,
        map_extent = extent_subsets[[min(used_subsets)]],
        addfun_matrix = addfun_mtx,
        addenv_matrix = addenv_mtx,
        path = dir_res2_q2a,
        ftag = paste(
          "Agreement_GISSM_vs_Shriver2018",
          default_subprj, cfun, id_eus,
          paste0(alpha_agree_limit, "perc-agreement"),
          vtag,
          sep = "_"
        ),
        pborders = pborders
      )
    }
  }
}

# Overview table:
if (do_table_overview) {
  N_geo <- sapply(def_subsets[used_subsets], sum)

  resp_dir <- c("decreases", "increases")
  p_agree <- c(0.75, 0.90, 0.95)
  resp_agree <- c("any", paste0(">", 100 * p_agree, "%"))

  resp_type <- expand.grid(
    Direction = resp_dir,
    Agreement = resp_agree,
    Model = vars_model,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  table_overview <- t(matrix(
    NA,
    nrow = length(ensemble_units),
    ncol = nrow(resp_type),
    dimnames = list(ensemble_units,
      apply(resp_type, 1, paste, collapse = ".")
    )
  ))

  for (k in seq_along(ensemble_units)) {
    n_gcms <- sum(grepl(
      rSW2analysis:::find_scen_part(0, ensemble_units[k]),
      unlist(SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]])
    ))

    for (iv in seq_along(vars_model)) {
      x_delta <- datEns_q2a[["delta_Current_abs"]][cfun, , vars[iv], ensemble_units[k], default_subprj]

      table_overview[paste0("decreases.any.", vars_model[iv]), k] <-
        sum(x_delta < 0, na.rm = TRUE) / N_geo[iv]
      table_overview[paste0("increases.any.", vars_model[iv]), k] <-
        sum(x_delta > 0, na.rm = TRUE) / N_geo[iv]

      x_agree <- 1 / n_gcms * datEns_q2a[["delta_Current_abs"]]["agreement", , vars[iv], ensemble_units[k], default_subprj]

      for (ip in seq_along(p_agree)) {
        ptag <- paste0(">", 100 * p_agree[ip], "%")
        table_overview[paste0("decreases.", ptag, ".", vars_model[iv]), k] <-
          sum(x_delta < 0 & x_agree > p_agree[ip], na.rm = TRUE) / N_geo[iv]
        table_overview[paste0("increases.", ptag, ".", vars_model[iv]), k] <-
          sum(x_delta > 0 & x_agree > p_agree[ip], na.rm = TRUE) / N_geo[iv]
      }
    }
  }

  dir.create(dir_res2_q2a, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    x = table_overview,
    file = file.path(dir_res2_q2a,
      paste0(
        paste(
          "Table_Agreement_GISSM_vs_Shriver2018_delta_Current_abs",
          default_subprj, tag_subset, cfun,
          sep = "_"
        ),
        ".csv"
      )
    ),
    row.names = TRUE
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
  for (k in seq_along(ensemble_units)) {
    id_eus <- ensemble_units[k]

    n_gcms <- sum(grepl(
      rSW2analysis:::find_scen_part(0, id_eus),
      unlist(SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]])
    ))

    tmp_data <- 1 / n_gcms * datEns_q2a[["delta_Current_abs"]]["agreement", , vars, id_eus, default_subprj]

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
      file = file.path(
        dir_res2_q2a,
        paste0(
          "Table_Agreement_GISSM_vs_Shriver2018_delta_Current_abs_",
          paste(default_subprj, cfun, id_eus, sep = "_"),
          "_EcoregionsIII_", vtag, ".csv"
        )
      ),
      row.names = FALSE
    )
  }
}


#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q2a, datEns_q2a)


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
