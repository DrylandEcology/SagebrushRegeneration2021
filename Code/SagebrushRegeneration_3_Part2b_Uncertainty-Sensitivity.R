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

pkgs <- c("reshape2", "visreg", "MuMIn", "car", "dominanceanalysis", "vioplot")
stopifnot(all(sapply(pkgs, requireNamespace)))


#------------------------------------------------------------------------------#

#--- Specify variables
defs_q2b_models <- c(
  which(vars_table[, "analysis_group"] == "GISSM"),
  which(vars_table[, "analysis_group"] == "Shriver2018")
)
model_names <- vars_table[defs_q2b_models, "analysis_group"]
models <- vars_table[defs_q2b_models, "name_original_rSFSW2"]

preds <- c("MAT_C_mean", "MAP_mm_mean")
defs_q2b_preds <- which(
  vars_table[, "name_original_rSFSW2"] %in% preds
)

defs_q2b <- c(defs_q2b_models, defs_q2b_preds)


#------------------------------------------------------------------------------#
# If sensitivity analysis is already complete, then we don't need to load
# simulation values
dats_q2b <- list()



#------------------------------------------------------------------------------#
#------ Analyse data: Determine most important factor


# Overall and for each future time period separately
sets_varpart <- c("overall", req_dTime[-1])

formula_sets_varpart <- list(
  # Up to two-way interactions
  #   - but not SoilExp:dTime (because soils didn't change with time),
  #   - plus dTime:RCP:GCM (biggest contributor among higher interaction terms)
  overall = as.formula(value ~ 0 + (SoilExp + dTime + RCP + GCM) ^ 2 - SoilExp:dTime + dTime:RCP:GCM),
  d40yrs = tmp <- as.formula(value ~ 0 + (SoilExp + RCP + GCM) ^ 2),
  d90yrs = tmp
)

do_parts_template <- c(SS = TRUE, R2 = FALSE, DA = FALSE)
have_parts_template <- do_parts_template
have_parts_template[] <- FALSE

lresps <- list(
  SS = c("varpart", "SS"),
  R2 = c("R2GLMM_tm", "R2LR_n91"),
  DA = "DA"
)

resps <- unlist(lresps)


#--- Calculate sensitivity to factors for each gridcell
ids_geo <- def_subsets[["BigSage"]]

dir_tmp <- file.path(dir_res_data, "tmp_varpart")
dir.create(dir_tmp, recursive = TRUE, showWarnings = FALSE)

for (iv in seq_along(defs_q2b)) {

  varname <- vars_table[defs_q2b[iv], "name_original_rSFSW2"]

  vartype_is_climate <- varname %in% c("MAT_C_mean", "MAP_mm_mean")

  for (kv in seq_along(sets_varpart)) {
    tag_vp <- sets_varpart[kv]
    tag_name <- paste0(tag_vp, "__", varname)

    # set up formula and terms
    formula_varpart <- formula_sets_varpart[[kv]]

    if (vartype_is_climate) {
      # Remove SoilExp from formula because they did not affect climate inputs
      # if included, then overhelm fitting functions due to floating point arithmetic
      formula_varpart <- update(formula_varpart, . ~ . - SoilExp:. + 0)
    } else {
      # ?binomial
      # replace response with a two-column integer matrix: the first column
      # gives the number of successes and the second the number of failures
      formula_varpart <- update(formula_varpart, cbind(npos, nneg) ~ .)
    }

    terms_varpart <- attr(terms(formula_varpart), "term.labels")

    # Track progress
    do_parts <- do_parts_template
    have_parts <- have_parts_template

    # Create output container
    dats_q2b[[tag_name]] <- array(
      dim = c(
        SFSW2_prj_meta[["sim_size"]][["runsN_sites"]],
        length(terms_varpart) + 1,
        length(resps)
      ),
      dimnames = list(NULL, c(terms_varpart, "Residuals"), resps)
    )

    # Check whether we have results on files
    for (k in seq_along(do_parts)) {
      if (do_parts[k]) {
        fname_sets_varpart_resps <- file.path(
          dir_tmp,
          paste0(
            "tmp_q3c_var-partition_",
            sets_varpart[kv], "-",
            vars_table[defs_q2b[iv], "name_manuscript"], "_",
            names(do_parts)[k], ".rds"
          )
        )

        if (file.exists(fname_sets_varpart_resps)) {
          tmp <- readRDS(fname_sets_varpart_resps)

          dats_q2b[[tag_name]][, , lresps[[names(do_parts)[k]]]] <- tmp

          do_parts[k] <- FALSE
          have_parts[k] <- TRUE
        }
      }
    }

    # Calculate those that were not on file
    if (any(do_parts)) {

      #--- load data
      print(paste(Sys.time(), "--", "Load simulation data"))

      ftmp <- file.path(dir_res_data, "tmp_dbOut", "dats_q2b.rds")

      if (file.exists(ftmp)) {
        dats_q2b[["vals"]] <- readRDS(ftmp)

      } else {
        dats_q2b[["vals"]] <- load_rSFSW2_BSR_data_from_netCDFs(
          path = dir_prj,
          locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
          vars_table = vars_table,
          ids_vars = defs_q2b,
          sim_scenarios = sim_scens_id,
          sim_subprojects = tag_subprj
        )

        dir.create(dirname(ftmp), recursive = TRUE, showWarnings = FALSE)
        saveRDS(dats_q2b[["vals"]], file = ftmp)
      }


      # subset to time slice
      sim_scen_subset <- if (isTRUE(tag_vp %in% req_dTime)) {
        grep(tag_vp, sim_scens_full, value = TRUE)
      } else {
        sim_scens_full
      }

      #for (k in seq_len(SFSW2_prj_meta[["sim_size"]][["runsN_sites"]])) {
        # print(paste0(Sys.time(), ": ", tag_name,
        #   " / site = ", k, " out of ",
        #   SFSW2_prj_meta[["sim_size"]][["runsN_sites"]]))

      ik <- 1
      ks <- sum(ids_geo)

      print(paste0(Sys.time(), ": processing ", shQuote(tag_name)))

      pb <- utils::txtProgressBar(max = ks, style = 3)

      for (k in which(ids_geo)) {
        utils::setTxtProgressBar(pb, ik)
        ik <- ik + 1

        #--- Prepare data for gridcell k (use GCMs with full RCP coverage)
        dtemp <- reshape2::melt(dats_q2b[["vals"]][k, varname, sim_scen_subset, ])

        dtemp[, "SoilExp"] <- as.character(dtemp[, "Var2"])
        dtemp[, "Var1"] <- as.character(dtemp[, "Var1"])
        ids_hist <- dtemp[, "Var1"] == sc_hist

        dtemp[!ids_hist, "dTime"] <- rSW2analysis:::find_reqDeltaYR(dtemp[, "Var1"])
        dtemp[!ids_hist, "RCP"] <- rSW2analysis:::find_reqCSs(dtemp[, "Var1"])
        dtemp[!ids_hist, "GCM"] <- rSW2analysis:::find_reqMs(dtemp[, "Var1"])

        if (any(ids_hist)) {
          # Make data balanced: copy historical for each RCP/GCM
          dtemp_hist <- as.data.frame(matrix(
            NA,
            nrow = length(tag_subprj) * length(reqCSs) * length(reqGCMs_full),
            ncol = ncol(dtemp),
            dimnames = list(NULL, colnames(dtemp))
          ))

          dtemp_hist[, "dTime"] <- req_dTime[1]

          ids <- seq_along(tag_subprj)
          i <- 1

          for (rcp in reqCSs) for (gcm in reqGCMs_full) {
            ids2 <- ids + (i - 1) * length(tag_subprj)

            dtemp_hist[ids2, c("value", "SoilExp")] <-
              dtemp[ids_hist, c("value", "SoilExp")]

            dtemp_hist[ids2, "RCP"] <- rcp
            dtemp_hist[ids2, "GCM"] <- gcm

            i <- i + 1
          }

          dtemp2 <- rbind(dtemp[!ids_hist, ], dtemp_hist)

        } else {
          dtemp2 <- dtemp
        }

        if (FALSE) {
          aggregate(dtemp2[, "value"],
            by = list(dtemp2[, "SoilExp"], dtemp2[, "RCP"]),
            FUN = function(x) c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))

          aggregate(dtemp2[, "value"],
            by = list(dtemp2[, "SoilExp"], dtemp2[, "GCM"]),
            FUN = function(x) c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))

         aggregate(dtemp2[, "value"],
            by = list(dtemp2[, "SoilExp"], dtemp2[, "dTime"]),
            FUN = function(x) c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
        }

        #--- Calculate ANOVA with binomial family:
        if (!all(is.na(dtemp2[, "value"]))) {
          if (vartype_is_climate) {
            mfit <- try(glm(
              formula_varpart,
              data = dtemp2,
              family = gaussian,
              na.action = na.fail
            ))

          } else {
            # As a two-column integer matrix: the first column gives the number
            # of successes and the second the number of failures
            dtemp2[, "npos"] <- as.integer(round(nyrs * dtemp2[, "value"]))
            dtemp2[, "nneg"] <- as.integer(round(nyrs * (1 - dtemp2[, "value"])))

            mfit <- try(glm(
              formula_varpart,
              data = dtemp2,
              family = binomial,
              na.action = na.fail
            ))
          }

          if (!inherits(mfit, "try-error")) {
            if (FALSE) {
              visreg::visreg(mfit, "GCM",
                by = "dTime",
                scale = "response",
                ylim = c(0, 1),
                overlay = TRUE
              )
            }

            if (do_parts["R2"]) {
              # Null model
              mnull <- update(mfit, . ~ 1)

              # Determine pseudo-R2:
              warning(
                "reconsider using Jaeger's R2beta instead;",
                "and/or use delta instead theoretical variance estimate if using Nakagawa's R2'"
              )
              temp <- try(MuMIn::r.squaredGLMM(mfit, null = mnull), silent = TRUE)
              if (!inherits(temp, "try-error")) {
                dats_q2b[[tag_name]][k, "Residuals", "R2GLMM_tm"] <-
                  if (vartype_is_climate) {
                    temp[1, "R2m"]
                  } else {
                    temp["theoretical", "R2m"]
                  }
              }

              temp <- try(MuMIn::r.squaredLR(mfit, null = mnull), silent = TRUE)
              if (!inherits(temp, "try-error")) {
                dats_q2b[[tag_name]][k, "Residuals", "R2LR_n91"] <-
                  attr(temp, "adj.r.squared")
              }
            }

            #--- Variable importance: likelihood-ratio chi-square statistics from anova-II
            if (do_parts["SS"]) {
              # Note: not affected by overdispersion --> don't fit with quasibinomial
              # likelihood-ratio chi-square statistics == F-test sum of squares
              antemp <- try(car::Anova(
                  mfit,
                  type = "II",
                  test.statistic = "F",
                  error.estimate = "deviance"
                ),
                silent = TRUE)

              if (!inherits(antemp, "try-error")) {
                vim_lrII <- antemp[, "Sum Sq"]
                names(vim_lrII) <- rownames(antemp)

                dats_q2b[[tag_name]][k, names(vim_lrII), "SS"] <- vim_lrII
                dats_q2b[[tag_name]][k, names(vim_lrII), "varpart"] <- vim_lrII / sum(vim_lrII)
              }
            }

            #--- Dominance analysis
            # Could be really great but very slow, particularly for many terms
            if (do_parts["DA"]) {
              # https://cran.r-project.org/web/packages/dominanceanalysis/vignettes/da-logistic-regression.html
              da_averageContribution <- function(x, fit.functions = "r2.m") {
                daRaw <- dominanceanalysis:::daRawResults(x)
                daAverageByLevel <- dominanceanalysis:::daAverageContributionByLevel(daRaw)

                colMeans(daAverageByLevel[[fit.functions]][, -1])
              }

              tmp <- da_averageContribution(mfit)
              daavgcontr <- c(tmp, Residuals = 1 - sum(tmp))

              dats_q2b[[tag_name]][k, names(daavgcontr), "DA"] <- daavgcontr
            }
          }
        }
      }

      close(pb)

      for (k in seq_along(do_parts)) {
        if (do_parts[k]) {
          fname_sets_varpart_resps <- file.path(
            dir_tmp,
            paste0(
              "tmp_q3c_var-partition_",
              sets_varpart[kv], "-",
              vars_table[defs_q2b[iv], "name_manuscript"], "_",
              names(do_parts)[k], ".rds"
            )
          )

          saveRDS(
            dats_q2b[[tag_name]][, , lresps[[names(do_parts)[k]]]],
            file = fname_sets_varpart_resps
          )

          do_parts[k] <- FALSE
          have_parts[k] <- TRUE
        }
      }
    }


    #--- Determine first and second most sensitive/important factor
    # for each gridcell (besides residuals)
    tmp_res <- data.frame(
      imp = rep(NA, length(ids_geo)),
      varpart = NA,
      imp_2nd = NA
    )
    tmpd <- dats_q2b[[tag_name]][ids_geo, terms_varpart, "varpart"]

    ids_tmpduse <- which(complete.cases(tmpd))
    tmpd2 <- tmpd[ids_tmpduse, , drop = FALSE]
    ids_geo2 <- which(ids_geo)[ids_tmpduse]

    # First most sensitive factor
    ids1 <- apply(tmpd2, 1, which.max)

    tmp_res[ids_geo2, "imp"] <- factor(
      terms_varpart[ids1],
      levels = terms_varpart
    )
    tmp_res[ids_geo2, "varpart"] <- sapply(
      seq_along(ids1),
      function(i) {
        dats_q2b[[tag_name]][ids_geo2[i], terms_varpart[ids1[i]], "varpart"]
      }
    )

    # Second most sensitive factor (set first most sensitive factor to NA)
    tmpd2[as.matrix(cbind(row = seq_along(ids1), col = ids1))] <- NA
    ids2 <- apply(tmpd2, 1, which.max)
    tmp_res[ids_geo2, "imp_2nd"] <- factor(
      terms_varpart[ids2],
      levels = terms_varpart
    )

    # Save results
    dats_q2b[[paste0(tag_name, "_imp")]] <- tmp_res
  }
}

rm(ids_geo)


do_sensitivity_maps <- FALSE
do_sensitivity_violins <- TRUE
do_violins_values_change <- c(only_values = FALSE, values_and_change = TRUE)



dir_res2_q2b <- file.path(
  dir_res2, "Output2_Uncertainty"
)
dir.create(dir_res2_q2b, recursive = TRUE, showWarnings = FALSE)



if (do_sensitivity_violins) {
  ids_fig_groups <- list(
    resps = models,
    climpred = c("MAT_C_mean", "MAP_mm_mean")
  )

  fall <- list(
    resps = c(
      sort(attr(terms(formula_sets_varpart[["overall"]]), "term.labels")),
      "Residuals"
    ),
    climpred = c(
      sort(attr(terms(
        update(
          old = formula_sets_varpart[["overall"]],
          new = . ~ . - SoilExp:. + 0
        )),
        "term.labels"
      )),
      "Residuals"
    )
  )


  xlim <- lapply(fall, function(x) c(1, length(x)) + 0.5 * c(-1, 1))


  for (ifg in seq_along(ids_fig_groups)) for (ifx in do_violins_values_change) {
    # Set up figure
    mfx <- if (ifx) 2L else 1L

    n_panels <- c(length(sets_varpart), mfx * length(ids_fig_groups[[ifg]]))
    N_prj <- prod(n_panels)

    w.panel <- 2.75
    w.edgeL <- 0.35; w.edgeI <- if (ifx) 0.35 else 0.05; w.edgeR <- 0.01 # Width of panel and left/interior/right edge
    h.panel <- 2.5
    h.edgeL <- 0.75; h.edgeI <- 0.25; h.edgeU <- 0.25 # Heigth of panel and lower/interior/upper edge

    fexp <- 1
    fexp_axis <- 0.75 * fexp

    # Figure layout
    lmat <- matrix(
      0L,
      nrow = 2 + 2 * n_panels[1] - 1L,
      ncol = 2 + 2 * n_panels[2] - 1,
      byrow = TRUE
    )

    stemp2 <- seq_len(n_panels[2L])
    for (k in seq_len(n_panels[1L])) { # byrow = TRUE
      #lmat[k + 1L, -c(1L, 2L + n_panels[2L])] <- (k - 1L) * n_panels[2L] + stemp2
      lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- (k - 1L) * n_panels[2L] + stemp2
    }

    temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
    layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
    temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
    layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

    # Figure device
    file <- file.path(
      dir_res2_q2b,
      paste0(
        "Fig_Sensitivity_ExperimentalFactors_",
        names(ids_fig_groups)[ifg],
        if (ifx) "_withValuesAndChange", "_", vtag,
        ".png"
      )
    )
    png(
      units = "in", res = 150,
      height = fexp * sum(layout_heights),
      width = fexp * sum(layout_widths),
      file = file
    )

    layout(lmat, heights = layout_heights, widths = layout_widths)

    par_prev <- par(mar = rep(0.25, 4), mgp = c(1, 0, 0), tcl = 0.3, cex = 1)

    # Panel count
    i <- 1
    id_ticks1 <- seq(N_prj - n_panels[2] + 1, N_prj)
    id_ticks2 <- if (ifx) {
        seq_len(N_prj)
      } else {
        (seq_len(n_panels[1]) - 1) * n_panels[2] + 1
      }

    # Data containers
    if (ifx) {
      tmp <- c("mean", "upper", "lower", "median", "q1", "q3")
      table_vp <- array(
        NA,
        dim = c(length(fall[[ifg]]), length(tmp), n_panels),
        dimnames = list(
          fall[[ifg]],
          tmp,
          sets_varpart,
          paste0(
            rep(ids_fig_groups[[ifg]], each = mfx),
            "__",
            c("ExplainedVariance", if (ifx) "Change")
          )
        )
      )
    }

    # Panels
    for (k1 in seq_len(n_panels[1])) {
      tag_vp <- sets_varpart[k1]

      for (k2 in seq_len(n_panels[2])) {
        is_k2_vals <- if (ifx) as.logical(k2 %% 2) else TRUE

        ik2 <- if (ifx) (1 + k2) %/% 2 else k2
        varname <- ids_fig_groups[[ifg]][ik2]
        iv <- which(varname == vars_table[defs_q2b, "name_original_rSFSW2"])
        vartag <- vars_table[defs_q2b[iv], "name_manuscript"]
        tag_name <- paste0(tag_vp, "__", varname)

        xall <- dats_q2b[[paste0("overall__", varname)]][, fall[[ifg]], "varpart"]
        x <- dats_q2b[[tag_name]][, , "varpart"]

        # Geographic subset
        if (
          grepl(
            "ArtemisiaTridentataNotSpecifiedSubspecies_Seedlings1stSeason_SuitableYears",
            varname
          )
        ) {
          tmp <- !def_subsets[[used_subsets[[1]]]]

        } else if (grepl("Shriver2018", varname)) {
          tmp <- !def_subsets[[used_subsets[[2]]]]

        } else {
          tmp <- !def_subsets[["BigSage"]]
        }

        if (!is.null(tmp)) {
          xall[tmp, ] <- NA
          x[tmp, ] <- NA
        }


        ats <- match(fall[[ifg]], colnames(x), nomatch = 0)
        id_fats <- which(ats > 0)
        x <- x[, ats]
        #ats <- which(fall[[ifg]] %in% colnames(x))

        if (k1 == 1 && !is_k2_vals) {
          plot.new()

        } else {

          if (is_k2_vals) {
            tmp <- x
            ylims <- c(0, 1)

          } else {
            tmp <- x - xall[, id_fats]
            ylims <- c(-1, 1)
          }

          xtmp <- as.data.frame(tmp)

          # Plot panel
          plot(
            NA,
            xlim = xlim[[ifg]],
            ylim = ylims,
            type = "n",
            axes = FALSE,
            ann = FALSE
          )

          resvp <- vioplot::vioplot(
            xtmp,
            ylim = ylims,
            at = id_fats,
            range = 0, # no whiskers, only interquartile range
            equalArea = TRUE, # no effect here because balanced experiment
            add = TRUE,
            xaxt = "n",
            yaxt = "n"
          )

          if (ifx) {
            stmp <- apply(xtmp, 2, sum, na.rm = TRUE)
            ttmp <- sum(apply(xtmp, 1, sum, na.rm = TRUE))
            table_vp[id_fats, "mean", k1, k2] <- stmp / ttmp

            table_vp[id_fats, -1, k1, k2] <- as.matrix(as.data.frame(resvp))
          }

          # Indicate constant time factor
          if (tag_vp != "overall") {
            itmp <- range(grep("dTime", fall[[ifg]]))

            rect(
              xleft = itmp[1],
              ybottom = 0,
              xright = itmp[2],
              ytop = 1,
              density = 5,
              angle = 45,
              col = "gray",
              border = TRUE
            )
          }

          # Annotate axes
          axis(
            1,
            at = seq_along(fall[[ifg]]),
            labels = FALSE,
            cex = fexp_axis
          )

          anny <- i %in% id_ticks2
          axis(2, labels = anny, cex = fexp_axis)
          if (anny) {
            mtext(side = 2, line = 1,
              text = if (is_k2_vals) "Explained variance" else "Change"
            )
          }

          # Beautify panel
          if (!is_k2_vals) {
            abline(h = 0, col = "red")
          }

          # Annotate panel
          if (i %in% id_ticks1) {
            text(
              x = seq_along(fall[[ifg]]),
              y = par("usr")[3],
              labels = fall[[ifg]],
              srt = 45, adj = c(1, 1.3), cex = fexp_axis, xpd = NA
            )
          }

          rSW2analysis::add_panel_identifier(
            i,
            add_label = TRUE,
            label = if (is_k2_vals) {
                paste0(
                  vars_table[defs_q2b[iv], "name_manuscript"], ": ",
                  if (k1 == 1) "overall" else req_CalendarPeriods[k1]
                )
              } else {
                ""
              },
            fexp = fexp
          )
        }

        i <- i + 1
      }
    }

    par(par_prev)
    dev.off()


    # Arrange table
    if (ifx) {
      tmp <- reshape2::melt(table_vp)
      ids <- as.character(tmp[, "Var2"]) %in% c("mean", "median", "q1", "q3")
      tmp2 <- tmp[ids, ]
      tmp2[, "value"] <- round(100 * tmp2[, "value"])

      res_vp <- reshape2::dcast(tmp2, Var2 + Var4 + Var3 ~ Var1)

      write.csv(
        res_vp,
        file = file.path(dir_res2_q2b,
          paste0(
            "Table_Sensitivity_ExperimentalFactors_",
            names(ids_fig_groups)[ifg], "_", vtag, ".csv"
          )
        )
      )
    }
  }
}





#------------------------------------------------------------------------------#
#--- Cleanup
rm(dats_q2b)



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
