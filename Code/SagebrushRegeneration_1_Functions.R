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


#------------------------------------------------------------------------------#
#--------- CUSTOM FUNCTIONS
#------------------------------------------------------------------------------#

if (!get0("SagebrushRegeneration_1_Functions", ifnotfound = FALSE)) {

pkgs <- c("rSW2st", "rSW2analysis", "raster", "sp")
stopifnot(all(sapply(pkgs, requireNamespace)))

if (!requireNamespace(
  "visreg",
  versionCheck = list(op = "==", version = "2.6-0")
)) {
  warning("Code works likely only with `visreg` v2.6-0!")
  # install.packages("https://cran.r-project.org/src/contrib/Archive/visreg/visreg_2.6-0.tar.gz", repos = NULL)
}

if (!requireNamespace("rSW2analysis")) {
  remotes::install_github("DrylandEcology/rSW2analysis")
}

if (!requireNamespace("rSW2st")) {
  remotes::install_github("DrylandEcology/rSW2st")
}



load_rSFSW2_BSR_data_from_netCDFs <- function(
  path,
  locations,
  vars_table,
  ids_vars,
  sim_scenarios,
  sim_subprojects
) {
  locations <- rSW2st::as_points(locations, to_class = "sf")

  res <- array(
    dim = c(
      nrow(locations),
      length(ids_vars),
      length(sim_scenarios),
      length(sim_subprojects)
    ),
    dimnames = list(
      NULL,
      vars_table[ids_vars, "name_original_rSFSW2"],
      sim_scenarios,
      sim_subprojects
    )
  )

  counter <- 0
  pb <- utils::txtProgressBar(max = prod(dim(res)[-1]), style = 3)

  for (k1 in seq_along(sim_subprojects)) {
    for (k2 in seq_along(sim_scenarios)) {

      if (sim_scenarios[k2] == "Current") {
        tag_mem <- "historical"
        tag_timerange <- "19800101-20101231"
        time_ids <- 1

      } else {
        tmp_scen <- strsplit(sim_scenarios[k2], split = ".", fixed = TRUE)[[1]]
        tag_mem <- paste(tmp_scen[4:3], collapse = "_")
        tag_timerange <- "20200101-20991231"
        time_ids <- if (tmp_scen[2] == "d90yrs") 2 else 1
      }


      for (k3 in seq_along(ids_vars)) {
        counter <- counter + 1
        utils::setTxtProgressBar(pb, counter)

        tmpv <- vars_table[ids_vars[k3], "name_original_rSFSW2"]

        if (tmpv %in% c("X_WGS84", "Y_WGS84")) {
          #--- Cases that require special treatment
          tmp <- sf::st_coordinates(
            sf::st_transform(locations, crs = 4326)
          )[, 1:2]

          res[, k3, k2, k1] <- tmp[, if (tmpv == "X_WGS84") 1 else 2]

        } else {

          #--- Locate netCDF
          ftmp <- sub(
            "TTT", tag_timerange,
            sub(
              "CCC", tag_mem,
              sub(
                "SoilSSS", sim_subprojects[k1],
                vars_table[ids_vars[k3], "filename_template_netCDF"]
              )
            )
          )

          dtmp <- file.path(path, vars_table[ids_vars[k3], "folder"])
          if (nchar(vars_table[ids_vars[k3], "subfolder"]) > 0) {
            dtmp <- file.path(dtmp, vars_table[ids_vars[k3], "subfolder"])
          }
          if (vars_table[ids_vars[k3], "function_of_soilsubproject"] == 1) {
            dtmp <- file.path(dtmp, paste0("SOILWAT2_", sim_subprojects[k1]))
          }

          fnc <- list.files(
            path = dtmp,
            pattern = ftmp,
            full.names = TRUE,
            recursive = TRUE
          )

          if (length(fnc) != 1) {
            browser()
          }
          #stopifnot(length(fnc) == 1)


          #--- Determine vertical index/count
          tmpz <- vars_table[ids_vars[k3], "z_level"]
          vertical_ids <- if (is.na(tmpz)) -1 else tmpz


          #--- Read netCDF & extract values at locations
          res[, k3, k2, k1] <- rSW2st::read_netCDF(
            fnc,
            var = vars_table[ids_vars[k3], "id_netCDF"],
            method = "xy_subset",
            xy_names = c("x", "y"),
            time_ids = time_ids,
            vertical_ids = vertical_ids,
            locations = locations,
          )
        }
      }
    }
  }

  close(pb)

  res
}



#' Calculate the probability of sagebrush establishment in response to
#' environmental variables in the year following seeding based on the
#' Shriver-2018 model
#'
#' @param Temp_mean A numerical vector; mean temperature from day 1 to 250
#'   in degree Celsius.
#' @param VWC_spring A numerical vector; mean soil moisture from day 70 to
#'   100 in the top soil layer (0–5 cm) in units of \var{m3 / m3}.
#'
#' @references Shriver, R. K., C. M. Andrews, D. S. Pilliod, R. S. Arkle,
#'   J. L. Welty, M. J. Germino, M. C. Duniway, D. A. Pyke, and
#'   J. B. Bradford. 2018. Adapting management to a changing world: Warm
#'   temperatures, dry soil, and interannual variability limit restoration
#'   success of a dominant woody shrub in temperate drylands.
#'   Global Change Biology 24:4972–4982.
p_Shriver2018 <- function(VWC_spring, Temp_mean) {
  logit_p <- 3.306 + 2.499 * VWC_spring - 0.289 * Temp_mean
  stats::plogis(logit_p) # inverse logit
}


find_siteIDs <- function(meta, xy, proj4string = NULL) {
  if (is.null(proj4string)) {
    proj4string <- raster::crs("+init=epsg:4326")
  }

  loc <- sp::SpatialPoints(coords = xy, proj4string = proj4string)
  loc <- sp::spTransform(loc, CRS = meta[["sim_space"]][["sim_crs"]])

  r <- raster::raster(meta[["sim_space"]][["sim_raster"]])
  r <- raster::init(r, fun = function(x) rep(0, x))
  ids <- raster::cellFromXY(r, sp::coordinates(loc))
  r[ids] <- 1

  ids <- raster::extract(r, meta[["sim_space"]][["run_sites"]])
  as.logical(ids)
}

find_siteIDs <- function(grid, locations, targets, proj4string = NULL) {
  if (is.null(proj4string)) {
    proj4string <- raster::crs("+init=epsg:4326")
  }

  loc <- sp::SpatialPoints(coords = targets, proj4string = proj4string)
  loc <- sp::spTransform(loc, CRS = raster::crs(grid))

  r <- raster::raster(grid)
  r <- raster::init(r, fun = function(x) rep(0, x))
  ids <- raster::cellFromXY(r, sp::coordinates(loc))
  r[ids] <- 1

  ids <- raster::extract(r, locations)
  as.logical(ids)
}



#' Intended for interactive use
check_models <- function(model, preds, data, datay = NULL,
  datax_scaled = NULL) {

  res <- list()

  #------ Collinearity
  #--- Variance inflation factors
  res[["vif_by_rms"]] <- try(rms::vif(model))
  res[["vif_by_car"]] <- car::vif(model)

  #--- Perturbation analysis to evaluate collinearity
  # If collinearity is a serious problem in the data, then the
  # estimates will be unstable and vary strongly.
  res[["perturb_cv"]] <- try({
    attach(data)
    check_m1_perturb <- perturb::perturb(model,
      pvars = preds,
      prange = rep(1, length(preds))
    )
    detach(data)
    tmp <- summary(check_m1_perturb, full = TRUE)
    tmp[["summ"]][, "s.d."] / tmp[["summ"]][, "mean"]
  })

  #--- Condition indexes and variance decomposition proportions
  # If the largest condition index (the condition number) is large
  # (Belsley et al suggest 30 or higher), then there may be collinearity
  # problems
  res[["condicolli"]] <- try({
    tmp <- perturb::colldiag(model)
    print(tmp, fuzz = .3)
    tmp
  })

  #--- Overall and individual tests for multicollinearity: Farrar-Glauber test
  # https://datascienceplus.com/multicollinearity-in-r/
  if (!is.null(datax_scaled) && !is.null(datay)) {
    res[["omcdiag"]] <- mctest::omcdiag(datax_scaled, datay)
    res[["imcdiag"]] <- mctest::imcdiag(datax_scaled, datay)
  } else {
    res[["omcdiag"]] <- res[["imcdiag"]] <- NULL
  }

  if (!is.null(datax_scaled)) {
    res[["pcors"]] <- ppcor::pcor(datax_scaled, method = "pearson")
  }

  res
}

get_DHARMA <- function(x, data_sp = NULL, data_ID = NULL) {
  rdh <- list()

  rdh[["DHARMa_sr"]] <- DHARMa::simulateResiduals(x)

  if (!is.null(data_ID)) {
    rdh[["DHARMa_sr_nongrouped"]] <- rdh[["DHARMa_sr"]]
    rdh[["DHARMa_sr"]] <- DHARMa::recalculateResiduals(
      simulationOutput = rdh[["DHARMa_sr"]],
      group = data_ID
    )
  }

  #--- Goodness-of-fit tests on the scaled residuals
  # Uniformity, dispersion, outliers, Zero-inflation
  rdh[["DHARMa_tests"]] <- list()

  rdh[["DHARMa_tests"]][["uniformity"]] <-
    DHARMa::testUniformity(rdh[["DHARMa_sr"]], plot = FALSE)

  rdh[["DHARMa_tests"]][["dispersion"]] <-
    DHARMa::testDispersion(rdh[["DHARMa_sr"]], plot = FALSE)

  rdh[["DHARMa_tests"]][["outliers"]] <-
    DHARMa::testOutliers(rdh[["DHARMa_sr"]], plot = FALSE)

  rdh[["DHARMa_tests"]][["zeroinflation"]] <-
    DHARMa::testZeroInflation(rdh[["DHARMa_sr"]], plot = FALSE)

  # Spatial autocorrelation
  rdh[["DHARMa_tests"]][["spatial"]] <- if (!is.null(data_sp)) {
    DHARMa::testSpatialAutocorrelation(
      rdh[["DHARMa_sr"]],
      x = data_sp[, "X_WGS84"],
      y = data_sp[, "Y_WGS84"],
      plot = FALSE
    )
  }

  rdh
}


#' Residual plots comparing two models
get_resid_plots <- function(m1, m2, m_names, data_sp = NULL, subsets = NULL,
  data_ID = NULL, path = ".", ftag = "", device = c("png", "pdf")) {

  device <- match.arg(device)

  if (is.null(subsets)) {
    if (is.null(data_sp)) {
      tmp <- NULL
    } else {
      tmp <- rep(TRUE, nrow(data_sp))
    }

    ids_data <- list(tmp, tmp)

  } else {
    if (inherits(subsets, "list") && length(subsets) == 2) {
      ids_data <- list(
        subsets[[1]],
        subsets[[2]]
      )
    } else {
      ids_data <- list(subsets, subsets)
    }
  }

  if (inherits(data_ID, "list")) {
    data_ID <- list(data_ID[[1]], data_ID[[2]])
  } else {
    data_ID <- list(data_ID, data_ID)
  }


  xdh <- list(
    get_DHARMA(m1, data_sp = data_sp[ids_data[[1]], ], data_ID = data_ID[[1]]),
    get_DHARMA(m2, data_sp = data_sp[ids_data[[2]], ], data_ID = data_ID[[2]])
  )
  names(xdh) <- m_names

  if (all(lengths(xdh) > 0)) {
    # QQ-plot and residuals versus predicted
    n_panels <- c(length(xdh), 2 + as.integer(!is.null(data_sp)))

    fname <- file.path(
      path,
      paste0("Fig_ScaledResidualPlots", ftag, ".", device)
    )

    if (device == "png") {
      grDevices::png(
        filename = fname,
        units = "in",
        res = 300,
        height = n_panels[1] * 3.5,
        width = n_panels[2] * 3.5
      )

    } else if (device == "pdf") {
      grDevices::pdf(
        file = fname,
        height = n_panels[1] * 3.5,
        width = n_panels[2] * 3.5
      )
    }


    par_prev <- graphics::par(
      mfrow = n_panels,
      mar = c(2.5, 2.5, 1.5, 0.5),
      mgp = c(1, 0, 0),
      tcl = 0.3,
      cex = 1
    )


    for (k1 in seq_along(xdh)) {
      # for a correctly specified model we would expect
      # (i) a uniform (flat) distribution of the overall residuals
      DHARMa::plotQQunif(xdh[[k1]][["DHARMa_sr"]], ci = TRUE)

      # (ii) uniformity in y direction if we plot against any predictor:
      #     lines should be straight, horizontal, and at y-values of
      #     0.25, 0.5 and 0.75
      DHARMa::plotResiduals(
        xdh[[k1]][["DHARMa_sr"]],
        rank = TRUE,
        quantreg = TRUE,
        xlab = "Predicted values (rank transformed)",
        ylab = "Standardized residuals"
      )
      graphics::mtext(side = 3, line = 0.5, text = m_names[k1], font = 2)
      graphics::mtext(side = 3, line = 0,
        text = "Expect red lines to be horizontal at 0.25, 0.5, and 0.75",
        cex = 0.65
      )

      # Spatial pattern
      if (!is.null(data_sp)) {
        ftmp <- grDevices::colorRamp(c("red", "white", "blue"))
        ctmp <- ftmp(xdh[[k1]][["DHARMa_sr"]][["scaledResiduals"]])

        graphics::plot(
          x = data_sp[ids_data[[k1]], "X_WGS84"],
          y = data_sp[ids_data[[k1]], "Y_WGS84"],
          xlab = "", ylab = "",
          col = grDevices::rgb(ctmp, maxColorValue = 255),
          pch = 15,
          cex = 0.5 * graphics::par("cex")
        )

        ttmp <- xdh[[k1]][["DHARMa_tests"]][["spatial"]]
        graphics::mtext(
          side = 3, line = 0.5, cex = 0.8,
          text = ttmp[["method"]]
        )
        graphics::mtext(
          side = 3, line = 0, cex = 0.5, adj = 0.01,
          text = paste0(
            paste(
              names(ttmp[["statistic"]]), "=", signif(ttmp[["statistic"]], 3),
              collapse = ", "
            ),
            ", p-value = ", signif(ttmp[["p.value"]], 3)
          )
        )
      }
    }

    graphics::par(par_prev)
    grDevices::dev.off()
  }
}



#' Print table of pseudo-R2s and coefficients
get_model_background <- function(res, vdim1, vdim2 = NULL,
  path = ".", ftag = "") {

  if (nchar(ftag) > 0) {
    ftag <- paste0("_", ftag)
  }

  vloop <- if (is.null(vdim2)) 1L else vdim2

  get_dims <- function(var, vdim2) {
    if (is.null(vdim2)) {
      c(length(var), length(vdim1))
    } else {
      c(length(var), length(vdim1), length(vdim2))
    }
  }

  get_dimnames <- function(var, vdim2) {
    if (is.null(vdim2)) list(var, vdim1) else list(var, vdim1, vdim2)
  }

  get_fun_value <- function(res, vdim1, vdim2 = NULL, fun, ...) {
    if (is.null(vdim2)) {
      lapply(
        X = vdim1,
        FUN = function(d1) do.call(fun, args = list(res[[d1]], ...))
      )

    } else {
      if (length(vdim2) > 1) {
        lapply(
          X = vdim1,
          FUN = function(d1) {
            lapply(
              X = vdim2,
              FUN = function(d2) do.call(fun, args = list(res[[d1]][[d2]], ...))
            )
          }
        )
      } else {
        lapply(
          X = vdim1,
          FUN = function(d1) do.call(fun, args = list(res[[d1]][[vdim2]], ...))
        )
      }
    }
  }

  #--- Tables with pseudo-R2
  var <- c(
    "pseudoR2_McFadden",
    "pseudoR2_Nakagawa_marginal", "pseudoR2_Nakagawa_conditional",
    "pseudoR2_Jaeger_sgv_marginal",
    "rmsd", "cv_rmsd"
  )
  table_pseudoR2 <- array(
    NA,
    dim = get_dims(var, vdim2),
    dimnames = get_dimnames(var, vdim2)
  )

  for (k in seq_along(var)) {
    tmp <- unlist(get_fun_value(
      res = res,
      vdim1 = vdim1,
      vdim2 = vdim2,
      fun = getElement,
      name = var[k]
    ))

    if (is.null(vdim2)) {
      table_pseudoR2[var[k], ] <- tmp
    } else {
      tmp2 <- array(
        tmp,
        dim = c(length(vdim2), length(vdim1)),
        dimnames = list(vdim2, vdim1)
      )
      table_pseudoR2[var[k], , ] <- t(tmp2)
    }
  }

  if (!is.null(vdim2)) {
    tmp <- reshape2::melt(table_pseudoR2)
    table_pseudoR2 <- reshape2::dcast(tmp, Var1 + Var3 ~ Var2)
  }

  write.csv(
    table_pseudoR2,
    row.names = TRUE,
    file = file.path(path, paste0("Table_PseudoR2", ftag, ".csv"))
  )


  #--- Tables with coefficients
  var <- c("coef_unstandard", "coef_standard")

  for (k1 in seq_along(var)) {
    ftag2a <- if (nchar(ftag) > 0) paste0(ftag, "_", var[k1]) else ftag

    for (k2 in seq_along(vloop)) {
      ftag2b <- if (nchar(ftag2a) > 0 && !is.null(vdim2)) {
        paste0(ftag2a, "_", vdim2[k2])
      } else {
        ftag2a
      }

      tmp <- get_fun_value(
        res,
        vdim1,
        vdim2[k2],
        fun = getElement,
        name = var[k1]
      )

      if (length(unlist(tmp)) > 0) {

        coefuns <- data.frame(Predictor = rownames(tmp[[1]]))

        for (k3 in seq_along(vdim1)) {
          coefuns[, paste0("MeanSE", k3)] <- paste0(
            signif(tmp[[k3]][, "Estimate"], 3),
            " ± ",
            signif(tmp[[k3]][, "Std. Error"], 3)
          )
        }

        write.csv(
          coefuns,
          row.names = FALSE,
          file = file.path(path, paste0("Table_Coefficients", ftag2b, ".csv"))
        )
      }
    }
  }

}

plot_panels_2way_predictors <- function(
  panel_id, model,
  preds, response,
  xlabs, ylab, data, datay, datax_scale, datax_center,
  xlim_probs = c(0, 1),
  ylim,
  id_ticks1, id_ticks2,
  include_0_on_xaxis = FALSE,
  last_plot = NULL,
  add_legend = FALSE,
  cex_legend = 0.65
) {

  response <- force(response)

  # Transform variables back (after scaling)
  funtrans <- function(iv) {
    function(x) x * datax_scale[iv] + datax_center[iv]
  }

  # Plot panels for different predictor variables
  for (iv in seq_along(preds)) {
    iv1 <- iv
    iv2 <- if (iv == 1) 2 else 1

    add_points <- cbind(
      x = funtrans(iv1)(data[, preds[iv1]]),
      y = datay
    )

    add_points <- add_points[complete.cases(add_points), ]

    vbreaks <- quantile(
      data[!is.na(datay), preds[iv2]],
      probs = c(0.025, 0.5, 0.975)
    )

    # `visreg` needs `xlim` on untransformed scale, but
    # passes `xlim` to `plot.visreg` which would
    # require `xlim` on the transformed scale
    obj_visreg <- visreg::visreg(
      model,
      preds[iv1],
      scale = "response",
      xtrans = funtrans(iv1),
      xlab = xlabs[iv1],
      alpha = 0.05,
      by = preds[iv2],
      breaks = vbreaks,
      ylab = ylab,
      xlim = quantile(
        data[, preds[iv1]],
        probs = xlim_probs,
        na.rm = TRUE
      ),
      ylim = ylim,
      strip.names = paste0(
        xlabs[iv2], ": ",
        round(funtrans(iv2)(vbreaks), 2)
      ),
      plot = FALSE
    )

    plot(
      obj_visreg,
      ylim = ylim,
      overlay = TRUE,
      partial = FALSE,
      rug = FALSE,
      legend = FALSE,
      axes = FALSE,
      ann = FALSE
    )

    # Add axes
    prev_xpd <- par(xpd = NA)
    tmp <- axTicks(side = 1)
    axis(
      side = 1,
      at = if (include_0_on_xaxis) sort(unique(c(0, tmp))) else tmp,
      labels = panel_id %in% id_ticks1
    )
    axis(side = 2, labels = panel_id %in% id_ticks2)

    if (panel_id %in% id_ticks1) {
      title(xlab = xlabs[iv1])
    }

    if (panel_id %in% id_ticks2) {
      title(ylab = ylab)
    }
    par(prev_xpd)

    # Add data points
    # based on smoothScatter and on grDevices:::.smoothScatterCalcDensity
    map <- KernSmooth::bkde2D(
      x = add_points,
      bandwidth = diff(
        apply(
          add_points,
          MARGIN = 2,
          stats::quantile,
          probs = c(0.05, 0.95),
          na.rm = TRUE,
          names = FALSE)
      ) / 25,
      gridsize = c(128, 128)
    )

    xlimt <- quantile(
      add_points[, "x"], # transformed
      probs = xlim_probs,
      na.rm = TRUE
    )

    ids_in1 <- map$x1 < xlimt[2] & map$x1 > xlimt[1]
    ids_in2 <- map$x2 < ylim[2] & map$x2 > ylim[1]
    graphics::image(
      x = map$x1[ids_in1],
      y = map$x2[ids_in2],
      z = map$fhat[ids_in1, ids_in2] ^ 0.25,
      col = {
        tmp <- c("white", blues9)
        tmp <- sapply(seq_along(tmp), function(k)
          adjustcolor(tmp[k],
            alpha.f = ((k - 1) / (2 * (length(tmp) - 1))) ^ 0.5)
        )
        colorRampPalette(colors = tmp,  alpha = TRUE)(256)
      },
      xlim = xlimt,
      ylim = ylim,
      add = TRUE
    )

    # Add post-panel plotting function
    if (!is.null(last_plot)) {
      eval(last_plot)
    }

    # Add legend
    if (add_legend) {
      loc_legend <- rSW2analysis::legend_location(
        add_points,
        xlim = xlimt,
        ylim = ylim
      )

      legend(
        x = loc_legend,
        bty = "n",
        legend = round(funtrans(iv2)(vbreaks), 2),
        title = sub(" (", "\n(", xlabs[iv2], fixed = TRUE),
        fill = visreg:::pal(n = length(vbreaks)),
        cex = cex_legend,
        inset = c(
          if (grepl("left", loc_legend)) 0.15 else 0.05,
          if (grepl("top", loc_legend)) 0.1 else 0.05
        )
      )
    }

    # Annotate
    prev_xpd <- par(xpd = NA)
    rSW2analysis::add_panel_identifier(panel_id, add_label = TRUE)
    par(prev_xpd)

    panel_id <- panel_id + 1
  }

  panel_id
}


plot_model_relationships <- function(
  res_fit, predictors,
  responses,
  xlim_probs = c(0, 1),
  ylim,
  xlabs, ylabs,
  include_0_on_xaxis = FALSE,
  panels_by_row = TRUE,
  last_plot = NULL,
  cex_legend = 0.65,
  path = ".",
  ftag = "",
  device = c("png", "pdf")
) {

  device <- match.arg(device)
  n_models <- length(responses)

  if (inherits(predictors, "list")) {
    predictors <- predictors[seq_len(n_models)]
    models_have_same_preds <- all(
      sapply(predictors, function(x) identical(predictors[[1]], x))
    )

  } else {
    predictors <- lapply(seq_len(n_models), function(k) predictors)
    models_have_same_preds <- TRUE
  }

  models_have_same_preds <- all(
    models_have_same_preds,
    identical(
      res_fit[[responses[1]]][["datax_scaled"]][!is.na(res_fit[[responses[1]]][["datay"]]), ],
      res_fit[[responses[2]]][["datax_scaled"]][!is.na(res_fit[[responses[2]]][["datay"]]), ]
    )
  )


  n_preds <- max(lengths(predictors))

  if (inherits(xlabs, "list")) {
    xlabs <- xlabs[seq_len(n_models)]

  } else {
    xlabs <- lapply(seq_len(n_models), function(k) xlabs)
  }


  # Prepare figure
  n_panels <- c(n_preds, n_models)
  N_prj <- prod(n_panels)

  w.panel <- 3
  w.edgeL <- 0.35; w.edgeI <- if (panels_by_row) 0.05 else 0.2; w.edgeR <- 0.01 # Width of panel and left/interior/right edge
  h.panel <- 2.5
  h.edgeL <- 0.35; h.edgeI <- if (panels_by_row && models_have_same_preds) 0.15 else 0.25; h.edgeU <- 0.15 # Heigth of panel and lower/interior/upper edge

  fexp <- 1

  # Figure layout
  lmat <- matrix(
    0L,
    nrow = 2 + 2 * n_panels[1] - 1L,
    ncol = 2 + 2 * n_panels[2] - 1,
    byrow = panels_by_row
  )

  stemp2 <- seq_len(n_panels[2L])
  for (k in seq_len(n_panels[1L])) {
    tmp <- if (panels_by_row) {
      # byrow = TRUE / par(mfrow):
      (k - 1L) * n_panels[2L] + stemp2
    } else {
      # byrow = FALSE / par(mfcol):
      (stemp2 - 1L) * n_panels[1L] + k
    }
    lmat[(k - 1) * 2 + 2, seq_len(ncol(lmat)) %% 2 == 0] <- tmp
  }

  temp <- rep(c(h.panel, h.edgeI), n_panels[1L])
  layout_heights <- c(h.edgeU, temp[-length(temp)], h.edgeL)
  temp <- rep(c(w.panel, w.edgeI), n_panels[2L])
  layout_widths <- c(w.edgeL, temp[-length(temp)], w.edgeR)

  # Figure device
  fname <- file.path(path, paste0("Fig_VisReg_", ftag, ".", device))
  if (device == "png") {
    grDevices::png(
      filename = fname,
      units = "in",
      res = 150,
      height = fexp * sum(layout_heights),
      width = fexp * sum(layout_widths)
    )

  } else if (device == "pdf") {
    grDevices::pdf(
      file = fname,
      height = fexp * sum(layout_heights),
      width = fexp * sum(layout_widths)
    )
  }

  layout(lmat, heights = layout_heights, widths = layout_widths)

  par_prev <- par(mar = rep(0.25, 4), mgp = c(1, 0, 0), tcl = 0.3, cex = 0.9)

  # Panel count
  i <- 1
  if (panels_by_row) {
    id_ticks1 <- if (models_have_same_preds) {
      seq(N_prj - n_panels[2] + 1, N_prj)
    } else {
      seq_len(N_prj)
    }
    id_ticks2 <- (seq_len(n_panels[1]) - 1) * n_panels[2] + 1

  } else {
    if (FALSE) {
      id_ticks1 <- seq_len(n_panels[2]) * n_panels[1]
      id_ticks2 <- seq_len(n_panels[1])
    } else {
      id_ticks1 <- id_ticks2 <- seq_len(N_prj)
    }
  }

  # Loop over model_names
  for (iy in seq_len(n_models)) {

    # Plot panels for predictor variables
    i <- plot_panels_2way_predictors(
      panel_id = i,
      model = res_fit[[responses[iy]]][["m"]],
      preds = predictors[[iy]],
      response = responses[iy],
      xlabs = xlabs[[iy]],
      ylab = ylabs[iy],
      data = res_fit[[responses[iy]]][["datax_scaled"]],
      datay = res_fit[[responses[iy]]][["datay"]],
      datax_scale = res_fit[[responses[iy]]][["datax_scale"]],
      datax_center = res_fit[[responses[iy]]][["datax_center"]],
      xlim_probs = xlim_probs,
      ylim = ylim,
      id_ticks1 = id_ticks1,
      id_ticks2 = id_ticks2,
      include_0_on_xaxis = include_0_on_xaxis,
      add_legend = if (models_have_same_preds) iy == 2 else TRUE,
      cex_legend = cex_legend,
      last_plot = last_plot
    )
  }

  par(par_prev)
  dev.off()
}



#' @section Notes: Distance matrices are calculated which use large amounts of
#'  memory and may kill the application
get_dcor_tests <- function(data, models, model_names,
  preds_all, preds_selected, R = 199, path = ".", ftag = "", verbose = FALSE) {

  out1 <- c("All_predictors", "Selected_predictors", preds_selected)
  out2 <- c("Statistic", "p.value")
  res <- matrix(NA,
    nrow = length(models) * length(out1),
    ncol = length(out2),
    dimnames = list(NULL, out2)
  )
  res <- data.frame(
    Model = rep(model_names, each = length(out1)),
    Predictor = out1,
    R = R,
    res
  )

  fname <- file.path(path,
    paste0("Table_DistanceCorrelations_", ftag, ".csv"))

  save_res_to_file <- function(x, fname) {
    write.csv(x, row.names = FALSE, file = fname)
  }

  for (iy in seq_along(models)) {
    k <- (iy - 1) * length(out1)
    i <- 1

    # All predictors
    if (length(preds_all) > 0) {
      if (verbose) {
        print(paste("dcor:", Sys.time(), model_names[iy], "all-predictors"))
      }

      model_dist <- stats::dist(data[, models[iy]])

      tmp <- energy::dcor.test(
        x = model_dist,
        y = data[, preds_all],
        R = R
      )

      res[k + i, out2] <- c(tmp[["statistic"]], tmp[["p.value"]])
      save_res_to_file(res, fname)
      gc()
    }
    i <- i + 1


    # All selected predictors
    if (verbose) {
      print(paste("dcor:", Sys.time(), model_names[iy], "selected-predictors"))
    }

    tmp <- energy::dcor.test(
      x = model_dist,
      y = data[, preds_selected],
      R = R
    )

    res[k + i, out2] <- c(tmp[["statistic"]], tmp[["p.value"]])
    save_res_to_file(res, fname)
    gc()
    i <- i + 1


    # Partial distance correlations for each selected predictor
    # removing all other selected predictors
    for (iv in seq_along(preds_selected)) {
      if (verbose) {
        print(paste("dcor:", Sys.time(), model_names[iy], preds_selected[iv]))
      }

      tmp <- energy::pdcor.test(
        x = model_dist,
        y = data[, preds_selected[iv]],
        z = data[, preds_selected[-iv]],
        R = R
      )

      res[k + i, out2] <- c(tmp[["statistic"]], tmp[["p.value"]])
      save_res_to_file(res, fname)
      gc()
      i <- i + 1
    }
  }

  res
}


#' Determine bivariate unbiased distance correlation
calc_model_dcor2d <- function(data, predictors, responses) {
  tmps <- c(
    "dcorU",
    "sqr_load", "cluster_id", "best",
    "max_dcorU_to_prev", "pred_sel"
  )
  res <- array(
    NA,
    dim = c(length(predictors), length(responses), length(tmps)),
    dimnames = list(predictors, responses, tmps)
  )

  for (iy in seq_along(responses)) {
    tmp <- sapply(predictors, function(var) {
      iuse <- complete.cases(data[, c(var, responses[iy])])
      if (any(iuse)) {
        energy::dcor2d(
          x = data[iuse, var],
          y = data[iuse, responses[iy]],
          type = "U"
        )
      } else {
        -Inf
      }
    })

    ids <- match(predictors, names(tmp), nomatch = 0)
    res[ids > 0, iy, "dcorU"] <- tmp[ids]
  }

  res
}

#' Ascendant hierarchical clustering of variables
calc_asc_hierarch_clust <- function(
  data, predictors,
  do_stability = FALSE,
  path = ".", ftag = "", device = c("png", "pdf")
) {

  device <- match.arg(device)

  #--- Ascendant hierarchical clustering of variables
  dats_preds_hclusts <- ClustOfVar::hclustvar(X.quanti = data[, predictors])

  fname <- file.path(
    path,
    paste0("Fig_PredictorSelection_TreeDiagram__hclust_", ftag, ".", device)
  )

  if (device == "png") {
    grDevices::png(
      filename = fname,
      units = "in",
      res = 600,
      height = 14,
      width = 14
    )

  } else if (device == "pdf") {
    grDevices::pdf(
      file = fname,
      height = 14,
      width = 14
    )
  }

  plot(dats_preds_hclusts)
  dev.off()

  #--- Stability of hclusts with bootstrap
  if (do_stability) {
    dats_stability_hclusts <- ClustOfVar::stability(dats_preds_hclusts)
    nclusters <- as.integer(1 + which.max(dats_stability_hclusts[["meanCR"]]))

    fname <- file.path(
      path,
      paste0("Fig_PredictorSelection_TreeStability__hclust_", ftag, ".", device)
    )
    if (device == "png") {
      grDevices::png(
        filename = fname,
        units = "in",
        res = 600,
        height = 10,
        width = 5
      )

    } else if (device == "pdf") {
      grDevices::pdf(
        file = fname,
        height = 10,
        width = 5
      )
    }

    par_prev <- par(mfrow = c(2, 1))
    plot(dats_stability_hclusts)
    boxplot(dats_stability_hclusts[["matCR"]])
    par(par_prev)
    dev.off()

  } else {
    nclusters <- NA
  }

  list(
    hclust = dats_preds_hclusts,
    nclusters = nclusters
  )
}


#' Locate best variable for each cluster
locate_preds_per_cluster <- function(hclust_sel, predictor_info, nbest = 3) {
  nclusters <- length(hclust_sel[["var"]])
  ntmp <- dimnames(predictor_info)
  predictors <- ntmp[[1]]
  responses <- ntmp[[2]]

  # Report cluster assignment for each variable (model independent)
  for (k in seq_len(nclusters)) {
    ids <- which(hclust_sel[["cluster"]] == k)
    predictor_info[names(ids), , "cluster_id"] <- k

    tmp <- hclust_sel[["var"]][[paste0("cluster", k)]]
    if (is.null(nrow(tmp))) {
      predictor_info[names(ids), , "sqr_load"] <- tmp["squared loading"]

    } else {
      predictor_info[rownames(tmp), , "sqr_load"] <- tmp[, "squared loading"]
    }
  }

  # Identify variable with highest `nbest` dcor for each cluster per response
  for (iy in seq_along(responses)) {
    for (k in seq_len(nclusters)) {
      ids <- which(predictor_info[, iy, "cluster_id"] == k)
      tmp <- predictor_info[ids, iy, "dcorU"]
      stmp <- sort(tmp, decreasing = TRUE)
      id <- tmp >= stmp[min(length(stmp), nbest)]
      predictor_info[ids[id], iy, "best"] <- 1
    }
  }

  predictor_info
}


#' Determine predictors that are relatively best across responses
consensus_best_preds <- function(predictor_info, data, nclusters,
  limit_cor = 0.5
) {

  clusters <- sort(unique(as.vector(predictor_info[, , "cluster_id"])))
  predictors <- dimnames(predictor_info)[[1]]


  # Sort by highest dcorU per cluster
  sid <- apply(
    predictor_info[, , c("dcorU", "cluster_id"), drop = FALSE],
    c(1, 3),
    median
  )

  asid <- aggregate(sid, by = list(sid[, "cluster_id"]), max)
  cloid <- asid[order(-asid[, "dcorU"]), "cluster_id"]

  max_dcorU <- apply(
    predictor_info[, , "dcorU", drop = FALSE],
    2,
    max,
    na.rm = TRUE
  )

  res <- as.data.frame(matrix(
    NA,
    nrow = length(clusters),
    ncol = 5,
    dimnames = list(
      NULL,
      c("cluster_id", "var", "median_dcorU", "range_dcorU", "max_dcorU_to_prev")
    )
  ))

  for (k in seq_along(clusters)) {
    res[k, "cluster_id"] <- cloid[k]

    ids <- apply(
      predictor_info[, , "cluster_id", drop = FALSE],
      1,
      function(x) any(x %in% res[k, "cluster_id"])
    )
    ids2 <- apply(
      predictor_info[ids, , "best", drop = FALSE],
      1,
      function(x) any(x %in% 1)
    )
    ids3 <- which(ids)[ids2]

    tmpx <- sweep(
      predictor_info[ids3, , "dcorU", drop = FALSE],
      MARGIN = 2,
      STATS = max_dcorU,
      FUN = "/"
    )

    tmp <- apply(
      tmpx,
      1,
      function(x) c(mean = median(x), range = diff(range(x)))
    )

    # Remove variables with high correlation with already selected variables
    # of at least the same dcorU
    if (k > 1) {
      ptmp <- predictors[ids3]

      itmp <-
        seq_along(clusters) %in% 1:(k - 1) &
        res[, "median_dcorU"] > max(tmp["mean", ])
      vtmp <- na.exclude(res[itmp, "var"])

      if (length(vtmp) > 0) {
        ctmp <- array(
          -Inf,
          dim = c(length(vtmp), length(ptmp)),
          dimnames = list(vtmp, ptmp)
        )

        for (v1 in vtmp) for (v2 in ptmp) {
          iuse <- complete.cases(data[, c(v1, v2)])
          if (any(iuse)) {
            ctmp[v1, v2] <- energy::dcor2d(
              x = data[iuse, v1],
              y = data[iuse, v2],
              type = "U"
            )
          }
        }

        ids_rmv <- apply(ctmp, 2, function(x) any(x > limit_cor))
        tmp <- tmp[, !ids_rmv, drop = FALSE]

        predictor_info[ids3, , "max_dcorU_to_prev"] <- apply(ctmp, 2, max)
        ids3 <- ids3[!ids_rmv]

      } else {
        tmp <- NULL
      }
    }

    if (!is.null(tmp) && NCOL(tmp) > 0) {
      # Find highest mean(dcorU)
      tmp2 <- which.max(tmp["mean", ])
      res[k, "var"] <- colnames(tmp)[tmp2]
      res[k, c("median_dcorU", "range_dcorU")] <- as.numeric(tmp[, tmp2])
      res[k, "max_dcorU_to_prev"] <- max(
        predictor_info[ids3[tmp2], , "max_dcorU_to_prev"]
      )
    }
  }

  res2 <- res[complete.cases(res[, c("var", "median_dcorU")]), ]
  res2 <- res2[!duplicated(res2[, "var"]), ]

  vtmp <- as.vector(res2[order(-res2[, "median_dcorU"]), "var"])

  predictor_info[vtmp, , "pred_sel"] <- seq_along(vtmp)

  list(
    predictor_info = predictor_info,
    summary = res2,
    preds_selected = vtmp,
    kpreds_selected = vtmp[seq_len(min(nclusters, length(vtmp)))]
  )
}


plot_correlation_checks <- function(data, predictors, responses,
  path = ".", ftag = "", device = c("png", "pdf")) {

  device <- match.arg(device)

  dats_cor <- cor(data[, predictors], use = "pairwise.complete.obs")

  no_cor <- apply(dats_cor, 1, function(x) unname(which(!is.finite(x))))
  tmp <- lengths(no_cor)
  ids_no_cor <- unique(unlist(no_cor[tmp < ncol(dats_cor)]))

  cor_tmp <- if (length(ids_no_cor) > 0) {
    dats_cor[-ids_no_cor, -ids_no_cor]
  } else {
    dats_cor
  }

  ns <- length(predictors)

  fname <- file.path(
    path,
    paste0("Fig_PairwiseCorrelations_", ftag, ".", device)
  )

  if (device == "png") {
    grDevices::png(
      height = max(7, ns),
      width = max(7, ns),
      units = "in",
      res = 600,
      filename = fname
    )

  } else if (device == "pdf") {
    grDevices::pdf(
      height = max(7, ns),
      width = max(7, ns),
      file = fname
    )
  }

  tmp <- corrplot::corrplot(
  cor_tmp,
    addCoef.col = "black", number.cex = 0.5,
    order = "hclust", hclust.method = "ward.D2", addrect = 2,
    tl.cex = 0.5
  )

  dev.off()


  #--- Relationship with model responses
  panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }

  ns <- length(predictors) + length(responses)
  fname <- file.path(
    path,
    paste0("Fig_PairwiseRelationships_", ftag, ".", device)
  )

  if (device == "png") {
    grDevices::png(
      height = 1.5 * ns,
      width = 1.5 * ns,
      units = "in",
      res = 600,
      filename = fname
    )

  } else if (device == "pdf") {
    grDevices::pdf(
      height = 1.5 * ns,
      width = 1.5 * ns,
      file = fname
    )
  }

  pairs(data[, c(responses, predictors)],
    lower.panel = panel.cor,
    diag.panel = panel.hist,
    upper.panel = panel.smooth
  )

  dev.off()
}


} # endif not exists "SagebrushRegeneration_1_Functions"

SagebrushRegeneration_1_Functions <- TRUE
