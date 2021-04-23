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


#------------------------------------------------------------------------------#
#--- Research questions
#
# (i) examine the geographic patterns of big sagebrush regeneration
# probabilities that the two different models project under historical
# conditions and future climate scenarios
#
# (ii) quantify the robustness of model projections, e.g., the consistency among
# climate models in projected changes in regeneration for future time periods
#
# (iii) identify how model predictions for regeneration potential relate to
# environmental site characteristics like climate, soil moisture, and soils.

#------------------------------------------------------------------------------#



if (
  !get0(
    "SagebrushRegeneration_2_DescribeExperimentAnalysis",
    ifnotfound = FALSE
  )
) {

#--- Check package dependencies
pkgs <- c("remotes", "raster", "sp", "sf")
stopifnot(all(sapply(pkgs, requireNamespace)))

if (!requireNamespace("rSW2st")) {
  remotes::install_github("DrylandEcology/rSW2st", build_vignettes = TRUE)
}

if (!requireNamespace("rSW2analysis")) {
  remotes::install_github("DrylandEcology/rSW2analysis", build_vignettes = TRUE)
}


#--- Load custom functions
source("SagebrushRegeneration_1_Functions.R")


#--- Define and create folder paths
# 1_Simulation_Setup: files describing/defining simulation experiments
# 2_Simulation_Output: simulation output values required for analysis
# 3_Prepared_Data: intermediate values
# 4_Analysis_Results: final output, e.g., for manuscript and talks

dir_here <- normalizePath(".") # this is `Code/`
dir_prj <- file.path(dir_here, "..")

dir_prj_in <- file.path(dir_prj, "1_Simulation_Setup")
dir.create(dir_prj_in, recursive = TRUE, showWarnings = FALSE)

dir_prj_in2 <- file.path(dir_prj, "1_Additional_Setup")
dir.create(dir_prj_in2, recursive = TRUE, showWarnings = FALSE)

dir_prj_out <- file.path(dir_prj, "2_Simulation_Output")
dir.create(dir_prj_out, recursive = TRUE, showWarnings = FALSE)

dir_res_data <- file.path(dir_prj, "3_Prepared_Data")
dir.create(dir_res_data, recursive = TRUE, showWarnings = FALSE)

dir_res2 <- file.path(dir_prj, "4_Analysis_Results")
dir.create(dir_res2, recursive = TRUE, showWarnings = FALSE)


#--- The experiment is organized in five parts each by a different
# soil specification
tag_subprj <- c(
  "SoilSiteSpecific",
  "SoilSBmedian",
  "SoilSBqClLo",
  "SoilSBqSaLo",
  "SoilSBqSiLo"
)
default_subprj <- tag_subprj[1]

# geographic extent: index to `def_subsets` for producing output
used_subsets <- c(GISSM = 1, Shriver2018 = 2)
cfun <- "median" # central tendency function for ensembles across GCMs

vtag <- "v4"


#--- Describe this analysis
nc_att_global <- list(
  title = "Understanding future big sagebrush regeneration potential",
  version = "v20181023",
  source = paste(
    "SOILWAT2 (v4.2.0);",
    "rSOILWAT2 (v2.3.2);",
    "rSFSW2 (v3.1.2)"
  ),
  source_id = "SOILWAT2",
  further_info_url = "https://github.com/DrylandEcology/",
  source_type = "LAND",
  realm = "land",
  product = "model-output",
  grid = "native Alberts projection grid with NAD83 datum (337x257 eastxnorth)",
  grid_label = "gn",
  nominal_resolution = "10 km",
  contact = "dschlaepfer@usgs.gov",
  institution = "Southwest Biological Science Center, U.S. Geological Survey"
)


# USA Contiguous Albers Equal Area Conic USGS version
# +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs
# http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#_albers_equal_area
nc_att_crs <- list(crs_wkt = sf::st_crs("EPSG:6350")$Wkt)

nc_att_xy <- list(name = c("x", "y"))


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#--- LOAD DESCRIPTION OF SIMULATION EXPERIMENT
# derived from rSFSW2 object with script "SagebrushRegeneration_0_DataSharing.R"

fvars <- file.path(dir_prj_in2, "SFSW2_project_variables_metadata.csv")
vars_table <- read.csv(fvars, stringsAsFactors = FALSE)


#---Create a rudimentary simulation experiment description
SFSW2_prj_meta <- list()

SFSW2_prj_meta[["fnames_out"]] <- list(
  dbOutput = ""
)

SFSW2_prj_meta[["sim_space"]] <- list()


# Simulation grid
fname_gridspec <- file.path(dir_prj_in, "gridspec_fx_SOILWAT2_BSR_gn.nc")

# The following `rSW2st::read_netCDF()` will generate messages/warnings; if
# they are as follows, then the code should still work correctly:

# [1] "vobjtovarid4: error #F: I could not find the requsted var (or dimvar) in the file!"
# [1] "var (or dimvar) name: crs: x y"
# [1] "file name: 1_Simulation_Setup/gridspec_fx_SOILWAT2_BSR_gn.nc"
# Warning messages:
#   1: In showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj = prefer_proj) :
#   Discarded datum NAD83 (National Spatial Reference System 2011) in Proj4 definition
# 2: In showSRID(uprojargs, format = "PROJ", multiline = "NO", prefer_proj = prefer_proj) :
#   Discarded datum Unknown based on GRS80 ellipsoid in Proj4 definition

SFSW2_prj_meta[["sim_space"]][["sim_raster"]] <- rSW2st::read_netCDF(
  x = fname_gridspec,
  var = "gridspec",
  method = "raster"
)
SFSW2_prj_meta[["sim_space"]][["sim_crs"]] <- raster::crs(
  x = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
)


# Simulated cells based on simulation grid
SFSW2_prj_meta[["sim_space"]][["crs_sites"]] <- raster::crs("EPSG:4326")

SFSW2_prj_meta[["sim_space"]][["run_sites"]] <- {
  tmp <- data.frame(
    sp::coordinates(SFSW2_prj_meta[["sim_space"]][["sim_raster"]]),
    raster::as.data.frame(SFSW2_prj_meta[["sim_space"]][["sim_raster"]])
  )
  tmp <- tmp[order(tmp[, 3], na.last = NA), 1:2]
  rownames(tmp) <- NULL
  colnames(tmp) <- c("X_WGS84", "Y_WGS84")
  tmp <- sp::SpatialPoints(
    coords = tmp,
    proj4string = SFSW2_prj_meta[["sim_space"]][["sim_crs"]]
  )

  sp::spTransform(
    x = tmp,
    CRSobj = SFSW2_prj_meta[["sim_space"]][["crs_sites"]]
  )
}


# Number of simulated units/runs
SFSW2_prj_meta[["sim_size"]] <- list(
  runsN_sites = length(SFSW2_prj_meta[["sim_space"]][["run_sites"]])
)

# Simulated scenarios and time periods
SFSW2_prj_meta[["sim_time"]] <- list(
  startyr = 1980,
  endyr = 2010,
  future_yrs = structure(
    c(40, 90, 2020, 2070, 2050, 2099),
    .Dim = 2:3,
    .Dimnames = list(
      c("d40yrs", "d90yrs"),
      c("delta", "DSfut_startyr", "DSfut_endyr")
    )
  )
)

SFSW2_prj_meta[["sim_scens"]] <- list(
  ambient = "Current",

  reqCSs = reqCSs <- c("RCP45", "RCP85"),

  reqCSsPerM = {
      temp <- rep(list(reqCSs), length = 37)
      temp[c(20, 22)] <- "RCP45" # GISS-E2-H-CC and GISS-E2-R-CC
      temp
    },

  reqMs = c(
    "ACCESS1-0", "ACCESS1-3", "bcc-csm1-1", "bcc-csm1-1-m", "BNU-ESM",
    "CanESM2", "CCSM4", "CESM1-BGC", "CESM1-CAM5", "CMCC-CM", "CNRM-CM5",
    "CSIRO-Mk3-6-0", "EC-EARTH", "FGOALS-g2", "FGOALS-s2", "FIO-ESM",
    "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H-CC", "GISS-E2-R",
    "GISS-E2-R-CC", "HadGEM2-AO", "HadGEM2-CC", "HadGEM2-ES", "inmcm4",
    "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC-ESM",
    "MIROC-ESM-CHEM", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3",
    "NorESM1-M", "NorESM1-ME"
  ),

  DeltaStr_yrs = c("d40yrs", "d60yrs"),

  method_DS = "hybrid-delta-3mod",

  id = c(
    "Current",
    "hybrid-delta-3mod.d40yrs.RCP45.ACCESS1-0", "hybrid-delta-3mod.d90yrs.RCP45.ACCESS1-0",
    "hybrid-delta-3mod.d40yrs.RCP45.ACCESS1-3", "hybrid-delta-3mod.d90yrs.RCP45.ACCESS1-3",
    "hybrid-delta-3mod.d40yrs.RCP45.bcc-csm1-1", "hybrid-delta-3mod.d90yrs.RCP45.bcc-csm1-1",
    "hybrid-delta-3mod.d40yrs.RCP45.bcc-csm1-1-m", "hybrid-delta-3mod.d90yrs.RCP45.bcc-csm1-1-m",
    "hybrid-delta-3mod.d40yrs.RCP45.BNU-ESM", "hybrid-delta-3mod.d90yrs.RCP45.BNU-ESM",
    "hybrid-delta-3mod.d40yrs.RCP45.CanESM2", "hybrid-delta-3mod.d90yrs.RCP45.CanESM2",
    "hybrid-delta-3mod.d40yrs.RCP45.CCSM4", "hybrid-delta-3mod.d90yrs.RCP45.CCSM4",
    "hybrid-delta-3mod.d40yrs.RCP45.CESM1-BGC", "hybrid-delta-3mod.d90yrs.RCP45.CESM1-BGC",
    "hybrid-delta-3mod.d40yrs.RCP45.CESM1-CAM5", "hybrid-delta-3mod.d90yrs.RCP45.CESM1-CAM5",
    "hybrid-delta-3mod.d40yrs.RCP45.CMCC-CM", "hybrid-delta-3mod.d90yrs.RCP45.CMCC-CM",
    "hybrid-delta-3mod.d40yrs.RCP45.CNRM-CM5", "hybrid-delta-3mod.d90yrs.RCP45.CNRM-CM5",
    "hybrid-delta-3mod.d40yrs.RCP45.CSIRO-Mk3-6-0", "hybrid-delta-3mod.d90yrs.RCP45.CSIRO-Mk3-6-0",
    "hybrid-delta-3mod.d40yrs.RCP45.EC-EARTH", "hybrid-delta-3mod.d90yrs.RCP45.EC-EARTH",
    "hybrid-delta-3mod.d40yrs.RCP45.FGOALS-g2", "hybrid-delta-3mod.d90yrs.RCP45.FGOALS-g2",
    "hybrid-delta-3mod.d40yrs.RCP45.FGOALS-s2", "hybrid-delta-3mod.d90yrs.RCP45.FGOALS-s2",
    "hybrid-delta-3mod.d40yrs.RCP45.FIO-ESM", "hybrid-delta-3mod.d90yrs.RCP45.FIO-ESM",
    "hybrid-delta-3mod.d40yrs.RCP45.GFDL-CM3", "hybrid-delta-3mod.d90yrs.RCP45.GFDL-CM3",
    "hybrid-delta-3mod.d40yrs.RCP45.GFDL-ESM2G", "hybrid-delta-3mod.d90yrs.RCP45.GFDL-ESM2G",
    "hybrid-delta-3mod.d40yrs.RCP45.GFDL-ESM2M", "hybrid-delta-3mod.d90yrs.RCP45.GFDL-ESM2M",
    "hybrid-delta-3mod.d40yrs.RCP45.GISS-E2-H-CC", "hybrid-delta-3mod.d90yrs.RCP45.GISS-E2-H-CC",
    "hybrid-delta-3mod.d40yrs.RCP45.GISS-E2-R", "hybrid-delta-3mod.d90yrs.RCP45.GISS-E2-R",
    "hybrid-delta-3mod.d40yrs.RCP45.GISS-E2-R-CC", "hybrid-delta-3mod.d90yrs.RCP45.GISS-E2-R-CC",
    "hybrid-delta-3mod.d40yrs.RCP45.HadGEM2-AO", "hybrid-delta-3mod.d90yrs.RCP45.HadGEM2-AO",
    "hybrid-delta-3mod.d40yrs.RCP45.HadGEM2-CC", "hybrid-delta-3mod.d90yrs.RCP45.HadGEM2-CC",
    "hybrid-delta-3mod.d40yrs.RCP45.HadGEM2-ES", "hybrid-delta-3mod.d90yrs.RCP45.HadGEM2-ES",
    "hybrid-delta-3mod.d40yrs.RCP45.inmcm4", "hybrid-delta-3mod.d90yrs.RCP45.inmcm4",
    "hybrid-delta-3mod.d40yrs.RCP45.IPSL-CM5A-LR", "hybrid-delta-3mod.d90yrs.RCP45.IPSL-CM5A-LR",
    "hybrid-delta-3mod.d40yrs.RCP45.IPSL-CM5A-MR", "hybrid-delta-3mod.d90yrs.RCP45.IPSL-CM5A-MR",
    "hybrid-delta-3mod.d40yrs.RCP45.IPSL-CM5B-LR", "hybrid-delta-3mod.d90yrs.RCP45.IPSL-CM5B-LR",
    "hybrid-delta-3mod.d40yrs.RCP45.MIROC-ESM", "hybrid-delta-3mod.d90yrs.RCP45.MIROC-ESM",
    "hybrid-delta-3mod.d40yrs.RCP45.MIROC-ESM-CHEM", "hybrid-delta-3mod.d90yrs.RCP45.MIROC-ESM-CHEM",
    "hybrid-delta-3mod.d40yrs.RCP45.MIROC5", "hybrid-delta-3mod.d90yrs.RCP45.MIROC5",
    "hybrid-delta-3mod.d40yrs.RCP45.MPI-ESM-LR", "hybrid-delta-3mod.d90yrs.RCP45.MPI-ESM-LR",
    "hybrid-delta-3mod.d40yrs.RCP45.MPI-ESM-MR", "hybrid-delta-3mod.d90yrs.RCP45.MPI-ESM-MR",
    "hybrid-delta-3mod.d40yrs.RCP45.MRI-CGCM3", "hybrid-delta-3mod.d90yrs.RCP45.MRI-CGCM3",
    "hybrid-delta-3mod.d40yrs.RCP45.NorESM1-M", "hybrid-delta-3mod.d90yrs.RCP45.NorESM1-M",
    "hybrid-delta-3mod.d40yrs.RCP45.NorESM1-ME", "hybrid-delta-3mod.d90yrs.RCP45.NorESM1-ME",
    "hybrid-delta-3mod.d40yrs.RCP85.ACCESS1-0", "hybrid-delta-3mod.d90yrs.RCP85.ACCESS1-0",
    "hybrid-delta-3mod.d40yrs.RCP85.ACCESS1-3", "hybrid-delta-3mod.d90yrs.RCP85.ACCESS1-3",
    "hybrid-delta-3mod.d40yrs.RCP85.bcc-csm1-1", "hybrid-delta-3mod.d90yrs.RCP85.bcc-csm1-1",
    "hybrid-delta-3mod.d40yrs.RCP85.bcc-csm1-1-m", "hybrid-delta-3mod.d90yrs.RCP85.bcc-csm1-1-m",
    "hybrid-delta-3mod.d40yrs.RCP85.BNU-ESM", "hybrid-delta-3mod.d90yrs.RCP85.BNU-ESM",
    "hybrid-delta-3mod.d40yrs.RCP85.CanESM2", "hybrid-delta-3mod.d90yrs.RCP85.CanESM2",
    "hybrid-delta-3mod.d40yrs.RCP85.CCSM4", "hybrid-delta-3mod.d90yrs.RCP85.CCSM4",
    "hybrid-delta-3mod.d40yrs.RCP85.CESM1-BGC", "hybrid-delta-3mod.d90yrs.RCP85.CESM1-BGC",
    "hybrid-delta-3mod.d40yrs.RCP85.CESM1-CAM5", "hybrid-delta-3mod.d90yrs.RCP85.CESM1-CAM5",
    "hybrid-delta-3mod.d40yrs.RCP85.CMCC-CM", "hybrid-delta-3mod.d90yrs.RCP85.CMCC-CM",
    "hybrid-delta-3mod.d40yrs.RCP85.CNRM-CM5", "hybrid-delta-3mod.d90yrs.RCP85.CNRM-CM5",
    "hybrid-delta-3mod.d40yrs.RCP85.CSIRO-Mk3-6-0", "hybrid-delta-3mod.d90yrs.RCP85.CSIRO-Mk3-6-0",
    "hybrid-delta-3mod.d40yrs.RCP85.EC-EARTH", "hybrid-delta-3mod.d90yrs.RCP85.EC-EARTH",
    "hybrid-delta-3mod.d40yrs.RCP85.FGOALS-g2", "hybrid-delta-3mod.d90yrs.RCP85.FGOALS-g2",
    "hybrid-delta-3mod.d40yrs.RCP85.FGOALS-s2", "hybrid-delta-3mod.d90yrs.RCP85.FGOALS-s2",
    "hybrid-delta-3mod.d40yrs.RCP85.FIO-ESM", "hybrid-delta-3mod.d90yrs.RCP85.FIO-ESM",
    "hybrid-delta-3mod.d40yrs.RCP85.GFDL-CM3", "hybrid-delta-3mod.d90yrs.RCP85.GFDL-CM3",
    "hybrid-delta-3mod.d40yrs.RCP85.GFDL-ESM2G", "hybrid-delta-3mod.d90yrs.RCP85.GFDL-ESM2G",
    "hybrid-delta-3mod.d40yrs.RCP85.GFDL-ESM2M", "hybrid-delta-3mod.d90yrs.RCP85.GFDL-ESM2M",
    "hybrid-delta-3mod.d40yrs.RCP85.GISS-E2-R", "hybrid-delta-3mod.d90yrs.RCP85.GISS-E2-R",
    "hybrid-delta-3mod.d40yrs.RCP85.HadGEM2-AO", "hybrid-delta-3mod.d90yrs.RCP85.HadGEM2-AO",
    "hybrid-delta-3mod.d40yrs.RCP85.HadGEM2-CC", "hybrid-delta-3mod.d90yrs.RCP85.HadGEM2-CC",
    "hybrid-delta-3mod.d40yrs.RCP85.HadGEM2-ES", "hybrid-delta-3mod.d90yrs.RCP85.HadGEM2-ES",
    "hybrid-delta-3mod.d40yrs.RCP85.inmcm4", "hybrid-delta-3mod.d90yrs.RCP85.inmcm4",
    "hybrid-delta-3mod.d40yrs.RCP85.IPSL-CM5A-LR", "hybrid-delta-3mod.d90yrs.RCP85.IPSL-CM5A-LR",
    "hybrid-delta-3mod.d40yrs.RCP85.IPSL-CM5A-MR", "hybrid-delta-3mod.d90yrs.RCP85.IPSL-CM5A-MR",
    "hybrid-delta-3mod.d40yrs.RCP85.IPSL-CM5B-LR", "hybrid-delta-3mod.d90yrs.RCP85.IPSL-CM5B-LR",
    "hybrid-delta-3mod.d40yrs.RCP85.MIROC-ESM", "hybrid-delta-3mod.d90yrs.RCP85.MIROC-ESM",
    "hybrid-delta-3mod.d40yrs.RCP85.MIROC-ESM-CHEM", "hybrid-delta-3mod.d90yrs.RCP85.MIROC-ESM-CHEM",
    "hybrid-delta-3mod.d40yrs.RCP85.MIROC5", "hybrid-delta-3mod.d90yrs.RCP85.MIROC5",
    "hybrid-delta-3mod.d40yrs.RCP85.MPI-ESM-LR", "hybrid-delta-3mod.d90yrs.RCP85.MPI-ESM-LR",
    "hybrid-delta-3mod.d40yrs.RCP85.MPI-ESM-MR", "hybrid-delta-3mod.d90yrs.RCP85.MPI-ESM-MR",
    "hybrid-delta-3mod.d40yrs.RCP85.MRI-CGCM3", "hybrid-delta-3mod.d90yrs.RCP85.MRI-CGCM3",
    "hybrid-delta-3mod.d40yrs.RCP85.NorESM1-M", "hybrid-delta-3mod.d90yrs.RCP85.NorESM1-M",
    "hybrid-delta-3mod.d40yrs.RCP85.NorESM1-ME", "hybrid-delta-3mod.d90yrs.RCP85.NorESM1-ME"
  )
)


# Transform sites = cell centers to same coordinate system as simulation raster
SFSW2_prj_meta[["sim_space"]][["run_sites_orig"]] <-
  SFSW2_prj_meta[["sim_space"]][["run_sites"]]
SFSW2_prj_meta[["sim_space"]][["run_sites"]] <- sf::st_transform(
  rSW2st::as_points(
    SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    to_class = "sf"
  ),
  crs = nc_att_crs[["crs_wkt"]]
)



#------------------------------------------------------------------------------#
#--- DEFINE SCENARIOS AND CLIMATE CONDITIONS OF SIMULATION EXPERIMENT

sc_hist <- SFSW2_prj_meta[["sim_scens"]][["ambient"]]
reqCSs <- SFSW2_prj_meta[["sim_scens"]][["reqCSs"]]
reqACSs <- c(sc_hist, reqCSs)

deltas_byCS <- {
  temp <- rep(sc_hist, length(reqCSs))
  names(temp) <- reqCSs
  temp
}

req_scens <- unname(unlist(sapply(
  X = reqCSs,
  FUN = function(x) grep(x, SFSW2_prj_meta[["sim_scens"]][["id"]], value = TRUE)
)))

req_Downs <- unique(sapply(
  X = strsplit(x = req_scens, split = ".", fixed = TRUE),
  FUN = function(x) x[[1]]
))
req_dTime <- c(
  "d0yrs",
  unique(sapply(
    X = strsplit(req_scens, split = ".", fixed = TRUE),
    FUN = function(x) x[[2]]
  ))
)

reqGCMs_perCS <- lapply(
  X = reqCSs,
  FUN = function(cs) {
    temp <- sapply(
      X = SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]],
      FUN = function(x) any(x == cs)
    )
    SFSW2_prj_meta[["sim_scens"]][["reqMs"]][temp]
  }
)
names(reqGCMs_perCS) <- reqCSs
reqMs <- sort(unique(unlist(reqGCMs_perCS)))

if (FALSE) {
  tmp <- sapply(
    X = reqCSs,
    FUN = function(cs) {
      tmp <- sapply(
        X = SFSW2_prj_meta[["sim_scens"]][["reqCSsPerM"]],
        FUN = function(x) any(x == cs)
      )
    }
  )
  colnames(tmp) <- reqCSs

  out <- data.frame(
    GCM = SFSW2_prj_meta[["sim_scens"]][["reqMs"]],
    tmp
  )
}

# GCMs with full RCP coverage
reqGCMs_full <- reqGCMs_perCS[["RCP45"]]
for (sc in reqCSs[-1]) {
  reqGCMs_full <- intersect(reqGCMs_full, reqGCMs_perCS[[sc]])
}

# Simulation scenarios
temp <- c(
  sc_hist,
  unlist(sapply(
    X = reqCSs,
    FUN = function(cs) sapply(
      X = reqGCMs_perCS[[cs]],
      FUN = function(x) grep(paste0(cs, ".", x, "$"), req_scens, value = TRUE)
    )
  ))
)
temp <- SFSW2_prj_meta[["sim_scens"]][["id"]] %in% temp
sim_scens_id <- SFSW2_prj_meta[["sim_scens"]][["id"]][temp]

temp <- c(
  sc_hist,
  unlist(sapply(
    X = reqCSs,
    FUN = function(cs) sapply(
      X = reqGCMs_full,
      FUN = function(x) grep(paste0(cs, ".", x, "$"), req_scens, value = TRUE)
    )
  ))
)
temp <- SFSW2_prj_meta[["sim_scens"]][["id"]] %in% temp
sim_scens_full <- SFSW2_prj_meta[["sim_scens"]][["id"]][temp]

# Ensemble units
ensemble_units <- apply(
  X = expand.grid(req_Downs, req_dTime[-1], reqCSs),
  MARGIN = 1,
  FUN = paste,
  collapse = "."
)

# Number of years in simulation period
nyrs <- SFSW2_prj_meta[["sim_time"]][["endyr"]] -
  SFSW2_prj_meta[["sim_time"]][["startyr"]] + 1

# Time periods in calendar years
req_CalendarPeriods <- c(
  paste0(
    SFSW2_prj_meta[["sim_time"]][["startyr"]],
    "-",
    SFSW2_prj_meta[["sim_time"]][["endyr"]]
  ),
  paste0(
    SFSW2_prj_meta[["sim_time"]][["future_yrs"]][, "DSfut_startyr"],
    "-",
    SFSW2_prj_meta[["sim_time"]][["future_yrs"]][, "DSfut_endyr"]
  )
)
names(req_CalendarPeriods) <- req_dTime


#------------------------------------------------------------------------------#

#--- AREAL EXTENT OF RASTER CELLS
fname_areacella <- file.path(dir_prj_in, "areacella_fx_SOILWAT2_BSR_gn.nc")
tmp <- rSW2st::read_netCDF(
  x = fname_areacella,
  var = "areacella",
  method = "xy_subset",
  xy_names = nc_att_xy[["name"]],
  locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]]
)
cells <- matrix(
  data = c(tmp / 1000, tmp / max(tmp)),
  ncol = 2,
  dimnames = list(NULL, c("km2", "rel"))
)




#------ STUDY AREA ------
#------ Mask of GISSM study area ------
fname_landmask <- file.path(dir_prj_in, "sftlf_fx_SOILWAT2_BSR_gn.nc")
tmp <- rSW2st::read_netCDF(
  x = fname_landmask,
  var = "maskGISSM",
  method = "xy_subset",
  xy_names = nc_att_xy[["name"]],
  locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]]
)
tmp[is.na(tmp)] <- 0
id_mask_GISSM <- as.logical(tmp)



#------ Mask of Shriver2018 study area ------
# Roughly a Great Basin extent (based on coverage in Shriver et al. 2018 GBC):
# 10.1.3 Northern Basin and Range
# 10.1.5 Central Basin and Range
# 10.1.8 Snake River Plain
tmp <- rSW2st::read_netCDF(
  x = fname_landmask,
  var = "maskShriver2018",
  method = "xy_subset",
  xy_names = nc_att_xy[["name"]],
  locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]]
)
tmp[is.na(tmp)] <- 0
id_mask_Shriver2018 <- as.logical(tmp)




#--- EPA/CEC Ecoregions Level III
fname_ecoregs <- file.path(
  "~",
  "BigData",
  "Ecoregions_CEC",
  "NA_CEC_Eco_Level3_WesternNA_AEA",
  "cec_ecoregions_iii_sp.rds"
)

if (!file.exists(fname_ecoregs)) {
  stop(
    "Download EPA ecoregions Level III ",
    "from https://www.epa.gov/eco-research/ecoregions-north-america, i.e., ",
    "ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip"
  )

} else {
  pecoregs <- readRDS(fname_ecoregs)
}



fname_ecoregs <- file.path(
  dir_prj_in,
  "EPAecoregionsL3_fx_SOILWAT2_BSR_gn.nc"
)
tmp <- rSW2st::read_netCDF(
  x = fname_ecoregs,
  var = "ecoregionsL3",
  method = "xy_subset",
  xy_names = nc_att_xy[["name"]],
  locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]]
)
id_ecoregs2 <- pecoregs[["NA_L3KEY"]][tmp]



# Create polygon shapefile of relevant ecoregions
ids_use <- pecoregs[["NA_L3KEY"]] %in% sort(unique(id_ecoregs2))
pecoregs2 <- sp::SpatialPolygonsDataFrame(
  Sr = rgeos::gSimplify(
    spgeom = pecoregs[ids_use, ],
    tol = 10,
    topologyPreserve = TRUE
  ),
  data = as(pecoregs, "data.frame")[ids_use, ]
)



# `spoly_mask_Shriver2018` is only used to add outline to figures
fname_spoly_mask_Shriver2018_rds <- file.path(
  dir_prj_in2,
  "spoly_mask_Shriver2018.rds"
)

if (file.exists(fname_spoly_mask_Shriver2018_rds)) {
  spoly_mask_Shriver2018 <- readRDS(fname_spoly_mask_Shriver2018_rds)

} else {
  if (exists("pecoregs")) {
    def_gbish <- c("10.1.3", "10.1.5", "10.1.8")
    gbish <- pecoregs[pecoregs[["NA_L3CODE"]] %in% def_gbish, ]
    gbish <- maptools::unionSpatialPolygons(gbish, rep(1, length(def_gbish)))

    spoly_mask_Shriver2018 <- sp::spTransform(
      x = gbish,
      raster::crs(SFSW2_prj_meta[["sim_space"]][["sim_raster"]])
    )

    saveRDS(spoly_mask_Shriver2018, file = fname_spoly_mask_Shriver2018_rds)

  } else {

    # Attempt to re-create `spoly_mask_Shriver2018`
    # from netCDF instead of non-existing `pecoregs` polygon-shapefile
    r <- rSW2st::read_netCDF(
      x = fname_landmask,
      var = "maskShriver2018",
      method = "raster"
    )

    spoly_mask_Shriver2018 <- raster::rasterToPolygons(
      r,
      fun = function(x) x == 1,
      dissolve = TRUE
    )


    if (requireNamespace("spatialEco")) {
      spoly_mask_Shriver2018 <- spatialEco::remove.holes(spoly_mask_Shriver2018)

    } else {
      message(
        "It would be desirable to remove holes from ",
        "`spoly_mask_Shriver2018`, but package `spatialEco` is not available."
      )
    }
  }
}



#--- Define different geographic subsets for analysis
def_subsets <- list(
  BigSage = id_mask_GISSM,
  GreatBasin = id_mask_Shriver2018
)


extent_subsets <- list()

for (ks in seq_along(def_subsets)) {
  tag_subset <- names(def_subsets)[ks]

  tmp <- def_subsets[[ks]]
  tmp[!tmp] <- NA
  r <- rSW2st::create_raster_from_variables(
    data = as.integer(tmp),
    site_locations = SFSW2_prj_meta[["sim_space"]][["run_sites"]],
    grid = SFSW2_prj_meta[["sim_space"]][["sim_raster"]]
  )

  r <- try(raster::trim(r), silent = TRUE)

  if (!inherits(r, "try-error")) {
    fextend <- switch(tag_subset, BigSage = 1.05, GreatBasin = 1.1, 1)
    extent_subsets[[tag_subset]] <- fextend * raster::extent(r)
  }
}





#--- Country borders, US states, Canadian provinces, and Mexican states
fname_borders <- file.path(
  "~",
  "BigData",
  "Maps",
  "Global",
  "NaturalEarth",
  "NaturalEarth_v410_WesternNA",
  "ne_borders_NA_sp_aea.rds"
)

if (file.exists(fname_borders)) {
  pborders <- readRDS(fname_borders)

} else {
  stop(
    "Download Natural Earth ",
    "from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/"
  )
}

} # endif not exists "SagebrushRegeneration_2_DescribeExperimentAnalysis"

SagebrushRegeneration_2_DescribeExperimentAnalysis <- TRUE


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
