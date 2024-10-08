## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  fig.show = "asis",
  fig.align = "centre",
  results = "markup",
  comment = "#>",
  fig.path = ""
)


## ----setup--------------------------------------------------------------------
library(poems)
library(raster)
library(sf)
library(scales)
library(stringi) # for randomly generating file names.

# function to round to any arbitrary value
round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}


## ----message = FALSE, warning=FALSE, results='markup', fig.show='asis',fig.align = "center", fig.width = 4, fig.height = 4----
# Region raster
data(tasmania_raster)
tasmania_raster

# Equal area projection
tasPrj <- 'PROJCS["Tasmania_Lambert_Azimuthal",
                 GEOGCS["GCS_WGS_1984",
                        DATUM["D_WGS_1984",
                              SPHEROID["WGS_1984",6378137.0,298.257223563]],
                        PRIMEM["Greenwich",0.0],
                        UNIT["Degree",0.0174532925199433]],
                 PROJECTION["Lambert_Azimuthal_Equal_Area"],
                 PARAMETER["False_Easting",0.0],
                 PARAMETER["False_Northing",0.0],
                 PARAMETER["Central_Meridian",147],
                 PARAMETER["Latitude_Of_Origin",-42.2],
                 UNIT["Meter",1.0]]'

# Template raster to project to
tempExt <- projectExtent(tasmania_raster, tasPrj)
res(tempExt) <- 10000 # 10 km resolution
tempExt

# Project the region
tasmania_raster <- projectRaster(tasmania_raster, tempExt,
  method = "ngb"
)
plot(tasmania_raster,
  main = "Tasmania raster",
  legend = FALSE,
  col = "#2E8B57", colNA = "grey75"
)


## -----------------------------------------------------------------------------
# Tasmania study region (735 non-NA cells stored in the order shown) #
region <- Region$new(template_raster = tasmania_raster)
region$region_raster

# Establish HS template and starting location #
# This will be our initial introduction point
int_ll <- sf_project(
  from = "EPSG:4326",
  to = tasPrj,
  pts = cbind(146.44, -41.18)
)
int_point <- region$region_indices[
  which(region$region_indices ==
    cellFromXY(tasmania_raster, xy = int_ll))
]

# row which corresponds to initial introduction site
int_index <- which(region$region_indices == int_point) # 114

# plot of region, and introduction locations
plot(region$region_raster,
  main = "Tasmanian study region (cell indices)",
  colNA = "grey75",
  addfun = function() {
    points(xyFromCell(region$region_raster, int_point), col = "red", pch = 16)
  }
)


## -----------------------------------------------------------------------------
# read in the land-use modifier
data(tasmania_modifier)

plot(tasmania_modifier,
  zlim = c(0, 1), colNA = "grey75",
  col = hcl.colors(100, "RdYlGn")
)

# Habitat suitability
data(thylacine_hs_raster)
hs_raster <- projectRaster(thylacine_hs_raster, region$region_raster, method = "bilinear")
hs_raster <- stretch(hs_raster, minv = 0, maxv = 1)
hs_raster

# initial_hs needed for generator
initial_hs <- hs_raster <- stack(replicate(n = nlayers(tasmania_modifier), hs_raster))


## -----------------------------------------------------------------------------
# static HS for the moment
## capacity generator will make it temporally dynamic
plot(hs_raster,
  zlim = c(0, 1), colNA = "grey75",
  col = hcl.colors(100, "RdYlGn"),
  addfun = function() {
    points(xyFromCell(region$region_raster, int_point), pch = 16, cex = 0.5)
  }
)


## ----message = FALSE----------------------------------------------------------
# Distance-based environmental correlation (via a compacted Cholesky decomposition)
env_corr <- SpatialCorrelation$new(
  region = region,
  amplitude = 0.496,
  breadth = 80,
  distance_scale = 1000
)
env_corr$calculate_compact_decomposition(decimals = 4)


## -----------------------------------------------------------------------------
# allow growth rates to vary by region using IBRA regions
# Tasmania study Interim Bioregionalisation of Australia (IBRA) bioregion cell distribution
data(tasmania_ibra_raster)
ibra_raster <- projectRaster(tasmania_ibra_raster, region$region_raster, method = "ngb")
plot(ibra_raster,
  colNA = "grey75",
  breaks = seq(1, 9, 1),
  main = "IBRA regions of Tasmania",
  col = hcl.colors(10, "Lajolla")
)

data(tasmania_ibra_data)
tasmania_ibra_data

# Calculate cell indices and counts for IBRA bioregions
ibra_indices <- lapply(
  as.list(tasmania_ibra_data$index),
  function(i) {
    which(ibra_raster[region$region_indices] == i)
  }
)
str(ibra_indices)

ibra_polygons <- rasterToPolygons(ibra_raster, dissolve = TRUE, na.rm = TRUE)
ibra_polygons@data <- merge(ibra_polygons@data, tasmania_ibra_data,
  by.x = "layer", by.y = "index"
)
ibra_polygons

plot(ibra_polygons, col = hcl.colors(9, "Lajolla"), border = "black")
text(ibra_polygons, labels = "abbr", cex = 1.2, halo = TRUE)

rmax_regional <- ibra_raster

# seed is set to keep example results constant
{
  set.seed(27)
  rmax <- round(rlnorm(9, 0.94, 0.3), 1)
}
for (val in 1:9) {
  rmax_regional[rmax_regional == val] <- rmax[val]
}
plot(rmax_regional,
  colNA = "grey75",
  legend = FALSE, main = "regional growth rates",
  zlim = range(rmax),
  addfun = function() {
    plot(ibra_polygons,
      border = hcl.colors(9, "Lajolla"),
      col = NA, add = TRUE
    )
    text(ibra_polygons, labels = rmax, halo = TRUE)
  },
  col = hcl.colors(100, "Zissou")
)

# set upper and lower growth rates per region
ibra_rmax <- cbind(tasmania_ibra_data,
  rmax_lower = round(rmax * 0.6, 2),
  rmax_mean = round(rmax, 2),
  rmax_upper = round(rmax / 0.75, 2)
)
ibra_rmax


## ----echo=TRUE----------------------------------------------------------------
# set up translocation locations and order
intro_trans_ll <- sf_project(
  from = "EPSG:4326",
  to = tasPrj,
  pts = cbind(
    c(148.01, 144.7, 147.9, 148.27, 145.24),
    c(-40.8, -40.7, -43.2, -42.02, -42.3)
  )
)
intro_trans_ll
intro_trans_point <- region$region_indices[which(region$region_indices %in%
  cellFromXY(region$region_raster,
    xy = intro_trans_ll
  ))]
intro_trans_point <- intro_trans_point[-1]
intro_cells <- intro_trans_point
intro_cells

intro_times <- c(2, 3, 6, 8)

# Introduction times and locations
cbind(intro_times, intro_cells)

plot(region$region_raster,
  main = "Introduction sites",
  col = hcl.colors(100, "Lajolla"),
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
    points(xyFromCell(region$region_raster, intro_cells),
      pch = 16,
      cex = 1.5, col = c("darkgreen", "blue2", "black", "goldenrod")
    )
    points(region$coordinates[which(region$region_indices %in% intro_cells), ],
      col = "firebrick", cex = 1.5, lwd = 2
    )
  }
)


## -----------------------------------------------------------------------------
# User-defined translocation function (list-nested) and alias ####

translocation <- list(

  # Function parameters (passed to function in params list)
  intro_cells = intro_cells, # cells where pops are introduced
  intro_timesteps = intro_times, # timesteps when introduced
  trans_n = 50, # translocated abundances. If not provided by LHS == 50
  region_indices = region$region_indices,

  # Function definition
  translocation_function = function(params) {
    # Unpack parameters (used at every time step)
    intro_cells <- params$intro_cells
    intro_timesteps <- params$intro_timesteps
    simulator <- params$simulator
    stages <- params$stages
    populations <- params$populations
    abundances <- params$abundance
    region_indices <- params$region_indices
    tm <- params$tm # timestep
    sa <- params$stage_abundance
    trans_n <- params$trans_n
    # if introduction at timestep, introduce pops
    if (tm %in% intro_timesteps) {
      # take stage abundance at timestep
      new_sa <- array(sa, c(stages, populations))
      # identifies location of introduction
      trans_loc <- which(region_indices == intro_cells[which(intro_timesteps == tm)])
      # add n individuals regardless of K
      new_sa[trans_loc] <- new_sa[trans_loc] + trans_n
      return(new_sa)
    } else {
      # else return pops as they are
      new_sa <- array(sa, c(stages, populations))
      return(new_sa)
    }
  }
)
translocation_aliases <- list(
  intro_cells = "translocation$intro_cells",
  intro_times = "translocation$intro_timesteps",
  trans_n = "translocation$trans_n",
  region_indices = "translocation$region_indices"
)


## -----------------------------------------------------------------------------
# Build a Rmax generator based on sampled IBRA Rmax range quantile
rmax_gen <- Generator$new(
  description = "Rmax",
  spatial_correlation = env_corr,
  generate_rasters = FALSE,
  ibra_data_rmax = ibra_rmax,
  ibra_indices = ibra_indices,
  region_cells = region$region_cells,
  inputs = c("rmax_quantile"),
  outputs = c("growth_rate_max"),
  generative_requirements = list(growth_rate_max = "function")
)

# growth_rate_max template
rmax_gen$add_function_template(
  "growth_rate_max",
  function_def = function(params) {
    growth_rate_max <- array(0, params$region_cells)
    for (i in 1:nrow(params$ibra_data_rmax)) {
      growth_rate_max[params$ibra_indices[[i]]] <-
        stats::qunif(params$rmax_quantile,
          min = params$ibra_data_rmax$rmax_lower[i],
          max = params$ibra_data_rmax$rmax_upper[i]
        )
    }
    return(growth_rate_max)
  },
  call_params = c("ibra_data_rmax", "ibra_indices", "region_cells", "rmax_quantile")
)

# test rmax generator at median values
rmax_gen_ex <- rmax_gen$generate(input_values = list(rmax_quantile = 0.5))
rmax_regional[region$region_indices] <- rmax_gen_ex$growth_rate_max
plot(rmax_regional,
  main = "median regional rmax",
  col = hcl.colors(100),
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
  }
)


## -----------------------------------------------------------------------------
# Test multiple quantiles
test_rmax <- lapply(seq(0, 1, 0.1), function(i) {
  region$raster_from_values(rmax_gen$generate(input_values = list(rmax_quantile = i))$growth_rate_max)
})
test_rmax <- stack(test_rmax)
names(test_rmax) <- paste0("Q", seq(0, 1, 0.1))

# plot
plot(test_rmax,
  colNA = "grey75",
  legend = TRUE,
  zlim = c(
    min(values(test_rmax), na.rm = TRUE),
    max(values(test_rmax), na.rm = TRUE)
  ),
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
  },
  col = hcl.colors(100)
)


## -----------------------------------------------------------------------------
# Dispersal generator ####
# Set for veriable mean distance, max hard-coded at 150
dispersal_gen <- DispersalGenerator$new(
  region = region,
  spatial_correlation = env_corr,
  generate_rasters = FALSE,
  dispersal_max_distance = 150,
  distance_classes = seq(5, 150, by = 10),
  distance_scale = 1000, # in km
  dispersal_friction = DispersalFriction$new(), # modify coastline distances
  inputs = c("dispersal_p", "dispersal_b"), # proportion and average distance
  decimals = 4
)
dispersal_gen$calculate_distance_data()
head(dispersal_gen$distance_data$base, 10)
table(dispersal_gen$distance_data$base$distance_class)


## -----------------------------------------------------------------------------
# plot dispersal curves for mean dispersal rates
disp_fun <- function(p, b, distance) {
  p * exp(-distance / b)
}

disp_mat <- data.frame(
  p = round(runif(1000, 5, 40) / 100, 2), # prop
  b = round(runif(1000, 5, 40)) # mean distance
)
head(disp_mat)

disp_test <- lapply(1:nrow(disp_mat), function(i) {
  p <- disp_mat[i, "p"]
  b <- disp_mat[i, "b"]
  disp_x <- disp_fun(p, b, seq(5, 150, 5))
  return(disp_x)
})

{
  par(mar = c(4, 4, 0.5, 0.5))
  matplot(
    x = seq(5, 150, 5), y = rep(NA, 30), type = "l", ylim = c(0, 0.4),
    xlab = "Disp. dist (km)", ylab = "Prop. disp.", yaxt = "n", xaxt = "n"
  )
  axis(1, at = seq(0, 150, 10))
  axis(2, at = seq(0, 40, 5) / 100, labels = seq(0, 40, 5))
  lapply(disp_test, function(i) {
    matplot(
      x = seq(5, 150, 5), y = unlist(i), type = "l", add = TRUE,
      col = c("#C9C9C944")
    )
  })
  lines(
    x = seq(5, 150, 5),
    y = apply(as.data.frame(disp_test), 1, mean), col = "firebrick"
  )
}
dev.off()


## -----------------------------------------------------------------------------
# Generate sampled dispersals for p = 0.35, b = 40 (km)
sample_dispersal_data <- dispersal_gen$generate(
  input_values = list(dispersal_p = 0.35, dispersal_b = 40)
)$dispersal_data
head(sample_dispersal_data[[1]], 10) # examine


## -----------------------------------------------------------------------------
capacity_gen <- Generator$new(
  description = "capacity",
  spatial_correlation = env_corr,
  generate_rasters = FALSE,
  time_steps = ncol(initial_hs),
  hs_raster = initial_hs[region$region_indices], # provide full stack of HS. Template attached
  hs_mod = tasmania_modifier[region$region_indices], # provide full stack of LULC modifier. Template attached
  int_index = int_point,
  trans_n = translocation$trans_n, # number of animals introduced
  region_indices = region$region_indices,
  inputs = c("max_dens", "q_thresh", "trans_n"),
  outputs = c("initial_abundance", "carrying_capacity"),
  generative_requirements = list(
    initial_abundance = "function",
    carrying_capacity = "function"
  )
)

capacity_gen$add_function_template(
  param = "initial_abundance",
  function_def = function(params) {
    distr_a <- params$hs_raster[, 1]
    ## 0 everywhere except the intro point at the first time step
    ## intro point to trans_n
    ## Could be above or below carrying capacity
    idx <- which(params$region_indices == params$int_index)
    distr_a[idx] <- params$trans_n
    distr_a[-idx] <- 0
    return(distr_a)
  },
  call_params = c("hs_raster", "int_index", "region_indices", "trans_n")
)

capacity_gen$add_function_template(
  "carrying_capacity",
  function_def = function(params) {
    idx <- which(params$region_indices == params$int_index)
    distr_k <- params$hs_raster
    distr_mod <- params$hs_mod
    stopifnot(
      "hs_raster and hs_mod have different number of layers" =
        dim(distr_k) == dim(distr_mod)
    )
    # stretch HS values based on q_thresh
    distr_k <- scales::rescale(distr_k, from = c(0, params$q_thresh), to = c(0, 1))
    distr_k[distr_k < 0] <- 0
    distr_k[distr_k > 1] <- 1
    # multiply thresholded HS by hs_modifier
    distr_k <- distr_k * distr_mod
    # rescale back to {0, 1}
    qMax <- max(distr_k, na.rm = TRUE)
    distr_k <- scales::rescale(distr_k, from = c(0, qMax), to = c(0, 1))
    distr_k[distr_k < 0] <- 0
    distr_k[distr_k > 1] <- 1
    # carrying capacity = (HS * maximum density)
    distr_k <- ceiling(distr_k * params$max_dens)
    distr_k[idx, 1] <- params$max_dens
    # distr_k[-idx, 1] <- 0
    return(distr_k)
  },
  call_params = c("hs_raster", "hs_mod", "int_index", "region_indices", "max_dens", "q_thresh")
)

# have all parameters been specified correctly
capacity_gen$generative_requirements_satisfied()


## -----------------------------------------------------------------------------
# Generate example initial abundance and declining carrying capacity time-series
generated_k <- capacity_gen$generate(input_values = list(
  max_dens = 100, q_thresh = 0.90,
  trans_n = 60
))
example_initial_abundance <- generated_k$initial_abundance
example_carrying_capacity <- generated_k$carrying_capacity

# Plot the example initial abundance
example_initial_n_raster <- region$raster_from_values(example_initial_abundance)
example_initial_n_raster
plot(example_initial_n_raster,
  main = "Example initial abundance",
  col = hcl.colors(100, "Lajolla", rev = TRUE), colNA = "grey75",
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
  }
)

# Plot the carrying capacity
## carrying capacity is forced to maximum theoretical value at first time step
example_k <- region$raster_from_values(example_carrying_capacity)
example_k[[c(1, 6, 11)]]
plot(example_k,
  col = hcl.colors(100, "RdYlGn", rev = TRUE), colNA = "grey75",
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
  },
  zlim = c(0, 100)
)


## -----------------------------------------------------------------------------
# Template model ####
model_template <- PopulationModel$new(
  region = region,
  time_steps = 11,
  years_per_step = 1,
  stage_matrix = 1, # single-stage
  populations = region$region_cells, # 735
  demographic_stochasticity = TRUE,
  standard_deviation = 0.18,
  density_dependence = "logistic", # Ricker
  harvest = FALSE, # No harvest
  dispersal = dispersal_gen,
  translocation = translocation,
  dispersal_source_n_k = list(threshold = 0.92, cutoff = 0),
  simulation_order = c("translocation", "results", "transition", "dispersal"),
  random_seed = 20230210,
  attribute_aliases = translocation_aliases,
  results_selection = c("abundance")
)

model <- model_template$clone()
model$set_attributes(
  initial_abundance = example_initial_abundance,
  carrying_capacity = example_carrying_capacity,
  growth_rate_max = rmax_gen_ex$growth_rate_max,
  translocation = translocation,
  trans_n = 75, # passed through to translocation function
  dispersal = sample_dispersal_data
)
# run poems simulator
results <- population_simulator(model)
results$all$abundance

# timeseries of total abundance
plot(1:11, results$all$abundance,
  type = "l",
  xlab = "timestep", ylab = "Total abundance"
)


## -----------------------------------------------------------------------------
abund_ras <- region$raster_from_values(results$abundance)
abund_ras[[c(1, 6, 11)]]
abd_max <- round_any(max(values(abund_ras), na.rm = TRUE), 20, f = ceiling)

# plot of abundances. log(x+1) transformed.
plot(log1p(abund_ras),
  col = hcl.colors(100),
  colNA = "grey75",
  addfun = function() {
    plot(ibra_polygons, border = "black", col = NA, add = TRUE)
  },
  zlim = c(0, log1p(abd_max))
)


## -----------------------------------------------------------------------------
model$set_attributes(
  initial_abundance = example_initial_abundance,
  carrying_capacity = example_carrying_capacity,
  growth_rate_max = rmax_gen_ex$growth_rate_max,
  translocation = NULL,
  dispersal = sample_dispersal_data
)
results_notransn <- population_simulator(model) # run poems simulator
results_notransn$all$abundance
results$all$abundance
abund_ras_notransn <- region$raster_from_values(results_notransn$abundance)
abund_ras_notransn[[c(1, 6, 11)]]
diff_ras <- abund_ras - abund_ras_notransn
diff_ras[[1:2]]


## -----------------------------------------------------------------------------
plotmax <- round_any(max(abs(values(diff_ras)), na.rm = TRUE), 10, ceiling)
plot(diff_ras[[c(1, 6, 11)]],
  zlim = c(-plotmax, plotmax),
  breaks = c(-plotmax, -100, -50, -20, 0, 20, 50, 100, plotmax),
  col = hcl.colors(9, "PuOr"),
  colNA = "grey75",
  addfun = function() {
    plot(ibra_polygons, col = NA, add = TRUE)
  }
)


## -----------------------------------------------------------------------------
# Latin-hypercube sampler ####
lhs_gen <- LatinHypercubeSampler$new()

# Habitat suitability threshold
lhs_gen$set_uniform_parameter("q_thresh", lower = 0.90, upper = 0.99, decimals = 2)

# Growth rate
lhs_gen$set_uniform_parameter("rmax_quantile", lower = 0, upper = 1, decimals = 2)
lhs_gen$set_uniform_parameter("standard_deviation", lower = 0.00, upper = 0.70, decimals = 2)

# Dispersal
lhs_gen$set_uniform_parameter("dispersal_p", lower = 0.05, upper = 0.40, decimals = 2)
## mean dispersal between 5 and 40 km
lhs_gen$set_uniform_parameter("dispersal_b", lower = 5, upper = 40, decimals = 0)
lhs_gen$set_uniform_parameter("dispersal_n_k_threshold", lower = 0.7, upper = 1.0, decimals = 2)

# Density max
## Density: animals/km2 needs to be scaled by grid size (10km x 10km)
## e.g. 1/km2 = (1 animal/km2 * (10*10) ) * frac_cell_used
## 1 km2 = 80 per grid cell = (1*(10*10))*0.8 # assuming 80% grid cell used
## Here I have assumed only 80% of cell is suitable. Upper/lower = 1/km - 6.25/km
lhs_gen$set_uniform_parameter("max_dens", lower = 80, upper = 500, decimals = 0)

# Translocation
lhs_gen$set_uniform_parameter("trans_n", lower = 10, upper = 100, decimals = 0)

sample_data <- lhs_gen$generate_samples(number = 10, random_seed = 42)
head(sample_data)

# Make unique row names for saving files
{
  set.seed(54612)
  sample_data$UniqueID <- paste0(
    stri_rand_strings(nrow(sample_data), 4, "[A-Z]"),
    stri_rand_strings(nrow(sample_data), 4, "[0-9]")
  )
}
sample_data <- sample_data[, c(9, 1:8)]
sample_data


## -----------------------------------------------------------------------------
OUTPUT_DIR <- tempdir()
model <- model_template$clone()
model$set_attributes(params = list(
  "standard_deviation" = NULL,
  "dispersal_source_n_k$threshold" = NULL,
  "dispersal_source_n_k$cutoff" = 0.00
))

# Build the simulation manager
sim_manager <- SimulationManager$new(
  sample_data = sample_data,
  model_template = model,
  # initial_hs = initial_hs,
  generators = list(dispersal_gen, capacity_gen, rmax_gen),
  parallel_cores = 1L,
  results_filename_attributes =
    c(NULL, "UniqueID", "results"),
  results_ext = ".RDS",
  results_dir = OUTPUT_DIR
)

# Takes <10 seconds to run 10 example sims on a single core.
system.time({
  run_output <- sim_manager$run()
})
run_output$summary


## -----------------------------------------------------------------------------
run_output$summary


## -----------------------------------------------------------------------------
# Extract timeseries of abundance from each of the sims
# Load our results (list) into a PopulationResults object
p_results <- PopulationResults$new(results = run_output)
res_manager <- ResultsManager$new(
  simulation_manager = sim_manager,
  simulation_results = p_results,
  generators = NULL,
  summary_matrices = c(
    "n",
    "distr_pop"
  ),
  summary_functions = list(
    # total pop abundance
    "n" = function(sim_results) {
      sim_results$all$abundance
    },
    # matrix of abundance
    ## can be made into raster
    "distr_pop" = function(sim_results) {
      sim_results$abundance
    }
  ),
  parallel_cores = 1L
)
gen_log <- res_manager$generate()
gen_log$summary

# matrix of total population abundances
## each row is a sim, each column a timestep
res_manager$summary_matrix_list$n

# plot
matplot(
  x = 1:ncol(res_manager$summary_matrix_list$n),
  y = t(res_manager$summary_matrix_list$n), type = "b",
  lty = 1, xlab = "timestep", ylab = "total abundance"
)


## -----------------------------------------------------------------------------
identical(
  unlist(res_manager$summary_matrix_list$n[, 1]),
  unlist(sample_data$trans_n)
)


## -----------------------------------------------------------------------------
best_sims <- c(2:4, 8)
dim(res_manager$summary_matrix_list$distr_pop[best_sims, ])

best_abund <- matrix(
  nrow = region$region_cells,
  ncol = 11, # 11 timesteps,
  data = round(colMeans(res_manager$summary_matrix_list$distr_pop[best_sims, ]))
)

best_abund <- region$raster_from_values(best_abund)
best_abund[[c(1, 6, 11)]]

abd_max <- round_any(max(values(best_abund), na.rm = TRUE),
  accuracy = 100, ceiling
)

# plot of log(x+1) abundances
plot(log1p(best_abund),
  col = hcl.colors(100, "Spectral", rev = TRUE),
  colNA = "grey75",
  addfun = function() {
    plot(ibra_polygons, border = "#000000", col = NA, add = TRUE)
  },
  zlim = c(0, log1p(abd_max))
)

