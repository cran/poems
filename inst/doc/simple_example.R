## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(poems)

## ----message = FALSE, fig.align = "center", fig.width = 4, fig.height = 4-----
# Demonstration example region (U Island)
coordinates <- data.frame(
  x = rep(seq(177.01, 177.05, 0.01), 5),
  y = rep(seq(-18.01, -18.05, -0.01), each = 5)
)
template_raster <- Region$new(coordinates = coordinates)$region_raster # full extent
template_raster[][-c(7, 9, 12, 14, 17:19)] <- NA # make U Island
region <- Region$new(template_raster = template_raster)
raster::plot(region$region_raster,
  main = "Example region (cell indices)",
  xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
  colNA = "blue"
)

## ----message = FALSE----------------------------------------------------------
# Distance-based environmental correlation (via a compacted Cholesky decomposition)
env_corr <- SpatialCorrelation$new(region = region, amplitude = 0.4, breadth = 500)
correlation <- env_corr$get_compact_decomposition(decimals = 2)
correlation # examine

## ----message = FALSE----------------------------------------------------------
# User-defined harvest function (list-nested) and alias
harvest <- list(
  rate = NA, # sample later
  function(params) round(params$stage_abundance * (1 - params$rate))
)
harvest_rate_alias <- list(harvest_rate = "harvest$rate")

## ----message = FALSE----------------------------------------------------------
# Population (simulation) model template for fixed parameters
stage_matrix <- matrix(
  c(
    0, 2.5, # Leslie/Lefkovitch matrix
    0.8, 0.5
  ),
  nrow = 2, ncol = 2, byrow = TRUE,
  dimnames = list(c("juv", "adult"), c("juv", "adult"))
)
stage_matrix # examine
model_template <- PopulationModel$new(
  region = region,
  time_steps = 10, # years
  populations = region$region_cells, # 7
  stages = 2,
  stage_matrix = stage_matrix,
  demographic_stochasticity = TRUE,
  standard_deviation = 0.05,
  correlation = correlation,
  density_dependence = "logistic",
  harvest = harvest,
  results_selection = c("abundance", "harvested"),
  attribute_aliases = harvest_rate_alias
)

## ----message = FALSE, fig.align = "center", fig.width = 4, fig.height = 4-----
# Example habitat suitability
example_hs <- c(0.8, 1, 0.7, 0.9, 0.6, 0.7, 0.8)
example_hs_raster <- region$region_raster
example_hs_raster[region$region_indices] <- example_hs
raster::plot(example_hs_raster,
  main = "Example habitat suitability",
  xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
  colNA = "blue"
)

## ----message = FALSE----------------------------------------------------------
# Initial abundance and carrying capacity generated via example habitat suitability
capacity_gen <- Generator$new(
  description = "Capacity generator",
  example_hs = example_hs, # template attached
  inputs = c("initial_n", "density_max"),
  outputs = c("initial_abundance", "carrying_capacity")
)
capacity_gen$add_generative_requirements(list(
  initial_abundance = "function",
  carrying_capacity = "function"
))
capacity_gen$add_function_template("initial_abundance",
  function_def = function(params) {
    stats::rmultinom(1,
      size = params$initial_n,
      prob = params$example_hs
    )[, 1]
  },
  call_params = c("initial_n", "example_hs")
)
capacity_gen$add_function_template("carrying_capacity",
  function_def = function(params) {
    round(params$density_max * params$example_hs)
  },
  call_params = c("density_max", "example_hs")
)
capacity_gen$generate(input_values = list(initial_n = 500, density_max = 100)) # test

## ----message = FALSE----------------------------------------------------------
# Distance-based dispersal generator
dispersal_gen <- DispersalGenerator$new(
  region = region,
  dispersal_max_distance = 3000, # in m
  dispersal_friction = DispersalFriction$new(),
  inputs = c("dispersal_p", "dispersal_b"),
  decimals = 5
)
dispersal_gen$calculate_distance_data() # pre-calculate
test_dispersal <- dispersal_gen$generate(input_values = list(
  dispersal_p = 0.5,
  dispersal_b = 700
))
head(test_dispersal$dispersal_data[[1]])

## ----message = FALSE----------------------------------------------------------
# Generate sampled values for variable model parameters via LHS
lhs_gen <- LatinHypercubeSampler$new()
lhs_gen$set_uniform_parameter("growth_rate_max", lower = 0.4, upper = 0.6, decimals = 2)
lhs_gen$set_uniform_parameter("harvest_rate", lower = 0.05, upper = 0.15, decimals = 2)
lhs_gen$set_uniform_parameter("initial_n", lower = 400, upper = 600, decimals = 0)
lhs_gen$set_uniform_parameter("density_max", lower = 80, upper = 120, decimals = 0)
lhs_gen$set_uniform_parameter("dispersal_p", lower = 0.2, upper = 0.5, decimals = 2)
lhs_gen$set_uniform_parameter("dispersal_b", lower = 400, upper = 1000, decimals = 0)
sample_data <- lhs_gen$generate_samples(number = 12, random_seed = 123)
sample_data # examine

## ----message = FALSE----------------------------------------------------------
# Create a simulation manager and run the sampled model simulations
OUTPUT_DIR <- tempdir()
sim_manager <- SimulationManager$new(
  sample_data = sample_data,
  model_template = model_template,
  generators = list(capacity_gen, dispersal_gen),
  parallel_cores = 2,
  results_dir = OUTPUT_DIR
)
run_output <- sim_manager$run()
run_output$summary
dir(OUTPUT_DIR, "*.RData") # includes 12 result files
dir(OUTPUT_DIR, "*.txt") # plus simulation log

## ----message = FALSE----------------------------------------------------------
results_manager <- ResultsManager$new(
  simulation_manager = sim_manager,
  simulation_results = PopulationResults$new(),
  summary_metrics = c("trend_n", "total_h"),
  summary_matrices = c("n", "h"),
  summary_functions = list(
    trend_n = function(results) {
      round(results$all$abundance_trend, 2)
    },
    total_h = function(results) {
      sum(results$harvested)
    },
    n = "all$abundance", # string
    h = "all$harvested"
  ),
  parallel_cores = 2
)
gen_output <- results_manager$generate()
gen_output$summary
dir(OUTPUT_DIR, "*.txt") # plus generation log
results_manager$summary_metric_data
results_manager$summary_matrix_list

## ----message = FALSE----------------------------------------------------------
# Create a validator for selecting the 'best' example models
validator <- Validator$new(
  simulation_parameters = sample_data,
  simulation_summary_metrics =
    results_manager$summary_metric_data[-1],
  observed_metric_targets = c(trend_n = 0, total_h = 600),
  output_dir = OUTPUT_DIR
)
suppressWarnings(validator$run(tolerance = 0.25, output_diagnostics = TRUE))
dir(OUTPUT_DIR, "*.pdf") # plus validation diagnostics (see abc library documentation)
validator$selected_simulations # top 3 models (stable abundance and high harvest)

## ----message = FALSE, fig.align = "center", fig.width = 6, fig.height = 5-----
# Plot the simulation, targets, and selected metrics
graphics::plot(
  x = results_manager$summary_metric_data$total_h,
  y = results_manager$summary_metric_data$trend_n,
  main = "Example model validation",
  xlab = "Total harvested", ylab = "Abundance trend"
)
graphics::points(x = 600, y = 0, col = "red", pch = 4)
selected_indices <- validator$selected_simulations$index
graphics::points(
  x = results_manager$summary_metric_data$total_h[selected_indices],
  y = results_manager$summary_metric_data$trend_n[selected_indices],
  col = "blue", pch = 3
)
graphics::legend("bottomleft",
  legend = c("Summary metrics", "Targets", "Selected"),
  col = c(1, "red", "blue"), pch = c(1, 4, 3), cex = 0.8
)

