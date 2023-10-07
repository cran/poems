## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(poems)
DEMONSTRATION <- TRUE # load pre-run data rather than running simulations
SAMPLES <- 20000
PARALLEL_CORES <- 2
OUTPUT_DIR <- tempdir()

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Raster of Tasmania (note: islands removed where there was no evidence of thylacine occupancy).
raster::plot(poems::tasmania_raster, main = "Tasmania raster",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Tasmania study region (795 cells stored in the order shown)
region <- Region$new(template_raster = tasmania_raster)
raster::plot(region$region_raster, main = "Tasmanian study region (cell indices)",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             colNA = "blue")

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Tasmania study Interim Bioregionalisation of Australia (IBRA) bioregion cell distribution
ibra_raster <- poems::tasmania_ibra_raster
raster::plot(ibra_raster, main = "Tasmania study IBRA bioregions", colNA = "blue",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)")
ibra_data <- poems::tasmania_ibra_data
ibra_data # examine

# Calculate cell indices and counts for IBRA bioregions
ibra_indices <- lapply(as.list(ibra_data$index),
                       function(i) {
                         which(ibra_raster[region$region_indices] == i)
                       })
ibra_indices[1:2] # examine
ibra_cells <- unlist(lapply(ibra_indices, length))
ibra_cells # examine

## ----message = FALSE----------------------------------------------------------
# Build the dispersal generator and calculate distance data
dispersal_gen <- DispersalGenerator$new(
  region = region,
  dispersal_max_distance = 50,
  distance_scale = 1000, # in km
  dispersal_friction = DispersalFriction$new(), # modify coastline distances
  inputs = c("dispersal_p", "dispersal_b"),
  decimals = 5)
dispersal_gen$calculate_distance_data()
head(dispersal_gen$distance_data$base) # examine

## ----message = FALSE----------------------------------------------------------
# Define neighborhoods (of up to 9 adjacent cells) based on a 14 km range from each 
# grid cell for density dependence calculations (using a dispersal generator)
distance_data <- dispersal_gen$distance_data[[1]]
nh_data <- distance_data[which(distance_data$distance_class <= 14), 2:1]
neighborhoods <- as.list(1:795)
for (i in 1:nrow(nh_data)) {
  neighborhoods[[nh_data$source_pop[i]]] <- c(neighborhoods[[nh_data$source_pop[i]]],
                                               nh_data$target_pop[i])
}
neighborhoods[1:3] # examine

## ----message = FALSE----------------------------------------------------------
# User-defined function for Ricker logistic density dependence via neighborhoods, with
# Allee effects; also remove fecundities if single thylacine in a neighborhood
density_dependence <- list(
  neighborhoods = neighborhoods,
  allee = 25, # Tasmania-wide Allee effect parameter
  function (params) {
    
    # Apply logistic density dependence using neighborhoods
    growth_rate_max <- params$growth_rate_max
    nh_density_abundance <- unlist(lapply(params$neighborhoods,
                                          function (nh_indices) {
                                            sum(params$density_abundance[nh_indices])
                                          }))
    nh_carrying_capacity <- unlist(lapply(params$neighborhoods,
                                          function (nh_indices) {
                                            sum(params$carrying_capacity[nh_indices])
                                          }))
    occupied_indices <- params$occupied_indices
    growth_rate <- growth_rate_max*(1 - (nh_density_abundance[occupied_indices]/
                                           nh_carrying_capacity[occupied_indices]))
    params$transition_array[, , occupied_indices] <-
      params$apply_multipliers(params$transition_array[, , occupied_indices],
                               params$calculate_multipliers(growth_rate))
    
    # Apply Tasmania-wide allee effect
    total_abundance <- sum(params$density_abundance)
    params$transition_array <-
      params$transition_array*total_abundance/(params$allee + total_abundance)
    
    # Remove fecundities for single thylacines
    single_indices <- which(nh_density_abundance == 1)
    params$transition_array[, , single_indices] <-
      (params$transition_array[, , single_indices]*as.vector(+(!params$fecundity_mask)))
    
    return(params$transition_array)
  }
)
density_aliases <- list(density_allee = "density_dependence$allee")

## ----message = FALSE----------------------------------------------------------
# Harvest bounty (economic) model user-defined function adapted from Bulte et al. (2003).
harvest <- list(

  # Function parameters (passed to function in params list)
  ohr = 0.05,        # opportunistic harvest rate
  t1 = 1888,         # first year of bounty
  tb = 1909,         # last year of bounty
  fb = 0.75,         # harvest fraction submitted for bounty
  B = c(1.6, 0.6),   # bounty/skin price in pounds, pre/post bounty
  w = 3.5,           # opportunity cost in pounds per year
  E0 = 25,           # effort in 1888 (no. hunters)
  q = 0.0025,        # catchability coefficient
  v1 = 0.02,         # entry rate
  v2 = 0.5,          # exit rate
  ibra_indices = ibra_indices, # bioregion cell (row) indices
  ibra_cells = ibra_cells,     # number of cells in bioregions

  # Function definition
  bounty_function = function(params) {

    # Unpack parameters (used at every time step)
    ohr <- params$ohr; t1 <- params$t1; tb <- params$tb; fb <- params$fb
    B <- params$B; w <- params$w; q <- params$q; v1 <- params$v1; v2 <- params$v2
    ibra_indices <- params$ibra_indices; ibra_cells <- params$ibra_cells
    ibra_number <- length(ibra_cells); stages <- params$stages
    populations <- params$populations; simulator <- params$simulator
    tm <- params$tm; x <- params$stage_abundance

    # Initialise (first time step only)
    if (tm == 1) { # attach variables and access results via simulator reference object
      simulator$attached$E <- params$E0 # current bounty effort
      simulator$attached$vi <- v1 # current bounty rate
      simulator$results$bounty <- array(0, c(ibra_number, params$time_steps))
    }

    # Access persistent parameters via simulator reference object
    E <- simulator$attached$E
    vi <- simulator$attached$vi

    # Next year's hunting effort and entry/exit rates based on this year's profit
    h <- max(0, round((ohr + q*E)*sum(x))) # harvest
    b <- round(h*fb*((tm + t1 - 1) <= tb)) # bounty submitted
    profit <- round(B[((tm + t1 - 1) > tb) + 1]*b + B[2]*(h - b) - w*E, 1)
    simulator$attached$E <-  max(0, round(E + vi*profit))
    simulator$attached$vi <- c(v1, v2)[(profit < 0) + 1]

    # Distribute harvest and bounty across bioregions based on each IBRA density
    staged_indices <- array(1:(stages*populations), c(stages, populations))
    rep_indices <- unlist(apply(matrix(staged_indices[, unlist(ibra_indices)]), 1,
                                function(i) rep(i, x[i])))
    distributed_h <- array(0, c(stages, populations))
    if (length(rep_indices) && h > 0) {
      ibra_x <- unlist(lapply(ibra_indices, function(indices) sum(x[, indices])))
      rep_ibra <- unlist(apply(matrix(1:ibra_number), 1, function(i) rep(i, ibra_x[i])))
      rep_prob <- 1/ibra_cells[rep_ibra]
      h_indices <- sample(1:length(rep_indices), min(h, sum(x)), prob = rep_prob)
      if (b > 0) {
        b_indices <- h_indices[sample(1:length(h_indices), b)]
        simulator$results$bounty[, tm] <- tabulate(rep_ibra[b_indices],
                                                   nbins = ibra_number)
      }
      for (i in rep_indices[h_indices]) distributed_h[i] <- distributed_h[i] + 1
    }

    # Return abundance
    return(x - distributed_h)
  }
)
harvest_aliases <- list(harvest_ohr = "harvest$ohr", harvest_fb = "harvest$fb",
                        harvest_w = "harvest$w", harvest_E0 = "harvest$E0",
                        harvest_q = "harvest$q", harvest_v1 = "harvest$v1",
                        harvest_v2 = "harvest$v2")

## ----message = FALSE----------------------------------------------------------
# Population (simulation) model template for fixed parameters
model_template <- PopulationModel$new(
  region = region,
  time_steps = 80, # years (1888-1967)
  populations = region$region_cells, # 795
  # initial_abundance : generated
  # stage_matrix: generated
  fecundity_max = 2,
  demographic_stochasticity = TRUE,
  # carrying_capacity : generated
  density_dependence = density_dependence, # user-defined
  harvest = harvest, # user-defined
  # dispersal : generated
  dispersal_target_k = 0.5,
  dispersal_target_n = list(threshold = 4, cutoff = 8),
  simulation_order = c("results", "harvest", "transition", "dispersal"),
  results_selection = c("abundance", "harvested"),
  attribute_aliases = c(density_aliases, harvest_aliases))

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Initial thylacine habitat suitability
hs_raster <- poems::thylacine_hs_raster
raster::plot(hs_raster, main = "Initial thylacine habitat suitability", colNA = "blue",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)")

## ----message = FALSE----------------------------------------------------------
# Build a carrying capacity generator model based on habitat suitability and sampled
# initial capacity, initial fraction (phi), & decline rate per year in selected bioregions.
capacity_gen <- Generator$new(
  description = "capacity",
  time_steps = 80, # years (1888-1967)
  initial_hs = hs_raster[region$region_indices],
  decline_indices = which(!ibra_raster[region$region_indices] %in% c(5, 8)),
  inputs = c("k_init", "k_decline", "k_phi"),
  outputs = c("initial_abundance", "carrying_capacity"),
  generative_requirements =  list(initial_abundance = "function",
                                  carrying_capacity = "function"))
capacity_gen$add_function_template(
  "initial_abundance",
  function_def = function (params) {
    distr_k <- round(params$initial_hs/sum(params$initial_hs)*params$k_init)
    a_init <- round(params$k_init*params$k_phi) # total initial
    distr_a <- array(0, length(distr_k))
    rep_indices <- unlist(apply(matrix(1:length(distr_k)), 1,
                                function(i) rep(i, distr_k[i])))
    sample_indices <- rep_indices[sample(1:length(rep_indices),
                                         min(a_init, length(rep_indices)))]
    for (i in sample_indices) distr_a[i] <- distr_a[i] + 1
    return(distr_a)
  },
  call_params = c("initial_hs", "k_init", "k_phi"))
capacity_gen$add_function_template(
  "carrying_capacity",
  function_def = function (params) {
    distr_k <- params$initial_hs/sum(params$initial_hs)*params$k_init
    decline_matrix <- array(1, c(length(distr_k), params$time_steps))
    decline_matrix[params$decline_indices,] <-
      matrix((1 - params$k_decline)^(0:(params$time_steps - 1)),
             nrow = length(params$decline_indices), ncol = params$time_steps,
             byrow = TRUE)
    return(distr_k*decline_matrix)
  },
  call_params = c("initial_hs", "time_steps", "decline_indices",
                  "k_init", "k_decline"))

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Generate example initial abundance and declining carrying capacity time-series
generated_k <- capacity_gen$generate(input_values = list(k_init = 2800,
                                                         k_decline = 0.04,
                                                         k_phi = 0.8))
example_initial_abundance <- generated_k$initial_abundance
example_carrying_capacity <- generated_k$carrying_capacity

# Plot the example initial abundance
example_initial_n_raster <- region$region_raster
example_initial_n_raster[region$region_indices] <- example_initial_abundance
raster::plot(example_initial_n_raster, main = "Example initial thylacines", 
             colNA = "blue", xlab = "Longitude (degrees)", ylab = "Latitude (degrees)")

# Plot the example final carrying capacity
example_final_raster <- region$region_raster
example_final_raster[region$region_indices] <- example_carrying_capacity[, 80]
raster::plot(example_final_raster, main = "Final thylacine carrying capacity", 
             colNA = "blue",  xlab = "Longitude (degrees)", ylab = "Latitude (degrees)", 
             zlim = c(0, 8))

## ----message = FALSE----------------------------------------------------------
# Build a stage matrix generator based on sampled growth rate
stage_matrix_gen <- Generator$new(
  description = "stage matrix",
  base_matrix = matrix(c(0.00, 0.57, 1.17,
                         0.50, 0.00, 0.00,
                         0.00, 0.80, 0.80), nrow = 3, ncol = 3, byrow = TRUE),
  inputs = c("growth_r"),
  outputs = c("stage_matrix"),
  generative_requirements = list(stage_matrix = "function"))
stage_matrix_gen$add_function_template(
  "stage_matrix",
  function_def = function (params) {
    return(params$base_matrix*(1 + params$growth_r)/
             Re((eigen(params$base_matrix)$values)[1]))
  },
  call_params = c("base_matrix", "growth_r"))

## ----message = FALSE----------------------------------------------------------
# Generate sampled stage matrix for growth rate: lambda = 1.25
gen_stage_m <- stage_matrix_gen$generate(input_values = list(growth_r = 0.25))
gen_stage_m # examine

## ----message = FALSE----------------------------------------------------------
# Generate sampled dispersals for p = 0.5, b = 7 (km)
sample_dispersal_data <- dispersal_gen$generate(
  input_values = list(dispersal_p = 0.5, dispersal_b = 7))$dispersal_data
head(sample_dispersal_data[[1]]) # examine

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Run the model with example parameters
model <- model_template$clone()
model$set_attributes(initial_abundance = example_initial_abundance,
                     carrying_capacity = example_carrying_capacity,
                     stage_matrix = gen_stage_m$stage_matrix,
                     dispersal = sample_dispersal_data)
results <- population_simulator(model) # run poems simulator

# Plot the total abundance and number harvested
plot(x = 1888:1967, y = results$all$abundance, xlab = "Year",
     ylab = "Number of thylacines", main = "Thylacine example model run",
     ylim = c(0, 2500), type = "l", col = "green", lwd = 2)
lines(x = 1888:1967, y = results$all$harvested, lty = 1, col = "blue", lwd = 2)
legend("topright", legend = c("Population size", "Simulated harvest"),
       col = c("green", "blue"), lty = c(1, 1), lwd = 2, cex = 0.8)

## ----message = FALSE----------------------------------------------------------
# Create a LHS object
lhs_gen <- LatinHypercubeSampler$new()

# Set capacity and growth parameters (as per Bulte et al., 2003)
lhs_gen$set_uniform_parameter("k_init", lower = 2100, upper = 3500, decimals = 0)
lhs_gen$set_uniform_parameter("k_decline", lower = 0.03, upper = 0.05, decimals = 3)
lhs_gen$set_uniform_parameter("k_phi", lower = 0.6, upper = 1.0, decimals = 2)
lhs_gen$set_uniform_parameter("growth_r", lower = 0.1875, upper = 0.3125, decimals = 2)

# Set density dependence allee effect parameter
lhs_gen$set_uniform_parameter("density_allee", lower = 0, upper = 50, decimals = 1)

# Set bio-economic harvest parameters (as per Bulte et al., 2003)
lhs_gen$set_uniform_parameter("harvest_ohr", lower = 0, upper = 0.1, decimals = 3)
lhs_gen$set_uniform_parameter("harvest_fb", lower = 0.5, upper = 1.0, decimals = 2)
lhs_gen$set_uniform_parameter("harvest_w", lower = 2.625, upper = 4.375, decimals = 1)
lhs_gen$set_uniform_parameter("harvest_E0", lower = 18.75, upper = 31.25, decimals = 0)
lhs_gen$set_uniform_parameter("harvest_q", lower = 0, upper = 0.005, decimals = 4)
lhs_gen$set_uniform_parameter("harvest_v1", lower = 0.015, upper = 0.025, decimals = 3)
lhs_gen$set_uniform_parameter("harvest_v2", lower = 0.375, upper = 0.625, decimals = 3)

# Set new spatial parameters for dispersal
lhs_gen$set_uniform_parameter("dispersal_p", lower = 0.3, upper = 0.7, decimals = 2)
lhs_gen$set_uniform_parameter("dispersal_b", lower = 4, upper = 10, decimals = 1)

# Generate samples
sample_data <- lhs_gen$generate_samples(number = SAMPLES, random_seed = 123)
head(sample_data) # examine
dim(sample_data) # dimensions

## ----message = FALSE----------------------------------------------------------
# Build the simulation manager
sim_manager <- SimulationManager$new(
  sample_data = sample_data,
  model_template = model_template,
  generators = list(capacity_gen, stage_matrix_gen, dispersal_gen),
  parallel_cores = PARALLEL_CORES,
  results_dir = OUTPUT_DIR)

# Run the simulations
if (DEMONSTRATION) {
  sim_manager$sample_data <- sample_data[1:2,]
}
run_output <- sim_manager$run()
run_output$summary
if (DEMONSTRATION) {
  dir(OUTPUT_DIR, "*.RData") # includes 2 result files
}
dir(OUTPUT_DIR, "*.txt") # plus simulation log

## ----message = FALSE----------------------------------------------------------
# Load our results (list) into a PopulationResults object
p_results <- PopulationResults$new(results = results,
                                   ibra_indices = ibra_indices)

# Summary metrics for IBRA bioregions and Tasmania-wide extinction
ibra_bounty <- p_results$get_attribute("bounty") # saved in harvest function
ibra_bounty_clone <- p_results$new_clone(results = list(abundance = ibra_bounty),
                                         trend_interval = (1888:1894) - 1887)
ibra_bounty_clone$all$abundance_trend # 1888-1894 total bounty slope
ibra_abundance <- t(array(unlist(lapply(p_results$get_attribute("ibra_indices"),
                                        function (indices) {
                                          colSums(p_results$abundance[indices,])
                                        })), c(80, 9)))
ibra_abundance_clone <- p_results$new_clone(results = list(abundance = ibra_abundance))
(1888:1967)[ibra_abundance_clone$extirpation] # IBRA extirpation
(1888:1967)[p_results$all$extirpation] # total extinction

## ----message = FALSE----------------------------------------------------------
# Set targets for our summary metrics (used to calculate combined metric errors)
slope_intervals <- c("1888-1894", "1895-1901", "1902-1909")
targets <- list(
  bounty_slope = array(c(2.36, 3.25, -17.71), dimnames = list(slope_intervals)),
  ibra_extirpation = array(c(c(NA, NA), c(1934, 1934), c(1912, 1919), c(1921, 1940),
                             c(1936, 1938), c(1935, 1935), c(1934, 1942), 
                             c(1934, 1934), c(1932, 1932)), c(2, 9), 
                           dimnames = list(c("lower", "upper"), ibra_data$abbr)),
  total_extinction = c(lower = 1936, upper = 1942) # CI
)

# Create a results manager for summary metrics and matrices
results_manager <- ResultsManager$new(
  simulation_manager = sim_manager,
  simulation_results = PopulationResults$new(ibra_indices = ibra_indices, # attachments
                                             targets = targets,
                                             extirp_NA_replace = 1968),
  result_attachment_functions = list( # attached for multiple use
    bounty_slope = function(results) { # via results object cloning
      bounty_slope <- array(NA, 3)
      ibra_bounty <- results$get_attribute("bounty") # saved in harvest function
      ibra_bounty_clone <- results$new_clone(results = list(abundance = ibra_bounty),
                                             trend_interval = (1888:1894) - 1887)
      bounty_slope[1] <- ibra_bounty_clone$all$abundance_trend
      ibra_bounty_clone <- results$new_clone(results = list(abundance = ibra_bounty),
                                             trend_interval = (1895:1901) - 1887)
      bounty_slope[2] <- ibra_bounty_clone$all$abundance_trend
      ibra_bounty_clone <- results$new_clone(results = list(abundance = ibra_bounty),
                                             trend_interval = (1902:1909) - 1887)
      bounty_slope[3] <- ibra_bounty_clone$all$abundance_trend
      bounty_slope
    },
    ibra_extirpation = function(results) { # via results object cloning
      ibra_abundance_clone <- results$new_clone(results = list(
        abundance = t(array(unlist(lapply(results$get_attribute("ibra_indices"),
                                          function (indices) {
                                            colSums(results$abundance[indices,])
                                          })), c(80, 9)))))
      (1888:1967)[ibra_abundance_clone$extirpation]
    }),
  summary_metrics = c("bounty_slope_error", "ibra_extirpation_error",
                      "total_extinction"),
  summary_matrices = c("extirpation", "total_bounty", "ibra_bounty",
                       "bounty_slope", "ibra_extirpation"),
  summary_functions = list(
    # Summary metrics
    bounty_slope_error = function(results) { # RMSE
      sqrt(mean((results$get_attribute("targets")$bounty_slope -
                   results$get_attribute("bounty_slope"))^2))
    },
    ibra_extirpation_error = function(results) { # RMSE with NAs replaced
      ibra_extirpation <- results$get_attribute("ibra_extirpation")
      ibra_extirpation[is.na(ibra_extirpation)] <-
        results$get_attribute("extirp_NA_replace")
      target_CI <- results$get_attribute("targets")$ibra_extirpation
      sqrt(mean(((ibra_extirpation < target_CI[1,])*(ibra_extirpation - target_CI[1,]) + 
           (ibra_extirpation > target_CI[2,])*(ibra_extirpation - target_CI[2,]))^2,
           na.rm = TRUE))
    },
    total_extinction = function(results) {
      (1888:1967)[results$all$extirpation]
    },
    # Summary matrices
    extirpation = function(results) { # for later use
      (1888:1967)[results$extirpation]
    },
    total_bounty = function(results) { # for later use
      colSums(results$get_attribute("bounty"))
    },
    ibra_bounty = function(results) { # for later use
      rowSums(results$get_attribute("bounty"))
    },
    bounty_slope = function(results) { # calculate RMSE later
      results$get_attribute("bounty_slope")
    },
    ibra_extirpation = function(results) { # calculate RMSE later
      results$get_attribute("ibra_extirpation")
    }),
  parallel_cores = PARALLEL_CORES)

# Generate the summary metrics and matrices
gen_output <- results_manager$generate()
gen_output$summary
dir(OUTPUT_DIR, "*.txt") # plus generation log
summary_metric_data <- results_manager$summary_metric_data
summary_matrix_list <- results_manager$summary_matrix_list
head(summary_metric_data) # examine
lapply(summary_matrix_list, dim) # dimensions
head(summary_matrix_list$bounty_slope) # examine
head(summary_matrix_list$ibra_extirpation) # examine

## ----message = FALSE----------------------------------------------------------
# Demonstrate calculating RMSE metrics from matrices 
if (DEMONSTRATION) { # Calculate RMSE for bounty slopes
  bounty_slope_error2 <- sqrt(rowMeans((summary_matrix_list$bounty_slope - 
                                          matrix(targets$bounty_slope, nrow = 2, 
                                                 ncol = 3, byrow = TRUE))^2))
  
  cbind(bounty_slope_error = summary_metric_data$bounty_slope_error, 
        bounty_slope_error2) # examine
}
if (DEMONSTRATION) { # Calculate RMSE for IBRA extirpation
  ibra_extirpation <- summary_matrix_list$ibra_extirpation
  ibra_extirpation[is.na(ibra_extirpation)] <- 1968
  target_CI <- array(targets$ibra_extirpation, c(dim(targets$ibra_extirpation), 2))
  ibra_extirpation_error2 <- sqrt(rowMeans(
    ((ibra_extirpation < t(target_CI[1,,]))*(ibra_extirpation - t(target_CI[1,,])) +
       (ibra_extirpation > t(target_CI[2,,]))*(ibra_extirpation - t(target_CI[2,,])))^2,
    na.rm = TRUE))
  cbind(ibra_extirpation_error = summary_metric_data$ibra_extirpation_error, 
        ibra_extirpation_error2) # examine
}

# Load full example metrics
if (DEMONSTRATION) { 
  summary_metric_data <- poems::thylacine_example_metrics
  dim(summary_metric_data) # dimensions
}

# Calculate the error from the CI of total extinction
extinct <- summary_metric_data$total_extinction
target_CI <- targets$total_extinction
summary_metric_data$total_extinction_error <-  
  ((extinct < target_CI[1])*(extinct - target_CI[1]) + 
     (extinct > target_CI[2])*(extinct - target_CI[2]))
head(summary_metric_data) # examine

## ----message = FALSE----------------------------------------------------------
# Create a validator for selecting the 'best' example models
validator <- Validator$new(
  simulation_parameters = sample_data,
  simulation_summary_metrics = summary_metric_data[c("bounty_slope_error",
                                                     "ibra_extirpation_error",
                                                     "total_extinction_error")],
  observed_metric_targets = c(bounty_slope_error = 0,
                              ibra_extirpation_error = 0,
                              total_extinction_error = 0),
  non_finite_replacements = list(total_extinction_error = function(x) {
    (1968 - targets$total_extinction[2])}),
  output_dir = OUTPUT_DIR)
suppressWarnings(validator$run(tolerance = 0.01, sizenet = 1, lambda = 0.0001,
                               output_diagnostics = TRUE))
dir(OUTPUT_DIR, "*.pdf") # plus validation diagnostics (see abc library documentation)
head(validator$selected_simulations) # examine
dim(validator$selected_simulations) # dimensions
selected_indices <- validator$selected_simulations$index
selected_weights <- validator$selected_simulations$weight

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Load pre-generated example matrices
if (DEMONSTRATION) {
  summary_matrix_list <- poems::thylacine_example_matrices
} else { # cell extirpation and total/ibra bounty for selected samples only
  summary_matrix_list$extirpation <- summary_matrix_list$extirpation[selected_indices,]
  summary_matrix_list$total_bounty <- summary_matrix_list$total_bounty[selected_indices,]
  summary_matrix_list$ibra_bounty <- summary_matrix_list$ibra_bounty[selected_indices,]
}
lapply(summary_matrix_list, dim) # dimensions

# Plot the simulation, targets, and selected metrics for bounty slopes
bounty_slope <- summary_matrix_list$bounty_slope
colnames(bounty_slope) <- slope_intervals
graphics::boxplot(bounty_slope, border = rep("gray", 3), outline = FALSE, range = 0, 
                  col = NULL, ylim = c(-35, 20), main = "Thylacine bounty slope", 
                  xlab = "Slope interval (years)", ylab = "Bounty regression slope")
graphics::boxplot(cbind(bounty_slope[selected_indices,], NA), # NA for relative width
                  border = rep("blue", 3), width = c(rep(0.5, 3), 1), outline = FALSE, 
                  range = 0, add = TRUE, col = NULL)
graphics::points(x = 1:3, y = targets$bounty_slope, col = "red", pch = 15)
legend("topright", c("All", "Target", "Selected"), fill = c("gray", "red", "blue"), 
       border = NA, cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 10----
# Plot the simulation, targets, and selected metrics for extirpation
extirpation <- cbind(summary_matrix_list$ibra_extirpation, 
                     summary_metric_data$total_extinction)
colnames(extirpation) <- c(as.character(ibra_data$abbr), "Total")
graphics::boxplot(extirpation[, 10:1], border = rep("gray", 3), horizontal = TRUE, 
                  outline = FALSE, range = 0, col = NULL, ylim = c(1888,1968), 
                  main = "Thylacine model IBRA extirpation and total extinction",
                  xlim = c(0.5, 11), xlab = "Extirpation/extinction time (year)", 
                  ylab = "IBRA bioregion/Tasmania-wide total")
graphics::boxplot(cbind(extirpation[selected_indices, 10:1], NA),  horizontal = TRUE, 
                  ylim = c(1888,1968),border = rep("blue", 10), 
                  width = c(rep(0.5, 10), 1), outline = FALSE, range = 0, 
                  add = TRUE, col = NULL)
for (i in 1:9) {
  graphics::points(x = targets$ibra_extirpation[, i], y = rep(11 - i, 2), 
                   col = "red", pch = 15)
  graphics::lines(x = targets$ibra_extirpation[, i], y = rep(11 - i, 2),
                  col = "red", lwd = 2)
}
graphics::points(x = targets$total_extinction, y = c(1, 1), col = "red", pch = 15)
graphics::lines(x = targets$total_extinction, y = c(1, 1), col = "red", lwd = 2)
legend("topright", c("All", "Target (CI)", "Selected"), fill = c("gray", "red", "blue"), 
       border = NA, cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Allee effect
hist(sample_data$density_allee, breaks = 30, main = "Model ensemble Allee effect", 
     xlab = "Allee effect", ylim = c(0, 1000), col = "gray", yaxt = "n")
hist(rep(sample_data$density_allee[selected_indices], 20), breaks = 20, col = "blue", add = TRUE)
legend("topright", c("All", "Selected"), fill = c("gray", "blue"), cex = 0.8)

# Harvest catchability
hist(sample_data$harvest_q, breaks = 30, main = "Model ensemble harvest catchability", 
     xlab = "Harvest catchability (q)", col = "gray") #, yaxt = "n")
hist(rep(sample_data$harvest_q[selected_indices], 20), breaks = 20, col = "blue", add = TRUE)
legend("topright", c("All", "Selected"), fill = c("gray", "blue"), cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 7-----
plot(x = sample_data$harvest_ohr, y = sample_data$harvest_q, ylim = c(0, 0.0055),
     xlab = "Opportunistic harvest rate", ylab = "Bio-economic harvest catchability",
     main = "Opportunistic harvest vs. catchability", col = "gray")
points(x = sample_data$harvest_ohr[selected_indices], 
     y = sample_data$harvest_q[selected_indices], col = "blue", pch = 3)
graphics::legend("topright", legend = c("All samples", "Selected"),
                 col = c("gray", "blue"), pch = c(1, 3), cex = 0.8)

## ----message = FALSE----------------------------------------------------------
# Run replicates of 10 for each selected model
sample_data_rerun <- cbind(sample = 1:nrow(sample_data), sample_data)
sample_data_rerun <- cbind(sample_data_rerun[rep(selected_indices, each = 10),],
                           rerun = rep(1:10, length(selected_indices)))
head(sample_data_rerun) # examine
if (DEMONSTRATION) {
  sim_manager$sample_data <- sample_data_rerun[1:2,]
} else {
  sim_manager$sample_data <- sample_data_rerun
}
sim_manager$results_filename_attributes <- c("sample", "rerun")
run_output <- sim_manager$run()
run_output$summary
if (DEMONSTRATION) {
  dir(OUTPUT_DIR, "*.RData") # includes 2 new result files
}

# Collate summary metrics for re-runs
if (DEMONSTRATION) {
  results_manager$sample_data <- sample_data_rerun[1:2,]
} else {
  results_manager$sample_data <- sample_data_rerun
}
results_manager$summary_matrices <- c("bounty_slope", "ibra_extirpation")
results_manager$results_filename_attributes <- c("sample", "rerun")
gen_output <- results_manager$generate()
gen_output$summary
if (DEMONSTRATION) {
  results_manager$summary_metric_data # examine demo
}
if (DEMONSTRATION) { # load full example metrics and (some) matrices
  summary_metric_data_rerun <- poems::thylacine_example_metrics_rerun
  summary_matrix_list_rerun <- poems::thylacine_example_matrices_rerun
} else {
  summary_metric_data_rerun <- results_manager$summary_metric_data
  summary_matrix_list_rerun <- results_manager$summary_matrix_list
}
head(summary_metric_data_rerun) # examine
dim(summary_metric_data_rerun) # dimensions
lapply(summary_matrix_list_rerun, dim) # dimensions

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Bounty slope error
bounty_slope_error <- summary_metric_data$bounty_slope_error
hist(bounty_slope_error, breaks = 140, main = "Thylacine bounty slope error",
     xlim = c(0, 50), xlab = "Bounty slope RMSE", col = "gray", yaxt = "n")
bounty_slope_error_r <- summary_metric_data_rerun$bounty_slope_error
hist(rep(bounty_slope_error_r, 2), breaks = 20, col = "gold3", add = TRUE)
hist(rep(bounty_slope_error[selected_indices], 5), breaks = 12,
     col = "blue", add = TRUE)
lines(x = rep(0, 2), y = c(0, 10000), col = "red", lwd = 2)
legend("topright", c("All", "Target", "Selected", "Replicates"),
       fill = c("gray", "red", "blue", "gold3"), cex = 0.8)

# IBRA extirpation error
ibra_extirpation_error <- summary_metric_data$ibra_extirpation_error
hist(bounty_slope_error, breaks = 140, main = "Thylacine IBRA extirpation error",
     xlim = c(0, 50), xlab = "IBRA extirpation RMSE", col = "grey", yaxt = "n")
ibra_extirpation_error_r <- summary_metric_data_rerun$ibra_extirpation_error
hist(rep(ibra_extirpation_error_r, 2), breaks = 50, col = "gold3", add = TRUE)
hist(rep(ibra_extirpation_error[selected_indices], 5), breaks = 20, col = "blue",
     add = TRUE)
lines(x = rep(0, 2), y = c(0, 10000), col = "red", lwd = 2)
legend("topright", c("All", "Target", "Selected", "Replicates"),
       fill = c("grey", "red", "blue", "gold3"), cex = 0.8)

# Extinction time
extinction_time <- summary_metric_data$total_extinction
persistent_number <- length(which(is.na(extinction_time)))
extinction_time_finite <- extinction_time[!is.na(extinction_time)]
extinction_time_modified <- 
hist(c(extinction_time_finite, rep(1968, round(persistent_number/10))), breaks = 81, 
     main = "Thylacine extinction", xlim = c(1888, 1968), 
     xlab = "Extinction time (year)", col = "gray40", yaxt = "n")
hist(extinction_time_finite, breaks = 81, col = "gray", add = TRUE)
extinction_time_rerun <- summary_metric_data_rerun$total_extinction
extinction_time_rerun[which(is.na(extinction_time_rerun))] <- 1968
hist(rep(extinction_time_rerun, 2), breaks = 50, col = "gold3", add = TRUE)
hist(rep(extinction_time[selected_indices], 5), breaks = 28, col = "blue", add = TRUE)
lines(x = rep(1931, 2), y = c(0, 10000), col = "red", lwd = 2)
lines(x = rep(1937, 2), y = c(0, 10000), col = "red", lwd = 2)
legend("topleft", c("All", "(persistent/10)", "Target (CI)", "Selected", "Replicates"),
       fill = c("gray", "gray40", "red", "blue", "gold3"), cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Compare weighted ensemble model to actual historical bounty time-series
historic_bounty <- poems::thylacine_bounty_record
selected_bounty <- summary_matrix_list$total_bounty
weighted_bounty <- colSums(selected_bounty*selected_weights/sum(selected_weights))

# Plot the simulated and actual bounty
plot(x = 1888:1909, y = weighted_bounty[1:22], xlab = "Year",
     ylab = "Number of thylacines", main = "Thylacine bounty submitted",
     ylim = c(0, 200), type = "l", col = "blue", lwd = 2)
lines(x = 1888:1909, y = historic_bounty$Total, lty = 1, col = "red", lwd = 2)
legend("topright", legend = c("Model ensemble bounty", "Actual bounty"),
       col = c("blue", "red"), lty = c(1, 1), lwd = 2, cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 7, fig.height = 5-----
# Compare weighted ensemble model to actual historical bioregion bounty values
selected_bounty <- cbind(summary_matrix_list$ibra_bounty, 
                         rowSums(summary_matrix_list$ibra_bounty))
weighted_bounty <- colSums(selected_bounty*selected_weights/sum(selected_weights))
combined_bounty <- rbind(weighted_bounty, 
                         actual_bounty <- colSums(historic_bounty[, c(3:11, 2)]))
combined_bounty[, 10] <- combined_bounty[, 10]/2

# Comparative plot of simulated and historic IBRA/total bounty
barplot(combined_bounty, xlab = "IBRA bioregion/Tasmania-wide total", beside = TRUE,
        ylab = "Number of thylacines", col = c("blue", "red"),
        main = "Thylacine bounty submitted by region", border = NA)
legend("topleft", c("Model ensemble bounty", "Actual bounty", "Note: Total/2"), 
       fill = c("blue", "red", NA), border = NA, cex = 0.8)

## ----message = FALSE, fig.align = "center", fig.width = 5, fig.height = 5-----
# Calculate the weighted cell extirpation dates
selected_extirp <- summary_matrix_list$extirpation
weighted_extirpation <- colSums(selected_extirp*selected_weights/sum(selected_weights))

# Plot the weighted cell extirpation dates
extirpation_raster <- region$raster_from_values(weighted_extirpation)
raster::plot(extirpation_raster, main = "Thylacine model ensemble extirpation",
             xlab = "Longitude (degrees)", ylab = "Latitude (degrees)",
             zlim = c(1888, 1940), col = grDevices::heat.colors(100), colNA = "blue")

