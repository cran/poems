#' Nested functions for a user-defined population abundance (and capacity) transformation.
#'
#' Modular functions for the population simulator for performing a transformation of
#' the stage abundance (and optionally carrying capacity) at a specified time step via
#' a user-defined function.
#'
#' @examples
#' simulator <- SimulatorReference$new()
#' # Example transformation: a random population is chosen for a severe disturbance event
#' # (only life cycle stage 3 individuals survive)
#' disturbance_function <- function(params) {
#'   params$simulator$attached$params <- params # attach to reference object
#'   random_population <- sample(params$occupied_indices, 1)
#'   new_stage_abundance <- params$stage_abundance
#'   new_stage_abundance[1:2, random_population] <- 0
#'   return(new_stage_abundance)
#' }
#' transformation_function <- population_transformation(
#'   replicates = 4, time_steps = 10, years_per_step = 1,
#'   populations = 7, demographic_stochasticity = TRUE,
#'   density_stages = c(0, 1, 1), transformation = disturbance_function,
#'   simulator
#' )
#' carrying_capacity <- rep(10, 7)
#' carrying_capacity <- rep(10, 7)
#' stage_abundance <- matrix(c(
#'   7, 13, 0, 26, 0, 39, 47,
#'   2, 0, 6, 8, 0, 12, 13,
#'   0, 3, 4, 6, 0, 9, 10
#' ), nrow = 3, ncol = 7, byrow = TRUE)
#' occupied_indices <- (1:7)[-5]
#' transformation_function(
#'   r = 2, tm = 6, carrying_capacity, stage_abundance,
#'   occupied_indices
#' )
#'
#' @param replicates Number of replicate simulation runs.
#' @param time_steps Number of simulation time steps.
#' @param years_per_step Number of years per time step.
#' @param populations Number of populations.
#' @param demographic_stochasticity Boolean for optionally choosing demographic stochasticity for the transformation.
#' @param density_stages Array of booleans or numeric (0,1) for each stage to indicate which stages are affected by density.
#' @param transformation A user-defined function (optionally nested in a list with additional attributes) for performing transformation: \code{function(params)}, where \emph{params} is a list passed to the function containing:
#'   \describe{
#'     \item{\code{replicates}}{Number of replicate simulation runs.}
#'     \item{\code{time_steps}}{Number of simulation time steps.}
#'     \item{\code{years_per_step}}{Number of years per time step.}
#'     \item{\code{populations}}{Number of populations.}
#'     \item{\code{stages}}{Number of life cycle stages.}
#'     \item{\code{demographic_stochasticity}}{Boolean for optionally choosing demographic stochasticity for the transformation.}
#'     \item{\code{density_stages}}{Array of booleans or numeric (0,1) for each stage to indicate which stages are affected by density.}
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{stage_abundance}}{Matrix of (current) abundance for each stage (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied at (current) time step.}
#'     \item{\code{simulator}}{\code{\link{SimulatorReference}} object with dynamically accessible \emph{attached} and \emph{results} lists.}
#'     \item{\code{additional attributes}}{Additional attributes when the transformation is optionally nested in a list.}
#'   }
#'   returns a transformed stage abundance matrix (or a list with stage abundance and carrying capacity)
#' @param simulator \code{\link{SimulatorReference}} object with dynamically accessible \emph{attached} and \emph{results} lists.
#' @param name Optional name for the transformation (default is "transformation").
#' @return Abundance (and capacity) transformation function: \code{function(r, tm, carrying_capacity, stage_abundance, occupied_indices)}, where:
#'   \describe{
#'     \item{\code{r}}{Simulation replicate.}
#'     \item{\code{tm}}{Simulation time step.}
#'     \item{\code{carrying_capacity}}{Array of carrying capacity values for each population at time step.}
#'     \item{\code{stage_abundance}}{Matrix of abundance for each stage (rows) and population (columns) at time step.}
#'     \item{\code{occupied_indices}}{Array of indices for populations occupied at time step.}
#'     \item{\code{returns}}{List with transformed stage abundance matrix (and optionally carrying capacity).}
#'   }
#' @export population_transformation

population_transformation <- function(
  replicates,
  time_steps,
  years_per_step,
  populations,
  demographic_stochasticity,
  density_stages,
  transformation,
  simulator,
  name = "transformation"
) {
  if (!is.null(transformation)) {
    # Derive stages
    stages <- length(density_stages)

    # Unpack transformation function and additional attributes from a list
    additional_attributes <- list()
    if (is.list(transformation)) {
      function_index <- which(unlist(lapply(transformation, is.function)))
      additional_attributes <- transformation[-function_index]
      transformation <- transformation[[function_index]]
    }

    # List of parameters to pass to the user-defined function
    params <- c(
      list(
        replicates = replicates,
        time_steps = time_steps,
        years_per_step = years_per_step,
        populations = populations,
        stages = stages,
        demographic_stochasticity = demographic_stochasticity,
        density_stages = density_stages,
        simulator = simulator
      ),
      additional_attributes
    )

    if (is.function(transformation)) {
      # user-defined function

      ## Create a nested function for applying user-defined transformation to stage abundance ##
      user_defined_function <- function(
        r,
        tm,
        carrying_capacity,
        stage_abundance,
        occupied_indices
      ) {
        # Add attributes to be made available to the user-defined function
        params$r <- r
        params$tm <- tm
        params$carrying_capacity <- carrying_capacity
        params$stage_abundance <- stage_abundance
        params$occupied_indices <- occupied_indices

        # Run user-defined transformation function
        tryCatch(
          {
            transformed <- transformation(params)
          },
          error = function(e) {
            stop(
              paste(
                "Error produced within user-defined",
                name,
                "function:",
                as.character(e)
              ),
              call. = FALSE
            )
          }
        )

        # Resolve and check returned transformed values
        if (is.list(transformed)) {
          if (
            !all(
              c("stage_abundance", "carrying_capacity") %in% names(transformed)
            )
          ) {
            stop(
              paste(
                "When returning a list, the user-defined",
                name,
                "function should contain stage_abundance and carrying_capacity"
              ),
              call. = FALSE
            )
          }
        } else {
          # assume stage abundance matrix and place a list
          transformed <- list(stage_abundance = transformed)
        }

        # Warn if any negative or non-finite
        if (any(!is.finite(transformed$stage_abundance))) {
          warning(
            paste(
              "Non-finite stage abundances returned by user-defined",
              name,
              "function"
            ),
            call. = FALSE
          )
        }
        if (
          "carrying_capacity" %in%
            names(transformed) &&
            any(!is.finite(transformed$carrying_capacity))
        ) {
          warning(
            paste(
              "Non-finite carrying capacities returned by user-defined",
              name,
              "function"
            ),
            call. = FALSE
          )
        }
        if (
          any(
            transformed$stage_abundance[which(is.finite(
              transformed$stage_abundance
            ))] <
              0
          )
        ) {
          warning(
            paste(
              "Negative stage abundances returned by user-defined",
              name,
              "function"
            ),
            call. = FALSE
          )
        }
        if (
          "carrying_capacity" %in%
            names(transformed) &&
            any(
              transformed$carrying_capacity[which(is.finite(
                transformed$carrying_capacity
              ))] <
                0
            )
        ) {
          warning(
            paste(
              "Negative carrying capacities returned by user-defined",
              name,
              "function"
            ),
            call. = FALSE
          )
        }

        return(transformed)
      }
      environment(user_defined_function)[["name"]] <- name

      return(user_defined_function)
    }
  }
}
