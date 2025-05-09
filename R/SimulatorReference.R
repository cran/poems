#' R6 class for a simulator reference
#'
#' @description
#' \code{\link[R6:R6Class]{R6}} class for dynamically attaching simulator attributes
#' and results (passed by reference).
#'
#' @examples
#' test_class <- SimulatorReference$new()
#' test_class$attached$attr1 <- "example1"
#' test_class$results$attr1 <- "example2"
#' str(test_class)
#'
#' @importFrom R6 R6Class
#' @export SimulatorReference

SimulatorReference <- R6Class(
  "SimulatorReference",
  public = list(
    ## Attributes ##

    #' @field attached A list of dynamically attached simulator attributes (name-value pairs).
    attached = list(),

    #' @field results A list of dynamically accessed simulator results (name-value pairs).
    results = list()
  ), # end public
)
