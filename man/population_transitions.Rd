% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population_transitions.R
\name{population_transitions}
\alias{population_transitions}
\title{Nested functions for stage-based population transitions.}
\usage{
population_transitions(
  populations,
  demographic_stochasticity,
  fecundity_matrix,
  fecundity_max,
  survival_matrix
)
}
\arguments{
\item{populations}{Number of populations.}

\item{demographic_stochasticity}{Boolean for choosing demographic stochasticity for transitions.}

\item{fecundity_matrix}{Matrix of transition fecundity rates (Leslie/Lefkovitch matrix with non-zero fecundities only).}

\item{fecundity_max}{Maximum transition fecundity rate (in Leslie/Lefkovitch matrix).}

\item{survival_matrix}{Matrix of transition survival rates (Leslie/Lefkovitch matrix with non-zero survivals only).}
}
\value{
Transition calculation function: \code{function(fecundity_array, survival_array, stage_abundance, occupied_indices)}, where:
  \describe{
    \item{\code{fecundity_array}}{3D array of fecundity rates (\emph{stages} rows by \emph{stages} columns by \emph{populations} deep).}
    \item{\code{survival_array}}{3D array of survival rates (\emph{stages} rows by \emph{stages} columns by \emph{populations} deep).}
    \item{\code{stage_abundance}}{Matrix of stage abundances for each population at time step (\emph{stages} rows by \emph{populations} columns).}
    \item{\code{occupied_indices}}{Array of indices for those populations occupied.}
    \item{\code{returns}}{Transitioned stage abundances.}
  }
}
\description{
Modular functions for the population simulator for performing staged-based
(Leslie/Lefkovitch matrix) transitions via 3D survival and fecundity arrays.
}
