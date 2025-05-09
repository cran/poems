% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenericManager.R
\name{GenericManager}
\alias{GenericManager}
\title{R6 class representing a generic manager.}
\description{
\code{\link[R6:R6Class]{R6}} class to represent a generic (abstract) manager for
generating or processing simulation results, as well as optionally generating values
via generators.
}
\examples{
generic_manager <- GenericManager$new(
  attr1 = 22:23,
  results_filename_attributes = c("attr1", "example")
)
generic_manager$get_results_filename(1)
generic_manager$get_results_filename(2)

}
\section{Super class}{
\code{\link[poems:GenericClass]{poems::GenericClass}} -> \code{GenericManager}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{attached}}{A list of dynamically attached attributes (name-value pairs).}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{sample_data}}{A data frame of sampled parameters for each simulation/result.}

\item{\code{generators}}{A list of generators (\code{\link{Generator}} or
inherited class) objects for generating simulation model values.}

\item{\code{parallel_cores}}{Number of cores for running the simulations in parallel.}

\item{\code{results_dir}}{Results directory path.}

\item{\code{results_ext}}{Result file extension (default is .RData).}

\item{\code{results_filename_attributes}}{A vector of: prefix (optional);
attribute names (from the sample data frame); postfix (optional);
utilized to construct results filenames.}

\item{\code{error_messages}}{A vector of error messages encountered.}

\item{\code{warning_messages}}{A vector of warning messages encountered.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-GenericManager-new}{\code{GenericManager$new()}}
\item \href{#method-GenericManager-get_attribute}{\code{GenericManager$get_attribute()}}
\item \href{#method-GenericManager-get_message_sample}{\code{GenericManager$get_message_sample()}}
\item \href{#method-GenericManager-get_results_filename}{\code{GenericManager$get_results_filename()}}
\item \href{#method-GenericManager-clone}{\code{GenericManager$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="poems" data-topic="GenericClass" data-id="new_clone"><a href='../../poems/html/GenericClass.html#method-GenericClass-new_clone'><code>poems::GenericClass$new_clone()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GenericManager-new"></a>}}
\if{latex}{\out{\hypertarget{method-GenericManager-new}{}}}
\subsection{Method \code{new()}}{
Initialization method sets any included attributes (\emph{sample_data},
\emph{generators}, \emph{parallel_cores}, \emph{results_dir},
\emph{results_filename_attributes}) and attaches other attributes
individually listed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GenericManager$new(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Parameters listed individually.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GenericManager-get_attribute"></a>}}
\if{latex}{\out{\hypertarget{method-GenericManager-get_attribute}{}}}
\subsection{Method \code{get_attribute()}}{
Returns a named manager or attached attribute.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GenericManager$get_attribute(param)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{param}}{Character string name of the attribute.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Selected attribute value.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GenericManager-get_message_sample"></a>}}
\if{latex}{\out{\hypertarget{method-GenericManager-get_message_sample}{}}}
\subsection{Method \code{get_message_sample()}}{
Substitutes the specified sample details into a status message (using sprintf) and returns the result.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GenericManager$get_message_sample(status_message, sample_index)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{status_message}}{Character string message with a placeholder for sample details.}

\item{\code{sample_index}}{Row index of sample data frame containing details of substitution parameters.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Status message with substituted sample details.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GenericManager-get_results_filename"></a>}}
\if{latex}{\out{\hypertarget{method-GenericManager-get_results_filename}{}}}
\subsection{Method \code{get_results_filename()}}{
Constructs and returns the results filename based on the sample data
frame index and results filename attributes.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GenericManager$get_results_filename(sample_index)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sample_index}}{Row index of sample data frame containing details
of substitution parameters.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Results filename with substituted sample details.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GenericManager-clone"></a>}}
\if{latex}{\out{\hypertarget{method-GenericManager-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GenericManager$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
