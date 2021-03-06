% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_state.R
\name{expr_freq_filter}
\alias{expr_freq_filter}
\title{Filter differential state analysis results by gene expression frequencies.}
\usage{
expr_freq_filter(results, sce, percent = 0.1)
}
\arguments{
\item{results}{Results returned by \code{\link[muscat]{pbDS}} or
\code{\link{pseudobulk_diff_state}}}

\item{sce}{A SingleCellExperiment object. The following
cell metadata (\code{colData}) columns have to be provided:
\describe{
\item{\code{sample_id}}{unique sample identifiers}
\item{\code{cluster_id}}{subpopulation (cluster) assignments}
\item{\code{group_id}}{experimental group/condition}
}}

\item{percent}{Minimum percentage of cells in which a gene expressed.}
}
\value{
Same as \code{tbl_list}
}
\description{
Only retain genes that are expressed in an average of given percentage of
cells in at least 1 group.
}
