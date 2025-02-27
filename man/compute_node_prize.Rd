% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrate_gene_network.R
\name{compute_node_prize}
\alias{compute_node_prize}
\title{Compute Node Prize}
\usage{
compute_node_prize(type_a, type_b, dir_out)
}
\arguments{
\item{type_a}{Cell type A}

\item{type_b}{Cell type B}

\item{dir_out}{Output directory}
}
\value{
None
}
\description{
For each cell type, compute overall node prize (integrate gene relevance and
PEM scores).
}
\examples{
{
type_a <- "BCells"
type_b <- "TCells"
dir_out <- "~/CytoTalk-output"

compute_node_prize(type_a, type_b, dir_out)
}

}
