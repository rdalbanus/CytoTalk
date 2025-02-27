% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{preprocess}
\alias{preprocess}
\title{Filter Lowly Expressed Genes}
\usage{
preprocess(proteins, type_a, type_b, cutoff_a, cutoff_b, dir_in, dir_out)
}
\arguments{
\item{proteins}{Character vector of protein-coding gene names}

\item{type_a}{Cell type A}

\item{type_b}{Cell type B}

\item{cutoff_a}{Minimum proportion of non-zero columns, per row, cell type A}

\item{cutoff_b}{Minimum proportion of non-zero columns, per row, cell type B}

\item{dir_in}{Input directory, contains scRNAseq files}

\item{dir_out}{Output directory}
}
\value{
None
}
\description{
Given cell types, cutoff values, and the names of protein-coding genes,
filter scRNAseq data so that only protein-coding genes with expression levels
over the cutoff remain.
}
\examples{
\dontrun{
proteins <- CytoTalk::pcg_mouse
type_a <- "BCells"
type_b <- "TCells"
cutoff_a <- 0.1
cutoff_b <- 0.1
dir_in <- "~/scRNAseq-data"
dir_out <- "~/CytoTalk-output"

preprocess(
    proteins, type_a, type_b, cutoff_a, cutoff_b, dir_in, dir_out
)
}

}
