% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cytotalk.R
\name{run_cytotalk}
\alias{run_cytotalk}
\title{Run CytoTalk Process}
\usage{
run_cytotalk(
  type_a,
  type_b,
  dir_in,
  dir_out = "CytoTalk-output",
  proteins = CytoTalk::pcg_human,
  ligands = CytoTalk::ligands_human,
  cutoff_a = 0.1,
  cutoff_b = 0.1,
  beta_max = 100,
  omega_min = 0.5,
  omega_max = 0.5,
  depth = 3
)
}
\arguments{
\item{type_a}{Name of cell type A that matches scRNA-seq file; for example,
\code{"BCells"}}

\item{type_b}{Name of cell type B that matches scRNA-seq file; for example,
\code{"TCells"}}

\item{dir_in}{Folder containing scRNA-seq data}

\item{dir_out}{Folder used for output; if not specified, a "CytoTalk-output"
folder will be generated}

\item{proteins}{A character vector, contains the names of protein coding
genes; by default, uses the \code{pcg_human} data. This package also includes
\code{pcg_mouse}, but you can also use your own data}

\item{ligands}{A dataframe or matrix object with two columns, ligands names
and the names of their receptors; by default, uses the \code{ligands_human} data.
This package also includes \code{ligands_mouse}, but you can also use your own
data}

\item{cutoff_a}{Proportional threshold for lowly expressed genes in cell
type A (range of [0-1]); for example, 0.1 means genes with some expression
in at least 10\% of cells are retained}

\item{cutoff_b}{Proportional expression threshold for cell type B (range of
[0-1])}

\item{beta_max}{Upper limit of the test values of the PCSF objective
function parameter $I^2$, which is inversely proportional to the total
number of genes in a given cell-type pair; suggested to be 100 (default) if
the total number of genes in a given cell-type pair is above 10,000; if the
total number of genes is below 5,000, increase to 500}

\item{omega_min}{Start point of omega range; omega represents the edge cost
of the artificial network, but has been found to be less significant than
beta. Recommended minimum of \code{0.5}.}

\item{omega_max}{End point of range between \code{omega_min} and \code{omega_max},
step size of \code{0.1}. Recommended maximum of \code{1.5}.}

\item{depth}{Starting at each ligand-receptor pair in the resultant network,
how many steps out from that pair should be taken to generate each
neighborhood?}
}
\value{
None
}
\description{
Runs all sub processes of CytoTalk. At minimum, the folder that contains
your scRNA-seq data, \code{dir_in}, and the two cell types corresponding to that
data, \code{type_a} and \code{type_b}, must be specified. It is also recommended you
set the output directory, \code{dir_out}. If you have data from a species other
that mice, please read the descriptions for the \code{proteins} and \code{ligands}
arguments. If you would like a more selective run (i.e. smaller, filters out
more lowly- expressed genes), please increase both cutoff values \code{cutoff_a}
and \code{cutoff_b}. Lastly, intermediate results will be reused, so if you want
a fresh run, either output to a new output folder or delete the files in a
particular output folder.
}
\examples{
\dontrun{
# set location of data and output folder
dir_in <- "~/scRNAseq-data"
dir_out <- "~/CytoTalk-output"

# examine if your filenames are correct
check_valid_names(dir_in)

# select which cell types to compare
type_a <- "BCells"
type_b <- "TCells"

# if desired, run a small test with highly-expressed genes
run_cytotalk(type_a, type_b, dir_in, dir_out = "my-test",
             cutoff_a = 0.75, cutoff_b = 0.75)

# finally, run the full process
run_cytotalk(type_a, type_b, dir_in, dir_out)
}

}
