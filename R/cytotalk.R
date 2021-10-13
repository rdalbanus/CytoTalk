#' CytoTalk: Construct signal transduction networks
#'
#' The CytoTalk package allows for *de novo* construction of a signaling
#' network (pathways from ligand-receptor pairs) between two cell types using
#' single- cell transcriptomics data (scRNA-seq).
#'
#' @section Running CytoTalk:
#' The central wrapper function, `run_cytotalk`, is
#' the main entry point to the package. Most users should start here. Check out
#' its MAN page by executing `?run_cytotalk` in your R console.
#'
#' @section Advanced Users:
#' If you're willing to do some digging, CytoTalk is
#' composed of many steps, some of which could be interchangeable with other
#' computational methods. For example, let's say that you have a different idea
#' for computing intra- cellular similarity (i.e. not a mutual information
#' matrix). You could skip steps 1-3 and attempt to run from step 4 onward
#' (simply execute `run_cytotalk` in your console to view the source code) with
#' differently integrated data. Currently, doing so is not user-friendly, so
#' this is recommended for developers only. However, it would be interesting to
#' make the different components of the overall process more modular to compare
#' sub processes.
#'
#' @docType package
#'
#' @name cytotalk
NULL

#' Run CytoTalk Process
#'
#' Runs all sub processes of CytoTalk. At minimum, the folder that contains
#' your scRNA-seq data, `dir_in`, and the two cell types corresponding to that
#' data, `type_a` and `type_b`, must be specified. It is also recommended you
#' set the output directory, `dir_out`. If you have data from a species other
#' that mice, please read the descriptions for the `proteins` and `ligands`
#' arguments. If you would like a more selective run (i.e. smaller, filters out
#' more lowly- expressed genes), please increase both cutoff values `cutoff_a`
#' and `cutoff_b`. Lastly, intermediate results will be reused, so if you want
#' a fresh run, either output to a new output folder or delete the files in a
#' particular output folder.
#'
#' @examples \dontrun{
#' # set location of data and output folder
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' # examine if your filenames are correct
#' check_valid_names(dir_in)
#'
#' # select which cell types to compare
#' type_a <- "BCells"
#' type_b <- "TCells"
#'
#' # if desired, run a small test with highly-expressed genes
#' run_cytotalk(type_a, type_b, dir_in, dir_out = "my-test",
#'              cutoff_a = 0.75, cutoff_b = 0.75)
#'
#' # finally, run the full process
#' run_cytotalk(type_a, type_b, dir_in, dir_out)
#' }
#'
#' @param type_a Name of cell type A that matches scRNA-seq file; for example,
#' `"BCells"`
#'
#' @param type_b Name of cell type B that matches scRNA-seq file; for example,
#' `"TCells"` 
#'
#' @param dir_in Folder containing scRNA-seq data
#'
#' @param dir_out Folder used for output; if not specified, a "CytoTalk-output"
#' folder will be generated 
#'
#' @param proteins A character vector, contains the names of protein coding
#' genes; by default, uses the `pcg_human` data. This package also includes
#' `pcg_mouse`, but you can also use your own data 
#'
#' @param ligands A dataframe or matrix object with two columns, ligands names
#' and the names of their receptors; by default, uses the `ligands_human` data.
#' This package also includes `ligands_mouse`, but you can also use your own
#' data
#'
#' @param cutoff_a Proportional threshold for lowly expressed genes in cell
#' type A (range of \[0-1\]); for example, 0.1 means genes with some expression
#' in at least 10% of cells are retained 
#'
#' @param cutoff_b Proportional expression threshold for cell type B (range of
#' \[0-1\]) 
#'
#' @param beta_max Upper limit of the test values of the PCSF objective
#' function parameter $I^2$, which is inversely proportional to the total
#' number of genes in a given cell-type pair; suggested to be 100 (default) if
#' the total number of genes in a given cell-type pair is above 10,000; if the
#' total number of genes is below 5,000, increase to 500 
#'
#' @param omega_min Start point of omega range; omega represents the edge cost
#' of the artificial network, but has been found to be less significant than
#' beta. Recommended minimum of `0.5`.
#'
#' @param omega_max End point of range between `omega_min` and `omega_max`,
#' step size of `0.1`. Recommended maximum of `1.5`. 
#'
#' @param depth Starting at each ligand-receptor pair in the resultant network,
#' how many steps out from that pair should be taken to generate each
#' neighborhood?
#'
#' @return None
#'
#' @export
run_cytotalk <- function(
    type_a, type_b, dir_in,
    dir_out="CytoTalk-output",
    proteins=CytoTalk::pcg_human,
    ligands=CytoTalk::ligands_human,
    cutoff_a=0.1, cutoff_b=0.1,
    beta_max=100, omega_min=0.5, omega_max=0.5,
    depth=3, cores=NULL) {

    # must have valid data directory
    type_names <- check_valid_names(dir_in)

    # check cell type names as well
    if (!all(c(type_a, type_b) %in% type_names)) {
        stop("one of the cell types cannot be found in the input directory")
    }

    # absolute paths to directories
    dir_in <- suppressWarnings(normalizePath(dir_in))
    dir_out <- suppressWarnings(normalizePath(dir_out))

    # make sure output directory exists
    if (!dir.exists(dir_out)) {
        dir.create(dir_out)
    }

    # TODO: could have a more robust check on input directory,
    # although errors will be determined in step 1 anyways...

    # filter out lowly expressed genes,
    # compute non-self-talk score (within each cell type),
    # compute preferential expression measure
    tick(1, "Preprocessing...")
    preprocess(proteins, type_a, type_b, cutoff_a, cutoff_b, dir_in, dir_out)
    compute_non_self_talk(ligands, type_a, type_b, dir_in, dir_out)
    compute_pem(dir_in, dir_out)

    # compute mutual information (within types)
    tick(2, "Mutual information matrix...")
    compute_mutual_information(dir_out)

    # use ARACNE.m to filter out indirect edges
    tick(3, "Indirect edge-filtered network...")
    filter_indirect_edges(dir_out)

    # integrate nodes and their prizes,
    # as well as edges and their costs
    tick(4, "Integrate network...")
    integrate_network(ligands, type_a, type_b, dir_in, dir_out)

    # run prize-collecting Steiner tree algorithm from Python
    tick(5, "PCSF...")
    run_pcst(dir_out, beta_max, omega_min, omega_max)

    # run Kolmogorov-Smirnov tests
    tick(6, "Determine best signaling network...")
    generate_signaling_network(dir_out)

    # generate SIF and SVG files
    tick(7, "Generate network output...")
    generate_pathways(type_a, type_b, dir_out, depth)

    # no output...
    NULL
}
