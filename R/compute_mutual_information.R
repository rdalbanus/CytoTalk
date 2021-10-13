#' @noRd
mi_matrix_fast <- function(mat) {
    # Miller Madow entropy
    entropy_mm <- function(x) {
        count <- table(x)
        freqs <- count / sum(count)
        shift <- (sum(count > 0) - 1) / (2 * sum(count))
        -sum(ifelse(freqs > 0, freqs * log(freqs), 0)) + shift
    }

    n <- ncol(mat)
    i <- NULL

    # calculate the entropy of each column first
    ent <- foreach::`%dopar%`(
        foreach::foreach(i = seq_len(n), .combine = "c"), {
        entropy_mm(mat[, i])
    })

    # use already caculated columns for mutual information
    res <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n)), {
        ent_i <- ent[i]
        vec <- vapply(i:n, function(j) {
            ent_i + ent[j] - entropy_mm(data.frame(mat[, i], mat[, j]))
        }, numeric(1))
        vec[vec < 0] <- 0
        vec
    })

    mat <- matrix(0, n, n)
    mat[lower.tri(mat, diag = TRUE)] <- unlist(res)
    mat + t(mat * lower.tri(mat))
}

#' @noRd
compute_mutual_information_type <- function(dir_out, type) {
    # format file path
    fpath_in <- file.path(dir_out, sprintf("Typ%sExp_rmRepGf_dm.txt", type))
    fpath_out <- file.path(dir_out, sprintf("MutualInfo_Typ%s_Para.txt", type))

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_in)) {
        stop(sprintf("cannot find input file: %s", fpath_in))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load in file
    mat <- suppressMessages(vroom::vroom(fpath_in, progress = FALSE))
    # set first column to rownames
    mat <- tibble::column_to_rownames(mat, names(mat)[1])

    # transpose matrix
    mat <- t(mat)

    # discretize continuous random variable
    mat <- infotheo::discretize(as.data.frame(mat))

    # compute mutual information matrix
    mat_out <- as.data.frame(mi_matrix_fast(mat))

    # write out
    rownames(mat_out) <- names(mat_out) <- names(mat)
    mat_out <- tibble::rownames_to_column(mat_out)
    vroom::vroom_write(mat_out, fpath_out, progress = FALSE)
    NULL
}


#' Compute Mutual Information Matrix
#'
#' Within each cell type, compute intracellular mutual information.
#'
#' @examples {
#' dir_out <- "~/CytoTalk-output"
#' compute_mutual_information(dir_out, cores = 2)
#' }
#'
#' @param dir_out Output directory
#' @param cores Sets the number of cores
#' @return None
#' @export
compute_mutual_information <- function(dir_out, cores=NULL) {
    # register cores for parallel computation
    if (is.null(cores)) {
        message("No number of cores were specified for MI - using all possible.")
        cores <- max(1, parallel::detectCores() - 2)
    } else {
        message(sprintf("Computing MI with %s cores", cores))
        cores <- max(1, cores)
    }
    doParallel::registerDoParallel(cores = cores)

    compute_mutual_information_type(dir_out, "A")
    compute_mutual_information_type(dir_out, "B")

    # unregister cores
    doParallel::stopImplicitCluster()
    NULL
}
