#' @noRd
noise <- function(n) {
    rng <- seq_len(n)
    (sample(rng, 1) == rng) * 1e-20
}

#' @noRd
compute_non_self_talk_type <- function(ligands, type, letter, dir_in, dir_out) {
    # format filename, join to full path
    fpath_scrna <- file.path(dir_in, sprintf("scRNAseq_%s.csv", type))
    fpath_out <- file.path(dir_out, sprintf("NonSelfTalkSco_Typ%s.txt", letter))

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_scrna)) {
        stop(sprintf("cannot find input file: %s", fpath_scrna))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load tumor expression matrix
    mat <- suppressMessages(vroom::vroom(fpath_scrna, progress = FALSE))
    mat <- tibble::column_to_rownames(mat, names(mat)[1])

    # casefold upper
    ligands <- apply(ligands, 2, toupper)
    rownames(mat) <- toupper(rownames(mat))

    # calculate number of bins
    n_cols <- ncol(mat)
    n_bins <- sqrt(n_cols)

    # add noise to rows with all zeros
    is_zero <- rowSums(mat) == 0
    mat[is_zero,] <- do.call(rbind, lapply(rep(n_cols, sum(is_zero)), noise))

    # match ligands to rownames
    index <- apply(ligands, 2, match, rownames(mat))

    # filter by rows without any NAs
    valid <- rowSums(is.na(index)) == 0
    index_valid <- index[valid,]
    ligands_valid <- ligands[valid,]

    # detect and register cores
    if (is.null(cores)) {
        message("No number of cores were specified - using all possible.")
        cores <- max(1, parallel::detectCores() - 2)
    } else {
        message(sprintf("Computing with %s cores"), cores)
    }
    
    doParallel::registerDoParallel(cores = cores)

    # parallel loop for MI distances
    i <- NULL
    score <- foreach::`%dopar%`(
        foreach::foreach(i = seq_len(nrow(index_valid)), .combine = rbind), {

        exp1 <- as.numeric(mat[index_valid[i, 1],])
        exp2 <- as.numeric(mat[index_valid[i, 2],])

        y2d <- entropy::discretize2d(exp1, exp2, n_bins, n_bins)
        H1 <- entropy::entropy(rowSums(y2d), method = "MM")
        H2 <- entropy::entropy(colSums(y2d), method = "MM")
        H12 <- entropy::entropy(y2d, method = "MM")
        mi <- H1 + H2 - H12

        norm <- min(c(H1, H2))
        mi_dist <- -log10(mi / norm)

        lig1 <- ligands_valid[i, 1]
        lig2 <- ligands_valid[i, 2]

        data.frame(lig1, lig2, mi_dist)
    })

    # write out
    vroom::vroom_write(score, fpath_out, progress = FALSE)

    # unregister cores
    doParallel::stopImplicitCluster()
    NULL
}

#' Compute Non-Self-Talk Score
#'
#' Given cell types and the names of ligand-receptor pairs, compute the
#' non-self-talk score (within each cell type).
#'
#' @examples {
#' ligands <- CytoTalk::ligands_mouse
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' compute_non_self_talk(ligands, type_a, type_b, dir_in, dir_out)
#' }
#'
#' @param ligands Character matrix of ligand-receptor pair names
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_in Input directory, contains scRNAseq files
#' @param dir_out Output directory
#' @return None
#' @export
compute_non_self_talk <- function(ligands, type_a, type_b, dir_in, dir_out) {
    compute_non_self_talk_type(ligands, type_a, "A", dir_in, dir_out)
    compute_non_self_talk_type(ligands, type_b, "B", dir_in, dir_out)
    NULL
}
