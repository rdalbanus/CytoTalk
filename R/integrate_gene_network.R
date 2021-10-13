# EXTRACT INTRACELLULAR EDGES
# (non-zero from ARACNE matrix)


#' @noRd
extract_intracell_edges_type <- function(dir_out, type) {
    # format filepaths
    fpath_in <- file.path(
        dir_out, sprintf("IntracellularNetwork_Type%s.txt", type))
    fpath_out <- file.path(dir_out, sprintf("MI_Typ%s.txt", type))

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
    mat <- tibble::column_to_rownames(mat, names(mat)[1])

    # make sure all nodes have an edge
    stopifnot(all(rowSums(mat) != 0))

    # create unique index
    bools <- (mat != 0) & lower.tri(mat)
    index <- which(bools, arr.ind = TRUE)
    value <- mat[index]

    # create new dataframe
    mat_out <- data.frame(
        row_ind = index[, 1],
        row = rownames(index),
        col_ind = index[, 2],
        col = colnames(mat)[index[, 2]],
        val = value
    )

    # write out
    vroom::vroom_write(mat_out, fpath_out, progress = FALSE)
    NULL
}


#' Extract Intracelluar Edges
#'
#' From the mutual information matrix, extract all non-zero edges and save them
#' in a data table.
#'
#' @examples {
#' dir_out <- "~/CytoTalk-output"
#' extract_intracell_edges(dir_out)
#' }
#'
#' @param dir_out Output directory
#' @return None
#' @export
extract_intracell_edges <- function(dir_out) {
    extract_intracell_edges_type(dir_out, "A")
    extract_intracell_edges_type(dir_out, "B")
    NULL
}


# PREFERENTIAL EXPRESSION MEASURE


#' @noRd
compute_gene_relevance_type <- function(ligands, dir_in, dir_out, type) {
    # create filepaths
    fpath_in <- file.path(
        dir_out, sprintf("IntracellularNetwork_Type%s.txt", type))
    fpath_out <- file.path(dir_out, sprintf("GeneRelevanceTyp%s.txt", type))

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_in)) {
        stop(sprintf("cannot find input file: %s", fpath_in))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load filtered intracellular network
    df_net <- suppressMessages(vroom::vroom(fpath_in, progress = FALSE))
    df_net <- tibble::column_to_rownames(df_net, names(df_net)[1])

    # load in species ligands
    df_lig <- apply(ligands, 2, toupper)

    # create a new matrix, same size
    mat <- as.matrix(df_net) * 0
    # set the diagonal as the network rowsums
    diag(mat) <- rowSums(df_net)

    # TODO: could use some work here!

    # negative square root
    mat_nsq <- corpcor::mpower(mat, -0.5)
    # symmetric matrix
    mat_wnorm <- mat_nsq %*% as.matrix(df_net) %*% mat_nsq

    names_lig <- unique(as.vector(df_lig))
    names_net <- toupper(names(df_net))
    vec_lig <- vapply(names_net, "%in%", integer(1), names_lig)

    # compute gene relevance value,
    # random walk with restart
    n <- 50
    alpha <- 0.9
    offset <- (1 - alpha) * vec_lig
    mat_relev <- as.matrix(vec_lig)

    vec_diff <- vector()
    for (i in seq_len(n)) {
        mat_old <- mat_relev
        mat_relev <- alpha * (mat_wnorm %*% mat_relev) + offset

        # compute difference between runs,
        # One-norm of a vector (the sum of absolute values)
        vec_diff[i] <- norm(mat_relev - mat_old) / norm(mat_old)
    }

    # matrix to dataframe
    df_out <- data.frame(relevance = mat_relev)
    df_out <- tibble::rownames_to_column(df_out)

    # write out
    vroom::vroom_write(df_out, fpath_out, progress = FALSE)
    NULL
}

#' Compute Gene Relevance
#'
#' Using a randomized walk with restart, compute gene relevance.
#'
#' @examples {
#' ligands <- CytoTalk::ligands_mouse
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' compute_gene_relevance(ligands, dir_in, dir_out)
#' }
#'
#' @param ligands Character matrix of ligand-receptor pair names
#' @param dir_in Input directory, contains scRNAseq files
#' @param dir_out Output directory
#' @return None
#' @export
compute_gene_relevance <- function(ligands, dir_in, dir_out) {
    compute_gene_relevance_type(ligands, dir_in, dir_out, "A")
    compute_gene_relevance_type(ligands, dir_in, dir_out, "B")
    NULL
}


# COMP NODE PRIZE COMBINE


#' @noRd
compute_node_prize_type <- function(type, dir_out, letter) {
    # format filepaths
    fpath_net <- file.path(
        dir_out, sprintf("IntracellularNetwork_Type%s.txt", letter))
    fpath_rel <- file.path(dir_out, sprintf("GeneRelevanceTyp%s.txt", letter))
    fpath_pem <- file.path(dir_out, "GeneCellTypeSpecific.txt")
    fpath_out <- file.path(dir_out, sprintf("GeneNodePrize%s.txt", letter))

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!all(file.exists(c(fpath_net, fpath_rel, fpath_pem)))) {
        stop("cannot find input file(s)")
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load in files
    df_net <- suppressMessages(vroom::vroom(fpath_net, progress = FALSE))
    df_rel <- suppressMessages(vroom::vroom(fpath_rel, progress = FALSE))
    df_pem <- suppressMessages(vroom::vroom(fpath_pem, progress = FALSE))

    # column to rownames
    df_net <- tibble::column_to_rownames(df_net, names(df_net)[1])
    df_rel <- tibble::column_to_rownames(df_rel, names(df_rel)[1])
    df_pem <- tibble::column_to_rownames(df_pem, names(df_pem)[1])

    # find cell type in PEM matrix
    vec_pem <- as.numeric(df_pem[, names(df_pem) == type])

    # match gene relevance names to cell type vector
    index <- match(rownames(df_rel), rownames(df_pem))
    vec_pem_match <- vec_pem[index]

    # relevance times cell specific,
    # if negative then zero out
    df_node_prize <- df_rel * ifelse(vec_pem_match < 0, 0, vec_pem_match)

    # format output dataframe
    df_out <- df_node_prize
    names(df_out) <- "node_prize"
    df_out <- tibble::rownames_to_column(df_out)

    # write out
    vroom::vroom_write(df_out, fpath_out, progress = FALSE)
    NULL
}


#' Compute Node Prize
#'
#' For each cell type, compute overall node prize (integrate gene relevance and
#' PEM scores).
#'
#' @examples {
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_out <- "~/CytoTalk-output"
#'
#' compute_node_prize(type_a, type_b, dir_out)
#' }
#'
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_out Output directory
#' @return None
#' @export
compute_node_prize <- function(type_a, type_b, dir_out) {
    compute_node_prize_type(type_a, dir_out, "A")
    compute_node_prize_type(type_b, dir_out, "B")
    NULL
}


# CROSS-TALK SCORE


#' Compute Cross-Talk Score
#'
#' Compute cross-talk scores between cell type ligand-receptor pairs.
#'
#' @examples {
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_out <- "~/CytoTalk-output"
#'
#' compute_crosstalk(type_a, type_b, dir_out)
#' }
#'
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_out Output directory
#' @return None
#' @export
compute_crosstalk <- function(type_a, type_b, dir_out) {
    # "nst": non-self talk

    # format filepaths
    fpath_pem <- file.path(dir_out, "GeneCellTypeSpecific.txt")
    fpath_nst_a <- file.path(dir_out, "NonSelfTalkSco_TypA.txt")
    fpath_nst_b <- file.path(dir_out, "NonSelfTalkSco_TypB.txt")
    fpath_out <- file.path(dir_out, "Crosstalk_TypATypB.txt")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!all(file.exists(c(fpath_pem, fpath_nst_a, fpath_nst_b)))) {
        stop("cannot find input file(s)")
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load in files
    df_pem <- suppressMessages(vroom::vroom(fpath_pem, progress = FALSE))
    df_nst_a <- suppressMessages(vroom::vroom(fpath_nst_a, progress = FALSE))
    df_nst_b <- suppressMessages(vroom::vroom(fpath_nst_b, progress = FALSE))

    # column to rownames
    df_pem <- tibble::column_to_rownames(df_pem, names(df_pem)[1])
    names_pem <- toupper(rownames(df_pem))

    # grab relevant PEM scores
    vec_pem_a <- as.numeric(df_pem[, names(df_pem) == type_a])
    vec_pem_b <- as.numeric(df_pem[, names(df_pem) == type_b])

    # zero out bad PEM scores
    vec_pem_a <- ifelse(vec_pem_a < 0 | is.na(vec_pem_a), 0, vec_pem_a)
    vec_pem_b <- ifelse(vec_pem_b < 0 | is.na(vec_pem_b), 0, vec_pem_b)

    # find negatives and NAs for non-self-talk
    index_a <- (df_nst_a$mi_dist < 0 | is.na(df_nst_a$mi_dist))
    index_b <- (df_nst_b$mi_dist < 0 | is.na(df_nst_b$mi_dist))

    # zero them out
    df_nst_a$mi_dist[index_a] <- 0
    df_nst_b$mi_dist[index_b] <- 0

    # merge the non-self talk scores together
    df_nst_a$index <- seq_len(nrow(df_nst_a))
    df_nst <- merge(
        df_nst_a, df_nst_b, by = c("lig1" = "lig1", "lig2" = "lig2")
    )
    df_nst <- df_nst[order(df_nst$index), ]
    df_nst$index <- NULL

    # initialize variables
    score_nst <- numeric(0)
    score_express <- numeric(0)
    mat_symbol <- matrix(NA, 0, 2)

    i <- 1
    for (i in seq_len(nrow(df_nst))) {
        lig_x <- df_nst[i, "lig1"]
        lig_y <- df_nst[i, "lig2"]

        # from type A to B (always run)
        pair <- c(sprintf("%s_TypA", lig_x), sprintf("%s_TypB", lig_y))
        mat_symbol <- rbind(mat_symbol, pair)

        # non-self talk score
        tmp <- sum(df_nst[i, c("mi_dist.x", "mi_dist.y")]) / 2
        score_nst <- c(score_nst, tmp)

        # find respective ligands
        index_x <- match(toupper(lig_x), names_pem)
        index_y <- match(toupper(lig_y), names_pem)

        # expressed score
        tmp <- (vec_pem_a[index_x] + vec_pem_b[index_y]) / 2
        score_express <- c(score_express, tmp)

        # from type B to A (sometimes run)
        if (!identical(lig_x, lig_y)) {
            pair <- c(sprintf("%s_TypB", lig_x), sprintf("%s_TypA", lig_y))
            mat_symbol <- rbind(mat_symbol, pair)

            # non-self talk score (same)
            score_nst <- c(score_nst, utils::tail(score_nst, 1))

            # expressed score (notice the difference!)
            tmp <- (vec_pem_b[index_x] + vec_pem_a[index_y]) / 2
            score_express <- c(score_express, tmp)
        }
    }

    # compute crosstalk score
    crosstalk <- minmax(score_express) * minmax(score_nst)

    # format output dataframe
    df_out <- as.data.frame(cbind(mat_symbol, crosstalk))
    names(df_out) <- c("lig1", "lig2", "score")

    # write out
    vroom::vroom_write(df_out, fpath_out, progress = FALSE)
    NULL
}


# NET COST


#' Write Out Integrated Network
#'
#' Uses node prize and edge cost data previously generated to write out a final
#' network configuration that can be processed by the PCST algorithm.
#'
#' @examples {
#' dir_out <- "~/CytoTalk-output"
#' write_integrated_net(dir_out)
#' }
#'
#' @param dir_out Output directory
#' @return None
#' @export
write_integrated_net <- function(dir_out) {
    # format filepaths
    fpath_node_a <- file.path(dir_out, "GeneNodePrizeA.txt")
    fpath_node_b <- file.path(dir_out, "GeneNodePrizeB.txt")
    fpath_edge_a <- file.path(dir_out, "MI_TypA.txt")
    fpath_edge_b <- file.path(dir_out, "MI_TypB.txt")
    fpath_edge_ct <- file.path(dir_out, "Crosstalk_TypATypB.txt")
    fpath_out <- file.path(dir_out, "IntegratedNetwork.cfg")

    # stop if input file doesn't exist,
    # skip if output already generated
    fpaths_in <- c(
        fpath_node_a, fpath_node_b,
        fpath_edge_a, fpath_edge_b, fpath_edge_ct
    )
    if (!all(file.exists(fpaths_in))) {
        stop("cannot find input file(s)")
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # load in files
    quietvroom <- function(fpath) {
        suppressMessages(vroom::vroom(fpath, progress = FALSE))
    }
    df_node_a <- quietvroom(fpath_node_a)
    df_node_b <- quietvroom(fpath_node_b)
    df_edge_a <- quietvroom(fpath_edge_a)
    df_edge_b <- quietvroom(fpath_edge_b)
    df_edge_ct <- quietvroom(fpath_edge_ct)

    # column to rownames
    df_node_a <- tibble::column_to_rownames(df_node_a, names(df_node_a)[1])
    df_node_b <- tibble::column_to_rownames(df_node_b, names(df_node_b)[1])

    # casefold names, determine node and edge types
    add_suffix <- function(v, type) paste0(toupper(v), sprintf("_Typ%s", type))
    node_names_a <- add_suffix(rownames(df_node_a), "A")
    node_names_b <- add_suffix(rownames(df_node_b), "B")
    edge_names_a <- apply(df_edge_a[, c("row", "col")], 2, add_suffix, "A")
    edge_names_b <- apply(df_edge_b[, c("row", "col")], 2, add_suffix, "B")

    # compile full nodes
    node_names_full <- unique(c(node_names_a, node_names_b))

    # validate crosstalk edges
    edge_names_ct <- as.matrix(df_edge_ct[, c(1, 2)])
    index <- rowSums(matrix(edge_names_ct %in% node_names_full, ncol = 2)) == 2
    edge_names_ct <- edge_names_ct[index, ]

    # compile full edges and costs
    edge_full <- rbind(edge_names_a, edge_names_b, edge_names_ct)
    cost <- list(df_edge_a$val, df_edge_b$val, df_edge_ct[index, 3])

    # form node prize section
    np_a <- paste(node_names_a, df_node_a[, 1], sep = "\t")
    np_b <- paste(node_names_b, df_node_b[, 1], sep = "\t")

    # normalize edge costs
    cost_norm <- unlist(lapply(cost, scale))
    cost_full <- 1 - minmax(cost_norm)

    # verify edge names belong to node names
    node_names_full <- unique(c(node_names_a, node_names_b))
    ec_full <- cbind(edge_full, cost_full)
    ec_done <- apply(ec_full, 1, paste0, collapse = "\t")

    # form lines
    lines <- c(
        "# Node names and prizes", np_a, np_b,
        "\n# Edges and costs", ec_done
    )

    # write out
    writeLines(lines, fpath_out)
    NULL
}


# INTEGRATE NETWORK


#' Integrate Cellular Network
#'
#' Given cell types and the names of ligand-receptor pairs, generate node prizes
#' from intracellular edges (gene relevance, mutual information), compute a
#' cross-talk score between cell types, and write out a configuration file to be
#' processed by the PCST algorithm.
#'
#' @examples {
#' ligands <- CytoTalk::ligands_mouse
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' integrate_network(ligands, type_a, type_b, dir_in, dir_out)
#' }
#'
#' @param ligands Character matrix of ligand-receptor pair names
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_in Input directory, contains scRNAseq files
#' @param dir_out Output directory
#' @return None
#' @export
integrate_network <- function(ligands, type_a, type_b, dir_in, dir_out) {
    # gene node prize
    extract_intracell_edges(dir_out)
    compute_gene_relevance(ligands, dir_in, dir_out)
    compute_node_prize(type_a, type_b, dir_out)

    # cross-talk score
    compute_crosstalk(type_a, type_b, dir_out)

    # network analysis
    write_integrated_net(dir_out)
    NULL
}
