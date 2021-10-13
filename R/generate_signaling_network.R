#' @noRd
generate_summary <- function(dir_out) {
    # format filepaths
    fpath_edge <- file.path(dir_out, "PCSF_EdgeOccurance.txt")
    fpath_summ <- file.path(dir_out, "PCSF_EdgeSummary.txt")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_edge)) {
        stop(sprintf("cannot find input file: %s", fpath_edge))
    } else if (file.exists(fpath_summ)) {
        message(sprintf("file already exists, continuing: %s", fpath_summ))
        return()
    }

    # load in file
    df_edge <- suppressMessages(vroom::vroom(fpath_edge, progress = FALSE))

    # remove artificial nodes
    df_edge <- df_edge[df_edge$node1 != "ARTI" & 0.5 <= df_edge$omega, ]

    # count occurances, filter to less than 2000
    tab <- table(df_edge$beta, df_edge$omega)
    tab_which <- which(tab < 2000, arr.ind = TRUE, useNames = FALSE)

    # summary dataframe
    df_summ <- data.frame(
        beta = rownames(tab)[tab_which[, 1]],
        omega = colnames(tab)[tab_which[, 2]],
        n_edges = tab[tab_which]
    )

    # write out
    vroom::vroom_write(df_summ, fpath_summ, progress = FALSE)
    NULL
}

#' @noRd
compute_kolmogorov_smirnov <- function(dir_out, cores=NULL) {
    # format filepaths
    fpath_edge <- file.path(dir_out, "PCSF_EdgeOccurance.txt")
    fpath_pval <- file.path(dir_out, "PCSF_EdgeTestValues.txt")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_edge)) {
        stop(sprintf("cannot find input file: %s", fpath_edge))
    } else if (file.exists(fpath_pval)) {
        message(sprintf("file already exists, continuing: %s", fpath_pval))
        return()
    }

    # load in file
    df_edge <- suppressMessages(vroom::vroom(fpath_edge, progress = FALSE))

    # remove artificial nodes
    df_edge <- df_edge[df_edge$node1 != "ARTI" & 0.5 <= df_edge$omega, ]

    vec_param <- paste(df_edge$beta, df_edge$omega)
    vec_edges <- paste(df_edge$node1, df_edge$node2)
    tab_edges <- table(vec_edges)

    vec_counts <- tab_edges[match(vec_edges, names(tab_edges))]
    lst_counts <- tapply(vec_counts, vec_param, c)

    # detect and register cores
    if (is.null(cores)) {
        message("No number of cores were specified - using all possible.")
        cores <- max(1, parallel::detectCores() - 2)
    } else {
        message(sprintf("Computing K-S tests with %s cores", cores))
    }
    
    doParallel::registerDoParallel(cores = cores)

    # parallel loop for Kolmogorov-Smirnov test
    i <- NULL
    vec_pval <- foreach::`%dopar%`(
        foreach::foreach(i = seq_len(length(lst_counts)), .combine = c), {
        ks <- stats::ks.test(
            lst_counts[[i]], unlist(lst_counts[-i]),
            alternative = "less"
        )
        ks$p.value
    })

    # order by low to hight p-values
    index <- order(vec_pval)
    mat_param <- do.call(rbind, strsplit(names(lst_counts)[index], " "))
    df_pval <- data.frame(mat_param, vec_pval[index])
    names(df_pval) <- c("beta", "omega", "pval")

    # write out
    vroom::vroom_write(df_pval, fpath_pval, progress = FALSE)

    # unregister cores
    doParallel::stopImplicitCluster()
    NULL
}

#' Generate Signaling Network
#'
#' From PCST runs, find the "best" choice for a network based on how unlikely
#' their distribution of edge counts is (higher the better).
#'
#' @examples {
#' dir_out <- "~/CytoTalk-output"
#' generate_signaling_network(dir_out)
#' }
#'
#' @param dir_out Output directory
#' @return None
#' @export
generate_signaling_network <- function(dir_out, cores=NULL) {
    generate_summary(dir_out)
    compute_kolmogorov_smirnov(dir_out, cores)
    NULL
}
