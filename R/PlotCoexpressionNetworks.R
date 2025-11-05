#' Plot a single community from one layer (simple)
#'
#' @param ml   mcx_multilayer from build_multilayer(); uses ml$layer_adj
#' @param comm community result (e.g., from detect_genlouvain); uses $state_membership
#' @param community_id integer community label to plot
#' @param layer character; which layer (must be in ml$layer_names)
#' @param seed Random seed for reproducibility.
#'
#' @return invisibly returns the igraph object used for plotting
#' @export
PlotCoexpressionNetworks <- function(ml, comm,
                                     community_id,
                                     layer,
                                     seed = 123) {
  set.seed(seed)
  stopifnot(inherits(ml, "BuildMultilayerNetworkResult"))
  if (!layer %in% ml$layer_names) stop("`layer` not found in ml$layer_names.")
  sm <- comm$state_membership
  req_cols <- c("gene", "layer", "community")
  if (!all(req_cols %in% names(sm))) {
    stop("comm$state_membership must have columns: gene, layer, community")
  }

  # 1) select genes in the requested (community, layer)
  genes_in_comm <- sm$gene[sm$layer == layer & sm$community == community_id]
  genes_in_comm <- unique(genes_in_comm)

  if (length(genes_in_comm) < 2L) {
    message("Nothing to plot: community ", community_id, " in layer '", layer,
            "' has < 2 genes (or none).")
    return(invisible(NULL))
  }

  # 2) slice the layer adjacency to that community’s genes
  A <- ml$layer_adj[[layer]]
  if (!all(genes_in_comm %in% rownames(A))) {
    genes_in_comm <- intersect(genes_in_comm, rownames(A))
  }
  A_sub <- A[genes_in_comm, genes_in_comm, drop = FALSE]

  # remove zero rows/cols (isolates)
  keep <- rowSums(abs(A_sub)) > 0
  A_sub <- A_sub[keep, keep, drop = FALSE]
  if (nrow(A_sub) < 2L) {
    message("Nothing to plot after removing isolates.")
    return(invisible(NULL))
  }

  # 3) build graph and plot (simple igraph)
  suppressPackageStartupMessages(requireNamespace("igraph", quietly = TRUE))
  g <- igraph::graph_from_adjacency_matrix(
    as.matrix(A_sub), mode = "undirected", weighted = TRUE, diag = FALSE
  )

  w <- abs(igraph::E(g)$weight)
  # basic width scaling (robust to range)
  if (length(w)) {
    w_scaled <- 0.5 + 2.5 * (w / max(w, na.rm = TRUE))
  } else {
    w_scaled <- 1
  }

  lay <- igraph::layout_with_fr(g)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mar = c(0,0,2,0))
  plot(
    g,
    layout = lay,
    vertex.label = igraph::V(g)$name,
    vertex.size  = 8,
    vertex.label.cex = 0.7,
    vertex.label.family = "sans",
    vertex.color = "white",
    vertex.frame.color = "grey40",
    edge.width = w_scaled,
    edge.color = "grey40",
    main = paste0("Community ", community_id, " (", layer, ")")
  )

  invisible(g)
}

# [END]
