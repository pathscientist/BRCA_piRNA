# 09_gene_network_analysis.R
# Gene network analysis module for signature piRNAs:
# - Build piRNA-gene and gene-gene network from correlation table
# - Compute centrality metrics
# - Detect communities/modules
# - Export Cytoscape-ready node/edge tables
# - Plot static network overview

source("pipeline/00_config.R")

load_correlated_table <- function(path = file.path(config$output$report_dir, "pearson_correlated_genes.csv")) {
  if (!file.exists(path)) stop("Correlation table not found: ", path)
  tbl <- data.table::fread(path, data.table = FALSE)
  req <- c("piRNA", "gene", "cor")
  miss <- setdiff(req, colnames(tbl))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
  tbl
}

build_network_edges <- function(cor_tbl, abs_cor_cutoff = 0.3, gene_gene_method = c("shared_piRNA", "none")) {
  gene_gene_method <- match.arg(gene_gene_method)
  sig <- cor_tbl[abs(cor_tbl$cor) >= abs_cor_cutoff, , drop = FALSE]

  edges_pg <- sig %>%
    dplyr::transmute(
      from = piRNA,
      to = gene,
      edge_type = "piRNA-gene",
      weight = abs(cor),
      sign = ifelse(cor >= 0, "positive", "negative")
    )

  if (gene_gene_method == "none") return(edges_pg)

  # Gene-gene links if two genes are connected to the same piRNA(s)
  gene_pairs <- sig %>%
    dplyr::select(piRNA, gene) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(., ., by = "piRNA", suffix = c("_a", "_b")) %>%
    dplyr::filter(gene_a < gene_b) %>%
    dplyr::count(gene_a, gene_b, name = "shared_pirna_n") %>%
    dplyr::filter(shared_pirna_n > 0)

  edges_gg <- gene_pairs %>%
    dplyr::transmute(
      from = gene_a,
      to = gene_b,
      edge_type = "gene-gene",
      weight = shared_pirna_n,
      sign = "undirected"
    )

  dplyr::bind_rows(edges_pg, edges_gg)
}

build_network_nodes <- function(edges, cor_tbl) {
  nodes <- tibble::tibble(name = unique(c(edges$from, edges$to))) %>%
    dplyr::mutate(node_type = ifelse(name %in% unique(cor_tbl$piRNA), "piRNA", "gene"))

  # Signed mean correlation for genes (from piRNA-gene edges)
  gene_cor <- cor_tbl %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(mean_cor = mean(cor, na.rm = TRUE), .groups = "drop")

  nodes <- nodes %>%
    dplyr::left_join(gene_cor, by = c("name" = "gene")) %>%
    dplyr::mutate(mean_cor = ifelse(is.na(mean_cor), 0, mean_cor))
  nodes
}

compute_network_stats <- function(nodes, edges) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph is required for network stats")
  g <- igraph::graph_from_data_frame(edges[, c("from", "to", "edge_type", "weight", "sign")], directed = FALSE, vertices = nodes)

  node_stats <- tibble::tibble(
    name = igraph::V(g)$name,
    degree = igraph::degree(g, mode = "all"),
    betweenness = igraph::betweenness(g, normalized = TRUE),
    closeness = igraph::closeness(g, normalized = TRUE),
    eigenvector = igraph::eigen_centrality(g)$vector
  )

  comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  node_stats$module <- igraph::membership(comm)[match(node_stats$name, names(igraph::membership(comm)))]

  graph_stats <- tibble::tibble(
    n_nodes = igraph::vcount(g),
    n_edges = igraph::ecount(g),
    density = igraph::edge_density(g),
    modularity = igraph::modularity(comm),
    n_modules = length(unique(node_stats$module))
  )

  list(graph = g, node_stats = node_stats, graph_stats = graph_stats)
}

plot_network_overview <- function(g, nodes, out_png = file.path(config$output$report_dir, "gene_network_overview.png")) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    warning("ggraph not installed; skipping network plot")
    return(invisible(NULL))
  }

  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(width = weight, color = edge_type), alpha = 0.25) +
    ggraph::geom_node_point(ggplot2::aes(color = node_type, shape = node_type), size = 2.8) +
    ggraph::geom_node_text(ggplot2::aes(label = ifelse(node_type == "piRNA", name, "")), repel = TRUE, size = 3) +
    ggplot2::theme_void() +
    ggplot2::labs(title = "piRNA-gene network overview")

  ggplot2::ggsave(out_png, p, width = 10, height = 8, dpi = 300)
}

export_network_outputs <- function(edges, nodes, stats) {
  data.table::fwrite(edges, file.path(config$output$report_dir, "gene_network_edges.csv"))
  data.table::fwrite(nodes, file.path(config$output$report_dir, "gene_network_nodes.csv"))
  data.table::fwrite(stats$node_stats, file.path(config$output$report_dir, "gene_network_node_centrality.csv"))
  data.table::fwrite(stats$graph_stats, file.path(config$output$report_dir, "gene_network_graph_summary.csv"))
}

run_gene_network_analysis <- function(cor_path = file.path(config$output$report_dir, "pearson_correlated_genes.csv"),
                                      abs_cor_cutoff = 0.3,
                                      gene_gene_method = "shared_piRNA") {
  cor_tbl <- load_correlated_table(cor_path)
  edges <- build_network_edges(cor_tbl, abs_cor_cutoff = abs_cor_cutoff, gene_gene_method = gene_gene_method)
  nodes <- build_network_nodes(edges, cor_tbl)
  stats <- compute_network_stats(nodes, edges)
  export_network_outputs(edges, nodes, stats)
  plot_network_overview(stats$graph, nodes)
  invisible(stats)
}

if (interactive()) {
  run_gene_network_analysis()
}
