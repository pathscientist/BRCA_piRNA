################################################################################
#                                                                              #
#   Breast Cancer piRNA — Enhanced Network Analysis                            #
#                                                                              #
#   Publication-quality piRNA ↔ Gene ↔ Pathway network visualizations          #
#                                                                              #
#   Produces 8 distinct network/interaction plots:                             #
#                                                                              #
#     N1. Per-piRNA radial networks (one per signature piRNA)                  #
#     N2. Combined bipartite network (piRNAs left ↔ genes right)              #
#     N3. Pathway-centric clustered network                                    #
#         (genes grouped/colored by KEGG pathway, piRNAs connected)            #
#     N4. Hierarchical layered network (piRNA → genes → pathways)             #
#     N5. piRNA-gene correlation heatmap (top genes × signature piRNAs)        #
#     N6. Shared-gene overlap network (piRNAs connected if they share genes)   #
#     N7. Circos-style chord diagram (piRNA ↔ pathway connections)            #
#     N8. Multi-panel figure combining key network views                       #
#                                                                              #
#   Requires objects from piRNA_functional_meta_analysis.R:                    #
#     - sig_cors, top_feats, combat_df_all, gene_cols                          #
#   OR runs its own correlation analysis if sig_cors doesn't exist             #
#                                                                              #
################################################################################

start_time_net <- Sys.time()

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
cran_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "igraph", "ggraph", "ggrepel",
  "RColorBrewer", "scales", "gridExtra", "cowplot",
  "pheatmap", "reshape2", "circlize", "ComplexHeatmap"
)

bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
  library(gridExtra)
  library(cowplot)
  library(pheatmap)
  library(reshape2)
  library(circlize)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

cat("All network analysis packages loaded.\n")

# ==============================================================================
# 0.1 SETTINGS
# ==============================================================================
SEED <- 2024
set.seed(SEED)

# piRNA color (red-family), gene color (blue-family), pathway color (green-family)
PIRNA_COLOR   <- "#D7191C"
GENE_COLOR    <- "#2B83BA"
PATHWAY_COLOR <- "#4DAF4A"

# Distinguishable palette for pathways (up to 12)
PATHWAY_PALETTE <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
)

# piRNA palette (individual piRNAs)
PIRNA_PALETTE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62"
)

pub_theme <- theme_bw() +
  theme(
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 11, color = "black"),
    axis.title        = element_text(size = 13, face = "bold"),
    plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey40"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

net_theme <- theme_void() +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.text   = element_text(size = 9),
    legend.title  = element_text(size = 10, face = "bold"),
    plot.margin   = margin(10, 10, 10, 10)
  )

dir.create("results/network", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0.2 LOAD REQUIRED DATA
# ==============================================================================
cat("\n========== LOADING DATA ==========\n")

# Load model + features
if (!exists("model")) {
  if (file.exists("results/models/final_model.rds"))
    model <- readRDS("results/models/final_model.rds")
}
if (!exists("top_feats")) {
  if (file.exists("results/models/final_features.rds")) {
    top_feats <- readRDS("results/models/final_features.rds")
  } else if (file.exists("results/feature_selection/final_features.txt")) {
    top_feats <- readLines("results/feature_selection/final_features.txt")
    top_feats <- top_feats[nchar(top_feats) > 0]
  }
}

if (!exists("combat_df_all"))
  stop("combat_df_all not found. Run piRNA_multicohort_pipeline.R first.")

gene_cols <- setdiff(colnames(combat_df_all),
                     c("Group", "Batch", "T_Score", "T_Score_binary", "T_Score_all",
                       "Age", "Stage", "Subtype", "Age_binary", "Stage_binary",
                       "Age_Group", "Outcome01", "OS_time", "OS_status"))

cat("Signature piRNAs:", paste(top_feats, collapse = ", "), "\n")

# ==============================================================================
# 0.3 COMPUTE CORRELATIONS (if not already available)
# ==============================================================================
if (!exists("sig_cors") || nrow(sig_cors) == 0) {
  cat("Computing Pearson correlations (piRNAs vs all features)...\n")

  non_sig_genes <- setdiff(gene_cols, top_feats)
  tumor_idx     <- combat_df_all$Group == "Tumor"
  expr_sig      <- as.matrix(combat_df_all[tumor_idx, top_feats])
  expr_others   <- as.matrix(combat_df_all[tumor_idx, non_sig_genes])

  cor_results <- list()
  for (feat in top_feats) {
    cors <- apply(expr_others, 2, function(col) {
      ct <- cor.test(expr_sig[, feat], col, method = "pearson")
      c(r = ct$estimate, p = ct$p.value)
    })
    cor_df <- data.frame(
      piRNA     = feat,
      Gene      = colnames(cors),
      Pearson_r = cors["r.cor", ],
      P_value   = cors["p", ],
      stringsAsFactors = FALSE
    )
    cor_df$P_adj <- p.adjust(cor_df$P_value, method = "BH")
    cor_df$Significant <- abs(cor_df$Pearson_r) > 0.3 & cor_df$P_adj < 0.05
    cor_results[[feat]] <- cor_df
  }

  all_cors <- do.call(rbind, cor_results)
  sig_cors <- all_cors[all_cors$Significant, ]
  cat("  Significant correlations:", nrow(sig_cors), "\n")
}

# --- Attempt KEGG enrichment for pathway mapping ---
kegg_sig <- NULL
correlated_genes <- unique(sig_cors$Gene)
gene_names_clean <- gsub("^hsa[-_]", "", correlated_genes)
gene_names_clean <- gsub("^piR[-_]", "", gene_names_clean)

entrez_map <- tryCatch(
  bitr(gene_names_clean, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db),
  error = function(e) tryCatch(
    bitr(gene_names_clean, fromType = "ALIAS", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e2) NULL
  )
)

if (!is.null(entrez_map) && nrow(entrez_map) >= 5) {
  entrez_ids <- unique(entrez_map$ENTREZID)
  kegg_res_net <- tryCatch(
    enrichKEGG(gene = entrez_ids, organism = "hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.1),
    error = function(e) NULL
  )
  if (!is.null(kegg_res_net)) {
    kegg_sig <- kegg_res_net@result[kegg_res_net@result$p.adjust < 0.05, ]
    if (nrow(kegg_sig) == 0) kegg_sig <- NULL
  }
}

# --- Build master lookup: gene → pathway ---
gene_to_pathway <- data.frame()
if (!is.null(kegg_sig) && nrow(kegg_sig) > 0) {
  for (i in seq_len(min(10, nrow(kegg_sig)))) {
    pw_genes <- unlist(strsplit(kegg_sig$geneID[i], "/"))
    gene_to_pathway <- rbind(gene_to_pathway, data.frame(
      Gene    = pw_genes,
      Pathway = kegg_sig$Description[i],
      stringsAsFactors = FALSE
    ))
  }
}

# --- Select top genes per piRNA for network building ---
TOP_N <- 8  # top correlated genes per piRNA for network plots

top_genes_per_pirna <- sig_cors %>%
  group_by(piRNA) %>%
  arrange(desc(abs(Pearson_r))) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

cat("Top genes for network:", nrow(top_genes_per_pirna),
    "edges across", length(top_feats), "piRNAs\n\n")


# ══════════════════════════════════════════════════════════════════════════════
#  N1. PER-piRNA RADIAL NETWORKS
#      Each piRNA at center, correlated genes radiate outward
#      Edge width = |r|, edge color = direction (red negative, blue positive)
# ══════════════════════════════════════════════════════════════════════════════
cat("========== N1: Per-piRNA Radial Networks ==========\n")

per_pirna_plots <- list()

for (i in seq_along(top_feats)) {
  feat <- top_feats[i]
  feat_edges <- top_genes_per_pirna[top_genes_per_pirna$piRNA == feat, ]
  if (nrow(feat_edges) < 2) {
    cat(sprintf("  %s: too few correlated genes, skipping.\n", feat))
    next
  }

  # Build graph
  edges <- data.frame(
    from      = feat_edges$piRNA,
    to        = feat_edges$Gene,
    r         = feat_edges$Pearson_r,
    abs_r     = abs(feat_edges$Pearson_r),
    direction = ifelse(feat_edges$Pearson_r > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )

  nodes <- data.frame(
    name = c(feat, feat_edges$Gene),
    type = c("piRNA", rep("Gene", nrow(feat_edges))),
    stringsAsFactors = FALSE
  )

  # Assign pathway if available
  nodes$pathway <- "Other"
  if (nrow(gene_to_pathway) > 0) {
    for (j in seq_len(nrow(nodes))) {
      pw <- gene_to_pathway$Pathway[gene_to_pathway$Gene == nodes$name[j]]
      if (length(pw) > 0) nodes$pathway[j] <- pw[1]
    }
  }
  nodes$pathway[nodes$type == "piRNA"] <- "piRNA"

  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

  set.seed(SEED + i)
  p <- ggraph(g, layout = "star") +
    geom_edge_link(
      aes(edge_width = abs_r, edge_color = direction),
      alpha = 0.7, show.legend = TRUE
    ) +
    scale_edge_width_continuous(range = c(0.5, 3), name = "|Pearson r|") +
    scale_edge_color_manual(
      values = c("Positive" = "#2166AC", "Negative" = "#B2182B"),
      name = "Correlation"
    ) +
    geom_node_point(
      aes(shape = type, fill = pathway),
      size = ifelse(V(g)$type == "piRNA", 10, 6),
      color = "black", stroke = 0.5
    ) +
    scale_shape_manual(
      values = c("piRNA" = 23, "Gene" = 21),
      name = "Node Type"
    ) +
    scale_fill_manual(
      values = c("piRNA" = PIRNA_COLOR,
                 setNames(PATHWAY_PALETTE[seq_along(setdiff(unique(nodes$pathway),
                   "piRNA"))], setdiff(unique(nodes$pathway), "piRNA"))),
      name = "Pathway"
    ) +
    geom_node_text(
      aes(label = name), repel = TRUE, size = 3.5,
      fontface = ifelse(V(g)$type == "piRNA", "bold", "plain"),
      max.overlaps = 20
    ) +
    labs(
      title = feat,
      subtitle = paste0("Top ", nrow(feat_edges), " correlated genes (|r|>0.3, adj.p<0.05)")
    ) +
    net_theme +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    )

  fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
  ggsave(paste0("results/network/N1_radial_", fname, ".png"),
         p, width = 10, height = 9, dpi = 300)

  per_pirna_plots[[feat]] <- p
  cat(sprintf("  %s: %d genes plotted.\n", feat, nrow(feat_edges)))
}


# ══════════════════════════════════════════════════════════════════════════════
#  N2. COMBINED BIPARTITE NETWORK
#      piRNAs on left, genes on right
#      Edge color = correlation direction, width = |r|
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N2: Combined Bipartite Network ==========\n")

# Build combined edge list
bipartite_edges <- top_genes_per_pirna %>%
  transmute(
    from      = piRNA,
    to        = Gene,
    r         = Pearson_r,
    abs_r     = abs(Pearson_r),
    direction = ifelse(Pearson_r > 0, "Positive", "Negative")
  )

all_genes_net <- unique(bipartite_edges$to)

bipartite_nodes <- data.frame(
  name = c(top_feats, all_genes_net),
  type = c(rep("piRNA", length(top_feats)), rep("Gene", length(all_genes_net))),
  stringsAsFactors = FALSE
)

# Assign pathway colors to genes
bipartite_nodes$pathway <- "Other"
if (nrow(gene_to_pathway) > 0) {
  for (j in seq_len(nrow(bipartite_nodes))) {
    pw <- gene_to_pathway$Pathway[gene_to_pathway$Gene == bipartite_nodes$name[j]]
    if (length(pw) > 0) bipartite_nodes$pathway[j] <- pw[1]
  }
}
bipartite_nodes$pathway[bipartite_nodes$type == "piRNA"] <- "piRNA"

g_bip <- graph_from_data_frame(bipartite_edges, directed = FALSE,
                               vertices = bipartite_nodes)

# Create bipartite layout: piRNAs left (x=0), genes right (x=1)
layout_bip <- data.frame(
  x = ifelse(V(g_bip)$type == "piRNA", 0, 1),
  y = NA
)
pirna_idx <- which(V(g_bip)$type == "piRNA")
gene_idx  <- which(V(g_bip)$type == "Gene")
layout_bip$y[pirna_idx] <- seq(0, 1, length.out = length(pirna_idx))
layout_bip$y[gene_idx]  <- seq(0, 1, length.out = length(gene_idx))

set.seed(SEED)
p_bip <- ggraph(g_bip, layout = "manual", x = layout_bip$x, y = layout_bip$y) +
  geom_edge_link(
    aes(edge_width = abs_r, edge_color = direction, edge_alpha = abs_r),
    show.legend = TRUE
  ) +
  scale_edge_width_continuous(range = c(0.3, 2.5), name = "|r|") +
  scale_edge_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
  scale_edge_color_manual(
    values = c("Positive" = "#2166AC", "Negative" = "#B2182B"),
    name = "Correlation"
  ) +
  geom_node_point(
    aes(shape = type, fill = pathway),
    size = ifelse(V(g_bip)$type == "piRNA", 8, 5),
    color = "black", stroke = 0.5
  ) +
  scale_shape_manual(values = c("piRNA" = 23, "Gene" = 21), name = "Node Type") +
  scale_fill_manual(
    values = c("piRNA" = PIRNA_COLOR,
               setNames(PATHWAY_PALETTE[seq_along(setdiff(unique(bipartite_nodes$pathway),
                 "piRNA"))], setdiff(unique(bipartite_nodes$pathway), "piRNA"))),
    name = "Pathway"
  ) +
  geom_node_text(
    aes(label = name),
    hjust = ifelse(V(g_bip)$type == "piRNA", 1.2, -0.2),
    size = ifelse(V(g_bip)$type == "piRNA", 3.5, 2.8),
    fontface = ifelse(V(g_bip)$type == "piRNA", "bold", "plain")
  ) +
  labs(
    title = "piRNA-Gene Bipartite Interaction Network",
    subtitle = paste0(length(top_feats), " signature piRNAs and their top correlated genes")
  ) +
  net_theme +
  theme(plot.margin = margin(10, 60, 10, 60)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)))

ggsave("results/network/N2_bipartite_network.png",
       p_bip, width = 14, height = max(8, length(all_genes_net) * 0.25 + 2), dpi = 300)
cat("  Bipartite network saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  N3. PATHWAY-CENTRIC CLUSTERED NETWORK
#      Genes grouped/colored by KEGG pathway
#      piRNAs connected to gene members
#      Layout emphasizes pathway clusters
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N3: Pathway-Centric Clustered Network ==========\n")

if (!is.null(kegg_sig) && nrow(kegg_sig) > 0) {
  # Build 3-layer edges: piRNA → gene + gene → pathway
  pw_edges_from_genes <- data.frame()
  pw_nodes <- data.frame()

  top_pw <- head(kegg_sig, 8)

  for (k in seq_len(nrow(top_pw))) {
    pw_gene_list <- unlist(strsplit(top_pw$geneID[k], "/"))
    matched_genes <- intersect(pw_gene_list, all_genes_net)
    if (length(matched_genes) == 0) next

    pw_edges_from_genes <- rbind(pw_edges_from_genes, data.frame(
      from      = matched_genes,
      to        = top_pw$Description[k],
      abs_r     = 0.5,
      direction = "pathway",
      stringsAsFactors = FALSE
    ))

    pw_nodes <- rbind(pw_nodes, data.frame(
      name    = top_pw$Description[k],
      type    = "Pathway",
      pathway = top_pw$Description[k],
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(pw_edges_from_genes) > 0) {
    # Combine edges
    combined_edges <- rbind(
      bipartite_edges %>% select(from, to, abs_r, direction),
      pw_edges_from_genes
    )

    # Combine nodes
    pw_nodes <- pw_nodes[!duplicated(pw_nodes$name), ]
    combined_nodes <- rbind(bipartite_nodes[, c("name", "type", "pathway")], pw_nodes)
    combined_nodes <- combined_nodes[!duplicated(combined_nodes$name), ]

    g_pw <- graph_from_data_frame(combined_edges, directed = FALSE,
                                  vertices = combined_nodes)

    # Assign node sizes
    node_sizes <- ifelse(V(g_pw)$type == "piRNA", 9,
                  ifelse(V(g_pw)$type == "Pathway", 7, 4.5))

    # Build unique pathway colors
    all_pathways <- unique(combined_nodes$pathway)
    pw_color_map <- c("piRNA" = PIRNA_COLOR, "Other" = GENE_COLOR)
    other_pws <- setdiff(all_pathways, c("piRNA", "Other"))
    if (length(other_pws) > 0) {
      pw_color_map <- c(pw_color_map,
                        setNames(PATHWAY_PALETTE[seq_along(other_pws)], other_pws))
    }

    set.seed(SEED)
    p_pw <- ggraph(g_pw, layout = "fr") +
      geom_edge_link(
        aes(edge_alpha = abs_r,
            edge_linetype = ifelse(direction == "pathway", "dotted",
                            ifelse(direction == "Negative", "dashed", "solid"))),
        edge_colour = "grey55", edge_width = 0.4, show.legend = FALSE
      ) +
      geom_node_point(
        aes(shape = type, fill = pathway, size = node_sizes),
        color = "black", stroke = 0.4
      ) +
      scale_shape_manual(
        values = c("piRNA" = 23, "Gene" = 21, "Pathway" = 22),
        name = "Node Type"
      ) +
      scale_fill_manual(values = pw_color_map, name = "Group") +
      scale_size_identity() +
      geom_node_text(
        aes(label = ifelse(type != "Gene", name, ""),
            fontface = ifelse(type == "piRNA", "bold", "plain")),
        repel = TRUE, size = 3.2, max.overlaps = 25
      ) +
      # Smaller labels for genes
      geom_node_text(
        aes(label = ifelse(type == "Gene", name, "")),
        repel = TRUE, size = 2.2, color = "grey30", max.overlaps = 30
      ) +
      labs(
        title = "Pathway-Centric Interaction Network",
        subtitle = paste0(length(top_feats), " piRNAs \u2192 correlated genes \u2192 KEGG pathways")
      ) +
      net_theme +
      guides(
        fill  = guide_legend(override.aes = list(shape = 21, size = 4)),
        shape = guide_legend(override.aes = list(size = 4, fill = "grey70"))
      )

    ggsave("results/network/N3_pathway_clustered_network.png",
           p_pw, width = 16, height = 12, dpi = 300)
    cat("  Pathway-centric network saved.\n")
  } else {
    cat("  No gene-pathway overlaps found for network N3.\n")
  }
} else {
  cat("  KEGG results unavailable; skipping pathway-centric network.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
#  N4. HIERARCHICAL LAYERED NETWORK (piRNA → Genes → Pathways)
#      Sugiyama-style top-to-bottom layout with 3 layers
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N4: Hierarchical Layered Network ==========\n")

if (!is.null(kegg_sig) && nrow(kegg_sig) > 0 && nrow(pw_edges_from_genes) > 0) {
  # Assign layer: piRNA=1, Gene=2, Pathway=3
  combined_nodes$layer <- ifelse(combined_nodes$type == "piRNA", 1,
                          ifelse(combined_nodes$type == "Gene", 2, 3))

  g_hier <- graph_from_data_frame(combined_edges, directed = TRUE,
                                  vertices = combined_nodes)

  # Manual 3-layer layout
  n_pirna   <- sum(combined_nodes$type == "piRNA")
  n_gene    <- sum(combined_nodes$type == "Gene")
  n_pathway <- sum(combined_nodes$type == "Pathway")

  layout_hier <- data.frame(x = NA, y = NA)
  row.names(layout_hier) <- combined_nodes$name

  pirna_names   <- combined_nodes$name[combined_nodes$type == "piRNA"]
  gene_names    <- combined_nodes$name[combined_nodes$type == "Gene"]
  pathway_names <- combined_nodes$name[combined_nodes$type == "Pathway"]

  layout_hier[pirna_names, "y"]   <- 3
  layout_hier[gene_names, "y"]    <- 2
  layout_hier[pathway_names, "y"] <- 1

  layout_hier[pirna_names, "x"]   <- seq(0, 1, length.out = max(n_pirna, 1))
  layout_hier[gene_names, "x"]    <- seq(0, 1, length.out = max(n_gene, 1))
  layout_hier[pathway_names, "x"] <- seq(0, 1, length.out = max(n_pathway, 1))

  # Match vertex order
  layout_mat <- as.matrix(layout_hier[V(g_hier)$name, c("x", "y")])

  set.seed(SEED)
  p_hier <- ggraph(g_hier, layout = "manual", x = layout_mat[, 1], y = layout_mat[, 2]) +
    # Layer labels
    annotate("text", x = -0.08, y = 3, label = "piRNAs", fontface = "bold",
             size = 4.5, hjust = 1, color = PIRNA_COLOR) +
    annotate("text", x = -0.08, y = 2, label = "Genes", fontface = "bold",
             size = 4.5, hjust = 1, color = GENE_COLOR) +
    annotate("text", x = -0.08, y = 1, label = "Pathways", fontface = "bold",
             size = 4.5, hjust = 1, color = PATHWAY_COLOR) +
    geom_edge_link(
      aes(edge_alpha = abs_r),
      edge_colour = "grey65", edge_width = 0.35, show.legend = FALSE
    ) +
    geom_node_point(
      aes(shape = type, fill = pathway),
      size = ifelse(V(g_hier)$type == "piRNA", 8,
              ifelse(V(g_hier)$type == "Pathway", 6, 4)),
      color = "black", stroke = 0.4
    ) +
    scale_shape_manual(values = c("piRNA" = 23, "Gene" = 21, "Pathway" = 22),
                       name = "Node Type") +
    scale_fill_manual(values = pw_color_map, name = "Group") +
    geom_node_text(
      aes(label = name),
      size = ifelse(V(g_hier)$type == "Gene", 2, 3),
      repel = TRUE, max.overlaps = 30,
      fontface = ifelse(V(g_hier)$type == "piRNA", "bold", "plain")
    ) +
    labs(
      title = "Hierarchical Network: piRNAs \u2192 Genes \u2192 Pathways",
      subtitle = "Three-layer functional architecture"
    ) +
    net_theme +
    coord_cartesian(xlim = c(-0.15, 1.1)) +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, size = 4)),
      shape = guide_legend(override.aes = list(size = 4, fill = "grey70"))
    )

  ggsave("results/network/N4_hierarchical_network.png",
         p_hier, width = 16, height = 10, dpi = 300)
  cat("  Hierarchical layered network saved.\n")
} else {
  cat("  Skipping hierarchical network (no pathway data).\n")
}


# ══════════════════════════════════════════════════════════════════════════════
#  N5. piRNA-GENE CORRELATION HEATMAP
#      Rows = top correlated genes, Cols = signature piRNAs
#      Values = Pearson r, clustered both ways
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N5: piRNA-Gene Correlation Heatmap ==========\n")

# Build matrix: genes (rows) × piRNAs (cols) with Pearson r
all_net_genes <- unique(top_genes_per_pirna$Gene)

cor_heatmap_mat <- matrix(0, nrow = length(all_net_genes), ncol = length(top_feats),
                          dimnames = list(all_net_genes, top_feats))

# Fill from all_cors (full correlation table with all piRNA-gene pairs)
if (exists("all_cors")) {
  for (feat in top_feats) {
    feat_data <- all_cors[all_cors$piRNA == feat & all_cors$Gene %in% all_net_genes, ]
    cor_heatmap_mat[feat_data$Gene, feat] <- feat_data$Pearson_r
  }
} else {
  # Fallback: use sig_cors
  for (feat in top_feats) {
    feat_data <- sig_cors[sig_cors$piRNA == feat & sig_cors$Gene %in% all_net_genes, ]
    cor_heatmap_mat[feat_data$Gene, feat] <- feat_data$Pearson_r
  }
}

# Annotation: which piRNA selected each gene
gene_source <- data.frame(row.names = all_net_genes)
gene_source$`Primary piRNA` <- sapply(all_net_genes, function(g) {
  sub <- top_genes_per_pirna[top_genes_per_pirna$Gene == g, ]
  sub$piRNA[which.max(abs(sub$Pearson_r))]
})

# Pathway annotation
gene_source$Pathway <- "Other"
if (nrow(gene_to_pathway) > 0) {
  for (g in all_net_genes) {
    pw <- gene_to_pathway$Pathway[gene_to_pathway$Gene == g]
    if (length(pw) > 0) gene_source[g, "Pathway"] <- pw[1]
  }
}

# Color scales
hm_colors <- colorRampPalette(c("#B2182B", "#F4A582", "#FDDBC7", "#F7F7F7",
                                 "#D1E5F0", "#92C5DE", "#2166AC"))(100)

# piRNA annotation colors
pirna_ann_colors <- setNames(PIRNA_PALETTE[seq_along(top_feats)], top_feats)
pathway_levels <- unique(gene_source$Pathway)
pw_ann_colors <- setNames(c(GENE_COLOR, PATHWAY_PALETTE)[seq_along(pathway_levels)],
                          pathway_levels)

ann_colors <- list(
  `Primary piRNA` = pirna_ann_colors,
  Pathway = pw_ann_colors
)

png("results/network/N5_pirna_gene_correlation_heatmap.png",
    width = max(8, length(top_feats) * 1.5 + 3),
    height = max(8, length(all_net_genes) * 0.3 + 3),
    units = "in", res = 300)

pheatmap(
  cor_heatmap_mat,
  color = hm_colors,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = (length(top_feats) > 2),
  annotation_row = gene_source,
  annotation_colors = ann_colors,
  main = "Pearson Correlation: Signature piRNAs vs Correlated Genes",
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 10,
  border_color = "grey85",
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  legend_breaks = c(-1, -0.5, 0, 0.5, 1)
)

dev.off()
cat("  piRNA-gene correlation heatmap saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  N6. SHARED-GENE OVERLAP NETWORK
#      piRNAs are nodes, connected if they share correlated genes
#      Edge width = number of shared genes, edge label = count
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N6: piRNA Shared-Gene Overlap Network ==========\n")

if (length(top_feats) >= 2) {
  # Compute gene overlap between each piRNA pair
  pirna_gene_sets <- lapply(top_feats, function(feat) {
    unique(sig_cors$Gene[sig_cors$piRNA == feat])
  })
  names(pirna_gene_sets) <- top_feats

  overlap_edges <- data.frame()
  for (i in seq_along(top_feats)) {
    for (j in seq_along(top_feats)) {
      if (j <= i) next
      shared <- intersect(pirna_gene_sets[[i]], pirna_gene_sets[[j]])
      if (length(shared) > 0) {
        overlap_edges <- rbind(overlap_edges, data.frame(
          from         = top_feats[i],
          to           = top_feats[j],
          shared_count = length(shared),
          shared_genes = paste(head(shared, 5), collapse = ", "),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  overlap_nodes <- data.frame(
    name   = top_feats,
    n_genes = sapply(pirna_gene_sets, length),
    stringsAsFactors = FALSE
  )

  if (nrow(overlap_edges) > 0) {
    g_overlap <- graph_from_data_frame(overlap_edges, directed = FALSE,
                                       vertices = overlap_nodes)

    set.seed(SEED)
    p_overlap <- ggraph(g_overlap, layout = "circle") +
      geom_edge_link(
        aes(edge_width = shared_count, label = shared_count),
        edge_colour = "#4393C3", alpha = 0.6,
        label_size = 3.5, label_colour = "#333333",
        angle_calc = "along", label_dodge = unit(3, "mm")
      ) +
      scale_edge_width_continuous(range = c(1, 6), name = "Shared Genes") +
      geom_node_point(
        aes(size = n_genes), shape = 23,
        fill = PIRNA_COLOR, color = "black", stroke = 0.8
      ) +
      scale_size_continuous(range = c(8, 16), name = "Total Correlated\nGenes") +
      geom_node_text(aes(label = name), size = 3.5, fontface = "bold",
                     repel = TRUE) +
      labs(
        title = "piRNA Shared-Gene Overlap Network",
        subtitle = "Edges = number of shared correlated genes between piRNA pairs"
      ) +
      net_theme

    ggsave("results/network/N6_shared_gene_overlap.png",
           p_overlap, width = 10, height = 9, dpi = 300)
    cat("  Shared-gene overlap network saved.\n")

    # Save overlap table
    write.csv(overlap_edges, "results/network/shared_gene_overlap_table.csv",
              row.names = FALSE)
  } else {
    cat("  No shared genes between piRNAs.\n")
  }
} else {
  cat("  Only 1 piRNA; skipping overlap network.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
#  N7. CIRCOS-STYLE CHORD DIAGRAM
#      Sectors = piRNAs + Pathways, chords = connections through shared genes
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N7: Circos Chord Diagram ==========\n")

if (!is.null(kegg_sig) && nrow(kegg_sig) > 0 && nrow(gene_to_pathway) > 0) {
  # Build piRNA → pathway adjacency matrix (count of genes connecting them)
  top_pw_names <- head(kegg_sig$Description, 8)
  chord_mat <- matrix(0, nrow = length(top_feats), ncol = length(top_pw_names),
                      dimnames = list(top_feats, top_pw_names))

  for (feat in top_feats) {
    feat_genes <- top_genes_per_pirna$Gene[top_genes_per_pirna$piRNA == feat]
    for (pw in top_pw_names) {
      pw_genes <- unlist(strsplit(
        kegg_sig$geneID[kegg_sig$Description == pw], "/"))
      chord_mat[feat, pw] <- length(intersect(feat_genes, pw_genes))
    }
  }

  # Remove empty rows/cols
  chord_mat <- chord_mat[rowSums(chord_mat) > 0, colSums(chord_mat) > 0, drop = FALSE]

  if (nrow(chord_mat) > 0 && ncol(chord_mat) > 0) {
    # Color for piRNAs and pathways
    pirna_cols <- setNames(PIRNA_PALETTE[seq_len(nrow(chord_mat))], rownames(chord_mat))
    pw_cols    <- setNames(PATHWAY_PALETTE[seq_len(ncol(chord_mat))], colnames(chord_mat))
    grid_cols  <- c(pirna_cols, pw_cols)

    png("results/network/N7_circos_chord_diagram.png",
        width = 10, height = 10, units = "in", res = 300)

    circos.clear()
    circos.par(start.degree = 90, gap.after = c(rep(2, nrow(chord_mat) - 1), 10,
                                                 rep(2, ncol(chord_mat) - 1), 10))

    chordDiagram(
      chord_mat,
      grid.col = grid_cols,
      annotationTrack = c("grid", "name"),
      annotationTrackHeight = mm_h(c(3, 4)),
      transparency = 0.35,
      big.gap = 15,
      small.gap = 2,
      preAllocateTracks = list(
        track.height = mm_h(3),
        track.margin = c(mm_h(1), 0)
      )
    )

    title("piRNA-Pathway Chord Diagram\n(connections through shared correlated genes)",
          cex.main = 1.2, font.main = 2, line = -2)

    circos.clear()
    dev.off()

    cat("  Circos chord diagram saved.\n")
  } else {
    cat("  Chord matrix is empty; skipping circos plot.\n")
  }
} else {
  cat("  No pathway data available; skipping circos plot.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
#  N8. MULTI-PANEL COMPOSITE FIGURE
#      Combine key network views into one publication figure
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========== N8: Multi-Panel Composite Figure ==========\n")

# Combine N2 (bipartite) + N5 (heatmap as ggplot) + N6 (overlap) if available
panel_plots <- list()

# Panel A: bipartite
if (exists("p_bip")) panel_plots[["A"]] <- p_bip + labs(title = NULL, subtitle = NULL)

# Panel B: overlap network
if (exists("p_overlap")) panel_plots[["B"]] <- p_overlap + labs(title = NULL, subtitle = NULL)

if (length(panel_plots) >= 2) {
  combined <- plot_grid(
    panel_plots[["A"]] + theme(legend.position = "bottom",
                               legend.text = element_text(size = 7)),
    panel_plots[["B"]] + theme(legend.position = "bottom",
                               legend.text = element_text(size = 7)),
    ncol = 2, labels = c("A", "B"), label_size = 16,
    rel_widths = c(1.2, 1)
  )

  title_grob <- ggdraw() +
    draw_label("piRNA-Gene Interaction Networks",
               fontface = "bold", size = 16, hjust = 0.5)

  final_panel <- plot_grid(title_grob, combined, ncol = 1, rel_heights = c(0.04, 1))

  ggsave("results/network/N8_composite_networks.png",
         final_panel, width = 22, height = 12, dpi = 300)
  cat("  Composite multi-panel figure saved.\n")
} else {
  cat("  Not enough panels for composite figure.\n")
}


# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  NETWORK ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Signature:", paste(top_feats, collapse = ", "), "\n")
cat("Correlated genes used:", length(all_net_genes), "\n\n")

output_files <- list.files("results/network", pattern = "\\.png$", full.names = TRUE)
cat("Generated plots (", length(output_files), "):\n")
for (f in output_files) cat("  ", basename(f), "\n")

cat("\nCSV tables:\n")
csv_files <- list.files("results/network", pattern = "\\.csv$", full.names = TRUE)
for (f in csv_files) cat("  ", basename(f), "\n")

end_time_net <- Sys.time()
cat("\nRuntime:", round(difftime(end_time_net, start_time_net, units = "mins"), 1), "min\n")
cat("\n*** NETWORK ANALYSIS COMPLETE ***\n")
