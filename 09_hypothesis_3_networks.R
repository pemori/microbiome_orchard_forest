# =============================================================================
# Script: 09_hypothesis_3_networks.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates hypothesis 3 by comparing co-occurrence patterns of
# bacterial-fungal communities across soils from forests and orchards with
# different surrounding forest cover. The workflow:
#   (i) prepares bacterial and fungal genus-level phyloseq objects,
#   (ii) adds niche-breadth categories to taxa,
#   (iii) merges bacterial and fungal datasets,
#   (iv) defines forest-cover classes,
#   (v) constructs co-occurrence networks with SPIEC-EASI / NetCoMi for
#       LFC, MFC, HFC and Forest samples,
#   (vi) visualizes each network separately,
#   (vii) quantifies positive/negative interactions among domains (B-B, B-F, F-F),
#   (viii) quantifies positive/negative interactions among niche-breadth
#         categories (G-G, O-O, S-S, G-O, G-S, O-S),
#   (ix) exports each panel separately as TIFF files.
#
# Required objects in the workspace:
#   - data_phylo_g_bac
#   - data_phylo_g_fun
#   - bacti_pal_graph_long
#   - fungi_pal_graph_long
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(NetCoMi)
  library(SpiecEasi)
  library(igraph)
  library(tidyverse)
  library(tidygraph)
  library(ggraph)
  library(graphlayouts)
  library(scales)
  library(ggClusterNet)
})

# =============================================================================
# 0. Output directories
# =============================================================================

outdir_h3 <- file.path(getwd(), "results_hyp3_networks")
figdir_h3 <- file.path(outdir_h3, "figures")
tabdir_h3 <- file.path(outdir_h3, "tables")
rdsdir_h3 <- file.path(outdir_h3, "rds")

dir.create(outdir_h3, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir_h3, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir_h3, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir_h3, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Prepare bacterial and fungal phyloseq objects
# =============================================================================

ps16 <- data_phylo_g_bac
psIT <- data_phylo_g_fun

# -----------------------------------------------------------------------------
# 1.1 Clean bacterial taxonomy
# -----------------------------------------------------------------------------

tax_bac <- tax_table(ps16)
colnames(tax_bac)[1] <- "Kingdom"

prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__")

tax16_clean2 <- apply(tax_bac, 2, function(col) {
  col <- gsub(paste0("^(", paste(prefixes, collapse = "|"), ")"), "", col)
  col[col %in% c("NA", "", " ")] <- NA
  return(col)
})

tax16_clean_df <- as.data.frame(tax_table(as.matrix(tax16_clean2)))
tax16_clean_df$Genus <- make.unique(as.character(tax16_clean_df$Genus), sep = "_")
rownames(tax16_clean_df) <- taxa_names(ps16)

tax_table(ps16) <- tax_table(as.matrix(tax16_clean_df))

# -----------------------------------------------------------------------------
# 1.2 Clean fungal taxonomy
# -----------------------------------------------------------------------------

tax_fun <- tax_table(psIT)

prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__")

taxIT_clean2 <- apply(tax_fun, 2, function(col) {
  col <- gsub(paste0("^(", paste(prefixes, collapse = "|"), ")"), "", col)
  col[col %in% c("NA", "", " ")] <- NA
  return(col)
})

taxIT_clean_df <- as.data.frame(tax_table(as.matrix(taxIT_clean2)))
taxIT_clean_df$Genus <- make.unique(as.character(taxIT_clean_df$Genus), sep = "_")
rownames(taxIT_clean_df) <- taxa_names(psIT)

tax_table(psIT) <- tax_table(as.matrix(taxIT_clean_df))

# =============================================================================
# 2. Add niche-breadth categories
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Bacteria
# -----------------------------------------------------------------------------

nb_cat_df <- bacti_pal_graph_long %>%
  select(Genus, NB_Category) %>%
  distinct()

tax16_clean_df <- as.data.frame(tax_table(ps16), stringsAsFactors = FALSE)
tax16_clean_df <- left_join(tax16_clean_df, nb_cat_df, by = "Genus")
rownames(tax16_clean_df) <- taxa_names(ps16)

tax_table(ps16) <- tax_table(as.matrix(tax16_clean_df))

# -----------------------------------------------------------------------------
# 2.2 Fungi
# -----------------------------------------------------------------------------

nb_cat_df <- fungi_pal_graph_long %>%
  select(Genus, NB_Category) %>%
  distinct()

taxIT_clean_df <- as.data.frame(tax_table(psIT), stringsAsFactors = FALSE)
taxIT_clean_df <- left_join(taxIT_clean_df, nb_cat_df, by = "Genus")
rownames(taxIT_clean_df) <- taxa_names(psIT)

tax_table(psIT) <- tax_table(as.matrix(taxIT_clean_df))

# =============================================================================
# 3. Harmonize samples and merge bacterial + fungal phyloseq objects
# =============================================================================

common_samples <- intersect(sample_names(ps16), sample_names(psIT))
ps16 <- prune_samples(common_samples, ps16)
psIT <- prune_samples(common_samples, psIT)

ps16 <- prune_taxa(taxa_sums(ps16) > 0, ps16)
psIT <- prune_taxa(taxa_sums(psIT) > 0, psIT)
ps16
psIT

ps <- ggClusterNet::merge16S_ITS(
  ps16s = ps16,
  psITS = psIT,
  N16s = 507,
  NITS = 589,
  scale = TRUE,
  onlygroup = FALSE,
  dat1.lab = "bac",
  dat2.lab = "fun"
)

# =============================================================================
# 4. Define forest-cover classes
# =============================================================================

sample_df <- data.frame(sample_data(ps), check.names = FALSE)
str(sample_df)
sample_df$FC <- ifelse(
  sample_df$Group == "Forest", "F",
  ifelse(
    sample_df$P_2500 <= 0.3584376, "LFC",
    ifelse(sample_df$P_2500 <= 0.5224324, "MFC", "HFC")
  )
)

sample_df$FC <- factor(sample_df$FC, levels = c("F", "LFC", "MFC", "HFC"))
sample_df
# preserve sample names as rownames
sample_df <- as.data.frame(sample_df, stringsAsFactors = FALSE)
sample_names_now <- rownames(sample_df)

# assign back to phyloseq as proper sample_data object
sample_data(ps) <- sample_data(sample_df)

# write plain data.frame to disk
sample_df_export <- as.data.frame(as.matrix(sample_data(ps)), stringsAsFactors = FALSE)
rownames(sample_df_export) <- sample_names_now

write.csv(
  sample_df_export,
  file.path(tabdir_h3, "sample_metadata_with_FC.csv"),
  row.names = TRUE
)

table(sample_df_export$FC)

# =============================================================================
# 5. Global plotting settings
# =============================================================================

fVar <- 1096

nb_colors <- c(
  "Generalist"  = "#1f78b4",
  "Opportunist" = "gold",
  "Specialist"  = "#33a02c",
  "Unassigned"  = "gray70"
)

assoc_colors <- c(
  "Positive" = "#56B4E9",
  "Negative" = "#D55E00"
)

# =============================================================================
# 6. LFC network
# =============================================================================

ps_LFC <- subset_samples(ps, FC == "LFC")
ps_LFC <- prune_taxa(taxa_sums(ps_LFC) > 0, ps_LFC)

# -----------------------------------------------------------------------------
# 6.1 OTU matrix preparation
# NOTE: preserved exactly as in your working workflow
# -----------------------------------------------------------------------------

otu_LFC <- as.data.frame(otu_table(ps_LFC))

dim(otu_LFC)
otu_LFC2 <- otu_LFC + 1e-6
colSums(otu_LFC2)
otu_LFC3 <- sweep(otu_LFC2, 2, colSums(otu_LFC2), FUN = "/")
colSums(otu_LFC3)
otu_LFC4 <- t(otu_LFC3)
dim(otu_LFC4)
rowSums(otu_LFC4)

# -----------------------------------------------------------------------------
# 6.2 NetCoMi / SPIEC-EASI
# -----------------------------------------------------------------------------

net_LFC <- netConstruct(
  otu_LFC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  measurePar  = list(
    method           = "mb",
    sel.criterion    = "stars",
    pulsar.params    = list(
      rep.num         = 30,
      thresh          = 0.05,
      subsample.ratio = 0.8,
      seed            = 2015
    ),
    nlambda          = 20,
    lambda.min.ratio = 0.05
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_LFC_anal <- netAnalyze(net_LFC)
summary(net_LFC_anal)

saveRDS(net_LFC,      file.path(rdsdir_h3, "net_LFC.rds"))
saveRDS(net_LFC_anal, file.path(rdsdir_h3, "net_LFC_anal.rds"))

edgelist_LFC <- net_LFC$edgelist1[
  order(net_LFC$edgelist1$adja, decreasing = TRUE), ]
write.csv(edgelist_LFC, file.path(tabdir_h3, "edgelist_LFC.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6.3 Build igraph object
# -----------------------------------------------------------------------------

edg_LFC <- net_LFC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g_LFC <- graph_from_data_frame(edg_LFC, directed = FALSE)

V(g_LFC)$type <- ifelse(grepl("^bac_", V(g_LFC)$name), "Bacteria",
                        ifelse(grepl("^fun_", V(g_LFC)$name), "Fungi", "Unknown"))
V(g_LFC)$deg  <- igraph::degree(g_LFC)
V(g_LFC)$eig  <- igraph::eigen_centrality(g_LFC)$vector
V(g_LFC)$comm <- cluster_fast_greedy(g_LFC)$membership

E(g_LFC)$weight <- edg_LFC$w_abs
E(g_LFC)$sign   <- edg_LFC$sign
E(g_LFC)$w_abs  <- edg_LFC$w_abs

# -----------------------------------------------------------------------------
# 6.4 Add NB categories
# -----------------------------------------------------------------------------

tax_df_LFC <- as.data.frame(tax_table(ps_LFC))
tax_df_LFC$node <- rownames(tax_df_LFC)

tax_g_df_LFC <- tax_df_LFC %>%
  filter(node %in% V(g_LFC)$name) %>%
  select(node, NB_Category)

tax_g_df_LFC$NB_Category[is.na(tax_g_df_LFC$NB_Category)] <- "Unassigned"
V(g_LFC)$NB_Category <- tax_g_df_LFC$NB_Category[match(V(g_LFC)$name, tax_g_df_LFC$node)]

table(V(g_LFC)$NB_Category)

# -----------------------------------------------------------------------------
# 6.5 Top positive / negative nodes
# -----------------------------------------------------------------------------

ed_df_LFC <- igraph::as_data_frame(g_LFC, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign_LFC <- ed_df_LFC %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0)

if (!"positive" %in% colnames(count_by_sign_LFC)) count_by_sign_LFC$positive <- 0
if (!"negative" %in% colnames(count_by_sign_LFC)) count_by_sign_LFC$negative <- 0

count_by_sign_LFC <- count_by_sign_LFC %>%
  rename(neg = negative, pos = positive)

top_pos_nodes_LFC <- count_by_sign_LFC %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes_LFC <- count_by_sign_LFC %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)

V(g_LFC)$label <- ifelse(V(g_LFC)$name %in% union(top_pos_nodes_LFC, top_neg_nodes_LFC), V(g_LFC)$name, NA)
V(g_LFC)$label_sign <- case_when(
  V(g_LFC)$name %in% top_pos_nodes_LFC & V(g_LFC)$name %in% top_neg_nodes_LFC ~ "both",
  V(g_LFC)$name %in% top_pos_nodes_LFC ~ "positive",
  V(g_LFC)$name %in% top_neg_nodes_LFC ~ "negative",
  TRUE ~ NA_character_
)

# -----------------------------------------------------------------------------
# 6.6 Filtered graph for plotting
# -----------------------------------------------------------------------------

g_LFC_filt <- igraph::delete_edges(
  g_LFC,
  E(g_LFC)[abs(weight) < 0.3]
)

E(g_LFC_filt)$sign <- as.character(E(g_LFC_filt)$sign)
E(g_LFC_filt)$sign[E(g_LFC_filt)$sign %in% c("Positive", "POSITIVE")] <- "positive"
E(g_LFC_filt)$sign[E(g_LFC_filt)$sign %in% c("Negative", "NEGATIVE")] <- "negative"
E(g_LFC_filt)$sign <- factor(E(g_LFC_filt)$sign, levels = c("positive", "negative"))

V(g_LFC_filt)$degree <- igraph::degree(g_LFC_filt)
V(g_LFC_filt)$deg <- igraph::degree(g_LFC_filt, mode = "all", loops = FALSE)
node_size_LFC <- scales::rescale(V(g_LFC_filt)$deg, to = c(2, 10))

V(g_LFC_filt)$NB_Category <- factor(V(g_LFC_filt)$NB_Category, levels = names(nb_colors))

set.seed(42)
lay_LFC <- create_layout(g_LFC_filt, layout = "circle")

lay_df <- as.data.frame(lay_LFC) %>% 
  filter(!is.na(label))

# -----------------------------------------------------------------------------
# 6.7 Network plot
# -----------------------------------------------------------------------------

p_net_LFC <- ggraph(lay_LFC) +
  geom_edge_link(aes(edge_colour = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_colour_manual(values = c(positive = "gray70",
                                      negative = "gray40"),
                           drop = FALSE,
                           name = "Association") +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.5, 0.8)) +
  guides(edge_width = "none", edge_alpha = "none") +
  geom_node_point(aes(fill = NB_Category),
                  shape = 21, size = node_size_LFC, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = nb_colors, name = "NB Category") +
  theme_void() +
  coord_equal()

p_net_LFC
# p_net_LFC # network_LFC 900*900

ggsave(file.path(figdir_h3, "network_LFC.tiff"),
       plot = p_net_LFC, width = 9, height = 9, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.8 Domain-level interactions
# -----------------------------------------------------------------------------

edge_df_LFC <- igraph::as_data_frame(g_LFC, what = "edges") %>%
  mutate(
    from_type = V(g_LFC)$type[match(from, V(g_LFC)$name)],
    to_type   = V(g_LFC)$type[match(to,   V(g_LFC)$name)],
    interaction = case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_df_LFC <- edge_df_LFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_LFC <- edge_df_LFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_LFC <- edge_counts_LFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_LFC <- left_join(edge_counts_LFC, total_counts_LFC, by = "interaction")

write.csv(edge_df_LFC,     file.path(tabdir_h3, "domain_edges_LFC.csv"), row.names = FALSE)
write.csv(prop_df_LFC,     file.path(tabdir_h3, "domain_edge_proportions_LFC.csv"), row.names = FALSE)
write.csv(edge_counts_LFC, file.path(tabdir_h3, "domain_edge_counts_LFC.csv"), row.names = FALSE)

negative_ff_LFC <- edge_df_LFC %>% filter(interaction == "F-F", association == "Negative")
negative_bf_LFC <- edge_df_LFC %>% filter(interaction == "B-F", association == "Negative")
positive_bf_LFC <- edge_df_LFC %>% filter(interaction == "B-F", association == "Positive")
negative_bb_LFC <- edge_df_LFC %>% filter(interaction == "B-B", association == "Negative")

write.csv(negative_ff_LFC, file.path(tabdir_h3, "negative_ff_LFC.csv"), row.names = FALSE)
write.csv(negative_bf_LFC, file.path(tabdir_h3, "negative_bf_LFC.csv"), row.names = FALSE)
write.csv(positive_bf_LFC, file.path(tabdir_h3, "positive_bf_LFC.csv"), row.names = FALSE)
write.csv(negative_bb_LFC, file.path(tabdir_h3, "negative_bb_LFC.csv"), row.names = FALSE)

y_max_LFC <- max(edge_counts_LFC$total, na.rm = TRUE)
neg_offset_LFC <- 0.05 * y_max_LFC

p_dom_LFC_lab <- ggplot(edge_counts_LFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_LFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_LFC, association == "Negative"),
    aes(y = total + neg_offset_LFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_dom_LFC_nl <- ggplot(edge_counts_LFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_dom_LFC_lab
p_dom_LFC_nl
# p_dom_LFC_lab # domain_assoc_LFC_labels 386*536
# p_dom_LFC_nl  # domain_assoc_LFC_nolabels 386*536

ggsave(file.path(figdir_h3, "domain_assoc_LFC_labels.tiff"),
       plot = p_dom_LFC_lab, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "domain_assoc_LFC_nolabels.tiff"),
       plot = p_dom_LFC_nl, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.9 Niche-breadth interactions
# -----------------------------------------------------------------------------

edge_nb_LFC <- igraph::as_data_frame(g_LFC, what = "edges") %>%
  mutate(
    from_cat = V(g_LFC)$NB_Category[match(from, V(g_LFC)$name)],
    to_cat   = V(g_LFC)$NB_Category[match(to,   V(g_LFC)$name)]
  ) %>%
  filter(!is.na(from_cat), !is.na(to_cat),
         from_cat != "Unassigned", to_cat != "Unassigned") %>%
  mutate(
    interaction = case_when(
      from_cat == "Generalist"   & to_cat == "Generalist"   ~ "G-G",
      from_cat == "Opportunist"  & to_cat == "Opportunist"  ~ "O-O",
      from_cat == "Specialist"   & to_cat == "Specialist"   ~ "S-S",
      (from_cat == "Generalist"  & to_cat == "Opportunist") |
        (from_cat == "Opportunist" & to_cat == "Generalist") ~ "G-O",
      (from_cat == "Generalist"  & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Generalist") ~ "G-S",
      (from_cat == "Opportunist" & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Opportunist") ~ "O-S",
      TRUE ~ "Other"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_nb_LFC <- edge_nb_LFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_nb_LFC <- edge_nb_LFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_nb_LFC <- edge_counts_nb_LFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_nb_LFC <- left_join(edge_counts_nb_LFC, total_counts_nb_LFC, by = "interaction")

write.csv(edge_nb_LFC,         file.path(tabdir_h3, "nb_edges_LFC.csv"), row.names = FALSE)
write.csv(prop_nb_LFC,         file.path(tabdir_h3, "nb_edge_proportions_LFC.csv"), row.names = FALSE)
write.csv(edge_counts_nb_LFC,  file.path(tabdir_h3, "nb_edge_counts_LFC.csv"), row.names = FALSE)

y_max_nb_LFC <- max(edge_counts_nb_LFC$total, na.rm = TRUE)
neg_offset_nb_LFC <- 0.05 * y_max_nb_LFC

p_nb_LFC_lab <- ggplot(edge_counts_nb_LFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_nb_LFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_nb_LFC, association == "Negative"),
    aes(y = total + neg_offset_nb_LFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_nb_LFC_nl <- ggplot(edge_counts_nb_LFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_nb_LFC_lab
p_nb_LFC_nl
# p_nb_LFC_lab # nb_assoc_LFC_labels 556*512
# p_nb_LFC_nl  # nb_assoc_LFC_nolabels 556*512

ggsave(file.path(figdir_h3, "nb_assoc_LFC_labels.tiff"),
       plot = p_nb_LFC_lab, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "nb_assoc_LFC_nolabels.tiff"),
       plot = p_nb_LFC_nl, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 7. MFC network
# =============================================================================

ps_MFC <- subset_samples(ps, FC == "MFC")
ps_MFC <- prune_taxa(taxa_sums(ps_MFC) > 0, ps_MFC)

otu_MFC <- as.data.frame(otu_table(ps_MFC))

dim(otu_MFC)
otu_MFC2 <- otu_MFC + 1e-6
colSums(otu_MFC2)
otu_MFC3 <- sweep(otu_MFC2, 2, colSums(otu_MFC2), FUN = "/")
colSums(otu_MFC3)
otu_MFC4 <- t(otu_MFC3)
dim(otu_MFC4)
rowSums(otu_MFC4)

net_MFC <- netConstruct(
  otu_MFC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  measurePar  = list(
    method           = "mb",
    sel.criterion    = "stars",
    pulsar.params    = list(
      rep.num         = 30,
      thresh          = 0.05,
      subsample.ratio = 0.8,
      seed            = 2015
    ),
    nlambda          = 20,
    lambda.min.ratio = 0.05
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_MFC_anal <- netAnalyze(net_MFC)
summary(net_MFC_anal)

saveRDS(net_MFC,      file.path(rdsdir_h3, "net_MFC.rds"))
saveRDS(net_MFC_anal, file.path(rdsdir_h3, "net_MFC_anal.rds"))

edgelist_MFC <- net_MFC$edgelist1[
  order(net_MFC$edgelist1$adja, decreasing = TRUE), ]
write.csv(edgelist_MFC, file.path(tabdir_h3, "edgelist_MFC.csv"), row.names = FALSE)

edg_MFC <- net_MFC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g_MFC <- graph_from_data_frame(edg_MFC, directed = FALSE)

V(g_MFC)$type <- ifelse(grepl("^bac_", V(g_MFC)$name), "Bacteria",
                        ifelse(grepl("^fun_", V(g_MFC)$name), "Fungi", "Unknown"))
V(g_MFC)$deg  <- igraph::degree(g_MFC)
V(g_MFC)$eig  <- igraph::eigen_centrality(g_MFC)$vector
V(g_MFC)$comm <- cluster_fast_greedy(g_MFC)$membership

E(g_MFC)$weight <- edg_MFC$w_abs
E(g_MFC)$sign   <- edg_MFC$sign
E(g_MFC)$w_abs  <- edg_MFC$w_abs

tax_df_MFC <- as.data.frame(tax_table(ps_MFC))
tax_df_MFC$node <- rownames(tax_df_MFC)

tax_g_df_MFC <- tax_df_MFC %>%
  filter(node %in% V(g_MFC)$name) %>%
  select(node, NB_Category)

tax_g_df_MFC$NB_Category[is.na(tax_g_df_MFC$NB_Category)] <- "Unassigned"
V(g_MFC)$NB_Category <- tax_g_df_MFC$NB_Category[match(V(g_MFC)$name, tax_g_df_MFC$node)]

ed_df_MFC <- igraph::as_data_frame(g_MFC, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign_MFC <- ed_df_MFC %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0)

if (!"positive" %in% colnames(count_by_sign_MFC)) count_by_sign_MFC$positive <- 0
if (!"negative" %in% colnames(count_by_sign_MFC)) count_by_sign_MFC$negative <- 0

count_by_sign_MFC <- count_by_sign_MFC %>%
  rename(neg = negative, pos = positive)

top_pos_nodes_MFC <- count_by_sign_MFC %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes_MFC <- count_by_sign_MFC %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)

V(g_MFC)$label <- ifelse(V(g_MFC)$name %in% union(top_pos_nodes_MFC, top_neg_nodes_MFC), V(g_MFC)$name, NA)
V(g_MFC)$label_sign <- case_when(
  V(g_MFC)$name %in% top_pos_nodes_MFC & V(g_MFC)$name %in% top_neg_nodes_MFC ~ "both",
  V(g_MFC)$name %in% top_pos_nodes_MFC ~ "positive",
  V(g_MFC)$name %in% top_neg_nodes_MFC ~ "negative",
  TRUE ~ NA_character_
)

g_MFC_filt <- igraph::delete_edges(
  g_MFC,
  E(g_MFC)[abs(weight) < 0.3]
)

E(g_MFC_filt)$sign <- as.character(E(g_MFC_filt)$sign)
E(g_MFC_filt)$sign[E(g_MFC_filt)$sign %in% c("Positive", "POSITIVE")] <- "positive"
E(g_MFC_filt)$sign[E(g_MFC_filt)$sign %in% c("Negative", "NEGATIVE")] <- "negative"
E(g_MFC_filt)$sign <- factor(E(g_MFC_filt)$sign, levels = c("positive", "negative"))

V(g_MFC_filt)$degree <- igraph::degree(g_MFC_filt)
V(g_MFC_filt)$deg <- igraph::degree(g_MFC_filt, mode = "all", loops = FALSE)
node_size_MFC <- scales::rescale(V(g_MFC_filt)$deg, to = c(2, 10))
V(g_MFC_filt)$NB_Category <- factor(V(g_MFC_filt)$NB_Category, levels = names(nb_colors))

set.seed(42)
lay_MFC <- create_layout(g_MFC_filt, layout = "circle")

lay_df <- as.data.frame(lay_MFC) %>% 
  filter(!is.na(label))

p_net_MFC <- ggraph(lay_MFC) +
  geom_edge_link(aes(edge_colour = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_colour_manual(values = c(positive = "gray70",
                                      negative = "gray40"),
                           drop = FALSE,
                           name = "Association") +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.5, 0.8)) +
  guides(edge_width = "none", edge_alpha = "none") +
  geom_node_point(aes(fill = NB_Category),
                  shape = 21, size = node_size_MFC, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = nb_colors, name = "NB Category") +
  theme_void() +
  coord_equal()

p_net_MFC
ggsave(file.path(figdir_h3, "network_MFC.tiff"),
       plot = p_net_MFC, width = 9, height = 9, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_df_MFC <- igraph::as_data_frame(g_MFC, what = "edges") %>%
  mutate(
    from_type = V(g_MFC)$type[match(from, V(g_MFC)$name)],
    to_type   = V(g_MFC)$type[match(to,   V(g_MFC)$name)],
    interaction = case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_df_MFC <- edge_df_MFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_MFC <- edge_df_MFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_MFC <- edge_counts_MFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_MFC <- left_join(edge_counts_MFC, total_counts_MFC, by = "interaction")

write.csv(edge_df_MFC,     file.path(tabdir_h3, "domain_edges_MFC.csv"), row.names = FALSE)
write.csv(prop_df_MFC,     file.path(tabdir_h3, "domain_edge_proportions_MFC.csv"), row.names = FALSE)
write.csv(edge_counts_MFC, file.path(tabdir_h3, "domain_edge_counts_MFC.csv"), row.names = FALSE)

y_max_MFC <- max(edge_counts_MFC$total, na.rm = TRUE)
neg_offset_MFC <- 0.05 * y_max_MFC

p_dom_MFC_lab <- ggplot(edge_counts_MFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_MFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_MFC, association == "Negative"),
    aes(y = total + neg_offset_MFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_dom_MFC_nl <- ggplot(edge_counts_MFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "domain_assoc_MFC_labels.tiff"),
       plot = p_dom_MFC_lab, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "domain_assoc_MFC_nolabels.tiff"),
       plot = p_dom_MFC_nl, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_nb_MFC <- igraph::as_data_frame(g_MFC, what = "edges") %>%
  mutate(
    from_cat = V(g_MFC)$NB_Category[match(from, V(g_MFC)$name)],
    to_cat   = V(g_MFC)$NB_Category[match(to,   V(g_MFC)$name)]
  ) %>%
  filter(!is.na(from_cat), !is.na(to_cat),
         from_cat != "Unassigned", to_cat != "Unassigned") %>%
  mutate(
    interaction = case_when(
      from_cat == "Generalist"   & to_cat == "Generalist"   ~ "G-G",
      from_cat == "Opportunist"  & to_cat == "Opportunist"  ~ "O-O",
      from_cat == "Specialist"   & to_cat == "Specialist"   ~ "S-S",
      (from_cat == "Generalist"  & to_cat == "Opportunist") |
        (from_cat == "Opportunist" & to_cat == "Generalist") ~ "G-O",
      (from_cat == "Generalist"  & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Generalist") ~ "G-S",
      (from_cat == "Opportunist" & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Opportunist") ~ "O-S",
      TRUE ~ "Other"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_nb_MFC <- edge_nb_MFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_nb_MFC <- edge_nb_MFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_nb_MFC <- edge_counts_nb_MFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_nb_MFC <- left_join(edge_counts_nb_MFC, total_counts_nb_MFC, by = "interaction")

write.csv(edge_nb_MFC,        file.path(tabdir_h3, "nb_edges_MFC.csv"), row.names = FALSE)
write.csv(prop_nb_MFC,        file.path(tabdir_h3, "nb_edge_proportions_MFC.csv"), row.names = FALSE)
write.csv(edge_counts_nb_MFC, file.path(tabdir_h3, "nb_edge_counts_MFC.csv"), row.names = FALSE)

y_max_nb_MFC <- max(edge_counts_nb_MFC$total, na.rm = TRUE)
neg_offset_nb_MFC <- 0.05 * y_max_nb_MFC

p_nb_MFC_lab <- ggplot(edge_counts_nb_MFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_nb_MFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_nb_MFC, association == "Negative"),
    aes(y = total + neg_offset_nb_MFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_nb_MFC_nl <- ggplot(edge_counts_nb_MFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "nb_assoc_MFC_labels.tiff"),
       plot = p_nb_MFC_lab, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "nb_assoc_MFC_nolabels.tiff"),
       plot = p_nb_MFC_nl, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 8. HFC network
# =============================================================================

ps_HFC <- subset_samples(ps, FC == "HFC")
ps_HFC <- prune_taxa(taxa_sums(ps_HFC) > 0, ps_HFC)

otu_HFC <- as.data.frame(otu_table(ps_HFC))

dim(otu_HFC)
otu_HFC2 <- otu_HFC + 1e-6
colSums(otu_HFC2)
otu_HFC3 <- sweep(otu_HFC2, 2, colSums(otu_HFC2), FUN = "/")
colSums(otu_HFC3)
otu_HFC4 <- t(otu_HFC3)
dim(otu_HFC4)
rowSums(otu_HFC4)

net_HFC <- netConstruct(
  otu_HFC4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  measurePar  = list(
    method           = "mb",
    sel.criterion    = "stars",
    pulsar.params    = list(
      rep.num         = 30,
      thresh          = 0.05,
      subsample.ratio = 0.8,
      seed            = 2015
    ),
    nlambda          = 20,
    lambda.min.ratio = 0.05
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_HFC_anal <- netAnalyze(net_HFC)
summary(net_HFC_anal)

saveRDS(net_HFC,      file.path(rdsdir_h3, "net_HFC.rds"))
saveRDS(net_HFC_anal, file.path(rdsdir_h3, "net_HFC_anal.rds"))

edgelist_HFC <- net_HFC$edgelist1[
  order(net_HFC$edgelist1$adja, decreasing = TRUE), ]
write.csv(edgelist_HFC, file.path(tabdir_h3, "edgelist_HFC.csv"), row.names = FALSE)

edg_HFC <- net_HFC$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g_HFC <- graph_from_data_frame(edg_HFC, directed = FALSE)

V(g_HFC)$type <- ifelse(grepl("^bac_", V(g_HFC)$name), "Bacteria",
                        ifelse(grepl("^fun_", V(g_HFC)$name), "Fungi", "Unknown"))
V(g_HFC)$deg  <- igraph::degree(g_HFC)
V(g_HFC)$eig  <- igraph::eigen_centrality(g_HFC)$vector
V(g_HFC)$comm <- cluster_fast_greedy(g_HFC)$membership

E(g_HFC)$weight <- edg_HFC$w_abs
E(g_HFC)$sign   <- edg_HFC$sign
E(g_HFC)$w_abs  <- edg_HFC$w_abs

tax_df_HFC <- as.data.frame(tax_table(ps_HFC))
tax_df_HFC$node <- rownames(tax_df_HFC)

tax_g_df_HFC <- tax_df_HFC %>%
  filter(node %in% V(g_HFC)$name) %>%
  select(node, NB_Category)

tax_g_df_HFC$NB_Category[is.na(tax_g_df_HFC$NB_Category)] <- "Unassigned"
V(g_HFC)$NB_Category <- tax_g_df_HFC$NB_Category[match(V(g_HFC)$name, tax_g_df_HFC$node)]

ed_df_HFC <- igraph::as_data_frame(g_HFC, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign_HFC <- ed_df_HFC %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0)

if (!"positive" %in% colnames(count_by_sign_HFC)) count_by_sign_HFC$positive <- 0
if (!"negative" %in% colnames(count_by_sign_HFC)) count_by_sign_HFC$negative <- 0

count_by_sign_HFC <- count_by_sign_HFC %>%
  rename(neg = negative, pos = positive)

top_pos_nodes_HFC <- count_by_sign_HFC %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes_HFC <- count_by_sign_HFC %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)

V(g_HFC)$label <- ifelse(V(g_HFC)$name %in% union(top_pos_nodes_HFC, top_neg_nodes_HFC), V(g_HFC)$name, NA)
V(g_HFC)$label_sign <- case_when(
  V(g_HFC)$name %in% top_pos_nodes_HFC & V(g_HFC)$name %in% top_neg_nodes_HFC ~ "both",
  V(g_HFC)$name %in% top_pos_nodes_HFC ~ "positive",
  V(g_HFC)$name %in% top_neg_nodes_HFC ~ "negative",
  TRUE ~ NA_character_
)

g_HFC_filt <- igraph::delete_edges(
  g_HFC,
  E(g_HFC)[abs(weight) < 0.3]
)

E(g_HFC_filt)$sign <- as.character(E(g_HFC_filt)$sign)
E(g_HFC_filt)$sign[E(g_HFC_filt)$sign %in% c("Positive", "POSITIVE")] <- "positive"
E(g_HFC_filt)$sign[E(g_HFC_filt)$sign %in% c("Negative", "NEGATIVE")] <- "negative"
E(g_HFC_filt)$sign <- factor(E(g_HFC_filt)$sign, levels = c("positive", "negative"))

V(g_HFC_filt)$degree <- igraph::degree(g_HFC_filt)
V(g_HFC_filt)$deg <- igraph::degree(g_HFC_filt, mode = "all", loops = FALSE)
node_size_HFC <- scales::rescale(V(g_HFC_filt)$deg, to = c(2, 10))
V(g_HFC_filt)$NB_Category <- factor(V(g_HFC_filt)$NB_Category, levels = names(nb_colors))

set.seed(42)
lay_HFC <- create_layout(g_HFC_filt, layout = "circle")

lay_df <- as.data.frame(lay_HFC) %>% 
  filter(!is.na(label))

p_net_HFC <- ggraph(lay_HFC) +
  geom_edge_link(aes(edge_colour = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_colour_manual(values = c(positive = "gray70",
                                      negative = "gray40"),
                           drop = FALSE,
                           name = "Association") +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.5, 0.8)) +
  guides(edge_width = "none", edge_alpha = "none") +
  geom_node_point(aes(fill = NB_Category),
                  shape = 21, size = node_size_HFC, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = nb_colors, name = "NB Category") +
  theme_void() +
  coord_equal()

p_net_HFC
ggsave(file.path(figdir_h3, "network_HFC.tiff"),
       plot = p_net_HFC, width = 9, height = 9, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_df_HFC <- igraph::as_data_frame(g_HFC, what = "edges") %>%
  mutate(
    from_type = V(g_HFC)$type[match(from, V(g_HFC)$name)],
    to_type   = V(g_HFC)$type[match(to,   V(g_HFC)$name)],
    interaction = case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_df_HFC <- edge_df_HFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_HFC <- edge_df_HFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_HFC <- edge_counts_HFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_HFC <- left_join(edge_counts_HFC, total_counts_HFC, by = "interaction")

write.csv(edge_df_HFC,     file.path(tabdir_h3, "domain_edges_HFC.csv"), row.names = FALSE)
write.csv(prop_df_HFC,     file.path(tabdir_h3, "domain_edge_proportions_HFC.csv"), row.names = FALSE)
write.csv(edge_counts_HFC, file.path(tabdir_h3, "domain_edge_counts_HFC.csv"), row.names = FALSE)

y_max_HFC <- max(edge_counts_HFC$total, na.rm = TRUE)
neg_offset_HFC <- 0.05 * y_max_HFC

p_dom_HFC_lab <- ggplot(edge_counts_HFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_HFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_HFC, association == "Negative"),
    aes(y = total + neg_offset_HFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_dom_HFC_nl <- ggplot(edge_counts_HFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "domain_assoc_HFC_labels.tiff"),
       plot = p_dom_HFC_lab, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "domain_assoc_HFC_nolabels.tiff"),
       plot = p_dom_HFC_nl, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_nb_HFC <- igraph::as_data_frame(g_HFC, what = "edges") %>%
  mutate(
    from_cat = V(g_HFC)$NB_Category[match(from, V(g_HFC)$name)],
    to_cat   = V(g_HFC)$NB_Category[match(to,   V(g_HFC)$name)]
  ) %>%
  filter(!is.na(from_cat), !is.na(to_cat),
         from_cat != "Unassigned", to_cat != "Unassigned") %>%
  mutate(
    interaction = case_when(
      from_cat == "Generalist"   & to_cat == "Generalist"   ~ "G-G",
      from_cat == "Opportunist"  & to_cat == "Opportunist"  ~ "O-O",
      from_cat == "Specialist"   & to_cat == "Specialist"   ~ "S-S",
      (from_cat == "Generalist"  & to_cat == "Opportunist") |
        (from_cat == "Opportunist" & to_cat == "Generalist") ~ "G-O",
      (from_cat == "Generalist"  & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Generalist") ~ "G-S",
      (from_cat == "Opportunist" & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Opportunist") ~ "O-S",
      TRUE ~ "Other"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_nb_HFC <- edge_nb_HFC %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_nb_HFC <- edge_nb_HFC %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_nb_HFC <- edge_counts_nb_HFC %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_nb_HFC <- left_join(edge_counts_nb_HFC, total_counts_nb_HFC, by = "interaction")

write.csv(edge_nb_HFC,        file.path(tabdir_h3, "nb_edges_HFC.csv"), row.names = FALSE)
write.csv(prop_nb_HFC,        file.path(tabdir_h3, "nb_edge_proportions_HFC.csv"), row.names = FALSE)
write.csv(edge_counts_nb_HFC, file.path(tabdir_h3, "nb_edge_counts_HFC.csv"), row.names = FALSE)

y_max_nb_HFC <- max(edge_counts_nb_HFC$total, na.rm = TRUE)
neg_offset_nb_HFC <- 0.05 * y_max_nb_HFC

p_nb_HFC_lab <- ggplot(edge_counts_nb_HFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_nb_HFC, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_nb_HFC, association == "Negative"),
    aes(y = total + neg_offset_nb_HFC, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_nb_HFC_nl <- ggplot(edge_counts_nb_HFC, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "nb_assoc_HFC_labels.tiff"),
       plot = p_nb_HFC_lab, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "nb_assoc_HFC_nolabels.tiff"),
       plot = p_nb_HFC_nl, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 9. Forest network
# =============================================================================

ps_F <- subset_samples(ps, FC == "F")
ps_F <- prune_taxa(taxa_sums(ps_F) > 0, ps_F)

otu_F <- as.data.frame(otu_table(ps_F))

dim(otu_F)
otu_F2 <- otu_F + 1e-6
colSums(otu_F2)
otu_F3 <- sweep(otu_F2, 2, colSums(otu_F2), FUN = "/")
colSums(otu_F3)
otu_F4 <- t(otu_F3)
dim(otu_F4)
rowSums(otu_F4)

net_F <- netConstruct(
  otu_F4,
  measure     = "spieceasi",
  normMethod  = "none",
  zeroMethod  = "pseudo",
  filtTax     = "highestVar",
  filtTaxPar  = list(highestVar = fVar),
  measurePar  = list(
    method           = "mb",
    sel.criterion    = "stars",
    pulsar.params    = list(
      rep.num         = 30,
      thresh          = 0.05,
      subsample.ratio = 0.8,
      seed            = 2015
    ),
    nlambda          = 20,
    lambda.min.ratio = 0.05
  ),
  sparsMethod = "none",
  verbose     = 3,
  seed        = 2015
)

net_F_anal <- netAnalyze(net_F)
summary(net_F_anal)

saveRDS(net_F,      file.path(rdsdir_h3, "net_Forest.rds"))
saveRDS(net_F_anal, file.path(rdsdir_h3, "net_Forest_anal.rds"))

edgelist_F <- net_F$edgelist1[
  order(net_F$edgelist1$adja, decreasing = TRUE), ]
write.csv(edgelist_F, file.path(tabdir_h3, "edgelist_Forest.csv"), row.names = FALSE)

edg_F <- net_F$edgelist1 %>%
  mutate(sign = ifelse(asso >= 0, "positive", "negative"),
         w_abs = abs(asso))

g_F <- graph_from_data_frame(edg_F, directed = FALSE)

V(g_F)$type <- ifelse(grepl("^bac_", V(g_F)$name), "Bacteria",
                      ifelse(grepl("^fun_", V(g_F)$name), "Fungi", "Unknown"))
V(g_F)$deg  <- igraph::degree(g_F)
V(g_F)$eig  <- igraph::eigen_centrality(g_F)$vector
V(g_F)$comm <- cluster_fast_greedy(g_F)$membership

E(g_F)$weight <- edg_F$w_abs
E(g_F)$sign   <- edg_F$sign
E(g_F)$w_abs  <- edg_F$w_abs

tax_df_F <- as.data.frame(tax_table(ps_F))
tax_df_F$node <- rownames(tax_df_F)

tax_g_df_F <- tax_df_F %>%
  filter(node %in% V(g_F)$name) %>%
  select(node, NB_Category)

tax_g_df_F$NB_Category[is.na(tax_g_df_F$NB_Category)] <- "Unassigned"
V(g_F)$NB_Category <- tax_g_df_F$NB_Category[match(V(g_F)$name, tax_g_df_F$node)]

ed_df_F <- igraph::as_data_frame(g_F, what = "edges") %>%
  transmute(from, to, sign = as.character(sign))

count_by_sign_F <- ed_df_F %>%
  pivot_longer(c(from, to), names_to = "end", values_to = "node") %>%
  count(node, sign) %>%
  pivot_wider(names_from = sign, values_from = n, values_fill = 0)

if (!"positive" %in% colnames(count_by_sign_F)) count_by_sign_F$positive <- 0
if (!"negative" %in% colnames(count_by_sign_F)) count_by_sign_F$negative <- 0

count_by_sign_F <- count_by_sign_F %>%
  rename(neg = negative, pos = positive)

top_pos_nodes_F <- count_by_sign_F %>% arrange(desc(pos)) %>% slice_head(n = 5) %>% pull(node)
top_neg_nodes_F <- count_by_sign_F %>% arrange(desc(neg)) %>% slice_head(n = 5) %>% pull(node)

V(g_F)$label <- ifelse(V(g_F)$name %in% union(top_pos_nodes_F, top_neg_nodes_F), V(g_F)$name, NA)
V(g_F)$label_sign <- case_when(
  V(g_F)$name %in% top_pos_nodes_F & V(g_F)$name %in% top_neg_nodes_F ~ "both",
  V(g_F)$name %in% top_pos_nodes_F ~ "positive",
  V(g_F)$name %in% top_neg_nodes_F ~ "negative",
  TRUE ~ NA_character_
)

g_F_filt <- igraph::delete_edges(
  g_F,
  E(g_F)[abs(weight) < 0.3]
)

E(g_F_filt)$sign <- as.character(E(g_F_filt)$sign)
E(g_F_filt)$sign[E(g_F_filt)$sign %in% c("Positive", "POSITIVE")] <- "positive"
E(g_F_filt)$sign[E(g_F_filt)$sign %in% c("Negative", "NEGATIVE")] <- "negative"
E(g_F_filt)$sign <- factor(E(g_F_filt)$sign, levels = c("positive", "negative"))

V(g_F_filt)$degree <- igraph::degree(g_F_filt)
V(g_F_filt)$deg <- igraph::degree(g_F_filt, mode = "all", loops = FALSE)
node_size_F <- scales::rescale(V(g_F_filt)$deg, to = c(2, 10))
V(g_F_filt)$NB_Category <- factor(V(g_F_filt)$NB_Category, levels = names(nb_colors))

set.seed(42)
lay_F <- create_layout(g_F_filt, layout = "circle")

lay_df <- as.data.frame(lay_F) %>% 
  filter(!is.na(label))

p_net_F <- ggraph(lay_F) +
  geom_edge_link(aes(edge_colour = sign,
                     edge_width = abs(weight),
                     edge_alpha = abs(weight)),
                 show.legend = TRUE) +
  scale_edge_colour_manual(values = c(positive = "gray70",
                                      negative = "gray40"),
                           drop = FALSE,
                           name = "Association") +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.5, 0.8)) +
  guides(edge_width = "none", edge_alpha = "none") +
  geom_node_point(aes(fill = NB_Category),
                  shape = 21, size = node_size_F, color = "grey25", stroke = 0.4) +
  scale_fill_manual(values = nb_colors, name = "NB Category") +
  theme_void() +
  coord_equal()

p_net_F
ggsave(file.path(figdir_h3, "network_Forest.tiff"),
       plot = p_net_F, width = 9, height = 9, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_df_F <- igraph::as_data_frame(g_F, what = "edges") %>%
  mutate(
    from_type = V(g_F)$type[match(from, V(g_F)$name)],
    to_type   = V(g_F)$type[match(to,   V(g_F)$name)],
    interaction = case_when(
      from_type == "Bacteria" & to_type == "Bacteria" ~ "B-B",
      from_type == "Fungi"    & to_type == "Fungi"    ~ "F-F",
      TRUE                                            ~ "B-F"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_df_F <- edge_df_F %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_F <- edge_df_F %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_F <- edge_counts_F %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_F <- left_join(edge_counts_F, total_counts_F, by = "interaction")

write.csv(edge_df_F,     file.path(tabdir_h3, "domain_edges_Forest.csv"), row.names = FALSE)
write.csv(prop_df_F,     file.path(tabdir_h3, "domain_edge_proportions_Forest.csv"), row.names = FALSE)
write.csv(edge_counts_F, file.path(tabdir_h3, "domain_edge_counts_Forest.csv"), row.names = FALSE)

y_max_F <- max(edge_counts_F$total, na.rm = TRUE)
neg_offset_F <- 0.05 * y_max_F

p_dom_F_lab <- ggplot(edge_counts_F, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_F, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_F, association == "Negative"),
    aes(y = total + neg_offset_F, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_dom_F_nl <- ggplot(edge_counts_F, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(x = "", y = "Number of Edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "domain_assoc_Forest_labels.tiff"),
       plot = p_dom_F_lab, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "domain_assoc_Forest_nolabels.tiff"),
       plot = p_dom_F_nl, width = 3.86, height = 5.36, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

edge_nb_F <- igraph::as_data_frame(g_F, what = "edges") %>%
  mutate(
    from_cat = V(g_F)$NB_Category[match(from, V(g_F)$name)],
    to_cat   = V(g_F)$NB_Category[match(to,   V(g_F)$name)]
  ) %>%
  filter(!is.na(from_cat), !is.na(to_cat),
         from_cat != "Unassigned", to_cat != "Unassigned") %>%
  mutate(
    interaction = case_when(
      from_cat == "Generalist"   & to_cat == "Generalist"   ~ "G-G",
      from_cat == "Opportunist"  & to_cat == "Opportunist"  ~ "O-O",
      from_cat == "Specialist"   & to_cat == "Specialist"   ~ "S-S",
      (from_cat == "Generalist"  & to_cat == "Opportunist") |
        (from_cat == "Opportunist" & to_cat == "Generalist") ~ "G-O",
      (from_cat == "Generalist"  & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Generalist") ~ "G-S",
      (from_cat == "Opportunist" & to_cat == "Specialist") |
        (from_cat == "Specialist" & to_cat == "Opportunist") ~ "O-S",
      TRUE ~ "Other"
    ),
    association = ifelse(as.character(sign) %in% c("positive", "Positive"), "Positive", "Negative")
  )

prop_nb_F <- edge_nb_F %>%
  group_by(interaction, association) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(interaction) %>%
  mutate(proportion = count / sum(count))

edge_counts_nb_F <- edge_nb_F %>%
  count(interaction, association, name = "n") %>%
  group_by(interaction) %>%
  mutate(
    percentage = n / sum(n) * 100,
    label = paste0(n, "\n(", round(percentage, 1), "%)")
  ) %>%
  ungroup()

total_counts_nb_F <- edge_counts_nb_F %>%
  group_by(interaction) %>%
  summarise(total = sum(n), .groups = "drop")

edge_counts_nb_F <- left_join(edge_counts_nb_F, total_counts_nb_F, by = "interaction")

write.csv(edge_nb_F,        file.path(tabdir_h3, "nb_edges_Forest.csv"), row.names = FALSE)
write.csv(prop_nb_F,        file.path(tabdir_h3, "nb_edge_proportions_Forest.csv"), row.names = FALSE)
write.csv(edge_counts_nb_F, file.path(tabdir_h3, "nb_edge_counts_Forest.csv"), row.names = FALSE)

y_max_nb_F <- max(edge_counts_nb_F$total, na.rm = TRUE)
neg_offset_nb_F <- 0.05 * y_max_nb_F

p_nb_F_lab <- ggplot(edge_counts_nb_F, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(
    data = dplyr::filter(edge_counts_nb_F, association == "Positive"),
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4
  ) +
  geom_text(
    data = dplyr::filter(edge_counts_nb_F, association == "Negative"),
    aes(y = total + neg_offset_nb_F, label = label),
    color = "black", size = 4
  ) +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

p_nb_F_nl <- ggplot(edge_counts_nb_F, aes(x = interaction, y = n, fill = association)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = assoc_colors) +
  labs(y = "Number of edges", fill = "Association") +
  scale_y_continuous(limits = c(0, 4000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black")
  )

ggsave(file.path(figdir_h3, "nb_assoc_Forest_labels.tiff"),
       plot = p_nb_F_lab, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

ggsave(file.path(figdir_h3, "nb_assoc_Forest_nolabels.tiff"),
       plot = p_nb_F_nl, width = 5.56, height = 5.12, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# Save graph objects
# =============================================================================

saveRDS(g_LFC, file.path(rdsdir_h3, "graph_LFC.rds"))
saveRDS(g_MFC, file.path(rdsdir_h3, "graph_MFC.rds"))
saveRDS(g_HFC, file.path(rdsdir_h3, "graph_HFC.rds"))
saveRDS(g_F,   file.path(rdsdir_h3, "graph_Forest.rds"))


# =============================================================================
# Table summary
# =============================================================================

get_network_summary <- function(g){
  
  # --- nodos ---
  n_nodes <- igraph::vcount(g)
  
  node_df <- data.frame(
    node = igraph::V(g)$name,
    type = igraph::V(g)$type,
    NB   = igraph::V(g)$NB_Category,
    stringsAsFactors = FALSE
  )
  
  n_bac <- sum(node_df$type == "Bacteria", na.rm = TRUE)
  n_fun <- sum(node_df$type == "Fungi", na.rm = TRUE)
  
  n_gen <- sum(node_df$NB == "Generalist", na.rm = TRUE)
  n_opp <- sum(node_df$NB == "Opportunist", na.rm = TRUE)
  n_spe <- sum(node_df$NB == "Specialist", na.rm = TRUE)
  n_un  <- sum(node_df$NB == "Unassigned", na.rm = TRUE)
  
  # --- edges ---
  edge_df <- igraph::as_data_frame(g, what = "edges") %>%
    dplyr::mutate(
      association = ifelse(sign == "positive", "Positive", "Negative"),
      from_type = igraph::V(g)$type[match(from, igraph::V(g)$name)],
      to_type   = igraph::V(g)$type[match(to,   igraph::V(g)$name)],
      from_nb   = igraph::V(g)$NB_Category[match(from, igraph::V(g)$name)],
      to_nb     = igraph::V(g)$NB_Category[match(to,   igraph::V(g)$name)]
    )
  
  # totales
  n_pos <- sum(edge_df$association == "Positive", na.rm = TRUE)
  n_neg <- sum(edge_df$association == "Negative", na.rm = TRUE)
  
  pos_perc <- n_pos / (n_pos + n_neg) * 100
  
  # --- topología ---
  clustering <- igraph::transitivity(g, type = "global")
  modularity_val <- igraph::modularity(igraph::cluster_fast_greedy(g))
  edge_density_val <- igraph::edge_density(g)
  
  avg_path <- suppressWarnings(
    igraph::average.path.length(g, directed = FALSE)
  )
  
  # natural connectivity
  A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  eig <- eigen(A, only.values = TRUE)$values
  nat_conn <- log(mean(exp(Re(eig))))
  
  # tipo interacción
  is_BB <- edge_df$from_type == "Bacteria" & edge_df$to_type == "Bacteria"
  is_FF <- edge_df$from_type == "Fungi"    & edge_df$to_type == "Fungi"
  is_BF <- !is_BB & !is_FF
  
  # NB interacción
  comb_nb <- paste(edge_df$from_nb, edge_df$to_nb, sep = "-")
  comb_nb <- ifelse(comb_nb == "Opportunist-Generalist", "Generalist-Opportunist", comb_nb)
  comb_nb <- ifelse(comb_nb == "Specialist-Generalist",  "Generalist-Specialist", comb_nb)
  comb_nb <- ifelse(comb_nb == "Specialist-Opportunist", "Opportunist-Specialist", comb_nb)
  
  edge_df$comb_nb <- comb_nb
  
  # --- conteos positivos ---
  pos_df <- edge_df %>% dplyr::filter(association == "Positive")
  
  pos_BB <- sum(is_BB & edge_df$association == "Positive", na.rm = TRUE)
  pos_BF <- sum(is_BF & edge_df$association == "Positive", na.rm = TRUE)
  pos_FF <- sum(is_FF & edge_df$association == "Positive", na.rm = TRUE)
  
  pos_GG <- sum(pos_df$comb_nb == "Generalist-Generalist", na.rm = TRUE)
  pos_GO <- sum(pos_df$comb_nb == "Generalist-Opportunist", na.rm = TRUE)
  pos_GS <- sum(pos_df$comb_nb == "Generalist-Specialist", na.rm = TRUE)
  pos_OO <- sum(pos_df$comb_nb == "Opportunist-Opportunist", na.rm = TRUE)
  pos_OS <- sum(pos_df$comb_nb == "Opportunist-Specialist", na.rm = TRUE)
  pos_SS <- sum(pos_df$comb_nb == "Specialist-Specialist", na.rm = TRUE)
  
  # --- conteos negativos ---
  neg_df <- edge_df %>% dplyr::filter(association == "Negative")
  
  neg_BB <- sum(is_BB & edge_df$association == "Negative", na.rm = TRUE)
  neg_BF <- sum(is_BF & edge_df$association == "Negative", na.rm = TRUE)
  neg_FF <- sum(is_FF & edge_df$association == "Negative", na.rm = TRUE)
  
  neg_GG <- sum(neg_df$comb_nb == "Generalist-Generalist", na.rm = TRUE)
  neg_GO <- sum(neg_df$comb_nb == "Generalist-Opportunist", na.rm = TRUE)
  neg_GS <- sum(neg_df$comb_nb == "Generalist-Specialist", na.rm = TRUE)
  neg_OO <- sum(neg_df$comb_nb == "Opportunist-Opportunist", na.rm = TRUE)
  neg_OS <- sum(neg_df$comb_nb == "Opportunist-Specialist", na.rm = TRUE)
  neg_SS <- sum(neg_df$comb_nb == "Specialist-Specialist", na.rm = TRUE)
  
  list(
    clustering = clustering,
    modularity = modularity_val,
    pos_perc   = pos_perc,
    density    = edge_density_val,
    nat_conn   = nat_conn,
    avg_path   = avg_path,
    nodes      = n_nodes,
    bac = n_bac,
    fun = n_fun,
    gen = n_gen,
    opp = n_opp,
    spe = n_spe,
    un  = n_un,
    pos = n_pos,
    neg = n_neg,
    pos_BB = pos_BB, pos_BF = pos_BF, pos_FF = pos_FF,
    pos_GG = pos_GG, pos_GO = pos_GO, pos_GS = pos_GS,
    pos_OO = pos_OO, pos_OS = pos_OS, pos_SS = pos_SS,
    neg_BB = neg_BB, neg_BF = neg_BF, neg_FF = neg_FF,
    neg_GG = neg_GG, neg_GO = neg_GO, neg_GS = neg_GS,
    neg_OO = neg_OO, neg_OS = neg_OS, neg_SS = neg_SS
  )
}

res_LFC <- get_network_summary(g_LFC)
res_MFC <- get_network_summary(g_MFC)
res_HFC <- get_network_summary(g_HFC)
res_F   <- get_network_summary(g_F)

# =============================================================================
# Build final topology table
# =============================================================================

format_node <- function(n, total){
  paste0(n, " (", round(100 * n / total, 1), ")")
}

format_edge <- function(n, total){
  paste0(n, " (", round(100 * n / total, 1), ")")
}

make_col <- function(res){
  
  total_nodes <- res$nodes
  total_pos   <- res$pos
  total_neg   <- res$neg
  
  c(
    round(res$clustering, 5),
    round(res$modularity, 5),
    round(res$pos_perc, 5),
    round(res$density, 5),
    round(res$nat_conn, 5),
    round(res$avg_path, 5),
    res$nodes,
    
    format_node(res$bac, total_nodes),
    format_node(res$fun, total_nodes),
    format_node(res$gen, total_nodes),
    format_node(res$opp, total_nodes),
    format_node(res$spe, total_nodes),
    format_node(res$un,  total_nodes),
    
    res$pos,
    format_edge(res$pos_BB, total_pos),
    format_edge(res$pos_BF, total_pos),
    format_edge(res$pos_FF, total_pos),
    format_edge(res$pos_GG, total_pos),
    format_edge(res$pos_GO, total_pos),
    format_edge(res$pos_GS, total_pos),
    format_edge(res$pos_OO, total_pos),
    format_edge(res$pos_OS, total_pos),
    format_edge(res$pos_SS, total_pos),
    
    res$neg,
    format_edge(res$neg_BB, total_neg),
    format_edge(res$neg_BF, total_neg),
    format_edge(res$neg_FF, total_neg),
    format_edge(res$neg_GG, total_neg),
    format_edge(res$neg_GO, total_neg),
    format_edge(res$neg_GS, total_neg),
    format_edge(res$neg_OO, total_neg),
    format_edge(res$neg_OS, total_neg),
    format_edge(res$neg_SS, total_neg)
  )
}

row_names <- c(
  "Clustering coefficient",
  "Modularity",
  "Positive edge percentage",
  "Edge density",
  "Natural connectivity",
  "Average path length",
  "Nodes",
  "N° of bacteria",
  "N° of fungi",
  "N° of generalist taxa",
  "N° of opportunist taxa",
  "N° of specialist taxa",
  "N° of unassigned taxa",
  "Positive edges",
  "Positive B-B edges",
  "Positive B-F edges",
  "Positive F-F edges",
  "Positive G-G edges",
  "Positive G-O edges",
  "Positive G-S edges",
  "Positive O-O edges",
  "Positive O-S edges",
  "Positive S-S edges",
  "Negative edges",
  "Negative B-B edges",
  "Negative B-F edges",
  "Negative F-F edges",
  "Negative G-G edges",
  "Negative G-O edges",
  "Negative G-S edges",
  "Negative O-O edges",
  "Negative O-S edges",
  "Negative S-S edges"
)

table_out <- data.frame(
  Network_property = row_names,
  LFC = make_col(res_LFC),
  MFC = make_col(res_MFC),
  HFC = make_col(res_HFC),
  F   = make_col(res_F),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

table_out

# =============================================================================
# Session info
# =============================================================================

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir_h3, "sessionInfo_networks.txt"))

