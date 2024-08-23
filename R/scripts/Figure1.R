
#load libraries-----
library(dplyr)
library(Rubrary)
library(data.table)
library(tibble)
library(ComplexHeatmap)
library(ggpubr)
library(xlsx)

setwd("/Users/Zero/Dropbox/Graeber/OVP/OV_NE")
source("./Figures/R/functions/Fig1_functions.R")
#
#load data------
ucla_ov_paired = read.delim(
  file = "ovarian_ucla12345_rsem_genes_upper_norm_counts_coding_log2_paired.txt",
  stringsAsFactors = F
)
uc5p.rbl <- fread("data/ucla5_patch_rbe2f_ssgsea.txt")
ov_anno_rc = read.delim("data/UCLAOvarian_CombinedBatches_Annotation_v13_rcp.txt", stringsAsFactors = F)
patch_anno_full <- read.delim(
  file = "data/patch.signature.scores_anno_full.txt",
  stringsAsFactors = F
)
patch.r = read.delim(
  file = "data/Patch_ovarian_rsem_genes_upper_norm_counts_coding_log2.txt",
  row.names = 1,
  stringsAsFactors = F
)
#
#GSEA Fig 1D & S1E-----

#(Bubble Plot - Fig 1D)
ucp_fgsea = read.delim(
  file = "./data/fGSEA_UCLAAllPatch_deseq_recur_fig.txt",
  stringsAsFactors = F
)
#select top & bottom 10 pathways based on NES
top_bottom_paths <- ucp_fgsea %>%
  arrange(desc(NES)) %>%
  top_n(10, NES) %>%
  bind_rows(
    ucp_fgsea %>%
      arrange(NES) %>%
      top_n(10, -NES)
  ) %>%
  arrange(desc(NES)) %>%
  mutate(PathRank = factor(pathway, levels = rev(pathway)))

# Splitting data into positive & negative
positive_paths <- top_bottom_paths %>%
  filter(NES >= 0)

negative_paths <- top_bottom_paths %>%
  filter(NES < 0) %>%
  mutate(PathRank = factor(pathway, levels = rev(levels(PathRank))))

# Create positive NES plot
positive_plot <- ggplot(positive_paths, aes(x = NES, y = PathRank, size = size, color = SignedLogpadj, fill = signedlogp)) +
  geom_point(shape = 21, stroke = 1.5, alpha = 0.8) +
  scale_size_area(max_size = 15) +
  scale_color_gradient(low = "#fcbad3", high = "#c94c4c") +
  scale_fill_gradient(low = adjustcolor("#fcbad3", alpha.f = 0.5), high = "#c94c4c", guide = F) +
  labs(
    x = "NES",
    y = "Gene Set",
    size = "Gene Set Size") +
  scale_x_continuous(limits = c(2, 2.8), breaks = seq(2, 2.8, by = 0.2)) +
  theme_classic2(base_size = 12) +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10)
  )
positive_plot

# Create negative NES plot with "flipped" aesthetics
negative_plot <- ggplot(negative_paths, aes(x = NES, y = PathRank, size = size, color = SignedLogpadj, fill = signedlogp)) +
  geom_point(shape = 21, stroke = 1.5, alpha = 0.8) +
  scale_size_area(max_size = 15) +
  scale_color_gradient(high = "#add8e6", low = "#4169e1") +
  scale_fill_gradient(high = adjustcolor("#add8e6", alpha.f = 0.5), low = "#4169e1", guide = F) +
  labs(
    x = "NES",
    y = "Gene Set",
    size = "Gene Set Size",
    color = "-Log10 p-adj") +
  scale_x_continuous(position = "top", limits = c(-2.5, -1.8), breaks = seq(-2.5, -1.8, by = 0.2))  +  
  scale_y_discrete(position = "right") +  
  theme_classic2(base_size = 12) +
  theme(legend.position = "left",
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10))
negative_plot 

# Convert plots to grobs
g1 <- ggplotGrob(positive_plot)
g2 <- ggplotGrob(negative_plot)

plot_width = 0.8   
plot_height = 0.5

# For the positive plot (top), bottom-left aligns at the midpoint horizontally
positive_viewport <- viewport(
  x = 0.15,
  y = 0.5,
  width = plot_width,
  height = plot_height,
  just = c("left", "bottom")
)

# For the negative plot (bottom), top-right aligns at the midpoint horizontally
negative_viewport <- viewport(
  x = 1.1,
  y = 0.5,
  width = plot_width * 1.3,
  height = plot_height,
  just = c("right", "top")
)

# Open a PDF device
pdf(
  "./Output/fGSEA_UCLAAllPatch_deseq_recur_bp.pdf", 
  width = 5, 
  height = 6.5
)

grid.newpage() 

pushViewport(positive_viewport)
grid.draw(g1)
upViewport()

pushViewport(negative_viewport)
grid.draw(g2)
upViewport()  

dev.off()

#Fig S1E
Rubrary::plot_GSEA_barplot(
  gsea_res = "data/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC_fig.txt",
  n_pws = 10,
  NES_cutoff = 2,
  sig_cutoff = c("padj", 0.05),
  pw_size = 4.5,
  colors = c("firebrick", "darkblue"),
  title = "",
  savename = "Output/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC.png",
  width = 10,
  height = 7
)

gsea_pws <- c(
  gsea_res[order(gsea_res$NES, decreasing = T), ]$pathway[c(2,4:7,9:18)], 
gsea_res[order(gsea_res$NES, decreasing = F), ]$pathway[1:n_pws])
#
#GSEA-squared Fig 1B, S1B, S1C-------
##Dotplot (S1B & C)
#UCLA + Patch
categories = c("Immune Response", "Cell Differentiation", "Lipid Metabolism", "Wound Healing", "EMT", 
               "JAK/STAT & LIF", "Neuronal Function", "Epigenetic Regulation", "Cell Cycle Regulation", "DNA Damage Repair", "RB-E2F Dysregulation")
cat_terms = c("INFLAM|IMMUN|INTERLUEKIN|LEUKOCYTE|TNF|MHC|CYTOKINE_|CHEMOKINE|ANTIGEN|LYMPH",
              "DIFFERENTIATION",
              "_COA_|LIPID|STEROL|FATTY|FAT",
              "WOUND",
              "EPITHELIAL_MESENCHYMAL_TRANSITION|EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
              "JAK_STAT|LEUKEMIA_INHIBITORY_FACTOR|FGFR",
              "NEURON|NEUROTRANSMITTER|SYNAP|VOLTAGE|AXON|CEREBRAL|CORTEX|NEURAL",
              "HISTONE|METHYL|ACETYL|SUMOYLATION|UBIQUITIN",
              "CELL_CYCLE|MITOTIC|DNA_REPLICATION|CHROMOSOME_SEGREGATION|SPINDLE|CELL_DIVISION|REPLICATION|STRAND|G2", 
              "REPAIR|STRESS|HDR|HRR|DAMAGE|FANCONI",
              "RB_1|RB1|RBL1|RBL2|RETINOBLASTOMA|E2F"
)

gseasq.UCLAAllPatch_deseq_recur_full_supp = Rubrary::run_GSEA_squared(
  df_GSEA = "data/fGSEA_UCLAAllPatch_deseq_recur.txt",
  categories = categories,
  cat_terms = cat_terms,
  title = "Chemonaive v. Recurrent",
  savename = "/Output/fGSEA_UCLAAllPatch_deseq_recur_full_supp",
  height = 7, width = 7
)

#UCLA + Starr 
gseasq.UCLAAllPatch_deseq_postnact_full_supp = Rubrary::run_GSEA_squared(
  df_GSEA = "data/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC.txt",
  categories = categories,
  cat_terms = cat_terms,
  title = "Chemonaive v. Post-NACT",
  savename = "Output/fGSEA_UCLAAllPatch_deseq_postnact_full_supp",
  height = 7, width = 7
)

##Density Plot (1B)
#UCLA & Patch
gseasq.UCLAAllPatch_deseq_recur_full = Rubrary::run_GSEA_squared(
  df_GSEA = "data/fGSEA_UCLAAllPatch_deseq_recur.txt",
  categories = categories[c(1:3,8:11)],
  cat_terms = cat_terms[c(1:3,8:11)],
  savename = "./data/fGSEA_UCLAAllPatch_deseq_recur_full",
  height = 4.5, width = 6.5
)
write.table(
  x = gseasq.UCLAAllPatch_deseq_recur_full[["pathways"]],
  file = "./Output/gseasq.fGSEA_UCLAAllPatch_deseq_recur_full.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#UCLA & Starr
gseasq.UCLAStarr_deseq_postnact_full = Rubrary::run_GSEA_squared(
  df_GSEA = "data/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC.txt",
  categories = categories[c(1:3,8:11)],
  cat_terms = cat_terms[c(1:3,8:11)],
  savename = "./data/fGSEA_UCLAStarr_GSEA_DESeq_PostNACT_full",
  height = 4.5, width = 6.5
)

Rubrary::plot_GSEAsq_density(
  GSEAsq_df1 = "./data/fGSEA_UCLAAllPatch_deseq_recur_full_GSEAsq_pathways.txt",
  GSEAsq_df2 = "./data/fGSEA_UCLAStarr_GSEA_DESeq_PostNACT_full_GSEAsq_pathways.txt",
  name1 = "Recurrent v. Chemonaive",
  name2 = "Post-NACT v. Chemonaive",
  plot_pval = F,
  cat_order = rev(categories[c(1:3,8:11)]),
  title = " ",
  colors = c("firebrick", "lightslateblue"),
  savename = "./Output/UCLAAllPatchRecur_UCLAStarrPostNACT",
  plot_fmt = "pdf",
  height = 5, width = 5
)

#
#DEG Heatmap - Recurrent (Fig 1C)----
# Generate Heatmap of top DEGs & their broader GSEA Category - UCLA & Patch (Fig 1C)

# Define file paths
deg_file <- "data/UCLAAllPatch_topgenes_de.le_dupgenefiltered_upd.txt"
ucla_expr_file <- "data/ovarian_ucla1234_rsem_genes_upper_norm_counts_coding_log2_paired_recurrent.txt"
display_genes_file <- "data/UCLAAllPatch_topgenes_de.le_dupgenefiltered_display.txt"


# UCLA: Merge DE genes with expression data, order, and prepare annotations
heatmap_data_rc_full <- fread(deg_file) %>%
  inner_join(fread(ucla_expr_file), by = "gene") %>%
  arrange(category)

row_annotation <- heatmap_data_rc_full %>%
  select(gene, category) %>%
  column_to_rownames(var = "gene")

heatmap_data_rc <- heatmap_data_rc_full %>%
  column_to_rownames(var = "gene") %>%
  select(-category, -co_rank_category, -deseq_padj, -gsea_leadingedge_frequency)
colnames(heatmap_data_rc) <- ov_anno_rc$Sample[match(colnames(heatmap_data_rc), ov_anno_rc$Sample_ID)]

col_annotation_rc <- colnames(heatmap_data_rc) %>%
  data.frame(Sample = ., stringsAsFactors = FALSE) %>%
  left_join(ov_anno_rc, by = "Sample") %>%
  left_join(uc5p.rbl, by = "Sample_ID") %>%
  column_to_rownames(var = "Sample") %>%
  select(Treatment_Status, Sample_Type, RB_E2F_WP)

col_color_rc <- list(
  Treatment_Status = c(Chemonaive = "#f9ceab", Recurrent = "#ae7d78"),
  `Sample Type` = c(Fluid = "#F0FFF0", Tumor = "#1b7837", `Lymph node` = "darkseagreen1"),
  `RB1-E2F Dysregulation` = circlize::colorRamp2(c(min(col_annotation_rc$RB_E2F_WP), max(col_annotation_rc$RB_E2F_WP)), 
                                                 c("white", "olivedrab"))
)

# Patch: Prepare data matrix and annotations
patch_col_anno <- patch_anno_full %>%
  rename(Sample_ID = Short_name) %>%
  left_join(uc5p.rbl, by = "Sample_ID") %>%
  column_to_rownames(var = "Sample") %>%
  select(Treatment_Status, Sample_Type, RB_E2F_WP)

heatmap_data.p <- patch.r %>%
  select(match(patch_anno_full$Short_name, colnames(.))) %>%
  setNames(patch_anno_full$Sample[match(colnames(.), patch_anno_full$Short_name)]) %>%
  filter(row.names(.) %in% rownames(row_annotation)) %>%
  slice(match(rownames(row_annotation), row.names(.)))

col_color_p <- list(
  Treatment_Status = c(Chemonaive = "#f9ceab", Recurrent = "#ae7d78"),
  `Sample Type` = c(Fluid = "#F0FFF0", Tumor = "#1b7837"),
  `RB1-E2F Dysregulation` = circlize::colorRamp2(c(min(patch_col_anno$RB_E2F_WP), max(patch_col_anno$RB_E2F_WP)), 
                                                 c("white", "olivedrab"))
)

# Define order of row and column splits
row_split_order <- rev(c("Immune Response", "Cell Differentiation", "Lipid Metabolism", "Neuronal Function", 
                         "Epigenetic Regulation", "Cell Cycle Regulation", "DNA Repair"))
col_split_order <- c("Chemonaive", "Recurrent")

# Read and filter genes for display on heatmap
f1hm_dispgs <- fread(display_genes_file) %>%
  mutate(mat_index = match(gene, rownames(heatmap_data.p))) %>%
  select(mat_index, gene) %>%
  na.omit()

# Create custom legend for row annotations
legend_annoblock_cat <- Legend(
  title = "Category",
  labels = row_split_order,
  legend_gp = gpar(fill = 2:8),
  title_gp = gpar(fontsize = 7.5, fontface = "bold"),
  nrow = 4, labels_gp = gpar(fontsize = 7.5)
)

# Generate UCLA and Patch heatmaps
cht.urc <- generate_deg_heatmap(
  heatmap_data = heatmap_data_rc,
  col_annotation = col_annotation_rc[, -which(names(col_annotation_rc) == "Treatment_Status")],
  col_color = col_color_rc,
  col_split = factor(col_annotation_rc$Treatment_Status, levels = col_split_order),
  col_split_order = col_split_order,
  column_title = "UCLA",
  colnames_fontsize = 5,
  row_annotation = row_annotation,
  row_split = factor(row_annotation$category, levels = row_split_order),
  row_split_labels = c("Repair", "Cycle", "Epigenetic", "Neuro", "Lipid", "Diff.", "Immune"),
  cluster_row_slices = FALSE,
  legend_title_fontsize = 7.5,
  legend_label_fontsize = 7.5,
  show_annotation_name = FALSE
)

cht.p <- generate_deg_heatmap(
  heatmap_data = heatmap_data.p,
  row_annotation = row_annotation,
  col_annotation = patch_col_anno[, -which(names(patch_col_anno) == "Treatment_Status")],
  col_color = col_color_p,
  col_split = factor(patch_col_anno$Treatment_Status, levels = col_split_order),
  col_split_order = col_split_order,
  column_title = "Patch et al.",
  colnames_fontsize = 5,
  show_row_names = FALSE,
  legend_title_fontsize = 7.5,
  legend_label_fontsize = 7.5,
  show_annotation_name = FALSE,
  display_gene_index = f1hm_dispgs$mat_index,
  display_gene_label = f1hm_dispgs$gene
)

# Save the combined heatmaps
pdf(
  "Output/ucp_deg_hm.pdf",
  height = 6.5,
  width = 8
)
draw(
  cht.urc + cht.p, 
  ht_gap = unit(0.2, "cm"), 
  annotation_legend_list = list(Category = legend_annoblock_cat), 
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()
#
#DEG Heatmap - PostNACT (Figure S1D)----
deg.le.ucs = fread("data/UCLAStarr_topgenes_de.le_dupgenefiltered.txt")
ucla.ov.paired.pn <- fread("data/ovarian_ucla12345_rsem_genes_upper_norm_counts_coding_log2_post_nact.txt")

heatmap_data_full.ucs <- merge(deg.le.ucs[, c("gene","category")], ucla.ov.paired.pn, by = "gene")

# Order genes by category for the heatmap, Remove rownames & set the col "gene" as rowname, Extract gene expression values
heatmap_data_full.ucs <- heatmap_data_full.ucs[order(heatmap_data_full.ucs$category),]
heatmap_data.ucs <- heatmap_data_full.ucs %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "gene") %>%
  select(-category)

row_annotation.ucs <- heatmap_data_full.ucs[, c(1,2)] %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "gene")

#Match anno names with those on heatmap matrix
colnames(heatmap_data.ucs) <- ov_anno_pn$Sample[match(colnames(heatmap_data.ucs), ov_anno_pn$Sample_ID)]

#UCLA Chemonaive
heatmap_data.pn <- heatmap_data.ucs %>% 
  select(all_of(
    ov_anno_pn %>%
      filter(Treatment_Status %in% c("Chemonaive", "Post-NACT") & Paired == T) %>%
      pull(Sample))) 

col_annotation_pn <- colnames(heatmap_data.pn) %>% 
  data.frame(Sample = ., stringsAsFactors = F) %>%
  left_join(ov_anno_pn, by = "Sample") %>%
  column_to_rownames(var = "Sample") %>%
  select(Treatment_Status, Sample_Type, RB_Loss_Malorni)

col_color.pn = list(
  Treatment_Status = c(Chemonaive = "#f9ceab", Post_NACT = "#f9b3ab"),
  Sample_Type = c(Fluid = "#F0FFF0", Tumor = "#1b7837",`Lymph node` ="darkseagreen1"),
  RB_Loss_Malorni =  circlize::colorRamp2(
    c(min(col_annotation_pn$RB_Loss_Malorni, na.rm = T), max(col_annotation_pn$RB_Loss_Malorni, na.rm = T)),
    c("white", "olivedrab")))

row_split_order.ucs = rev(c("Immune Response", "Cell Differentiation", "Lipid Metabolism", "Epigenetic Regulation", "Cell Cycle Regulation", "DNA Repair"))
col_split_order.ucs = c("Chemonaive", "Post-NACT")
cht.upn = generate_deg_heatmap(
  heatmap_data = heatmap_data.pn,
  col_annotation = col_annotation_pn[, -which(names(col_annotation_pn) == "Treatment_Status")],
  col_color = col_color.pn,
  col_split = factor(col_annotation_pn$Treatment_Status, levels = col_split_order.ucs),
  col_split_order = col_split_order.ucs,
  column_title = "UCLA",
  row_annotation = row_annotation,
  row_split = factor(row_annotation.ucs$category, levels = row_split_order.ucs),
  row_split_labels = c("Repair", "Cycle", "Epigenetic", "Lipid", "Diff.", "Immune"),
  cluster_row_slices = F,
  show_annotation_name = F
)
cht.upn

#Starr
starr_col_anno <- starr_anno_full %>%
  select(Treatment_Status, Sample_Type, RB1_Loss)

heatmap_data.s <- starr.r %>%
  select(match(starr_anno_full$Sample_ID, colnames(.))) %>%
  `colnames<-`(starr_anno_full$Sample[match(colnames(.), starr_anno_full$Sample_ID)]) %>%
  filter(row.names(.) %in% rownames(row_annotation.ucs)) %>%
  slice(match(rownames(row_annotation.ucs), rownames(.)))

col_color.s = list(
  Treatment_Status = c(Chemonaive = "#f9ceab", Post_NACT = "#f9b3ab"),
  Sample_Type = c(Fluid = "#F0FFF0", Tumor = "#1b7837"),
  RB1_Loss = circlize::colorRamp2(
    c(min(starr_col_anno$RB1_Loss), max(starr_col_anno$RB1_Loss)), 
    c("white", "olivedrab"))
)
cht.s = generate_deg_heatmap(
  heatmap_data = heatmap_data.s,
  row_annotation = row_annotation.ucs,
  row_names_fontsize = 3.5,
  col_annotation = starr_col_anno[, -which(names(starr_col_anno) == "Treatment_Status")],
  col_color = col_color.s,
  col_split = factor(starr_col_anno$Treatment_Status, levels = col_split_order.ucs),
  col_split_order = col_split_order.ucs,
  column_title = "Starr"
)
cht.s

cht.upn + cht.s
pdf(
  "~/Dropbox/Graeber/OVP/OV_NE/Figures/Output/Fig1/ucs_deg_hm.pdf", 
  height = 6.5,
  width = 9
)
draw(cht.upn + cht.s, ht_gap = unit(0.2, "cm"))
dev.off()


#
