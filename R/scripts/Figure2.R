
#Load library, function, data-----
library(dplyr)
library(ggpubr)
library(Rubrary)
library(cowplot)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(patchwork)
library(tibble)
library(data.table)

source("/RFunctions/FGSEA_functions.R")
source("R/functions/bulkRNAseq_functions.R")

ucp.deseq.recur.allsamps = read.delim(
  file = "data/OV_UCLA&Patch_recur_prevpost_slogp_DESeq_pconly_S13b.txt",
  stringsAsFactors = F
)
ucla.ov.paired.rc <- read.delim(
  file = "data/ovarian_ucla1234_rsem_genes_upper_norm_counts_coding_log2_paired_recurrent.txt",
  stringsAsFactors = F
)
scn.sig <- read.delim("data/SCN_Signature.txt") %>%
  arrange(desc(PC1.v)) 
scn.f1a = read.delim(
  file = "~/Dropbox/Graeber/OVP/Ovarian Project/log2_coding_expression_datasets/PRAD.norm_Beltran_LUAD.norm_SCLC_LUAD.subset_rsem_genes_upper_norm_counts_coding_log2.txt",
  stringsAsFactors = F
)
human.info.scn.anno <- read.delim(
  "data/human.info.scn.anno.txt",
  row.names = 1,
  stringsAsFactors = F
)
ov_anno_rc = read.delim("data/UCLAOvarian_CombinedBatches_Annotation_v13_rcp.txt", stringsAsFactors = F)
disp_gs_hm = data.table::fread("data/scn_genes_hm20.txt")
paletteTS = toupper(c(
  "#f9ceab", 
  "#f9b3ab", 
  "#ae7d78", 
  "#f9abd8", 
  "#ae7897" 
))
#
#Enrichment of Resistance Signatures in HGSOC Data (2A & S2A)-----
# geneset prep
resis.genesets <- list(
  "SCN_Balanis" = head(scn.sig$Gene, 500),
  "RBLoss_Chen" = read.delim("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/RB/RB1Loss_Signature_Chen.txt", stringsAsFactors = F)[,1],
  "RBLoss_Malorni" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/RB/Malorni_RBsig_genes.csv")[,1],
  "NE_Alshalalfa" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/NE/AlshalalfaNEsig.csv", header = F)[,1],
  "NE_Zhang" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/NE/Minna_NE.txt")[,1],
  "ASC_Smith" = read.delim("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/StemCell/ASC_Signature/adult_stem_cell_sig.txt")[,1],
  "StemCell_Palmer" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/StemCell/Palmer_StemCell.csv")[,1],
  "ImmuneEvasion_Zhang" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/Immune/ZhangImmuneEvasion.csv")[,1],
  "ImmuneRepressed_Jerby-Arnon" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/Immune/arnon_immune_repressed.csv")[,1],
  "ImmuneInduced_Jerby-Arnon" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/Immune/arnon_immune_induced.csv")[,1],
  "EMT_Hallmark" = read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/EMT/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt")[,1]
)
#debugonce(run_fgsea)
resis_fgsea_rc = run_fgsea(
  pathways = resis.genesets, 
  deg.df = data.table::fread("data/OV_UCLA&Patch_recur_prevpost_slogp_DESeq_pconly_S13b.txt"), 
  deg.df_slpval = "sign_log_p", 
  maxSize = 500, 
  nPermSimple = 10000,
  title = "",
  save = F
)
ggsave(
  "Output/resis_sigs_fgsea_bp.pdf", 
  resis_fgsea_rc[["plot"]], 
  device = "pdf",
  #dpi = 600, 
  width = 7.5, #6.3
  height = 5.8
)

resis_fgsea_pn = run_fgsea(
  pathways = resis.genesets, 
  deg.df = fread("data/UCLAStarr_DESeq_PostNACT_pconly.txt"), 
  deg.df_slpval = "sign_log_p", 
  maxSize = 500, 
  nPermSimple = 10000,
  title = "Resistance Signatures Enrichment: Post-NACT v. Chemonaive",
  save  = F
)
ggsave(
  "Output/resis_sigs_fgsea_pn.pdf", 
  resis_fgsea_pn[["plot"]], 
  device = "pdf", 
  width = 7.5, 
  height = 5.8
)
#
#boxplots of SCN & RB1 Loss Scores (2B, C)----
#combine UCLA & Patch Annotations
ovp_anno = rbind(
  ov_anno_rc[,c("Sample","Treatment_Status","Dataset", "SCN_Score", "RB_Loss_Malorni", "RB_Loss_Chen")], #RB1_Loss
  patch_anno_full[,c("Sample","Treatment_Status","Dataset", "SCN_Score", "RB_Loss_Malorni", "RB_Loss_Chen")]) %>%
  arrange(Sample) %>%
  arrange(Treatment_Status) %>%
  mutate(Dataset = factor(Dataset, levels = c("UCLA", "Patch")))

#combine UCLA, Starr; JS Annotations (JS values are tiny so plotted separately)
ovs_anno = rbind(
  ov_anno_pn[,c("Sample","Treatment_Status","Dataset", "SCN_Score", "RB_Loss_Malorni")], 
  starr_anno_full %>%
    left_join(s.rbl, by = "Sample_ID") %>%
    mutate(Sample = Sample_ID) %>%
    select(Sample, Treatment_Status, Dataset, SCN_Score, RB_Loss_Malorni)
) %>%
  arrange(Sample) %>%
  arrange(Treatment_Status) %>%
  mutate(Dataset = factor(Dataset, levels = c("UCLA", "Starr")))

js_anno_bp = js_anno %>%
  left_join(js.rbl, by = join_by(Sample == Sample_ID)) %>%
  select(Sample, Treatment_Status, Dataset, SCN_Score, RB_Loss_Malorni) %>%
  arrange(Sample) %>%
  arrange(Treatment_Status)

##SCN
bp.ucp.scn <- create_paired_boxplot(
  data = ovp_anno,
  x = "Treatment_Status", 
  y = "SCN_Score",
  ylab = "SCN Score",
  xlab = "Treatment Status",
  title = "Recurrent Pairs",
  fill = "Treatment_Status",
  legend.position = "none",
  #legend.justification = "left", 
  palette = paletteTS[c(1,3)],
  y_limits = c(-60, 75),
  y_breaks = seq(-60, 75, by = 30),
  p_method = "combined",
  p_anno_n = 1,
  p_anno_x = 1.2
)
print(bp.ucp.scn)

bp.ucs.scn = create_paired_boxplot(
  data = ovs_anno,
  x = "Treatment_Status", 
  y = "SCN_Score",
  ylab = "SCN Score",
  xlab = "Treatment Status",
  title = "Post-NACT Pairs",
  fill = "Treatment_Status",
  legend.position = "none",
  #legend.justification = "left",
  palette = paletteTS[c(1,2)],
  y_limits = c(-60, 75),
  y_breaks = seq(-60, 75, by = 30),
  p_method = "combined",
  p_anno_lab = factor(c("UCLA", "Starr"), levels = c("UCLA", "Starr")),
  p_anno_x = 1.2
)
print(bp.ucs.scn)

bp.js.scn = create_paired_boxplot(
  data = js_anno_bp,
  x = "Treatment_Status", 
  y = "SCN_Score",
  ylab = "SCN Score",
  xlab = "Treatment Status",
  title = "",
  fill = "Treatment_Status",
  legend.position = "none",
  palette = paletteTS[c(1,2)],
  y_limits = c(30, 52),
  y_breaks = seq(30, 52, by = 5),
  p_method = "wilcox.test",
  p_label = "p.format",
  p_label_x = 0.6,
  p_label_y = 51
)
print(bp.js.scn)

#combine scn plots
bp.rpn = bp.ucp.scn + 
  (bp.ucs.scn + bp.js.scn + patchwork::plot_layout(
    axis_titles = 'collect_y',
    widths = c(2,1)
  )) +
  patchwork::plot_layout(
    widths = c(2,3))
print(bp.rpn)

ggsave(
  "Output/ucpsj_delta.scn_bp2.png", 
  bp.rpn, device = "png", width = 15, height = 3.5, dpi = 600)

##RB1 Loss
bp.ucp.rb <- create_paired_boxplot(
  data = ovp_anno,
  x = "Treatment_Status", 
  y = "RB_Loss_Malorni",
  xlab = "Treatment Status",
  ylab = "RB Loss Score",
  title = "",
  fill = "Treatment_Status",
  legend.position = "none",
  palette = paletteTS[c(1,3)],
  y_limits = c(-0.45, 0.5),
  y_breaks = seq(-0.40, 0.5, by = 0.2),
  p_method = "combined",
  p_anno_n = 1,
  p_anno_y = 0.5,
  p_anno_x = 1.2
)
print(bp.ucp.rb)

bp.ucs.rb = create_paired_boxplot(
  data = ovs_anno,
  x = "Treatment_Status", 
  y = "RB_Loss_Malorni",
  xlab = "Treatment Status",
  ylab = "RB Loss Score",
  title = "",
  fill = "Treatment_Status",
  legend.position = "none",
  #legend.justification = "left",
  palette = paletteTS[c(1,2)],
  y_limits = c(-0.45, 0.5),
  y_breaks = seq(-0.40, 0.5, by = 0.2),
  p_method = "combined",
  p_anno_y = 0.5,
  p_anno_x = 1.2
)
print(bp.ucs.rb)

bp.js.rb = create_paired_boxplot(
  data = js_anno_bp,
  x = "Treatment_Status", 
  y = "RB_Loss_Malorni",
  xlab = "Treatment Status",
  ylab = "RB Loss Score",
  title = "",
  fill = "Treatment_Status",
  legend.position = "none",
  palette = paletteTS[c(1,2)],
  y_limits = c(-0.45, 0.5),
  y_breaks = seq(-0.40, 0.5, by = 0.2),
  p_method = "wilcox.test",
  p_label = "p.format",
  p_label_x = 0.6,
  p_label_y = 0.45
)
print(bp.js.rb)

#combine rbl plots
bp.rpn.rb = bp.ucp.rb + 
  (bp.ucs.rb + bp.js.rb + patchwork::plot_layout(
    axis_titles = 'collect_y',
    widths = c(2,1)
  )) +
  patchwork::plot_layout(
    widths = c(2,3))
print(bp.rpn.rb)

ggsave(
  "Output/ucpsj_delta.rb_bp2.png", 
  bp.rpn.rb, device = "png", width = 15, height = 3.5, dpi = 600)

#
#Patients with all three treatment timepoints (S2D)-----
tt3 <- ov_anno_full %>%
  filter(Cohort %in% c("Post-NACT & Recurrent"), Paired == TRUE) %>%
  select(Sample) %>%
  mutate(Sample = substr(Sample, 1, nchar(Sample) - 4)) %>%
  pull(Sample)

ov_anno_cpr <- ov_anno_full %>%
  mutate(Patient = substr(Sample, 1, nchar(Sample) - 4)) %>%
  filter(Patient %in% tt3) %>%
  arrange(Patient, Treatment_Status) %>%
  group_by(Patient) %>%
  mutate(
    Change_in_SCN_Score = SCN_Score - first(SCN_Score),
    Change_in_RB1_Loss_Score = RB1_Loss - first(RB1_Loss)
  )

# SCN
tt3.s <- ggplot(
  ov_anno_cpr, 
  aes(x = Treatment_Status, 
      y = Change_in_SCN_Score, 
      color = Patient, 
      group = Patient)) +
  geom_line() +
  geom_point() +
  theme_light() +
  labs(x = "", y = "Change in SCN Score") +
  ggsci::scale_color_jco()
tt3.s

# RB loss
tt3.r <- ggplot(
  ov_anno_cpr, 
  aes(x = Treatment_Status, 
      y = Change_in_RB1_Loss_Score, 
      color = Patient, 
      group = Patient)) +
  geom_line() +
  geom_point() +
  theme_light() +
  labs(x = "Treatment Status", y = "Change in RB1 Loss Score") +
  ggsci::scale_color_jco(guide = FALSE) 
tt3.r

# Combine plots, add common y-axis title + central legend
combined_plot <- (tt3.s / tt3.r) +
  plot_layout(guides = 'collect', axes = "collect_x") + 
  plot_annotation(title = "Patients with 3 Treatment Timepoints") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

combined_plot

ggsave(
  filename = "Output/delta.sig.score_3tt.png",
  plot = combined_plot,
  device = "png",
  height = 5,
  width = 5
)
#Dotplot - Relevant Genes Connected to SCN (2D)-----
ucp.deseq.recur.allsamps_hmg <- ucp.deseq.recur.allsamps %>%
  mutate(signed_logpadj = -log10(padj) * sign(log2FoldChange)) %>%
  mutate(percentile_rank = percent_rank(signed_logpadj)) %>%
  filter(gene %in% data.table::fread("data/hgsoc-scn_hmgenes.txt")$gene)

#pos
ucp.deseq.recur.allsamps_hmgpos <- ucp.deseq.recur.allsamps_hmg %>%
  filter(log2FoldChange > 0 | gene %in% c("SYP", "FOXA2", "NCAM1"))

dot_plotA <- ggplot(ucp.deseq.recur.allsamps_hmgpos, aes(x = gene, y = percentile_rank)) +
  geom_point(aes(size = abs(log2FoldChange), color = signed_logpadj), alpha = 0.8) +  
  scale_color_gradientn(
    colors = c("#6A0DAD", "#D8BFD8", "#B03A5B"), 
    values = scales::rescale(c(-5, -1.3, 1.3, 5)),
    limits = c(-5, 5),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(3, 8)) + 
  theme_minimal() +
  labs(title = "Expression Change in Recurrent Pairs",
       x = "Up in SCN",
       y = "Percentile Rank in HGSOC", 
       color = "Significance \n", 
       size = "Log Fold \nChange"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5), 
        axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 10, face = "bold"),
        legend.direction = "horizontal",
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 8.5),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "lightgray", linewidth = 0.1)
  ) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, by = 0.25)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "lightgray", linewidth = 0.35) 
dot_plotA

#neg
ucp.deseq.recur.allsamps_hmgneg <- ucp.deseq.recur.allsamps_hmg %>%
  filter(log2FoldChange < 0 & !gene %in% c("SYP", "FOXA2", "NCAM1"))
dot_plotB <- ggplot(ucp.deseq.recur.allsamps_hmgneg, aes(x = gene, y = percentile_rank)) +
  geom_point(aes(size = abs(log2FoldChange), color = signed_logpadj), alpha = 0.8) +  
  scale_color_gradientn(
    colors = c("#6A0DAD", "#D8BFD8", "#B03A5B"), 
    values = scales::rescale(c(-5, -1.3, 1.3, 5)),
    limits = c(-5, 5),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(3, 8)) + 
  theme_minimal() +
  labs(title = " ",
       x = "Down in SCN",
       y = "Percentile Rank in HGSOC",
       color = "Significance \n",
       size = "Log Fold \nChange"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5), 
        axis.text.y = element_text(size = 7.5),
        plot.title = element_text(size = 10, face = "bold"),
        legend.direction = "horizontal",
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 8.5),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "lightgray", linewidth = 0.1)
  ) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, by = 0.25)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "lightgray", linewidth = 0.35) 
dot_plotB

dpab = dot_plotA + dot_plotB + patchwork::plot_layout(
  axis_titles = 'collect_y',
  axes = 'collect_y',
  guides = 'collect', 
  widths = c(3,1)
)
print(dpab)

ggsave(
  "Output/scn_relevant_genes_dotplot.pdf", 
  dpab, 
  device = "pdf", 
  width = 12, 
  height = 3
)
#
#Heatmap of top 1000 SCN Genes - Combine UCLA & SCN data & and anno (2E)-----

#filter to top 1000 SCN genes
scn.sig.filt = head(scn.sig$Gene, 1000)

#filter scn.f1a data to filtered scn genes & scn samples
scn.f1a.filt <- scn.f1a[scn.f1a$gene %in% scn.sig.filt,] %>%
  `row.names<-`(., NULL) %>% 
  tibble::column_to_rownames(var = "gene") %>%
  select(all_of(rownames(human.info.scn.anno))) %>%
  as.matrix()

##ucla
#filter ov samps to just recurrent
ov_anno_rc.filt = ov_anno_rc[ov_anno_rc$Treatment_Status %in% "Recurrent",]  %>%
  left_join(
    data.table::fread("data/ovarian_ucla12345_HGSOCSubtypes_consensusOV.txt") %>%
      dplyr::select("V1", "Subtype"), 
    by = c("Sample_ID" = "V1")
  ) %>%
  `row.names<-`(., NULL) %>% 
  tibble::column_to_rownames(var = "Sample_ID") %>%
  dplyr::select(Treatment_Status, SCN_Score, RB_Loss_Malorni, Subtype) %>%
  dplyr::rename(RB1_Loss_Score = RB_Loss_Malorni)

#filter ucla.ov.paired.rc to top & bottom scn genes, and samples in ov_anno_rc.filt
ucla.ov.paired.rc.filt = ucla.ov.paired.rc[ucla.ov.paired.rc$gene %in% scn.sig.filt,] %>%
  `row.names<-`(., NULL) %>% 
  tibble::column_to_rownames(var = "gene") %>%
  select(all_of(rownames(ov_anno_rc.filt))) %>%
  as.matrix()

#merging scn & ucla data & anno
scn.ucla.paired.rc = merge(
  scn.f1a.filt, 
  ucla.ov.paired.rc.filt, 
  "row.names") %>%
  column_to_rownames(var = "Row.names") %>%
  as.matrix()

scn.uc.anno <- human.info.scn.anno %>%
  bind_rows(ov_anno_rc.filt)

scn.uc.col = list(
  Treatment_Status = c("#ae7d78", "#f9abd8", "#8f7dae"),
  Subtype = c(
    "None" = "white",
    "SCLC-A" = "#D55E00",          
    "SCLC-P" = "#E69F00",          
    "SCLC-N" = "#009E73",          
    "SCLC-Y" = "#0072B2",         
    "Differentiated" = "#F0E442", 
    "Immunoreactive" = "#56B4E9",  
    "Mesenchymal" = "#CC79A7",     
    "Proliferative" = "#000000"    
  ),
  `SCN Score` =  circlize::colorRamp2(
    c(min(scn.uc.anno$SCN_Score), max(scn.uc.anno$SCN_Score)), 
    c("white", "#4d004b")
  ),
  `RB1 Loss Score` = circlize::colorRamp2(
    c(min(scn.uc.anno$RB1_Loss_Score), max(scn.uc.anno$RB1_Loss_Score)), 
    c("white", "olivedrab")
  )
)

# Manually create legend for the Treatment_Status anno_block
treatment_status_legend <- Legend(
  title = "Treatment Status",
  labels = c("Recurrent", "SCLC", "NEPC"),
  legend_gp = gpar(fill = scn.uc.col$Treatment_Status),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 12)
)

# Read and filter genes for display on heatmap
f2hm_dispgs <- disp_gs_hm %>%
  mutate(mat_index = match(gene, rownames(scn.ucla.paired.rc))) %>%
  select(mat_index, gene) %>%
  na.omit()

# Generate heatmap
scn.ovr.chm = Heatmap(
  matrix = scn.ucla.paired.rc,
  name = "gex",
  border = T,
  col = viridis::magma(500)[50:450],
  top_annotation = HeatmapAnnotation(
    Treatment_Status = anno_block(
      labels = c("Recurrent", "SCLC", "NEPC"),
      labels_gp = gpar(fontsize = 11, col = "white", fontface = "bold"),
      gp = gpar(fill = scn.uc.col$Treatment_Status)
    ),
    `SCN Score` = scn.uc.anno$SCN_Score,
    `RB1 Loss Score` = scn.uc.anno$RB1_Loss_Score,
    Subtype = scn.uc.anno$Subtype,
    which = "col", 
    border = T,
    col = scn.uc.col,
    show_annotation_name = F,
    annotation_legend_param = list(
      legend_direction = "horizontal",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12), nrow = 5
    )
  ),
  clustering_method_rows = "ward.D2",
  column_split = factor(scn.uc.anno$Treatment_Status, levels = c("Recurrent", "SCLC", "NEPC")),
  cluster_column_slices = F,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  show_row_names = F,
  show_column_names = F,
  column_names_rot = 45,
  row_title = "Top 1000 SCN Genes",
  column_title = "Patient Samples",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "Gene Expression (log2 UQ)", 
    color_bar = "continuous",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12))
)
scn.ovr.chm

pdf(
  "Output/top1kscn_lg.pdf",
  width = 10, 
  height = 4.5
)
draw(
  scn.ovr.chm,
  annotation_legend_list = list(Treatment_Status = treatment_status_legend),
  merge_legend = T,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()
#
#Correlation Matrix Heatmap (Fig S3E)-----
scn.uc.col2 = list(
  `Sample Type` = c(
    "Recurrent" = "#ae7d78", 
    "SCLC" = "#f9abd8", 
    "NEPC" = "#8f7dae"), 
  Subtype = c(
    "None" = "white",
    "SCLC-A" = "#D55E00",          # Strong Orange
    "SCLC-P" = "#E69F00",          # Golden Yellow
    "SCLC-N" = "#009E73",          # Bluish Green
    "SCLC-Y" = "#0072B2",          # Sky Blue
    "Differentiated" = "#000000",  # Yellow
    "Immunoreactive" = "#000000",  # Blue
    "Mesenchymal" = "#000000",     # Red-Purple
    "Proliferative" = "#000000"    # Black
  ),
  `SCN Score` =  circlize::colorRamp2(
    c(min(scn.uc.anno$SCN_Score), max(scn.uc.anno$SCN_Score)), 
    c("white", "#4d004b")
  ),
  `RB1 Loss Score` = circlize::colorRamp2(
    c(min(scn.uc.anno$RB1_Loss_Score), max(scn.uc.anno$RB1_Loss_Score)), 
    c("white", "olivedrab")
  )
)

scn.ovr.chm_cor = ComplexHeatmap::Heatmap(
  matrix = cor(scn.ucla.paired.rc),
  name = "gex",
  border = T,
  col = viridis::magma(500)[50:450],
  top_annotation = HeatmapAnnotation(
    `Sample Type` = scn.uc.anno$Treatment_Status,
    Subtype = scn.uc.anno$Subtype,
    which = "col", 
    border = T,
    col = scn.uc.col2,
    show_annotation_name = T,
    annotation_legend_param = list(
      legend_direction = "horizontal",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 12), nrow = 5
    )
  ),
  right_annotation = rowAnnotation(
    `Sample Type` = scn.uc.anno$Treatment_Status,
    Subtype = scn.uc.anno$Subtype,
    border = T,
    col = scn.uc.col2,
    show_annotation_name = F,
    show_legend = F
  ),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  show_row_names = F,
  show_column_names = F,
  column_names_rot = 45,
  row_title = "Top 1000 SCN Genes",
  column_title = "Patient Samples",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "Correlation", 
    color_bar = "continuous",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12))
)
scn.ovr.chm_cor
pdf(
  "Output/top1kscn_lg_cormat.pdf",
  width = 8.5, 
  height = 8.5
)
draw(
  scn.ovr.chm_cor,
  annotation_legend_list = list(Treatment_Status = Legend(
    title = "Sample Type",
    labels = c("Recurrent", "SCLC", "NEPC"),
    legend_gp = gpar(fill = scn.uc.col$Treatment_Status),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12)
  )),
  merge_legend = T,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()
#
