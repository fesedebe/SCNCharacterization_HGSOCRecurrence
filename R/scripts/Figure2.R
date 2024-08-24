
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

source("/RFunctions/FGSEA_functions.R")
source("R/functions/Fig1_functions.R")

scn.sig <- read.delim("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/SCN_Signature/SCN_Signature.txt") %>%
  arrange(desc(PC1.v)) 
#
#Enrichment of Resistance Signatures in HGSOC Data-----
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
  deg.df = fread("data/OV_UCLA&Patch_recur_prevpost_slogp_DESeq_pconly_S13b.txt"), 
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
#boxplots of SCN & RB1 Loss Scores----
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
#Patients with all three treatment timepoints-----
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
#Heatmap of top 1000 NE Genes - Combine UCLA & SCN data & and anno-----
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

pdf(
  "/Output/top1kscn_lg.pdf",
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
