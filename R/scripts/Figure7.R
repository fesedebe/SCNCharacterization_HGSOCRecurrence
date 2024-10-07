
library(dplyr)
library(ggplot2)
library(gtable)
source("R/functions/Fig1_functions.R")

#Tables (7A & E)-----
library(gt)

#short term model
st_tab <- data.frame(
  TreatmentStatus = c(
    "PDX8 (Before Injection)",
    "PDX8 In Vivo Pre Therapy",
    "PDX8 In Vivo Post Therapy (1st cycle)",
    "PDX8 In Vivo Post Therapy (2nd Cycle)",
    "PDX8 In Vivo Untreated"
  ),
  ID = c(
    "T0",
    "T1a, T1b, T1c",
    "T2a, T2b, T2c",
    "T3a, T3b, T3c",
    "T4a, T4b, T4c"
  )
)

# Create gt table
gt_table <- gt(st_tab) %>%
  tab_header(
    title = "Short Term Carbocycling"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_title(groups = "title")
  ) %>%
  cols_label(
    TreatmentStatus = "Treatment Status",
    ID = "Sample ID"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = c(TreatmentStatus, ID))
  ) %>%
  tab_options(
    table.font.size = 20,
    heading.title.font.size = 21,
    heading.align = "center",
    table.width = pct(100),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    table_body.border.top.color = "black",
    table.border.left.color = "black",
    table.border.right.color = "black"
    #table_body.border.color = "black"
  ) 

gt_table

gtsave(gt_table, "Output/short_term_carbocycling_table.png")
gtsave(gt_table, "/Output/short_term_carbocycling_table.pdf")

#long term model
lt_tab <- data.frame(
  TreatmentStatus = c(
    "PDX8 In Vivo Pre Therapy",
    "PDX8 In Vivo Post Therapy (2nd Cycle)",
    "PDX8 In Vivo Tumor Relapse"
  ),
  ID = c(
    "N1a, N1b, N1c",
    "N2a, N2b, N2c",
    "N3a, N3b, N3c"
  )
)

lt_gt <- gt(lt_tab) %>%
  tab_header(
    title = "Long Term Carbocycling"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_title(groups = "title")
  ) %>%
  cols_label(
    TreatmentStatus = "Treatment Status",
    ID = "Sample ID"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = c(TreatmentStatus, ID))
  ) %>%
  tab_options(
    table.font.size = 20,
    heading.title.font.size = 21,
    heading.align = "center",
    table.width = pct(100),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black",
    table_body.border.top.color = "black",
    table.border.left.color = "black",
    table.border.right.color = "black"
    #table_body.border.color = "black"
  ) 

lt_gt
gtsave(lt_gt, "/Output/long_term_carbocycling_table.pdf")
#
#Short Term Model (7C & D): SCN & RB Scores-----
PDXB2_anno <- readxl::read_xlsx("data/OV_PDXB2_CarboCyc.xlsx")
PDXB2_log2UQ <- Rubrary::rread("data/OV_PDXB2_CarboCyc_rsem_genes_upper_norm_counts_coding_log2.txt", row.names = 1) %>%
  select(PDXB2_anno$ID)

#RB1 Loss Score
rb_lst_ssgsea <- list(
  "RB_Loss_Chen" = data.table::fread('signatures/RB/RB1Loss_Signature_Chen.txt')$`Highly_expressed_in_tumors_with 2+_RB1_mutations`,
  "RB_Loss_Malorni" = read.csv('signatures/RB/Malorni_RBsig_genes.csv')[,1]
)

#normalized ssgsea scores
pdxst.rbl = calc_ssgsea(
  exp.df = PDXB2_log2UQ,
  gset = rb_lst_ssgsea,
  col2 = "RB_Loss_Chen"
) 

#merge annos and scores
pdxst_rblm = pdxst.rbl %>%
  merge(PDXB2_anno, by.x = 'Sample_ID', by.y = "ID")
pdxst_rblm$Stage <- factor(
  pdxst_rblm$Stage,
  levels = c("Before injection", "Pre-therapy", "1st recurrence", "2nd recurrence", "Vehicle-treated"))
pos <- position_jitter(width = 0.25, seed = 13)

#RB Loss violin plot
st_rb = ggplot(pdxst_rblm, aes(x = Stage, y = RB_Loss_Malorni, fill = Stage)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.1) +
  xlab("Treatment Status") +
  ylab("Scores") +
  ggrepel::geom_label_repel(aes(label = Sample_ID), color = "black", position = pos, size = 2) +
  labs(title = "RB Loss", 
       fill = "Treatment Status")   +
  theme_classic() +
  ggpubr::stat_compare_means(
    method = "anova",
    label.y = 0.8, label.x = 0.8,
    label = "p.format"
  ) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(margin = margin(l = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold", family = "Arial"),
    plot.tag.location = "margin"
  ) 
st_rb

ggsave(
  "Output/PDXB2_ST_RB1_violin.png",
  plot = last_plot(),
  width = 6.5,
  height = 5,
  dpi = 600
)

#SCN Score violin plot
st_scn = ggplot(pdxst_rblm, aes(x = Stage, y = SCN_PCAProjPC1, fill = Stage)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.1) +
  xlab("Treatment Status") +
  ylab("Scores") +
  ggrepel::geom_label_repel(aes(label = Sample_ID), color = "black", position = pos, size = 2) +
  labs(title = "SCN", 
       fill = "Treatment Status")   +
  theme_classic() +
  ggpubr::stat_compare_means(
    method = "anova",
    label.y = 80, label.x = 0.8,
    label = "p.format"
  ) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(), 
    axis.title.y = element_text(margin = margin(l = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold", family = "Arial")
  )
st_scn

ggsave(
  "Output/PDXB2_ST_SCN_violin.png",
  plot = st_scn,
  width = 5.5,
  height = 5,
  dpi = 600
)

#combine both plots
st_scnrb = st_scn + st_rb + 
  patchwork::plot_layout(
    ncol = 2,
    nrow = 1,
    guides = 'collect', 
    axis_titles = 'collect'#,
    #widths = c(1,2)
  ) #+
  # patchwork::plot_annotation(
  #   title = "Short-Term Model",
  #   theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  # )
st_scnrb
ggsave(
  "Output/PDXB3_ST_SCN_RBL_violin_anovap.png",
  plot = st_scnrb,
  width = 5.5,
  height = 8.5,
  dpi = 600
)
ggsave(
  "Output/PDXB3_ST_SCN_RBL_violin_anovap.pdf",
  plot = st_scnrb,
  height = 3.5,
  width = 8
)

#
#Long Term Model (7I & J): SCN & RB Scores-----
PDXB3_anno <- readxl::read_xlsx("~/Dropbox/Graeber/OVP/Ovarian Project/PDX_CarboCycling/Batch3/PDXB3_anno.xlsx")
PDXB3_N_anno <- PDXB3_anno %>%
  filter(grepl("N", Sample_ID)) %>%
  mutate(Timepoint = factor(Timepoint, levels = unique(Timepoint)))
PDXB3_log2UQ <- Rubrary::rread("~/Dropbox/Graeber/OVP/Ovarian Project/PDX_CarboCycling/Batch3/PDXB3_rsem_genes_upper_norm_coding_log2.txt", row.names = 1)

PDXB3_N_log2UQ <- PDXB3_log2UQ %>%
  select(PDXB3_N_anno$Sample_ID)

#SCN Score
lt_scn = ggplot(PDXB3_N_anno, aes(x = Timepoint, y = SCN_Proj, fill = Timepoint)) +
  geom_violin(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.1) +
  xlab("Treatment Status") +
  ylab("Scores") +
  ggrepel::geom_label_repel(aes(label = Sample_ID), color = "black", position = pos, size = 2) +
  labs(title = "SCN",
       fill = "Treatment Status") +
  theme_classic() +
  ggpubr::stat_compare_means(
    method = "anova",
    label.y = 71, label.x = 1,
    label = "p.format"
  ) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(margin = margin(l = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold", family = "Arial")
  )
lt_scn
ggsave(
  "/Users/Zero/Dropbox/Graeber/OVP/PDX/Output/PDXB3_LT_SCN_violin.png",
  plot = lt_scn,
  width = 5.5, #6
  height = 5, #4.5
  dpi = 300
)

#RB1 Loss Score
#normalized ssgsea scores
pdxlt.rbl = calc_ssgsea(
  exp.df = PDXB3_N_log2UQ,
  gset = rb_lst_ssgsea,
  col2 = "RB_Loss_Chen"
) 
pdxlt.rblf <- merge(PDXB3_N_anno,pdxlt.rbl,by = 'Sample_ID')

#violin plot
lt_rb = ggplot(pdxlt.rblf, aes(x = Timepoint, y = RB_Loss_Malorni, fill = Timepoint)) +
  geom_violin(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.1) +
  xlab("Treatment Status") +
  ylab("Scores") +
  ggrepel::geom_label_repel(aes(label = Sample_ID), color = "black", position = pos, size = 2) +
  labs(title = "RB Loss", #tag = "J",
       fill = "Treatment Status") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(), # Remove x-axis text
    axis.title.y = element_text(margin = margin(l = 15)),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold", family = "Arial")
  ) +
  ggpubr::stat_compare_means(
    method = "anova",
    label.y = 1.3, label.x = 1,
    label = "p.format"
  ) 
lt_rb
ggsave(
  "Output/PDXB3_LT_RB1_violin.png",
  plot = last_plot(),
  width = 6,
  height = 4.5,
  dpi = 300
)

#combine both plots
lt_scnrb = lt_scn + lt_rb + 
  patchwork::plot_layout(
    ncol = 2,
    nrow = 1,
    guides = 'collect', 
    axis_titles = 'collect',#
    #widths = c(1,2)
  ) #+
  # patchwork::plot_annotation(
  #   title = "Long-Term Model",
  #   theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  # )
lt_scnrb

ggsave(
  "Output/PDXB3_LT_SCN_RBL_violin.pdf",
  plot = lt_scnrb,
  height = 3.5,
  width = 8
)
#
#Volcano Plots (7E & K)----

debugonce(create_volcano_plot)
#short term model
st_vc = prep_volcano_plot_data(
  deseq_data = data.table::fread("data/PDXB2_DESeq_R2vPT_slogp_PC.txt"),
  save_files = F
) %>% 
  create_volcano_plot(
    prepped_deseq_data = .[["deseq_data"]],
    top_genes_vp = .[["top_genes_vp"]],
    color_palette = c("Up" = scales::muted("#50C878"), "Down" = scales::muted("#673AB7"), "No Change" = "#b0b0b0"),
    point_size = 1,
    save = F
  ) +
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Arial"))

ggsave(
  "Output/PDXB2_ST_deseq_volc.pdf",
  plot = st_vc,
  width = 3.5,
  height = 5
)

#long term model
lt_vc = prep_volcano_plot_data(
  deseq_data = data.table::fread("data/PDXB3N_DESeq_PreRecur_results_pcgenes.txt"),
  save_files = F
) %>% 
  create_volcano_plot(
    prepped_deseq_data = .[["deseq_data"]],
    top_genes_vp = .[["top_genes_vp"]],
    color_palette = c("Upregulated" = "#50C878", "Downregulated" = "#673AB7", "No Change" = "#b0b0b0"),
    save = F
  ) +
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Arial")) 

ggsave(
  "Output/PDXB3_LT_deseq_volc.pdf",
  plot = lt_vc,
  width = 3.5,
  height = 5,
)
#
#GSEA-Squared Density (7K)------
categories = c("Immune Response", "Cell Differentiation", "Lipid Metabolism", "Epigenetic Regulation", "Cell Cycle Regulation", "DNA Damage Repair", "RB-E2F Dysregulation")
cat_terms = c("INFLAM|IMMUN|INTERLUEKIN|LEUKOCYTE|TNF|MHC|CYTOKINE_|CHEMOKINE|ANTIGEN|LYMPH",
              "DIFFERENTIATION",
              "_COA_|LIPID|STEROL|FATTY|FAT",
              "HISTONE|METHYL|ACETYL|SUMOYLATION|UBIQUITIN",
              "CELL_CYCLE|MITOTIC|DNA_REPLICATION|CHROMOSOME_SEGREGATION|SPINDLE|CELL_DIVISION|REPLICATION|STRAND|G2", 
              "REPAIR|STRESS|HDR|HRR|DAMAGE|FANCONI",
              "RB_1|RB1|RBL1|RBL2|RETINOBLASTOMA|E2F" 
)
#GSEA-squared
gseasq.lt = Rubrary::run_GSEA_squared(
  df_GSEA = "data/GSEA_PDXB3N_PreRecur_slogp_PC.txt",
  categories = categories,
  cat_terms = cat_terms,
  savename = "data/GSEA_PDXB3N_PreRecur_slogp_PC",
  height = 4.5, width = 6.5
)

gseasq.st = Rubrary::run_GSEA_squared(
  df_GSEA = "data/GSEA_PDXB3N_PrePost_slogp_PC.txt",
  categories = categories,
  cat_terms = cat_terms,,
  savename = "data/GSEA_PDXB3N_PrePost_slogp_PC",
  height = 4.5, width = 6.5
)

ts_cols <- setNames(
  nm = c("Chemonaive", "Post-NACT", "Recurrent", "NA (CellLine)", "NA (Control)"),
  c(scales::hue_pal()(3), "grey30", "grey60"))
scales::show_col(ts_cols)

#PDX Short & Long Term
Rubrary::plot_GSEAsq_density(
  GSEAsq_df1 = "data/GSEA_PDXB3N_PreRecur_slogp_PC_GSEAsq_pathways.txt",
  GSEAsq_df2 = "data/GSEA_PDXB3N_PrePost_slogp_PC_GSEAsq_pathways.txt",
  name1 = "Long-Term Model",
  name2 = "Short-Term Model",
  plot_pval = F,
  cat_order = rev(categories),
  title = " ",
  colors = unname(ts_cols[c(3, 2)]),
  savename = "Output/PDX_ST_LT_GSEAsq_dens.pdf",
  plot_fmt = "pdf",
  height = 5, width = 5
)

#PDX Short vs. Patient CPN
Rubrary::plot_GSEAsq_density(
  GSEAsq_df1 = "data/GSEA_PDXB3N_PrePost_slogp_PC_GSEAsq_pathways.txt",
  GSEAsq_df2 = "Output/Fig1/fGSEA_UCLAStarr_GSEA_DESeq_PostNACT_full_GSEAsq_pathways.txt",
  name1 = "PDX Short-Term",
  name2 = "Patient Post-NACT v. Chemonaive",
  plot_pval = F,
  cat_order = rev(categories),
  title = " ",
  colors = c(unname(ts_cols[2]), "lightslateblue"),
  savename = "Output/PDX_ST_PT_CPN_GSEAsq_dens.png",
  plot_fmt = "png",
  height = 5, width = 6
)

#PDX Long vs. Patient CRC
Rubrary::plot_GSEAsq_density(
  GSEAsq_df1 = "data/GSEA_PDXB3N_PreRecur_slogp_PC_GSEAsq_pathways.txt",
  GSEAsq_df2 = "data/fGSEA_UCLAAllPatch_deseq_recur_full_GSEAsq_pathways.txt",
  name1 = "PDX Long-Term",
  name2 = "Patient Recurrent v. Chemonaive",
  plot_pval = F,
  cat_order = rev(categories),
  title = " ",
  colors = c(unname(ts_cols[3]), "firebrick"),
  savename = "/Output/PDX_ST_PT_CRC_GSEAsq_dens.png",
  plot_fmt = "png",
  height = 5, width = 6
)

#



