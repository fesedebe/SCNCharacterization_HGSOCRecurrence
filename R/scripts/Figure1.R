
#load libraries-----
library(dplyr)
library(Rubrary)
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
  df_GSEA = "~/Library/CloudStorage/Box-Box/OVProject/NE_Features_in_HGSOC_Manuscript/Manuscript/Data/DEGenes&Pathways_FullFiles/PatientSamples_BulkRNASeq/txt_files/UCLA_Patch_Recurrent/fGSEA_UCLAAllPatch_deseq_recur.txt",
  categories = categories,
  cat_terms = cat_terms,
  title = "Chemonaive v. Recurrent",
  savename = "~/Dropbox/Graeber/OVP/OV_NE/Figures/Output/Fig1/fGSEA_UCLAAllPatch_deseq_recur_full_supp",
  height = 7, width = 7
)

#UCLA + Starr 
gseasq.UCLAAllPatch_deseq_postnact_full_supp = Rubrary::run_GSEA_squared(
  df_GSEA = "~/Dropbox/Graeber/OVP/Ovarian Project/signatures/DESeq_PostNACT_UCLAStarr/UCLAStarr/GSEA/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC.txt",
  categories = categories,
  cat_terms = cat_terms,
  title = "Chemonaive v. Post-NACT",
  savename = "~/Dropbox/Graeber/OVP/OV_NE/Figures/Output/Fig1/fGSEA_UCLAAllPatch_deseq_postnact_full_supp",
  height = 7, width = 7
)

##Density Plot (1B)
#UCLA & Patch
gseasq.UCLAAllPatch_deseq_recur_full = Rubrary::run_GSEA_squared(
  df_GSEA = "~/Library/CloudStorage/Box-Box/OVProject/NE_Features_in_HGSOC_Manuscript/Manuscript/Data/DEGenes&Pathways_FullFiles/PatientSamples_BulkRNASeq/txt_files/UCLA_Patch_Recurrent/fGSEA_UCLAAllPatch_deseq_recur.txt",
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
  df_GSEA = "~/Dropbox/Graeber/OVP/Ovarian Project/signatures/DESeq_PostNACT_UCLAStarr/UCLAStarr/GSEA/UCLAStarr_GSEA_DESeq_PostNACT_slogp_PC.txt",
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