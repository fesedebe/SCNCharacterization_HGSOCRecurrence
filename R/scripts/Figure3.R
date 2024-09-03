#load libraries----
library(dplyr)
library(ggplot2)
library(Seurat)

source("./R/functions/scRNAseq_Integration.R")
source("clean_UCLA_sobj.R")
#
#load data----
sobj_rpca_sct = readRDS("data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpca.rds") 
sobj_rpca_sct_clust = readRDS("data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpcaSCT_leidenclust_sigscores.rds") 

#
#run normalization (SCT) on filtered data while regressing out cell cycle genes & integration (SRPCA) ----

options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB

#Score cells with cell cycle markers from Tirosh et al, 2015 (loaded with Seurat)
sobj_rpca_sct_clust <- Seurat::CellCycleScoring(
  object = readRDS("data/ucla.sobj.dx.filt.epi.pat.rds"),
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
) %>%
  
  #SCT Normalized + Seurat RPCA Integration — regularized negative binomial regression
  run_integration_with_norm(
    sobj = .,
    norm.method = "SCT",
    vars.to.regress = c("S.Score", "G2M.Score"),
    integration.method = "rpca",
    save = T,
    savename = "data/ucla.sobj.dx.filt.epi.pat.sct_ccr" 
  )
#run clustering (leiden) on integrated data----
#leiden normalization (might require reticulate/python calling)
sobj_rpca_sct_clust = run_leiden_clustering(
  seurat_obj = readRDS("data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpca.rds"), # sobj_rpca_sct_clust,
  resolution = 0.4, 
  reduction = "pca",
  savename = "data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpca"
)
DimPlot(
  sobj_rpca_sct_clust, 
  group.by = c("seurat_clusters", "sample", "TreatmentStatus", "SampleType"),
  ncol = 2
)
#add signature scores + umap scores-----
sobj_rpca_sct_clust <- sobj_rpca_sct_clust %>% 
  clean_UCLA_sobj(.) %>%
  AddModuleScore(
    .,
    features = list(read.csv("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/RB/Malorni_RBsig_genes.csv")[,1]),
    search = T,
    ctrl = 75, 
    name = "RBLoss_Malorni"
  ) %>%
  AddModuleScore(
    .,
    features = list(read.delim("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/RB/RB1Loss_Signature_Chen.txt", stringsAsFactors = F)[,1]),
    search = T,
    ctrl = 75, 
    name = "RBLoss_Chen"
  ) %>% 
  AddModuleScore(
    .,
    features = list(
      read.delim("~/Dropbox/Graeber/OVP/Ovarian Project/signatures/SCN_Signature/SCN_Signature.txt") %>%
        arrange(desc(PC1.v)) %>%
        slice_head(n = 500) %>%
        pull(Gene)),
    search = T,
    ctrl = 75,
    name = "SCN_Balanis"
  )

#UMAP Plot
dimset <- data.frame(
  group = c("seurat_clusters", "NEScore", "RBLoss_Malorni1", "SampleType", "TreatmentStatus", "SCN_Balanis1", "RBLoss_Chen1", "Phase"),
  name = c("Harmony Integration", "SCN Score", "RB Loss Score", "Sample Type", "Treatment Status", "SCN Score (ams)", "RB Loss Chen", "Cell Cycle Phase"),
  label = c(T, F, F, F, F, F, F, F)
)
Rubrary::plot_dimplot_grid(
  sobj = sobj_rpca_sct_clust,
  reduction = "umap",
  set = dimset,
  ncol = 4, 
  width = 15.5, height = 6.5,
  savename = "output/UMAP_dimset_rpca.sct.lei_all.png"
)
###RPCA
#SCN-hi cluster - 9
sobj_rpca_sct_clust@meta.data$ClusterGroup <- ifelse(
  sobj_rpca_sct_clust@meta.data$seurat_clusters %in% c(9), "HighSCN", "Other"
)

#ID of SCN Subpopulation----
##Boxplot of clusters & their signature scores
bp_scnclust_sig_final = ggplot(
  sobj_rpca_sct_clust@meta.data %>%
    select(seurat_clusters, RBLoss_Malorni1, SCN_Balanis1, ClusterGroup) %>%
    tidyr::pivot_longer(cols = c(RBLoss_Malorni1, SCN_Balanis1),
                        names_to = "ScoreType", 
                        values_to = "Score")%>%
    mutate(ScoreType = factor(ScoreType, levels = c("SCN_Balanis1", "RBLoss_Malorni1"))), 
  aes(x = as.factor(seurat_clusters), y = Score, color = ClusterGroup)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.5) +
  facet_wrap(~ ScoreType, scales = "free_y", ncol = 1, labeller = as_labeller(c(
    RBLoss_Malorni1 = "RB Loss",
    SCN_Balanis1 = "SCN"))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.25),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(x = "Seurat Clusters", y = "Signature Score", title = "Signature Scores by Seurat Cluster") +
  scale_color_manual(values = c("Other" = "#AEC6CF", "HighSCN" = "#FFB347")) +
  guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
print(bp_scnclust_sig_final)

ggsave(
  filename = "output/SCNPC1_RBlossAMS_Scores_byCluster_rpca.sct.lei.png",
  dpi = 700,
  plot = bp_scnclust_sig_final,
  height = 6,
  width = 9
)
#find (scn) markers-----
#SCN
Idents(sobj_rpca_sct_clust) <- "ClusterGroup"
mrks_mast_sobj_rsl_scn <- FindMarkers(
  object = PrepSCTFindMarkers(sobj_rpca_sct_clust), 
  ident.1 = "HighSCN", 
  test.use = "MAST",
  assay = "SCT")
mrks_mast_sobj_rsl_scn_full = mrks_mast_sobj_rsl_scn %>% 
  mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin, p_val),
         p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
         gene = rownames(.)) %>%
  mutate(signed_log_p_val = sign(avg_log2FC) * -log10(p_val),
         signed_log_p_val_adj = sign(avg_log2FC) * -log10(p_val_adj))
write.table(
  x = mrks_mast_sobj_rsl_scn_full,
  file = "./data/scnmarkers_MAST_padjfull_rsl.txt",
  sep = "\t", quote = F, row.names = F
)

mrks_mast_sobj_rsl_scn_filt <- mrks_mast_sobj_rsl_scn_full %>%
  filter(p_val_adj < 0.01,
         abs(avg_log2FC) >= 1)
write.table(
  x = mrks_mast_sobj_rsl_scn_filt,
  file = "./data/scnmarkers_MAST_padjfilt_rsl.txt",
  sep = "\t", quote = F, row.names = F
)

#saveRDS
saveRDS(
  sobj_rpca_sct_clust,
  file = "data/ucla.sobj.dx.filt.epi.pat.sct_ccr_integ.srpca_leidenclust_sigscores.rds"
)

#find all markers - Seurat
Idents(sobj_rpca_sct_clust) <- "seurat_clusters"

mrks_mast_sobj_rsl <- FindAllMarkers(
  object = PrepSCTFindMarkers(sobj_rpca_sct_clust), 
  test.use = "MAST",
  min.pct = 0.25,
  assay = "SCT"
)
mrks_mast_sobj_rsl_full = mrks_mast_sobj_rsl %>% 
  mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin, p_val),
         p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
         gene = rownames(.)) %>%
  mutate(signed_log_p_val = sign(avg_log2FC) * -log10(p_val),
         signed_log_p_val_adj = sign(avg_log2FC) * -log10(p_val_adj))
write.table(
  x = mrks_mast_sobj_rsl_full,
  file = "./data/topmarkers_MAST_padjfilt_rsl.txt",
  sep = "\t", quote = F, row.names = F
)

library(msigdbr)

#GSEA-----
#SCN9
debugonce(run_fgsea)
fgsea_rsl_scn9 <- run_fgsea(
  pathways = get.MSigDB.genesets(
    msig_df = rbind(
      msigdbr(species = "Homo sapiens", category = "C2"),
      msigdbr(species = "Homo sapiens", category = "C5"),
      msigdbr(species = "Homo sapiens", category = "H")
    ), #restart R session if error
    genesets = c("CP", "GO", "H$")
  ),
  deg.df = mrks_mast_sobj_rsl_scn_full,
  #deg.df = data.table::fread("data/UCLA_srpca_noPDX_rsl_allmarkers_scn_MAST_padjfilt.txt"),
  deg.df_slpval = "signed_log_p_val_adj",
  deg.df_gene = "gene", type = "paste0",
  topn = 15, bottomn = 15, #minSize = 10,
  savename = paste0("rsl_scn9_fgsea.txt")
)

gseasq_SCN9 = Rubrary::run_GSEA_squared(
  df_GSEA = "data/rsl_scn9_fgsea.txt_fgsea.txt",
  categories = c("Immune Response", "Cell Differentiation", "Lipid Metabolism", "Epigenetic Regulation", "Cell Cycle Regulation", "DNA Damage Repair", "RB-E2F Dysregulation"),
  cat_terms = c("INFLAM|IMMUN|INTERLUEKIN|LEUKOCYTE|TNF|MHC|CYTOKINE_|CHEMOKINE|ANTIGEN|LYMPH",
                "DIFFERENTIATION",
                "_COA_|LIPID|STEROL|FATTY|FAT",
                "HISTONE|METHYL|ACETYL|SUMOYLATION|UBIQUITIN",
                "CELL_CYCLE|MITOTIC|DNA_REPLICATION|CHROMOSOME_SEGREGATION|SPINDLE|CELL_DIVISION|REPLICATION|STRAND|G2",
                "REPAIR|STRESS|HDR|HRR|DAMAGE|FANCONI",
                "RB_1|RB1|RBL1|RBL2|RETINOBLASTOMA|E2F" 
  ),
  savename = "./output/rsl_scn9_fGSEAsq",
  height = 4, width = 5
)
#
#Heatmap of SCN Markers-----



#
