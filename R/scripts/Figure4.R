
library(dplyr)
library(ggpubr)
library(Rubrary)
library(ggsignif)
library('GSVA')
library(readxl)

#SCN Histology Score Violin Plot - Chemonaive, Post-NACT, Recurrent (S4A)----------------
hist3b <- data.table::fread("./data/Hist/Fig3B.csv") %>%
  filter(SSN != '.') %>%
  mutate(
    SSN = as.numeric(SSN) * 100,
    Group = ifelse(Group == "Post", "Post-NACT", Group) 
  )

# Create the violin plot
hist3b_vp = ggplot(hist3b, aes(x = Group, y = SSN, fill = Group)) +
  geom_violin(scale = "width", adjust = 1, width = 0.8, trim = TRUE) +
  geom_jitter(width = 0.05, size = 1, color = "gray30", alpha = 0.7) + 
  scale_fill_manual(values = c('#F9CEAB','#F9B3AB','#AE7D78')) +  
  labs(y = "SCN Histology Score (%)", x = NULL) +  
  scale_y_continuous(breaks = seq(0, 105, by = 20), limits = c(0, 105)) +  
  theme_classic() +  
  theme(
    axis.line = element_line(linewidth = 0.5),  
    text = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), 
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5) 
  ) + 
  geom_segment(aes(x = 1, xend = 2, y = 93, yend = 93), color = "gray40") +
  geom_point(x = 1, y = 93, color = "gray40") + geom_point(x = 2, y = 93, color = "gray40")  + 
  annotate("text", x = 1.5, y = 96, label = "p = 0.018", size = 4) +
  geom_segment(aes(x = 2, xend = 3, y = 95, yend = 95), color = "gray40") + 
  geom_point(x = 2, y = 95, color = "gray40") + geom_point(x = 3, y = 95, color = "gray40") +
  annotate("text", x = 2.5, y = 98, label = "p < 0.001", size = 4) +
  geom_segment(aes(x = 1, xend = 3, y = 102, yend = 102) , color = "gray40") +
  geom_point(x = 1, y = 102, color = "gray40") + geom_point(x = 3, y = 102, color = "gray40") +
  annotate("text", x = 2, y = 105, label = "p = 0.021", size = 4)

hist3b_vp
ggsave(
  "Output/scnlikecells_cnpnrc_vp.png", 
  plot = hist3b_vp, 
  height = 5.5, 
  width = 6, 
  dpi = 600
) 
#
#SCN Histology Score Boxplot - Paired Chemonaive to Recurrent (4A)-------
paired_histology <- data.table::fread('./data/Hist/Fig3D.csv') %>%
  mutate(avg_NE = avg_NE*100)

rows_cn <- paired_histology %>%
  filter(TreatmentStatus == "Chemonaive") %>%
  select(Patient, TreatmentStatus, Tissue_ID, avg_NE) %>%
  rename(TS_CN = TreatmentStatus, ID_CN = Tissue_ID, Patho_Avg_CN = avg_NE)

segs_RC <- paired_histology %>%
  filter(TreatmentStatus == "Recurrent") %>%
  select(Patient, TreatmentStatus, Tissue_ID, avg_NE) %>%
  rename(TS_X = TreatmentStatus, ID_X = Tissue_ID, Patho_Avg_X = avg_NE) %>%
  left_join(., rows_cn, by = "Patient")

n_rc <- paired_histology %>%
  filter(TreatmentStatus == "Recurrent") %>%
  pull(Patient) %>%
  table() %>%
  `names<-`(as.character(names(.)))

pt_df <- function(pt){
  cn_row <- paired_histology %>%
    filter(Patient == pt) %>%
    filter(TreatmentStatus == "Chemonaive") %>%
    slice(rep(1:n(), each = n_rc[pt])) %>%
    mutate(TS2 = paste0(pt, "_", TS2, 1:n_rc[pt]))
  
  rc_row <- paired_histology %>%
    filter(Patient == pt) %>%
    filter(TreatmentStatus == "Recurrent") %>%
    mutate(TS2 = paste0(pt, "_", TS2))
  
  pt_rows <- bind_rows(cn_row, rc_row) %>%
    select(Patient, TreatmentStatus, Tissue_ID, avg_NE, TS2) %>%
    mutate(pair = gsub("[[:alpha:]]", "", TS2)) %>%
    arrange(pair)
}

hist3d_vp = ggplot(paired_histology, aes(x = TreatmentStatus, y = avg_NE)) +
  geom_boxplot(
    outlier.shape = NA, 
    aes(fill = TreatmentStatus), 
    staplewidth = 0.2, 
    varwidth = T,
    width = 0.65  
  ) +
  geom_segment(
    data = segs_RC,
    aes(x = TS_CN, xend = TS_X, y = Patho_Avg_CN, yend = Patho_Avg_X),
    color = "darkgray",
    linewidth = 0.5,  
    alpha = 0.5 
  ) +
  geom_point(size = 1) +
  scale_fill_manual(values = c('#F9CEAB','#AE7D78')) +
  labs(y = "SCN Histology Score (%)", x = NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.5),
    text = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) + 
  geom_segment(aes(x = 1, xend = 2, y = 93, yend = 93), color = "gray20", linewidth = 0.25) +
  geom_point(x = 1, y = 93, color = "gray20") + geom_point(x = 2, y = 93, color = "gray20")  + 
  annotate("text", x = 1.5, y = 95, label = "p = 0.047", size = 3.5) #+
hist3d_vp

ggsave(
  "Output/scnlikecells_paired_cnrc_vp_alt.png", 
  plot = hist3d_vp, 
  height = 5, 
  width = 5.5,
  dpi = 600
)
#Optimal Debulking vs. NACT VP (4D)-------

#change Post-NACT to NACT in Treatment column of treatment1
treatment1 <- data.table::fread('./data/Hist/fig3F.csv') %>% 
  mutate(Treatment = ifelse(Treatment %in% "Post-NACT", "NACT", Treatment))

ggplot(treatment1, aes(x = Treatment, y = SCN_Score_perc,fill= Treatment)) +
  geom_boxplot()  +
  geom_jitter(width = 0) + labs(x = NULL, y = "% SCN-like cells") + 
  scale_fill_manual(values = c('#F9CEAB','#AE7D78'),guide= 'none') +
  theme_classic() 

opt_debulk = ggplot(treatment1, aes(x = Treatment, y = SCN_Score_perc, fill = Treatment)) +
  geom_violin(scale = "width", adjust = 1, width = 0.8, trim = TRUE) +
  geom_jitter(width = 0.05, size = 1, color = "gray30", alpha = 0.7) +
  scale_fill_manual(values = c("#f9abd8", "#8f7dae")) + 
  labs(y = "SCN Histology Score (%)", x = NULL) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.5), 
    text = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5) 
  ) +
  stat_compare_means(
    comparisons = list(c("NACT", "Optimal Debulking")),
    method = "wilcox.test",
    label = "p.format"
  )
opt_debulk
ggsave(
  "Output/opt_debulk_vp.png", 
  plot = opt_debulk, 
  height = 5, 
  width = 5.5, 
  dpi = 600
)

#
#Survival Analyses (4E)-----
library(survival)
surv <- data.table::fread("data/Hist/fig3H.csv") %>% 
  mutate(
    SCN_status = ifelse(`SCN Score` >= 0.1, "SCN score >= 0.1", "SCN score < 0.1"),
    `PFS (months)` = as.numeric(`PFS (months)`)
  )

km_trt_fit <- survival::survfit(survival::Surv(surv$`PFS (months)`, surv$Recurrence) ~ surv$SCN_status, data=surv)
survival::survdiff(survival::Surv(`PFS (months)`, Recurrence) ~ SCN_status, data = surv)

plot(km_trt_fit, xlab = "Time in months", ylab = "Survival probability", main = "Survival Curve", xaxt = 'n')

# Adding custom x-axis
axis(1, at = seq(from = 0, to = 120, by = 6), labels = seq(from = 0, to = 120, by = 6))

p <- survminer::ggsurvplot(km_trt_fit, data = surv, xlab = "Progression Free Survival (months)", ylab = "Probability", ggtheme = theme_classic(), palette = 'Dark2',pval = T) 
p$plot + scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + scale_x_continuous(breaks = seq(0, 120, by = 6)) +
  theme(
    text = element_text(size = 14, color = "black"),      
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"), 
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black") 
  )

# Create the survival plot with enhancements
p <- survminer::ggsurvplot(
  km_trt_fit, 
  data = surv, 
  xlab = "Progression Free Survival (months)", 
  ylab = "Probability", 
  ggtheme = theme_classic(), 
  palette = 'Dark2', 
  conf.int = TRUE,
  legend.title = "SCN Score",
  legend.labs = c("Low", "High")
)

# Enhance the plot aesthetics
p$plot <- p$plot +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 120, by = 12)) +
  theme(
    text = element_text(size = 15, color = "black"),      
    axis.title = element_text(size = 15, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"), 
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"), 
    panel.grid.minor = element_line(color = "gray90", size = 0.25),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = 16, hjust = 0.5) 
  ) + 
  annotate("text", x = 60, y = 1, label = "p = 0.081", size = 3.5)
p
ggsave(
  "Output/kpm_survcurv.png", 
  plot = p$plot, 
  height = 5, 
  width = 6.5, 
  dpi = 600
)


