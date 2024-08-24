
library(dplyr)
library(ggpubr)
library(Rubrary)
library(ggsignif)
library('GSVA')
library(readxl)

#SCN Histology Score Violin Plot - Chemonaive, Post-NACT, Recurrent----------------
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
    axis.line = element_line(size = 0.5),  
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
  annotate("text", x = 2.5, y = 98, label = "p <= 0.001", size = 4) +
  geom_segment(aes(x = 1, xend = 3, y = 102, yend = 102) , color = "gray40") +
  geom_point(x = 1, y = 102, color = "gray40") + geom_point(x = 3, y = 102, color = "gray40") +
  annotate("text", x = 2, y = 105, label = "p <= 0.021", size = 4)

hist3b_vp
ggsave(
  "output/HIST/scnlikecells_cnpnrc_vp.png", 
  plot = hist3b_vp, 
  height = 5.5, 
  width = 6, 
  dpi = 600
)
#
