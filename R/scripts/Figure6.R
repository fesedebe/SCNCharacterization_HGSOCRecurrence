
library(ggplot2)

#GSEA - Normal vs FTE
gsea_filtered = data.table::fread("data/ducie_fte_gsea_fig.txt") %>%
  mutate(signed_log_padj = sign(NES) * -log10(padj))

ggplot(gsea_filtered, aes(x = reorder(pathway, NES), y = NES, fill = signed_log_padj)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradient2(low = "steelblue", mid = "lightgray", high = "darkred", midpoint = 0) +
  labs(
    title = " ",
    x = "Geneset",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(padj)"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) 

ggsave("Output/ducie_fte_gsea_fig.png", width = 6.3, height = 4, dpi = 600)
