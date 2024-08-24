
#run_fgsea----
#' run_fgsea
#' #http://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
#'
#' @param pathways list, MSigDB gene sets; uses Nick's script
#' @param deg.df df, deseq2 output
#' @param deg.df_slpval string, column name containing signed log p-values
#' @param deg.df_gene string, column name containing gene
#' @param topn int, number of top pathways to plot 
#' @param bottomn int, number of bottom pathways to plot 
#' @param collapsedPathways bool, select only independent pathways
#' @param save bool, save fgsea df
#' @param savename string
#'
#' @return df, fgsea
#' @export
#'
#' @examples
run_fgsea = function(cluster_id, pathways = NULL, subset = F, deg.df, deg.df_slpval, deg.df_gene = "gene", type = "get", maxSize = 500, nPermSimple = 1000, topn = 25, 
                     bottomn = 25, title = "GSEA", pw_size = 5, save = T, savename){
  
  require(msigdbr)
  require(fgsea)
  
  #get msigdb pathways
  if(is.null(pathways)){
    pathways = get.MSigDB.genesets(
      msig_df = rbind(
        msigdbr(species = "Homo sapiens", category = "C2"),
        msigdbr(species = "Homo sapiens", category = "C5"),
        msigdbr(species = "Homo sapiens", category = "H")
      ), #restart R session if error
      genesets = c("CP", "GO", "H$")
    )
  }
  
  #subset if needed
  if(subset){
    deg.df_full <- deg.df
    deg.df <- deg.df_full %>%
      filter(cluster == cluster_id)
  }
  
  #creates a named vector of DE metrics/rank
  if(type == "get")
  stats = setNames(
    deg.df[,get(deg.df_slpval)],
    deg.df[,get(deg.df_gene)]
  ) else{
  stats = setNames(
    deg.df[,paste0(deg.df_slpval)], #paste0
    deg.df[,paste0(deg.df_gene)]
  )
  }
  
  #run fgsea
  set.seed(999)
  fgseaRes <- fgsea(
    pathways = pathways,
    stats    = stats,
    eps      = 0.0,
    minSize  = 15,
    maxSize  = maxSize,
    nPermSimple = nPermSimple
  ) %>%
    arrange(desc(NES)) %>%
    mutate(signed_logp = sign(NES) * -log10(pval)) %>%
    na.omit(.)
  
  #plot genesets
  fgseaRes.top = fgseaRes[c(tail(order(NES), n = topn), head(order(NES), n = bottomn)), ]
  ggpubr::ggbarplot(
    fgseaRes.top,
    x = "pathway", 
    y = "NES",
    # fill = "NES",          
    # color = "NES",          
    fill = "signed_logp",          
    color = "signed_logp",          
    #palette = "jco",  
    #sort.val = sort.val,          
    sort.by.groups = FALSE,     
    x.text.angle = 90,   
    title = "GSEA",
    ylab = "NES",
    #legend.title = legend.title,
    #xlab = paste0(keyword),
    rotate = TRUE
    #ggtheme = theme_minimal()
  )
  
  # Plot
  col_gradient <- scale_fill_gradient2(
    low = scales::muted("red"), mid = "white", high = scales::muted("blue"),
    midpoint = 0, limits = range(fgseaRes$signed_logp)
  )
  ggplot(
    fgseaRes.top, 
    aes(x = NES, y = reorder(pathway, NES), fill = signed_logp)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +
    geom_vline(xintercept = -1.5, linetype = "dashed", color = "black") +
    labs(x = "NES", y = "Pathways", title = "GSEA") +
    col_gradient +
    theme_minimal() 
  # 
  
  #Ruby's portion
  g = Rubrary::plot_GSEA_barplot(
    gsea_res = fgseaRes,
    #gsea_pws = topPathways,
    gsea_pws = fgseaRes.top[, pathway],
    NES_cutoff = 2,
    sig_cutoff = c("pval", 0.05),
    #pw_format = TRUE,
    pw_size = pw_size,
    colors = c("firebrick", "darkblue"),
    title = title
  )
  
  fgseaRes$signedlogp = -log10(fgseaRes$pval) * sign(fgseaRes$NES)
  fgseaRes = fgseaRes[order(fgseaRes$signedlogp, decreasing = T),]
  fgseaRes$leadingEdge = apply(fgseaRes, 1, function(x){
    paste(
      unlist(x$leadingEdge), 
      collapse = ","
    )
  })
  fgseaRes$NAME = fgseaRes$pathway
  
  #save/output fgsea file
  if(save){
    filename = ifelse(
      is.null(cluster_id),
      paste0(savename, "_fgsea.txt"),
      paste0(savename, "_", cluster_id, "_fgsea.txt")
    )
    write.table(
      x = fgseaRes,
      file = filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  print(g)
  #return(g)
  #return(fgseaRes)
  return(list(df = fgseaRes, plot = g))
}




#
#Retrieve Top DEGs from GSEA categories----
#' retrieve_top_DEGs_gseacat
#'#' Retrieve top differentially expressed genes (DESeq2) & leading edge genes (GSEA) for each gsea-squared category
#'
#' @param gseasq.file file, gsea-squared output file
#' @param keyword.labels chr vector, keyword categories
#' @param gsea.n int, number of gene-sets to plot
#' @param pos.NES bool, do keyword categories have pos or neg NES?
#' @param sort.val str, asc or desc barplot
#' @param leading.edge.col int/str, col name or number containing leading edge genes
#' @param deseq.df df/mat, deseq output
#' @param deseq.pCutoff int, deseq p-val cutoff
#' @param deseq.n int, number of top genes to plot
#' @param colors chr vector, colrs corresponding to keyword categories
#' @param save bool, save barplot, volcano plot & write top category genes to file
#' @param h int; pathway NES waterfall height
#' @param w int; pathway NES waterfall width
#' @param len_pw int; max length of pathway name if it's too long
#' @param vw int; volcano plot width
#' @param vh int; volcano plot height
#' @param voltitle string; volcano plot title
#'
#' @return
#' @export
#'
#' @examples
retrieve_top_DEGs_gseacat = function(gseasq.file, keyword.labels, gsea.n = 50, pos.NES = T, sort.val = "desc", leading.edge.col, len_pw = 100,
                                     deseq.file, deseq.pCutoff = 5e-2, deseq.n = 20, colors = NULL, save = T, h = 11, w = 16,
                                     vw = 15, vh = 10, voltitle = "DE Genes"){
  
  df = read.delim(
    file = gseasq.file,
    stringsAsFactors = F
  )
  deseq.df = read.delim(
    file = deseq.file,
    stringsAsFactors = F
  )
  
  #retrieve genesets in selected category
  for(i in 1:length(keyword.labels)){
    keyword = keyword.labels[i]
    print(keyword)
    df.keyword = df[which(df$Category %in% keyword),]
    if(pos.NES){
      df.keyword = df.keyword[order(df.keyword$NES, decreasing = T),]
      df.keyword = df.keyword[df.keyword$NES > 0,]
      type = "pos"
    } else {
      df.keyword = df.keyword[order(df.keyword$NES),]
      df.keyword = df.keyword[df.keyword$NES < 0,]
      type = "neg"
    }
    title = paste0(tools::toTitleCase(keyword), " Genesets")
    gsea.n <- ifelse(nrow(df.keyword) < gsea.n, nrow(df.keyword), gsea.n)
    df.keyword$NAME <- substr(df.keyword$NAME, 1, len_pw)
    
    #plot top genesets in catagory
    gb = ggpubr::ggbarplot(
      df.keyword[1:gsea.n,],
      x = "NAME", 
      y = "NES",
      fill = "signedlogp",
      color = "signedlogp",
      sort.val = sort.val,          
      sort.by.groups = FALSE,     
      x.text.angle = 90,   
      title = paste0(title, " ", type),
      ylab = FALSE,
      xlab = "NES",
      legend.title = "signed log p-values",
      # xlab = paste0(keyword),
      rotate = TRUE,
      ggtheme = theme_minimal()
    )
    print(gb)
    
    #retrieve & find frequency of top leading-edge genes from significant genesets in GSEA-sq category
    df.keyword.sig <- df.keyword[which(df.keyword$padj <= deseq.pCutoff),]
    table.top.genes = as.data.frame(table(unlist(strsplit(df.keyword.sig[,leading.edge.col], ","))), 
                                    stringsAsFactors = F)
    
    #merge with significant DEGs from deseq2 list 
    table.top.genes$padj = deseq.df$padj[match(table.top.genes$Var1, deseq.df$gene)]
    table.top.genes.p = table.top.genes[which(table.top.genes$padj <= deseq.pCutoff),]
    table.top.genes.p = table.top.genes.p[order(table.top.genes.p$Freq, decreasing = T),]
    colnames(table.top.genes.p) = c("gene", "gsea_leadingedge_frequency", "deseq_padj")
    
    #co-rank of frequency & deseq p-val
    table.top.genes.p = table.top.genes.p %>% 
      arrange(desc(gsea_leadingedge_frequency), deseq_padj) %>%
      mutate(freq_rank = row_number()) %>%
      arrange(deseq_padj, desc(gsea_leadingedge_frequency)) %>%
      mutate(p_rank = row_number(), co_rank_category = (freq_rank + p_rank)/2) %>%
      select(-freq_rank, -p_rank) %>%
      arrange(co_rank_category) %>%
      mutate(category = keyword)
    
    #assign volcano plot colors to catagory
    color = ifelse(is.null(colors), "red", paste0(colors[i]))
    deseq.df$genetype = ifelse(
      deseq.df$gene %in% table.top.genes.p$gene, 
      paste0(color), 
      "black"
    )
    keyvals = deseq.df$genetype
    names(keyvals)[keyvals == 'black'] <- "other"
    names(keyvals)[keyvals == paste0(color)] <- paste0(keyword, " genes")
    
    #generate volcano plot
    gv = EnhancedVolcano::EnhancedVolcano(
      toptable = deseq.df,
      lab = deseq.df$gene,
      selectLab = table.top.genes.p$gene[1:deseq.n],
      x = 'log2FoldChange',
      y = 'padj',
      pCutoff = deseq.pCutoff,
      labSize = 2,
      pointSize = 1.5,
      colCustom = keyvals,
      boxedLabels = TRUE,
      drawConnectors = TRUE,
      title = voltitle,
      subtitle = paste0(tools::toTitleCase(keyword), " Genes"),
      caption = paste0("p-value cutoff = ", deseq.pCutoff, ", logFC cutoff = 1") 
    )
    print(gv)
    
    if(save){
      ggsave(
        filename = paste0(gsub(" ", "_", title), "_", type, ".png"),
        plot = gb,
        height = h,
        width = w
      )
      
      ggsave(
        filename = paste0(keyword, "_genes_volpt.png"),
        plot = gv,
        height = vh,
        width = vw
      )
      
      write.table(
        x = table.top.genes.p,
        file = paste0(gsub(" ", "_", title), "_leadingedge_genes.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    
    #combine top ~20 genes for each category
    if(i == 1){
      table.top.genes.full = table.top.genes.p[c(1:deseq.n),]
      next()
    }
    table.top.genes.full = rbind(table.top.genes.full, table.top.genes.p[c(1:deseq.n),])
  }
  return(table.top.genes.full)
}

#
#calc_ssgsea----
calc_ssgsea <- function(exp.df, gset, gsva.met = "ssgsea", ssgsea.norm = T, col1 = "Sample_ID", col2 = "Score", min.sz = 5, max.sz = 500){
  res = GSVA::gsva(
    expr = as.matrix(exp.df),
    gset.idx.list = gset,
    method = gsva.met,
    min.sz = min.sz, max.sz = max.sz,
    tau = 0.25, # default for ssGSEA by Barbie et al. 2009
    verbose = TRUE, ssgsea.norm = ssgsea.norm
    #BPPARAM = BiocParallel::MulticoreParam(workers = multicoreWorkers())
  ) %>%
    t(.) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(col1)
  colnames(res)[2] = col2
  return(res)
}
#

#Heatmap of Top DEGenes-----
generate_deg_heatmap <- function(heatmap_data, scale = "row", heatmap_color = viridis::viridis(200), col_annotation, col_color, col_split, col_split_order = NULL, show_heatmap_legend = TRUE,
                                 cluster_column_slices = F, row_annotation, row_split = NULL, row_split_labels = NULL, row_names_fontsize = 4.5, colnames_fontsize = 4.5, left_annotation = NULL, right_annotation = NULL,
                                 show_annotation_name = TRUE, display_gene_index = NULL, display_gene_label, cluster_row_slices = T, row_title = " ", column_title = " ", show_row_names = TRUE, 
                                 heatmap_legend = "Standardized Expression", top_anno_fontsize = 10, left_anno_fontsize = 5.5, legend_title_fontsize = 10, legend_label_fontsize = 10){
  
  set.seed(666)
  
  #Prep data 
  heatmap_matrix <- as.matrix(heatmap_data)
  
  if(!is.null(scale)){
    if(scale == "row"){
      heatmap_matrix <- t(scale(t(heatmap_matrix)))
    }
    if(scale == "col"){
      heatmap_matrix <- scale(heatmap_matrix)
    }
  }
  
  #Create annotations
  top_annotation <- HeatmapAnnotation(
    #df = col_annotation, 
    Treatment_Status = anno_block(
      labels = col_split_order,
      gp = gpar(fill = col_color$Treatment_Status),
      labels_gp = gpar(fontsize = top_anno_fontsize, fontface = "bold")),
    `Sample Type` = col_annotation$Sample_Type,
    `RB1-E2F Dysregulation` = col_annotation$RB_E2F_WP,
    which = "col", 
    border = T,
    col = col_color,
    show_annotation_name = show_annotation_name,
    annotation_legend_param = list(
      legend_direction = "horizontal",
      title_gp = gpar(fontsize = legend_title_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_label_fontsize)
    )
  )
  
  if(!is.null(row_split_labels)){
    left_annotation = rowAnnotation(
      category = anno_block(
        labels = row_split_labels,
        gp = gpar(fill = 2:8), #customize colors
        labels_gp = gpar(fontsize = left_anno_fontsize, 
                         fontface = "bold"))
    )
  }
  
  if(!is.null(display_gene_index)){
    right_annotation = rowAnnotation(
      foo = anno_mark(
        at = display_gene_index, 
        labels = display_gene_label,
        labels_gp = gpar(fontsize = 5)
      )
    ) 
  }
  
  #Generate Heatmap
  cht = ComplexHeatmap::Heatmap(
    matrix = heatmap_matrix, 
    col = heatmap_color,
    name = heatmap_legend, 
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    right_annotation = right_annotation,
    column_split = col_split, 
    row_split = row_split,
    cluster_row_slices = cluster_row_slices,
    cluster_column_slices = cluster_column_slices,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    column_dend_reorder = T,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    border = T,
    #row_gap = unit(3, "mm"),
    row_names_gp = gpar(fontsize = row_names_fontsize),
    column_names_gp = gpar(fontsize = colnames_fontsize),
    column_names_rot = 45,
    row_title = row_title, 
    column_title = column_title,
    heatmap_legend_param = list(title = heatmap_legend, 
                                color_bar = "continuous",
                                legend_direction = "horizontal",
                                title_gp = gpar(fontsize = legend_title_fontsize, fontface = "bold"),
                                show = show_heatmap_legend, labels_gp = gpar(fontsize = legend_label_fontsize))
  )
  cht
  
  #draw(cht, heatmap_legend_side = "bottom", annotation_legend_side = "right", main = "Paired UCLA samples")
  return(cht)
}

#process_heatmap_data-----
process_heatmap_data <- function(deg_data, expression_data, annotation_data, col_annotation,
                                 gene_col = "gene", category_col = "category", 
                                 sample_id_col= "Sample_ID", sample_col = "Sample") {
  
  # Merge DEGs with expression data
  heatmap_data_full <- merge(deg_data[, c(gene_col, category_col)], expression_data, by = gene_col) %>%
    arrange(!!sym(category_col))
  
  # Order genes by category and prepare the heatmap data
  heatmap_data <- heatmap_data_full %>%
    column_to_rownames(var = gene_col) %>%
    select(-all_of(category_col)) 
  
  # Match annotation names with those on heatmap matrix
  colnames(heatmap_data) <- annotation_data$Sample[match(colnames(heatmap_data), annotation_data[[sample_id_col]])] 
  heatmap_matrix <- as.matrix(
    heatmap_data[, match(rownames(col_annotation), colnames(heatmap_data))]
  )
  
  # Prepare row annotation
  row_annotation <- heatmap_data_full %>%
    select(all_of(c(gene_col, category_col))) %>%
    rename(category = category_col) %>%
    column_to_rownames(var = gene_col)
  
  # Return both heatmap data and row annotations
  output = list(heatmap_matrix = heatmap_matrix, row_annotation = row_annotation)
  return(output)
}

#

#Paired Boxplot-----
library(ggplot2)
library(ggpubr)
library(cowplot)

# Function to combine p-values using Fisher's method
combine_pvalues_fisher <- function(p_values) {
  chisq_stat <- -2 * sum(log(p_values))
  combined_p_value <- pchisq(chisq_stat, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p_value)
}

# Define the customized paired plot function
create_paired_boxplot <- function(data, x, y, xlab = NULL, ylab = "Signature Score", title, facet_by = "Dataset", color = "black", fill, 
                                  line_color = "gray", line_size = 0.4, text_size = 14.5, legend.position = "top", legend.justification = "center", 
                                  p_anno_lab = NULL, p_anno_n = 1, p_anno_x = 1.1, p_anno_y = 75, p_label_y = NULL,
                                  ggtheme = theme_cowplot(), palette = paletteTS[c(1,3)], y_limits = NULL, y_breaks = NULL,
                                  p_method = "wilcox.test", p_label = "p.signif", p_label_x = 1.3, p_paired = TRUE, plot_tag = NULL) {
  
  # Create the paired plot
  p <- ggpaired(
    data = data,
    x = x, 
    y = y,
    xlab = xlab,
    ylab = ylab,
    title = title,
    facet.by = facet_by,
    color = color, 
    fill = fill,
    line.color = line_color, 
    line.size = line_size,
    ggtheme = ggtheme +  # Start with the given theme
      theme(
        text = element_text(size = text_size),  # Increase base text size
        #axis.text = element_text(size = text_size), # Increase axis text size
        axis.title = element_text(size = 16, face = "bold"),  # Increase axis title size
        plot.title = element_text(size = 18),  # Increase and center plot title size
        strip.text = element_text(size = 16),  # Increase facet label size
        legend.position = legend.position,
        legend.justification = legend.justification,
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.title.x = element_blank() 
      ),
    palette = palette
  ) + labs(fill = "Treatment Status")
  
  #p-value type
  if(p_method != "combined"){
    p <- p + stat_compare_means(
      paired = p_paired, 
      label = p_label, 
      label.x = p_label_x, 
      label.y = p_label_y,
      method = p_method
    )
  } else {
    # Calculate p-values for each dataset and combine them using Fisher's method
    p_values <- data %>%
      group_by(across(all_of(facet_by))) %>%
      summarize(p_value = wilcox.test(get(y) ~ get(x), paired = p_paired)$p.value) %>%
      pull(p_value)
    
    print(p_values)
    
    combined_p_value <- format(
      combine_pvalues_fisher(p_values),
      scientific = TRUE, digits = 3)
    
    # Annotate the combined p-value on the plot
    if(is.null(p_anno_lab)){
      p_anno_lab <- unique(data$Dataset)
    }
    
    p <- p +
      geom_text(
        data = data.frame(
          Dataset = p_anno_lab, 
          label = c(
            paste0("Combined p = ", combined_p_value), 
            rep(" ", p_anno_n))
          ) ,
        x = p_anno_x,
        y = p_anno_y, 
        aes(label = label),
        size = 4
      )
    # data.frame(
    #   Dataset = unique(ovp_anno$Dataset), 
    #   label = c(paste0("Combined p = ", combined_p_value), " "))
  }
  
  # Set y-axis limits and breaks if specified
  if (!is.null(y_limits) || !is.null(y_breaks)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks)
  }
  
  # Add plot tag if specified
  if (!is.null(plot_tag)) {
    p <- p + 
      labs(tag = plot_tag) + 
      theme(plot.tag = element_text(size = 14, face = "bold"), 
            plot.tag.position = c(0.01, 0.98))
  }
  
  return(p)
}

#
#prep_volcano_plot_data----
prep_volcano_plot_data <- function(deseq_data, lfc = 1, pval = 0.01, pval_column = "padj", num_up_genes = 12, num_down_genes = 10, pca = TRUE, pca_columns = c("log2FoldChange", "signed_logpvalue"), rank_column = "signed_logpvalue",
                                   save_files = TRUE, include_genes = NULL, exclude_genes = NULL, savename = "df") 
{
  
  library(dplyr) 
  
  # Prep data for volcano plot
  deseq_data <- deseq_data %>%
    filter(!is.na(padj)) %>%
    mutate(color = case_when(
      ((log2FoldChange >= lfc) & (get(pval_column) <= pval & sign(log2FoldChange) > 0)) ~ "Upregulated",
      ((log2FoldChange <= -lfc) & (get(pval_column) <= pval & sign(log2FoldChange) < 0)) ~ "Downregulated",
      TRUE ~ "No Change"
    ))
  
  # Combine LFC & SignedLogP with PCA to rank genes; filter based on threshold
  top_genes <- deseq_data %>%
    mutate(
      # Normalize and perform PCA on the two variables
      PC1 = {
        pca_res <- stats::prcomp(
          .[, ..pca_columns],
          scale. = TRUE, 
          center = TRUE
        )$x[,"PC1"]  # Extracting the first principal component
      },
      PC1 = ifelse(sign(log2FoldChange) != sign(PC1), -PC1, PC1) # Adjust PCA sign to match log2FoldChange sign
    ) %>%
    filter(get(pval_column) < pval, abs(log2FoldChange) > lfc)
  
  # Rank genes based on PCA or other column (pval)
  if (pca) {
    top_genes <- top_genes %>% arrange(desc(abs(PC1)))
  } else {
    top_genes <- top_genes %>% arrange(desc(abs(rank_column)))
  }
  
  # Select top genes shown on volcano plot
  top_genes_vp <- rbind(
    top_genes %>% 
      filter(log2FoldChange > 0) %>%
      slice_head(n = num_up_genes),
    top_genes %>%
      filter(log2FoldChange < 0) %>%
      slice_head(n = num_down_genes)
  )
  
  # Add specific genes to top_genes_vp
  if (!is.null(include_genes)) {
    include_genes_df <- top_genes %>% 
      filter(gene %in% include_genes)
    top_genes_vp <- rbind(top_genes_vp, include_genes_df)
  }
  
  # Optionally filter out specific genes
  if (!is.null(exclude_genes)) {
    top_genes_vp <- top_genes_vp %>%
      filter(!gene %in% exclude_genes)
  }
  
  # Save files for final plot
  if (save_files) {
    write.table(
      x = deseq_data,
      file = paste0(savename, "_deseq_volc_data.txt"),
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )
    write.table(
      x = top_genes_vp,
      file = paste0(savename, "_deseq_volc_data_top.txt"),
      row.names = FALSE,
      quote = FALSE,
      sep = '\t'
    )
  }
  
  return(list(deseq_data = deseq_data, top_genes_vp = top_genes_vp))
}
#
#generate_volcano_plot----
create_volcano_plot <- function(prepped_deseq_data, top_genes_vp, data, lfc = 1, pval_column = "padj", pval_threshold = 0.01, point_shape = 21, point_size = 3, point_alpha = 0.5, 
                                color_palette = c("Upregulated" = "#e5534b", "Downregulated" = "#4f81c7", "No Change" = "#b0b0b0"), xlab = "Log2 Fold Change", ylab = "-Log10(padj)", 
                                legend_position = "top", text_size = 13, label_size = 3, x_limits = c(-4, 4), plot_height = 5.5, plot_width = 6, plot_dpi = 600, save = T, savename = "plot") {
  
  library(ggplot2)
  volcano_plot <- ggplot(prepped_deseq_data, aes(x = log2FoldChange, y = -log10(get(pval_column)), color = color)) +
    geom_point(aes(fill = color), shape = point_shape, size = point_size, stroke = lfc, alpha = point_alpha) +  # Using some transparency to manage overplotting
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = xlab, 
      y = ylab,
      color = " ",
      fill = " ") +
    theme_classic() +  # Minimal theme to keep it clean and professional
    theme(
      legend.position = legend_position,  # Set legend to top
      text = element_text(size = text_size, face = "bold")) +  # Center the title
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "darkgrey") +  # Threshold lines for LFC
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkgrey") +  # Threshold line for p-value
    scale_x_continuous(limits = x_limits) +
    ggrepel::geom_label_repel(
      data = top_genes_vp, 
      aes(label = gene), 
      size = label_size, 
      box.padding = unit(0.35, "lines"),  # Distance between text and points
      point.padding = unit(0.5, "lines")  # Avoid text overlapping with points
    )
  print(volcano_plot)
  
  if(save){
    ggsave(
      filename = paste0(savename, "_deseq_volcplt.png"),
      plot = volcano_plot, 
      height = plot_height, 
      width = plot_width, 
      dpi = plot_dpi
    )
  }
  
  return(volcano_plot)
}
#
#create_rrho_density_plot------

# Define the RRHO density plot function
create_rrho_density_plot <- function(file1, file2, pathway_col, filt_pathways = NULL, nes_col, output_file, plot_title = "RRHO Density Plot", plot_tag = "",
                                     x_labels = c("Label1", "Label2"), y_labels = c("Label3", "Label4"), axis.text.size = 15, legend.position = "right",
                                     x_labels_breaks = c(1650, 5200), y_labels_breaks = c(1650, 5200), axis.title.size = 17, x_axis_label = "", y_axis_label = "", width = 7.25, height = 6, save = F) {
  
  # Read the data
  Rank_1 <- data.table::fread(file1)
  Rank_2 <- data.table::fread(file2)
  
  # Find shared pathways
  Shared_Pathways <- intersect(Rank_1[[pathway_col]], Rank_2[[pathway_col]])
  
  #Filter shared pathways
  if(!is.null(filt_pathways)){
    Shared_Pathways <- Shared_Pathways[Shared_Pathways %in% filt_pathways]
  }
  
  # Filter and rank data
  Rank_1_filt <- Rank_1 %>%
    filter(Rank_1[[pathway_col]] %in% Shared_Pathways) %>%
    arrange(desc(.data[[nes_col]])) %>%
    mutate(rank_1 = row_number()) %>%
    select(pathway = .data[[pathway_col]], rank_1)
  
  Rank_2_filt <- Rank_2 %>%
    filter(Rank_2[[pathway_col]] %in% Shared_Pathways) %>%
    arrange(desc(.data[[nes_col]])) %>%
    mutate(rank_2 = row_number()) %>%
    select(pathway = .data[[pathway_col]], rank_2)
  
  # Merge ranks
  rank_merged <- merge(Rank_1_filt, Rank_2_filt, by = "pathway")
  
  # Calculate correlation and slope
  correlation <- cor(rank_merged$rank_1, rank_merged$rank_2, method = "spearman")
  lm_model <- lm(rank_1 ~ rank_2, data = rank_merged)
  slope <- coef(lm_model)[2]
  
  # Create the plot
  gp = ggplot(rank_merged, aes(x = rank_1, y = rank_2)) +
    theme(axis.ticks = element_blank()) +
    stat_density_2d(
      aes(fill = ..density..),
      geom = "raster", 
      contour = FALSE
    ) +
    scale_fill_distiller(
      palette = "Spectral",
      name = "Density"
    ) +
    scale_x_continuous(
      breaks = x_labels_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = y_labels_breaks,
      labels = y_labels,
      expand = c(0, 0)
    ) +
    theme(
      legend.position = legend.position, 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = .75), 
      #text = element_text(size = 20),
      plot.tag = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5, size = axis.text.size),
      axis.text.x = element_text(size = axis.text.size),
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = axis.title.size)
    ) +
    xlab(x_axis_label) + ylab(y_axis_label) +
    geom_point(size = .25, alpha = 0.2) +
    ggtitle(plot_title) +
    scale_fill_gradientn(colors = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100), name = "Density") 
    #scale_fill_viridis_c(name = "Density", option = "inferno") +
    labs(tag = plot_tag)
  
  #print(gp)
  
  if(save){  
    ggsave(
      output_file, 
      plot = gp, 
      width = width, 
      height = height, 
      dpi = 600
    )
  }
  
  return(gp)
}


#
