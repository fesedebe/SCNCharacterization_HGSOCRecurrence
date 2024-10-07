library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr) 

#run_fgsea----
#' Run FGSEA on DESeq2 Results and Plot Pathways
#'
#' This function performs FGSEA (Fast Gene Set Enrichment Analysis) on DESeq2 differential expression results.
#' It ranks genes based on signed log p-values, runs the FGSEA algorithm, and generates visualizations for the top pathways.
#' The function supports MSigDB gene sets and allows pathway filtering, plotting of top and bottom pathways, and optional saving of results.
#'
#' @param pathways A list of gene sets (e.g., from MSigDB). If NULL, default gene sets from MSigDB are used.
#' @param deg.df A data frame containing DESeq2 output (with at least gene names and signed log p-values).
#' @param deg.df_slpval The name of the column in `deg.df` containing the signed log p-values.
#' @param deg.df_gene The name of the column in `deg.df` containing gene names. Defaults to "gene".
#' @param topn Integer, number of top pathways to plot. Defaults to 25.
#' @param bottomn Integer, number of bottom pathways to plot. Defaults to 25.
#' @param collapsedPathways Logical, whether to select only independent pathways. Defaults to FALSE.
#' @param save Logical, whether to save the FGSEA results. Defaults to TRUE.
#' @param savename The file name prefix for saving FGSEA results (without extension).
#'
#' @return A list containing the FGSEA results data frame and two plots (`df`, `plot`, `plot2`).
#' @export
#'
#' @examples
#' # Example usage:
#' fgsea_results <- run_fgsea(
#'   deg.df = deg_df, 
#'   deg.df_slpval = "signed_log_p", 
#'   savename = "my_fgsea_results"
#' )
run_fgsea = function(cluster_id, pathways = NULL, subset = F, deg.df, deg.df_slpval, deg.df_gene = "gene", type = "get", 
                     minSize = 15, maxSize = 500, nPermSimple = 1000, topn = 25, sort.val = NULL,
                     bottomn = 25, title = "GSEA", pw_size = 5, save = T, savename){
  
  require(msigdbr)
  require(fgsea)
  
  if(is.null(pathways)){
    pathways = get.MSigDB.genesets(
      msig_df = rbind(
        msigdbr(species = "Homo sapiens", category = "C2"),
        msigdbr(species = "Homo sapiens", category = "C5"),
        msigdbr(species = "Homo sapiens", category = "H")
      ),
      genesets = c("CP", "GO", "H$")
    )
  }
  
  #subset if needed
  if(subset){
    deg.df_full <- deg.df
    deg.df <- deg.df_full %>%
      filter(cluster == cluster_id)
  }
  
  #create a named vector of DE metrics/rank
  if(type == "get")
    stats = setNames(
      deg.df[,get(deg.df_slpval)],
      deg.df[,get(deg.df_gene)]
    ) else{
      stats = setNames(
        deg.df[,paste0(deg.df_slpval)], 
        deg.df[,paste0(deg.df_gene)]
      )
    }
  
  #run fgsea
  set.seed(999)
  fgseaRes <- fgsea(
    pathways = pathways,
    stats    = stats,
    eps      = 0.0,
    minSize  = minSize,
    maxSize  = maxSize,
    nPermSimple = nPermSimple
  ) %>%
    arrange(desc(NES)) %>%
    mutate(signed_logp = sign(NES) * -log10(pval)) %>%
    na.omit(.)
  
  #plot genesets
  fgseaRes.top = fgseaRes[c(tail(order(NES), n = topn), head(order(NES), n = bottomn)), ]
  fgseaRes.top$logp = abs(fgseaRes.top$signed_logp)
  g1 = ggpubr::ggbarplot(
    fgseaRes.top,
    x = "pathway", 
    y = "NES",
    fill = "logp", 
    color = NA,
    sort.by.groups = FALSE,    
    sort.val = sort.val,
    title = "GSEA: Enrichment of MHC-Related Pathways in Patch Recurrent Samples",
    ylab = "Normalized Enrichment Score (NES)",
    xlab = "Pathway",
    rotate = TRUE
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
      axis.text.y = element_text(size = 12, angle = 0),                            
      axis.title.x = element_text(size = 14, face = "bold"),            
      axis.title.y = element_text(size = 14, face = "bold"),            
      legend.title = element_text(size = 14, face = "bold"),            
      legend.text = element_text(size = 12) ,                            
      legend.position = "right"
    ) + 
    labs(fill = "-log10(pval)")
  print(g1)
  
  #Plot
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
  
  #Rubrary Plot
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
  return(list(df = fgseaRes, plot = g, plot2 = g1))
}

#
#Retrieve Top DEGs from GSEA categories----
#' Retrieve Top DEGs and Leading Edge Genes from GSEA Categories
#'
#' This function retrieves the top differentially expressed genes (DEGs) from a DESeq2 dataset and leading edge genes 
#' from GSEA results for each GSEA category. It creates barplots of top GSEA categories and volcano plots of the 
#' top DEGs, optionally saving the results to files.
#'
#' @param gseasq.file A file containing GSEA-squared output data.
#' @param keyword.labels A character vector of keyword categories to search for in the GSEA output.
#' @param gsea.n Integer, the number of gene sets to plot. Defaults to 50.
#' @param pos.NES Logical, indicating whether the keyword categories have positive or negative NES. Defaults to TRUE.
#' @param sort.val A string indicating the sort order for the bar plot, either "asc" or "desc". Defaults to "desc".
#' @param leading.edge.col The column name or index containing the leading edge genes in the GSEA output.
#' @param deseq.df A data frame or matrix containing DESeq2 output.
#' @param deseq.pCutoff A numeric p-value cutoff for selecting significant DEGs from the DESeq2 output. Defaults to 0.05.
#' @param deseq.n Integer, the number of top DEGs to plot for each category. Defaults to 20.
#' @param colors A character vector of colors corresponding to the keyword categories.
#' @param save Logical, whether to save the bar plot, volcano plot, and top genes to file. Defaults to TRUE.
#' @param h Integer, the height of the GSEA bar plot. Defaults to 11.
#' @param w Integer, the width of the GSEA bar plot. Defaults to 16.
#' @param len_pw Integer, the maximum length of pathway names in the GSEA plot. Defaults to 100.
#' @param vw Integer, the width of the volcano plot. Defaults to 15.
#' @param vh Integer, the height of the volcano plot. Defaults to 10.
#' @param voltitle String, the title of the volcano plot. Defaults to "DE Genes".
#'
#' @return A data frame containing the top leading edge genes and their DESeq2 p-values.
#' @export
#'
#' @examples
#' # Retrieve top DEGs and leading edge genes, and generate plots
#' top_genes <- retrieve_top_DEGs_gseacat(
#'   gseasq.file = "gsea_results.txt", 
#'   keyword.labels = c("inflammation", "immune response"),
#'   deseq.file = "deseq_results.txt"
#' )
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
#' Calculate ssGSEA Scores for Expression Data
#'
#' This function calculates single-sample GSEA (ssGSEA) or GSVA scores for a given gene expression data frame.
#' The gene set enrichment scores are computed using the GSVA package and can be normalized if required.
#'
#' @param exp.df A data frame of expression values (genes as rows, samples as columns).
#' @param gset A list of gene sets to be used for enrichment scoring (e.g., MSigDB gene sets).
#' @param gsva.met The method to use for GSVA. Defaults to `"ssgsea"`. Other options include `"gsva"`, `"zscore"`, etc.
#' @param ssgsea.norm Logical, whether to normalize the ssGSEA scores. Defaults to `TRUE`.
#' @param col1 A string specifying the column name for sample identifiers in the output. Defaults to `"Sample_ID"`.
#' @param col2 A string specifying the column name for the ssGSEA or GSVA score in the output. Defaults to `"Score"`.
#' @param min.sz The minimum size of a gene set for scoring. Defaults to 5.
#' @param max.sz The maximum size of a gene set for scoring. Defaults to 500.
#'
#' @return A data frame with sample IDs and corresponding ssGSEA or GSVA scores.
#' @export
#'
#' @examples
#' # Calculate ssGSEA scores for a gene expression matrix
#' ssGSEA_scores <- calc_ssgsea(
#'   exp.df = expression_matrix, 
#'   gset = msigdb_gene_sets, 
#'   gsva.met = "ssgsea", 
#'   ssgsea.norm = TRUE
#' )
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
#' Generate Heatmap of Top Differentially Expressed Genes
#'
#' This function generates a customizable heatmap of top differentially expressed genes (DEGs) using the ComplexHeatmap package.
#' It supports row and column annotations, clustering, and customized color palettes for high-impact visualizations.
#'
#' @param heatmap_data A matrix or data frame of gene expression values (genes as rows, samples as columns).
#' @param scale Scaling method for the data. Can be "row", "col", or `NULL`. Defaults to "row".
#' @param heatmap_color A color palette for the heatmap values. Defaults to `viridis::viridis(200)`.
#' @param col_annotation A data frame containing annotations for columns (samples) in the heatmap.
#' @param col_color A list of colors to use for column annotations.
#' @param col_split A vector or factor to split the columns into groups.
#' @param col_split_order A vector defining the order of column splits.
#' @param show_heatmap_legend Logical, whether to display the heatmap legend. Defaults to TRUE.
#' @param cluster_column_slices Logical, whether to cluster column slices. Defaults to FALSE.
#' @param row_annotation A data frame or list for annotating rows.
#' @param row_split A vector or factor to split the rows into groups.
#' @param row_split_labels Labels for the row splits. Defaults to `NULL`.
#' @param row_names_fontsize Font size for row names. Defaults to 4.5.
#' @param colnames_fontsize Font size for column names. Defaults to 4.5.
#' @param left_annotation Annotations to display on the left side of the heatmap.
#' @param right_annotation Annotations to display on the right side of the heatmap.
#' @param show_annotation_name Logical, whether to show annotation names. Defaults to TRUE.
#' @param display_gene_index A vector of indices for genes to highlight with labels.
#' @param display_gene_label A vector of gene labels to display, corresponding to `display_gene_index`.
#' @param cluster_row_slices Logical, whether to cluster row slices. Defaults to TRUE.
#' @param row_title The title for the rows in the heatmap. Defaults to a space.
#' @param column_title The title for the columns in the heatmap. Defaults to a space.
#' @param show_row_names Logical, whether to display row names. Defaults to TRUE.
#' @param heatmap_legend Title for the heatmap legend. Defaults to "Standardized Expression".
#' @param top_anno_fontsize Font size for the top annotation text. Defaults to 10.
#' @param left_anno_fontsize Font size for the left annotation text. Defaults to 5.5.
#' @param legend_title_fontsize Font size for the heatmap legend title. Defaults to 10.
#' @param legend_label_fontsize Font size for the heatmap legend labels. Defaults to 10.
#'
#' @return A heatmap object generated using ComplexHeatmap.
#' @export
#'
#' @examples
#' # Generate a DEG heatmap with custom annotations and color palette
#' heatmap_obj <- generate_deg_heatmap(
#'   heatmap_data = deg_expression_matrix, 
#'   col_annotation = sample_annotations, 
#'   col_color = annotation_colors, 
#'   col_split = factor(sample_group), 
#'   row_split = factor(gene_category)
#' )
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
#' Process Data for Heatmap Visualization
#'
#' This function processes differential expression data, expression matrix, and annotation data 
#' to prepare a matrix for heatmap visualization. It merges DEGs with expression data, reorders 
#' the genes by category, and returns both the heatmap matrix and row annotations.
#'
#' @param deg_data A data frame containing differential expression data (e.g., gene and category columns).
#' @param expression_data A data frame or matrix of expression values (genes as rows, samples as columns).
#' @param annotation_data A data frame containing sample annotations, including sample identifiers.
#' @param col_annotation A data frame containing column annotations (samples) to match with the heatmap.
#' @param gene_col A string specifying the gene column in the `deg_data`. Defaults to `"gene"`.
#' @param category_col A string specifying the category column in the `deg_data`. Defaults to `"category"`.
#' @param sample_id_col A string specifying the column in the `annotation_data` that contains sample IDs. Defaults to `"Sample_ID"`.
#' @param sample_col A string specifying the sample column in `annotation_data` that corresponds to the expression matrix. Defaults to `"Sample"`.
#'
#' @return A list containing the processed heatmap matrix and row annotations.
#' @export
#'
#' @examples
#' # Example usage:
#' heatmap_data <- process_heatmap_data(
#'   deg_data = deg_df, 
#'   expression_data = expr_matrix, 
#'   annotation_data = sample_annotations, 
#'   col_annotation = col_anno
#' )
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
#' Create a Paired Boxplot with Statistical Comparison
#'
#' This function generates a paired boxplot with customizable aesthetics and statistical comparison. 
#' It supports faceting, custom themes, p-value annotation, and Fisher's combined p-value method.
#'
#' @param data A data frame containing the data to plot.
#' @param x A string specifying the x-axis variable.
#' @param y A string specifying the y-axis variable (e.g., signature score).
#' @param xlab Label for the x-axis. Defaults to NULL (no label).
#' @param ylab Label for the y-axis. Defaults to `"Signature Score"`.
#' @param title The title of the plot.
#' @param facet_by A string specifying the column to facet by. Defaults to `"Dataset"`.
#' @param color Color of the boxplot outline. Defaults to `"black"`.
#' @param fill The variable or color to use for filling the boxplot.
#' @param line_color Color of the lines connecting paired points. Defaults to `"gray"`.
#' @param line_size Size of the lines connecting paired points. Defaults to 0.4.
#' @param text_size Size of the base text elements. Defaults to 14.5.
#' @param legend.position Position of the legend. Defaults to `"top"`.
#' @param legend.justification Justification of the legend. Defaults to `"center"`.
#' @param p_anno_lab Labels for p-value annotations. Defaults to NULL (uses dataset labels).
#' @param p_anno_n Number of lines in the p-value annotation. Defaults to 1.
#' @param p_anno_x The x-axis position for p-value annotations. Defaults to 1.1.
#' @param p_anno_y The y-axis position for p-value annotations. Defaults to 75.
#' @param p_label_y The y-axis position for the p-value label (used with `stat_compare_means`). Defaults to NULL.
#' @param ggtheme A custom ggplot theme. Defaults to `theme_cowplot()`.
#' @param palette A color palette for the plot. Defaults to `paletteTS[c(1,3)]`.
#' @param y_limits Limits for the y-axis. Defaults to NULL (automatic).
#' @param y_breaks Breaks for the y-axis. Defaults to NULL (automatic).
#' @param p_method The method for statistical comparison (e.g., `"wilcox.test"`). Use `"combined"` for Fisher's method. Defaults to `"wilcox.test"`.
#' @param p_label The label format for the p-values (e.g., `"p.signif"`). Defaults to `"p.signif"`.
#' @param p_label_x The x-axis position for p-value labels. Defaults to 1.3.
#' @param p_paired Logical, whether the comparison is paired. Defaults to `TRUE`.
#' @param plot_tag A tag to add to the plot (e.g., "A", "B"). Defaults to NULL (no tag).
#'
#' @return Returns a ggplot object representing the paired boxplot.
#' @export
#'
#' @examples
#' # Create a paired boxplot with statistical comparison
#' create_paired_boxplot(
#'   data = my_data, 
#'   x = "Treatment", 
#'   y = "Score", 
#'   title = "Paired Boxplot", 
#'   facet_by = "Dataset", 
#'   p_method = "combined"
#' )
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
    ggtheme = ggtheme +  
      theme(
        text = element_text(size = text_size),  
        axis.title = element_text(size = 16, face = "bold"),  
        plot.title = element_text(size = 18),  
        strip.text = element_text(size = 16),  
        legend.position = legend.position,
        legend.justification = legend.justification,
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = text_size),
        axis.text.y = element_text(size = text_size),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
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

# Function to combine p-values using Fisher's method
combine_pvalues_fisher <- function(p_values) {
  chisq_stat <- -2 * sum(log(p_values))
  combined_p_value <- pchisq(chisq_stat, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p_value)
}


#
#prep_volcano_plot_data----
#' Prepare Data for Volcano Plot from DESeq2 Results
#'
#' This function processes DESeq2 results to prepare data for generating a volcano plot. 
#' It filters significant upregulated and downregulated genes based on log2 fold change (LFC) 
#' and p-value thresholds, and ranks genes based on PCA or a custom ranking column.
#'
#' @param deseq_data A data frame containing DESeq2 results.
#' @param lfc Numeric, the log2 fold change threshold for significance. Defaults to 1.
#' @param pval Numeric, the p-value threshold for significance. Defaults to 0.01.
#' @param pval_column A string specifying the column name for p-values. Defaults to `"padj"`.
#' @param num_up_genes Integer, the number of top upregulated genes to display. Defaults to 12.
#' @param num_down_genes Integer, the number of top downregulated genes to display. Defaults to 10.
#' @param pca Logical, whether to rank genes using PCA. Defaults to `TRUE`.
#' @param pca_columns A character vector specifying columns to use for PCA. Defaults to `c("log2FoldChange", "signed_logpvalue")`.
#' @param rank_column A string specifying the column to rank genes by if not using PCA. Defaults to `"signed_logpvalue"`.
#' @param save_files Logical, whether to save the filtered DESeq2 data and top genes to files. Defaults to `TRUE`.
#' @param include_genes A vector of specific genes to include in the top genes list. Defaults to `NULL`.
#' @param exclude_genes A vector of specific genes to exclude from the top genes list. Defaults to `NULL`.
#' @param savename A string prefix for the output file names. Defaults to `"df"`.
#'
#' @return A list containing the full DESeq2 data (`deseq_data`) and the top selected genes for the volcano plot (`top_genes_vp`).
#' @export
#'
#' @examples
#' # Example usage
#' result <- prep_volcano_plot_data(
#'   deseq_data = deseq_df, 
#'   lfc = 1, 
#'   pval = 0.01, 
#'   savename = "volcano_data"
#' )
prep_volcano_plot_data <- function(deseq_data, lfc = 1, pval = 0.01, pval_column = "padj", num_up_genes = 12, num_down_genes = 10, pca = TRUE, pca_columns = c("log2FoldChange", "signed_logpvalue"), rank_column = "signed_logpvalue",
                                   save_files = TRUE, include_genes = NULL, exclude_genes = NULL, savename = "df") 
{
  
  library(dplyr) 
  
  # Prep data for volcano plot
  deseq_data <- deseq_data %>%
    filter(!is.na(padj)) %>%
    mutate(color = case_when(
      ((log2FoldChange >= lfc) & (get(pval_column) <= pval & sign(log2FoldChange) > 0)) ~ "Up",
      ((log2FoldChange <= -lfc) & (get(pval_column) <= pval & sign(log2FoldChange) < 0)) ~ "Down",
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
#' Generate Volcano Plot for DESeq2 Results
#'
#' This function creates a volcano plot from preprocessed DESeq2 data, highlighting upregulated, downregulated, 
#' and non-significant genes based on log2 fold change (LFC) and p-value thresholds. It allows for flexible customization 
#' of the plot appearance and can optionally save the plot as a PNG file.
#'
#' @param prepped_deseq_data A data frame containing DESeq2 results, including columns for log2 fold change and p-values.
#' @param top_genes_vp A data frame containing the top upregulated and downregulated genes to label on the plot.
#' @param data A data frame containing the raw DESeq2 data.
#' @param lfc Numeric, the log2 fold change threshold for significance. Defaults to 1.
#' @param pval_column A string specifying the column name for p-values in `prepped_deseq_data`. Defaults to `"padj"`.
#' @param pval_threshold Numeric, the p-value threshold for significance. Defaults to 0.01.
#' @param point_shape Numeric, the shape of points in the plot. Defaults to 21.
#' @param point_size Numeric, the size of points in the plot. Defaults to 3.
#' @param point_alpha Numeric, the transparency level of points in the plot. Defaults to 0.5.
#' @param color_palette A named character vector specifying colors for upregulated, downregulated, and non-significant genes. 
#' Defaults to `c("Upregulated" = "#e5534b", "Downregulated" = "#4f81c7", "No Change" = "#b0b0b0")`.
#' @param xlab A string specifying the label for the x-axis. Defaults to `"Log2 Fold Change"`.
#' @param ylab A string specifying the label for the y-axis. Defaults to `"-Log10(padj)"`.
#' @param legend_position The position of the legend on the plot. Defaults to `"top"`.
#' @param text_size Numeric, the size of the text elements in the plot. Defaults to 13.
#' @param label_size Numeric, the size of the labels for top genes. Defaults to 3.
#' @param x_limits A numeric vector specifying the limits of the x-axis. Defaults to `c(-4, 4)`.
#' @param plot_height Numeric, the height of the saved plot in inches. Defaults to 5.5.
#' @param plot_width Numeric, the width of the saved plot in inches. Defaults to 6.
#' @param plot_dpi Numeric, the resolution of the saved plot in dots per inch (DPI). Defaults to 600.
#' @param save Logical, whether to save the plot as a PNG file. Defaults to `TRUE`.
#' @param savename A string specifying the file name prefix for the saved plot. Defaults to `"plot"`.
#'
#' @return A ggplot object representing the volcano plot.
#' @export
#'
#' @examples
#' # Generate and save a volcano plot from DESeq2 data
#' volcano_plot <- create_volcano_plot(
#'   prepped_deseq_data = prepped_data, 
#'   top_genes_vp = top_genes, 
#'   lfc = 1, 
#'   pval_threshold = 0.01, 
#'   savename = "my_volcano_plot"
#' )
create_volcano_plot <- function(prepped_deseq_data, top_genes_vp, data, lfc = 1, pval_column = "padj", pval_threshold = 0.01, point_shape = 21, point_size = 3, point_alpha = 0.5, 
                                color_palette = c("Upregulated" = "#e5534b", "Downregulated" = "#4f81c7", "No Change" = "#b0b0b0"), xlab = "Log2 Fold Change", ylab = "-Log10(padj)", 
                                legend_position = "top", text_size = 13, label_size = 3, x_limits = c(-4, 4), plot_height = 5.5, plot_width = 6, plot_dpi = 600, save = T, savename = "plot") {
  
  library(ggplot2)
  volcano_plot <- ggplot(prepped_deseq_data, aes(x = log2FoldChange, y = -log10(get(pval_column)), color = color)) +
    geom_point(aes(fill = color), shape = point_shape, size = point_size, stroke = lfc, alpha = point_alpha) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = xlab, 
      y = ylab,
      color = " ",
      fill = " ") +
    theme_classic() + 
    theme(
      legend.position = legend_position,  
      text = element_text(size = text_size, face = "bold")) +  
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "darkgrey") + 
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "darkgrey") + 
    scale_x_continuous(limits = x_limits) +
    ggrepel::geom_label_repel(
      data = top_genes_vp, 
      aes(label = gene), 
      size = label_size, 
      box.padding = unit(0.35, "lines"),  
      point.padding = unit(0.5, "lines")  
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
#' Generate RRHO Density Plot for Pathway Rankings
#'
#' This function generates a Rank-Rank Hypergeometric Overlap (RRHO) density plot comparing pathway rankings between two datasets.
#' It computes and plots the density of shared pathways based on normalized enrichment scores (NES) and their rankings.
#'
#' @param file1 A string specifying the file path of the first dataset containing pathway rankings.
#' @param file2 A string specifying the file path of the second dataset containing pathway rankings.
#' @param pathway_col A string specifying the column name for pathways in both datasets.
#' @param filt_pathways A vector of pathways to filter the data. Defaults to NULL (no filtering).
#' @param nes_col A string specifying the column name for normalized enrichment scores (NES).
#' @param output_file A string specifying the file path to save the plot. Defaults to NULL.
#' @param plot_title A string specifying the title of the plot. Defaults to `"RRHO Density Plot"`.
#' @param plot_tag A string specifying a tag to display on the plot. Defaults to an empty string.
#' @param x_labels A vector specifying labels for the x-axis. Defaults to `c("Label1", "Label2")`.
#' @param y_labels A vector specifying labels for the y-axis. Defaults to `c("Label3", "Label4")`.
#' @param axis.text.size Numeric, the size of the axis text. Defaults to 15.
#' @param legend.position The position of the legend on the plot. Defaults to `"right"`.
#' @param x_labels_breaks A numeric vector specifying the breaks for the x-axis labels. Defaults to `c(1650, 5200)`.
#' @param y_labels_breaks A numeric vector specifying the breaks for the y-axis labels. Defaults to `c(1650, 5200)`.
#' @param axis.title.size Numeric, the size of the axis titles. Defaults to 17.
#' @param x_axis_label A string specifying the label for the x-axis. Defaults to an empty string.
#' @param y_axis_label A string specifying the label for the y-axis. Defaults to an empty string.
#' @param width Numeric, the width of the saved plot. Defaults to 7.25.
#' @param height Numeric, the height of the saved plot. Defaults to 6.
#' @param save Logical, whether to save the plot as a PNG file. Defaults to `FALSE`.
#'
#' @return A ggplot object representing the RRHO density plot.
#' @export
#'
#' @examples
#' # Generate an RRHO density plot between two pathway ranking datasets
#' rrho_plot <- create_rrho_density_plot(
#'   file1 = "data/rankings_1.txt", 
#'   file2 = "data/rankings_2.txt", 
#'   pathway_col = "pathway", 
#'   nes_col = "NES",
#'   output_file = "output/rrho_density_plot.png",
#'   plot_title = "RRHO Density Plot",
#'   x_labels = c("Upregulated in Condition A", "Downregulated in Condition A"),
#'   y_labels = c("Upregulated in Condition B", "Downregulated in Condition B")
#' )
create_rrho_density_plot <- function(file1, file2, pathway_col, filt_pathways = NULL, nes_col, output_file, plot_title = "RRHO Density Plot", plot_tag = "",
                                     x_labels = c("Label1", "Label2"), y_labels = c("Label3", "Label4"), axis.text.size = 15, legend.position = "right",
                                     x_labels_breaks = c(1650, 5200), y_labels_breaks = c(1650, 5200), axis.title.size = 17, x_axis_label = "", y_axis_label = "", 
                                     width = 7.25, height = 6, save = F) {
  
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
  labs(tag = plot_tag)
  
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
