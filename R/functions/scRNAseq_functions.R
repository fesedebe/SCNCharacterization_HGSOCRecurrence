library(Seurat)
library(ggplot2)
library(dplyr)

#run_integration_with_norm-----
#' Integrate Datasets with Normalization and Optional Elbow Plot
#'
#' This function splits a Seurat object by a specified variable (e.g., sample or condition), 
#' normalizes the datasets using either SCT or LogNormalize, optionally performs integration 
#' using RPCA or Harmony, and returns the integrated Seurat object. The function supports 
#' customized variable feature selection, regression of unwanted variables, and generation of 
#' an elbow plot to help determine the number of principal components to retain.
#'
#' @param sobj A Seurat object to be split, normalized, and integrated.
#' @param split.by A column in the Seurat object metadata to use for splitting the object. Defaults to "sample".
#' @param norm.method The normalization method to use. Can be either "SCT" or "LogNormalize". Defaults to "SCT".
#' @param vars.to.regress Variables to regress out during normalization (e.g., "percent.mt", "S.Score", "G2M.Score"). Defaults to NULL.
#' @param selection.method Method for variable feature selection when using LogNormalize. Defaults to "vst".
#' @param vst.flavor The flavor of the SCTransform algorithm to use. Defaults to "v2".
#' @param return.only.var.genes Logical flag indicating whether to return only variable genes during SCT normalization. Defaults to FALSE.
#' @param norm.only Logical flag. If TRUE, the function performs only normalization and skips the integration step. Defaults to FALSE.
#' @param integration.method The method to use for integration, either "rpca" or "harmony". Defaults to "rpca".
#' @param nfeatures The number of features to select for integration. Defaults to 2000.
#' @param elbow_plot Logical flag. If TRUE, the function generates an elbow plot for each dataset to help select the number of principal components to retain. Defaults to FALSE.
#' @param save Logical flag. If TRUE, the function saves the integrated object as an RDS file. Defaults to TRUE.
#' @param savename The file name (without extension) to use when saving the integrated Seurat object. Only relevant if `save = TRUE`.
#'
#' @return Returns an integrated Seurat object.
#' @export
#'
#' @examples
#' # Example usage
#' integrated_obj <- run_integration_with_norm(sobj, split.by = "sample", norm.method = "SCT", integration.method = "rpca", nfeatures = 2000, elbow_plot = TRUE, save = TRUE, savename = "my_integration")
#' 
run_integration_with_norm <- function(sobj, split.by = "sample", norm.method = c("SCT", "LogNormalize"), vars.to.regress = NULL, 
                                      selection.method = "vst", vst.flavor = "v2", return.only.var.genes = F,
                                      norm.only = F, integration.method = c("rpca", "harmony"), nfeatures = 2000, 
                                      elbow_plot = F, save = T, savename){
  
  #split dataset into a list (of independent datasets or samples)
  sobj.list = SplitObject(
    object = sobj,
    split.by = split.by
  )
  
  # normalize and identify variable features for each dataset independently
  if(norm.method == "SCT"){
    sobj.list <- lapply(
      X = sobj.list,
      FUN = function(x) {
        x <- SCTransform(
          x, 
          vars.to.regress = vars.to.regress, #c("percent.mt", "S.Score", "G2M.Score"),
          vst.flavor = vst.flavor,
          return.only.var.genes = return.only.var.genes,
          verbose = FALSE
        )})}
  #Note: SCTransform inherently/automatically selects highly variable features as part of its processing, so there's no need for a separate FindVariableFeatures step.
  
  if(norm.method == "LogNormalize"){
    sobj.list <- lapply(
      X = sobj.list, 
      FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(
          x, 
          selection.method = selection.method, 
          nfeatures = nfeatures
        )})}
  
  # select features that are repeatedly variable across datasets for integration 
  features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = nfeatures)
  
  if(norm.only | integration.method == "harmony"){
    merged_sobj <- merge(
      x = sobj.list[[1]],
      y = sobj.list[2:length(sobj.list)],
      merge.data = TRUE)
    
    VariableFeatures(merged_sobj) <- features
    if(norm.only) {return(merged_sobj)}
    
  }
  
  if(integration.method == "rpca"){
    # perform scaling and run PCA on each dataset using selected variable features
    if(norm.method == "SCT"){
      sobj.list <- PrepSCTIntegration(object.list = sobj.list, anchor.features = features)
      normalization.method = "SCT"
    }
    if(norm.method == "LogNormalize"){
      sobj.list <- lapply(X = sobj.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
      })
      normalization.method = "LogNormalize"
    }
    sobj.list <- lapply(X = sobj.list, FUN = function(x) {
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    
    if(elbow_plot){
      elbow_plots <- list()
      for(i in seq_along(sobj.list)) {
        elbow_plot <- ElbowPlot(sobj.list[[i]], ndims = 50)
        elbow_plots[[i]] <- elbow_plot
      }
      combined_plot <- patchwork::wrap_plots(elbow_plots, ncol = 2)
      combined_plot
    }
    
    # Identify & use these anchors to integrate the datasets together
    uc.anchors <- FindIntegrationAnchors(
      object.list = sobj.list, 
      anchor.features = features, 
      normalization.method = normalization.method,
      reduction = integration.method
    )
    
    # this command creates an 'integrated' data assay
    sobj.integ.srpca <- IntegrateData(
      anchorset = uc.anchors,
      normalization.method = normalization.method
    )
    
    # specify that we will perform downstream analysis on the corrected data 
    # note that the original unmodified data still resides in the 'RNA' assay
    DefaultAssay(sobj.integ.srpca) <- "integrated"
    
    if(save){
      saveRDS(
        sobj.integ.srpca, 
        file = paste0(savename, "_integ.srpca.rds")
      )
    }
    return(sobj.integ.srpca)
  }
  
  if(integration.method == "harmony"){
    merged_sobj <- merged_sobj %>%
      RunPCA(., verbose = F) %>%
      harmony::RunHarmony(., group.by.vars = split.by)
    
    if(save){
      saveRDS(
        merged_sobj, 
        file = paste0(savename, "_integ.harm.rds")
      )
    }
    return(merged_sobj)
  }
}
#
#run_leiden_clustering-----
#' Perform Leiden Clustering on an Integrated Seurat Object
#'
#' This function performs Leiden clustering on an integrated Seurat object. It can normalize the data, 
#' run dimensionality reduction (PCA, UMAP, t-SNE), find neighbors, and identify clusters using the 
#' Leiden algorithm. The function also provides options to save the processed object.
#'
#' @param seurat_obj A Seurat object on which clustering will be performed.
#' @param norm_method The normalization method used in the Seurat object. Can be "SCT" or "LogNormalize". Defaults to "SCT".
#' @param resolution The resolution parameter for the clustering. Higher values lead to more clusters. Defaults to 0.8.
#' @param run_pca Logical, whether to run PCA if it has not already been done. Defaults to TRUE.
#' @param reduction The reduction method to use. Can be either "pca" or "harmony". Defaults to "pca".
#' @param n_pcs The number of principal components to use for dimensionality reduction. Defaults to 30.
#' @param algorithm The clustering algorithm to use. Set to 4 for the Leiden algorithm. Defaults to 4.
#' @param run_umap Logical, whether to run UMAP for visualization. Defaults to TRUE.
#' @param run_tsne Logical, whether to run t-SNE for visualization. Defaults to FALSE.
#' @param k_param The number of neighbors to consider when constructing the graph. Defaults to 20.
#' @param verbose Logical, whether to print progress and diagnostic messages. Defaults to TRUE.
#' @param save Logical, whether to save the Seurat object after clustering. Defaults to TRUE.
#' @param savename The base name of the file when saving the Seurat object (without file extension). Defaults to "obj".
#'
#' @return Returns a Seurat object with Leiden clustering results and optional dimensionality reduction (PCA, UMAP, t-SNE).
#' @export
#'
#' @examples
#' # Run Leiden clustering on a Seurat object with default parameters
#' seurat_obj <- run_leiden_clustering(seurat_obj, norm_method = "SCT", resolution = 0.8, run_pca = TRUE)
#'
#' # Run Leiden clustering and save the object with a custom name
#' seurat_obj <- run_leiden_clustering(seurat_obj, resolution = 1.0, savename = "custom_clustering")
run_leiden_clustering <- function(seurat_obj, norm_method = "SCT", resolution = 0.8, run_pca = TRUE, reduction = c("pca", "harmony"),
                                  n_pcs = 30, algorithm = 4, run_umap = TRUE, run_tsne = FALSE, 
                                  k_param = 20, verbose = TRUE, save = T, savename = "obj") {
  
  # Normalize the data if Log Normalization was used
  if(norm_method == "LogNormalize") {
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  }
  
  # Run PCA if not already done
  if (run_pca) {
    if (is.null(seurat_obj@reductions$pca)) {
      if (verbose) cat("PCA not found, running PCA...\n")
      seurat_obj <- RunPCA(
        seurat_obj, 
        npcs = n_pcs, 
        verbose = verbose
      )
    } else {
      if (verbose) cat("PCA already exists, skipping...\n")
    }
  }
  
  # Find Neighbors & Perform (Leiden) Clustering
  if (verbose) cat("Finding neighbors & clusters...\n")
  seurat_obj <- FindNeighbors(
    seurat_obj, 
    reduction = reduction, 
    dims = 1:n_pcs, 
    k.param = k_param, 
    verbose = verbose
  )
  seurat_obj <- FindClusters(
    seurat_obj, 
    resolution = resolution, 
    algorithm = algorithm, # Algorithm 4 is Leiden in Seurat
    verbose = verbose,
    method = "igraph"
  )  
  
  # Run UMAP & TSNE if requested
  if (run_umap) {
    seurat_obj <- RunUMAP(
      seurat_obj, 
      reduction = reduction, 
      dims = 1:n_pcs, 
      verbose = verbose
    )
  }
  
  if (run_tsne) {
    seurat_obj <- RunTSNE(seurat_obj, reduction = reduction, dims = 1:n_pcs, verbose = verbose)
  }
  
  if(save){
    saveRDS(
      seurat_obj, 
      file = paste0(savename, "_leidenclust.rds")
    )
  }
  
  return(seurat_obj)
}
#
#runSeuratWorkflow------
#' Run Standard Seurat Workflow for Clustering and Visualization
#'
#' This function runs the standard Seurat workflow for data normalization, dimensionality reduction, 
#' clustering, and visualization. It includes optional saving of the processed Seurat object.
#'
#' @param sobj A Seurat object to be processed.
#' @param norm.method The normalization method to use. Can be "SCT" or "LogNormalize". Defaults to "SCT".
#' @param npcs The number of principal components to compute. Defaults to 30.
#' @param resolution The resolution parameter for clustering. Higher values yield more clusters. Defaults to 0.5.
#' @param dims The dimensions to use for clustering and visualization. Defaults to 1:30.
#' @param save Logical, whether to save the processed Seurat object as an RDS file. Defaults to TRUE.
#' @param savename The file name (without extension) to use when saving the Seurat object. Only relevant if `save = TRUE`.
#'
#' @return Returns the processed Seurat object with PCA, UMAP, t-SNE, neighbors, and clusters.
#' @export
#'
#' @examples
#' # Run the standard workflow on a Seurat object with default parameters
#' processed_obj <- runSeuratWorkflow(sobj, norm.method = "SCT", npcs = 30, resolution = 0.5)
#'
#' # Run the workflow and save the processed object with a custom name
#' processed_obj <- runSeuratWorkflow(sobj, norm.method = "LogNormalize", npcs = 20, savename = "my_seurat")
runSeuratWorkflow <- function(sobj, norm.method = "SCT", npcs = 30, resolution = 0.5, dims = 1:30, save = T, savename) {
  # Normalize the data if Log Normalization was used
  if(norm.method == "LogNormalize") {
    sobj <- ScaleData(sobj, verbose = FALSE)
  }
  
  # Run PCA
  sobj <- RunPCA(sobj, npcs = npcs, verbose = FALSE)
  
  # Run UMAP
  sobj <- RunUMAP(sobj, reduction = "pca", dims = dims)
  
  # Run t-SNE
  sobj <- RunTSNE(sobj, reduction = "pca", dims = dims)
  
  # Find neighbors
  sobj <- FindNeighbors(sobj, dims = dims)
  
  # Find clusters
  sobj <- FindClusters(sobj, resolution = resolution)
  
  if(save){
    saveRDS(
      sobj, 
      file = paste0(savename, "_clust.rds")
    )
  }
  
  # Return the Seurat object with clustering and visualization results
  return(sobj)
}

#
#plot_signature_scores------
#' Plot Signature Scores by Seurat Clusters with Customizable Options
#'
#' This function generates a jitter plot of signature scores for specified clusters in a Seurat object. 
#' It allows customizable coloring, faceting, and crossbar summaries. The function can visualize specific 
#' clusters of interest and compare them against other clusters.
#'
#' @param seurat_obj A Seurat object containing the metadata and expression data.
#' @param clusters_of_interest A numeric vector of clusters to be labeled as a high group. Defaults to cluster 9.
#' @param high_group_label A string label for the high group. Defaults to "HighSCN".
#' @param other_group_label A string label for the other group. Defaults to "Other".
#' @param score_columns A vector of column names representing signature scores to plot. Defaults to `c("SCN_Balanis1", "RBLossMalorni1")`.
#' @param cluster_column The metadata column that contains the cluster identifiers. Defaults to `"seurat_clusters"`.
#' @param color_column The metadata column to use for coloring points in the plot. Defaults to `"ClusterGroup"`.
#' @param facet_labels A named vector for customizing the facet labels in the plot. Defaults to `c(RBLoss_Malorni1 = "RB Loss", SCN_Balanis1 = "SCN")`.
#' @param palette A named vector for custom color values. Defaults to `c("Other" = "#AEC6CF", "HighSCN" = "#FFB347")`.
#' @param title The title of the plot. Defaults to `"Signature Scores by Seurat Cluster"`.
#' @param jitter_width The width of the jitter points for better separation. Defaults to 0.2.
#' @param jitter_alpha Transparency level of the jitter points. Defaults to 0.4.
#' @param jitter_size The size of the jitter points. Defaults to 3.
#' @param crossbar_color The color of the median crossbars. Defaults to `"black"`.
#' @param crossbar_width The width of the crossbars. Defaults to 0.5.
#' @param crossbar_linewidth The linewidth of the crossbars. Defaults to 0.5.
#' @param plot Logical flag to control whether the plot is printed. Defaults to TRUE.
#'
#' @return Returns the modified Seurat object with the new cluster group labels, if applicable.
#' @export
#'
#' @examples
#' # Example usage: Plot signature scores for cluster 9
#' plot_signature_scores(seurat_obj = seurat_obj, clusters_of_interest = 9)
#'
#' # Customize palette and plot signature scores
#' plot_signature_scores(seurat_obj = seurat_obj, palette = c("Other" = "#B3E2CD", "HighSCN" = "#FDCDAC"))
plot_signature_scores <- function(seurat_obj, clusters_of_interest = 9, high_group_label = "HighSCN", 
                                  other_group_label = "Other", score_columns = c("SCN_Balanis1", "RBLossMalorni1"), cluster_column = "seurat_clusters", color_column = "ClusterGroup",
                                  facet_labels = c(RBLoss_Malorni1 = "RB Loss", SCN_Balanis1 = "SCN"), palette = c("Other" = "#AEC6CF", "HighSCN" = "#FFB347"),
                                  title = "Signature Scores by Seurat Cluster", jitter_width = 0.2, jitter_alpha = 0.4, 
                                  jitter_size = 3, crossbar_color = "black", crossbar_width = 0.5, crossbar_linewidth = 0.5, plot = T) {
  
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(ClusterGroup = ifelse(!!sym(cluster_column) %in% clusters_of_interest, high_group_label, other_group_label))
  
  # Prepare data for plotting
  plot_data <- seurat_obj@meta.data %>%
    select(all_of(cluster_column), all_of(score_columns), all_of(color_column)) %>%
    tidyr::pivot_longer(cols = all_of(score_columns), names_to = "ScoreType", values_to = "Score") %>%
    mutate(ScoreType = factor(ScoreType, levels = score_columns))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = as.factor(!!sym(cluster_column)), y = Score, color = !!sym(color_column))) +
    geom_jitter(width = jitter_width, alpha = jitter_alpha, size = jitter_size) +
    stat_summary(fun = median, geom = "crossbar", width = crossbar_width, color = crossbar_color, linewidth = crossbar_linewidth) +
    facet_wrap(~ ScoreType, scales = "free_y", ncol = 1, labeller = as_labeller(facet_labels)) +
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
    labs(x = "Seurat Clusters", y = "Signature Score", title = title) +
    scale_color_manual(values = palette) +
    guides(color = guide_legend(override.aes = list(size = 9, alpha = 1)))
  
  if (plot){
    print(p)
  }
  
  return(seurat_obj)
}

#
#find_markers------
#' Find Differential Markers in a Seurat Object with Custom Options
#'
#' This function identifies differential markers in a Seurat object based on a specified identity column, 
#' using various testing methods and assays. It can also filter the results based on log fold change 
#' and adjusted p-value, and optionally save the results to a file.
#'
#' @param seurat_obj A Seurat object from which to identify markers.
#' @param ident_col The identity column in the Seurat object metadata to use for finding markers. Defaults to "ClusterGroup".
#' @param ident_1 The identity class to compare (e.g., cluster or group). Defaults to "HighSCN".
#' @param test_use The statistical test to use for marker identification. Defaults to "MAST".
#' @param assay_use The assay to use for finding markers (e.g., "SCT" for SCTransform). Defaults to "SCT".
#' @param output_dir The directory where results will be saved. Defaults to "./data/".
#' @param output_prefix The prefix for the output file names. Defaults to "markers".
#' @param min_logfc Minimum log fold change to use for filtering markers. Defaults to 1.
#' @param p_val_cutoff Adjusted p-value cutoff for filtering markers. Defaults to 0.01.
#' @param find_all Logical, whether to find markers for all identity classes. Defaults to FALSE.
#' @param min_pct Minimum percentage of cells expressing the feature (used only when `find_all = TRUE`). Defaults to 0.25.
#' @param write Logical, whether to write the full and filtered marker results to a file. Defaults to TRUE.
#'
#' @return Returns a dataframe of all identified markers.
#' @export
#'
#' @examples
#' # Find markers for "HighSCN" using MAST test and save the results
#' markers <- find_markers(seurat_obj, ident_col = "ClusterGroup", ident_1 = "HighSCN", test_use = "MAST")
#'
#' # Find markers for all clusters and save them
#' all_markers <- find_markers(seurat_obj, ident_col = "seurat_clusters", find_all = TRUE)
find_markers <- function(seurat_obj, ident_col = "ClusterGroup",  ident_1 = "HighSCN",  test_use = "MAST", assay_use = "SCT", 
                         output_dir = "./data/", output_prefix = "markers", min_logfc = 1, p_val_cutoff = 0.01, find_all = FALSE, min_pct = 0.25, write = T) {
  
  # Set identity class
  Idents(seurat_obj) <- ident_col
  
  # Prepare markers based on the condition
  if (!find_all) {
    markers <- FindMarkers(
      object = PrepSCTFindMarkers(seurat_obj), 
      ident.1 = ident_1, 
      test.use = test_use,
      assay = assay_use
    )
  } else {
    markers <- FindAllMarkers(
      object = PrepSCTFindMarkers(seurat_obj), 
      test.use = test_use,
      min.pct = min_pct,
      assay = assay_use
    )
  }
  
  # Modify markers dataframe with log-p values
  markers_full <- markers %>%
    mutate(p_val = ifelse(p_val == 0, .Machine$double.xmin, p_val),
           p_val_adj = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
           gene = rownames(.)) %>%
    mutate(signed_log_p_val = sign(avg_log2FC) * -log10(p_val),
           signed_log_p_val_adj = sign(avg_log2FC) * -log10(p_val_adj))
  
  # Save the full markers data
  if (write){
    full_output_path <- paste0(output_dir, output_prefix, "_MAST_padjfull.txt")
    write.table(
      x = markers_full,
      file = full_output_path,
      sep = "\t" , 
      quote = FALSE, row.names = FALSE
    )}
  
  # Apply filters to markers
  markers_filt <- markers_full %>%
    filter(p_val_adj < p_val_cutoff,
           abs(avg_log2FC) >= min_logfc)
  
  if (write){
    # Save the filtered markers data
    filt_output_path <- paste0(output_dir, output_prefix, "_MAST_padjfilt.txt")
    write.table(
      x = markers_filt,
      file = filt_output_path,
      sep = "\t", 
      quote = FALSE, row.names = FALSE
    )}
  
  return(markers_full)
}
#


