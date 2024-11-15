#' Custom Ensemble Visualization Step
#'
#' @description This function generates customized UMAP visualizations and combines them into a grid layout.
#' It includes options for coloring by ensemble output, ScType classifications, and other metadata.
#'
#' @param seurat_object A Seurat object containing the data to visualize.
#' @return A combined ggplot object representing the ensemble visualization.
#' @export
#' @name custom_visualize_ensemble
library(Seurat)
custom_visualize_ensemble <- function(seurat_object) {
  # Define plot colors
  plot_cols <- c("malignant" = "#5560AC", "healthy" = "#E68D3D", "NA" = "#808080")
  sctype_colors <- c(
    "Basophils" = "#E68D3D",
    "CD8+ NKT-like cells" = "#FDD379",
    "Classical Monocytes" = "#5560AC",
    "Erythroid-like and erythroid precursor cells" = "#C6DDED",
    "Memory CD4+ T cells" = "#378000",
    "Naive B cells" = "#B2DF8A",
    "Neutrophils" = "#f36569",
    "Platelets" = "#ff9d9f",
    "Pre-B cells" = "#b48fef",
    "Progenitor cells" = "#e88ae8",
    "Unknown" = "#e07eac"
  )

  # Generate individual plots with custom colors
  p1 <- DimPlot(seurat_object, group.by = "ensemble_output", cols = plot_cols, label.size = 3) +
    ggtitle("ensemble output") +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Check if "sctype_classification" exists and use custom colors if it does
  if ("sctype_classification" %in% colnames(seurat_object@meta.data)) {
    p2 <- DimPlot(seurat_object, group.by = "sctype_classification", label = TRUE, label.size = 4, cols = sctype_colors) +
      ggtitle("scType classification") +
      NoLegend() +
      theme(plot.title = element_text(size = 10, face = "bold"))
  } else {
    print("sctype_classification does not exist. For an automated cell type annotation: 'run_scType'")
    p2 <- NULL
  }

  # Other plots without legends, using the defined colors
  p3 <- DimPlot(seurat_object, group.by = "sctype_malignant_healthy", cols = plot_cols) +
    ggtitle("malignant vs. healthy by scType") +
    NoLegend() +
    theme(plot.title = element_text(size = 10, face = "bold"))

  p4 <- DimPlot(seurat_object, group.by = "copyKat_output", cols = plot_cols) +
    ggtitle("copyKat output") +
    NoLegend() +
    theme(plot.title = element_text(size = 10, face = "bold"))

  p5 <- DimPlot(seurat_object, group.by = "SCEVAN_output", cols = plot_cols) +
    ggtitle("SCEVAN output") +
    NoLegend() +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Combine plots into a grid layout
  combined_plot <- plot_grid(
    p1, plot_grid(p2, p3, nrow = 1), plot_grid(p4, p5, nrow = 1),
    nrow = 3, rel_heights = c(1, 0.8, 0.8),
    labels = c("A", "B", "C"), label_size = 15
  )

  # Return the plot object if you want to view it in R
  return(combined_plot)
}
