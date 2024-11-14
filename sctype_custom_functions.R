# 自定义 gene_sets_prepare 函数，记录非标准基因符号的来源
gene_sets_prepare_custom <- function(path_to_db_file, cell_type) {
  
  cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  cell_markers <- cell_markers[cell_markers$tissueType == cell_type,]
  cell_markers$geneSymbolmore1 <- gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub(" ", "", cell_markers$geneSymbolmore2)
  
  # Initialize a list to collect warning messages
  warning_messages <- list()
  
  # Helper function to standardize gene symbols and capture non-standard symbols
  standardize_gene_symbols <- function(gene_column) {
    sapply(1:nrow(cell_markers), function(i) {
      markers_all <- gsub(" ", "", unlist(strsplit(cell_markers[[gene_column]][i], ",")))
      markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all <- sort(markers_all)
      
      if (length(markers_all) > 0) {
        # 使用 checkGeneSymbols 更正基因符号并捕获未批准的符号
        symbols_checked <- checkGeneSymbols(markers_all)
        
        # 记录非标准符号信息
        non_approved_symbols <- symbols_checked$x[is.na(symbols_checked$Suggested.Symbol)]
        if (length(non_approved_symbols) > 0) {
          warning_messages[[length(warning_messages) + 1]] <- paste(
            "Non-approved gene symbols found in", gene_column, 
            "row", i, ":", paste(non_approved_symbols, collapse = ", ")
          )
        }
        
        # 标准化符号
        standardized_symbols <- ifelse(is.na(symbols_checked$Suggested.Symbol), 
                                       symbols_checked$x,   # 保留原符号
                                       symbols_checked$Suggested.Symbol)  # 替换为标准符号
        paste0(unique(standardized_symbols), collapse = ",")
      } else {
        ""
      }
    })
  }
  
  # 应用基因符号标准化到上调和下调基因集
  cell_markers$geneSymbolmore1 <- standardize_gene_symbols("geneSymbolmore1")
  cell_markers$geneSymbolmore2 <- standardize_gene_symbols("geneSymbolmore2")
  
  cell_markers$geneSymbolmore1 <- gsub("///", ",", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub("///", ",", cell_markers$geneSymbolmore2)
  
  # 输出警告信息
  if (length(warning_messages) > 0) {
    message("Non-standard gene symbols identified:\n", paste(warning_messages, collapse = "\n"))
  }
  
  # 创建基因集列表
  gs <- lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]), ","))))
  names(gs) <- cell_markers$cellName
  gs2 <- lapply(1:nrow(cell_markers), function(j) gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]), ","))))
  names(gs2) <- cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}


# 自定义 run_sctype 函数
run_sctype_custom <- function(seurat_object, known_tissue_type = NULL, custom_marker_file = NULL, plot = FALSE, name = "sctype_classification") {
  db_ = sctype_source()
  if (is.null(seurat_object)) {
    stop("Argument 'seurat_object' is missing")
  }
  if (!inherits(seurat_object, "Seurat")) {
    stop("Argument 'seurat_object' must be a Seurat object")
  }
  if (is.null(custom_marker_file)) {
    custom_marker_file = db_
  }
  if (is.null(known_tissue_type)) {
    tissue_type = auto_detect_tissue_type(path_to_db_file = custom_marker_file, seuratObject = seurat_object, scaled = TRUE, assay = "RNA")
    rownames(tissue_type) = NULL
    tissue_type = tissue_type$tissue[1]
  } else {
    tissue_type = known_tissue_type
  }
  
  # 使用自定义 gene_sets_prepare
  gs_list = gene_sets_prepare_custom(custom_marker_file, tissue_type)
  
  package_type = substr(packageVersion("Seurat"), 1, 1)
  
  if (package_type == 5) {
    print("Running with Seuratv5")
    es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]$scale.data,
                          scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  } else {
    print("Running with Seuratv4")
    es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]@scale.data,
                          scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  }
  
  cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl) {
    es.max.cl = sort(rowSums(es.max[, rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters == cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] = "Unknown"
  seurat_object_res = seurat_object
  seurat_object_res@meta.data[name] = ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type = sctype_scores[sctype_scores$cluster == j, ]
    seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j, name] = as.character(cl_type$type[1])
  }
  if (plot) {
    plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = name)
    print(plot_)
  }
  text_ = paste("New metadata added: ", name)
  print(text_)
  return(seurat_object_res)
}

