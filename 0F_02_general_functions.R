##########################################################
# Conversion of gene symbols to Ensembl IDs and vice versa
##########################################################
Symbol_to_ENSEMBL <- function(gi){
  genelist <- NULL
  for(i in 1:length(gi)){
    genelist <- c(genelist, gene_symbol[gene_symbol$gene_name %in% gi[i],]$gene_id)
  }
  return(genelist)
}

ENSEMBL_to_Symbol <- function(gi){
  genelist <- NULL
  for(i in 1:length(gi)){
    genelist <- c(genelist,gene_symbol[gene_symbol$gene_id %in% gi[i],]$gene_name)#Export ENSEMBL ids. If there is more genes with same symbol, it exports more ENSEMBL ids. If there is no, it exports nothing. You need to be carefull. If your gene symbol or ensembl id is in another column, or columns have different names, you need to change this [genenames$V2 %in% gi[i],1]
  }
  return(genelist) #Return created list
}

########################
# Themes for gg plotting
########################
theme_mk <- theme_bw() + theme(
  text = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.text = element_text(size = 6, margin = margin(t = 0, r = 0, b = -1, l = 1)),  
  strip.text = element_text(size = 6),
  plot.title = element_blank(),
  legend.key.size = unit(2.5, "mm"),
  line = element_line(linewidth = 0.3),
  legend.background = element_blank(),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Margins around the legend
  #legend.key.spacing.y = unit(0.0, "mm"),
  #legend.key.spacing.x = unit(0.1, "mm"),
  legend.title = element_text(size = 6),
  panel.border = element_blank(),  # Remove all borders first
  axis.line = element_line(color = "black", linewidth = 0.3)  # Add bottom and left borders
)

theme_mk_title <- theme_bw() + theme(
  text = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.text = element_text(size = 6, margin = margin(t = 0, r = 0, b = -1, l = 1)),  
  strip.text = element_text(size = 6),
  plot.title = element_text(hjust = 0.5, size = 6, vjust = 1),
  legend.key.size = unit(1.5, "mm"),
  line = element_line(linewidth = 0.3),
  legend.background = element_blank(),
  legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Margins around the legend
  #legend.key.spacing.y = unit(0.0, "mm"),
  #legend.key.spacing.x = unit(0.1, "mm"),
  legend.title = element_text(size = 6),
  panel.border = element_blank(),  # Remove all borders first
  axis.line = element_line(color = "black", linewidth = 0.3)  # Add bottom and left borders
)

# remove grid from ggplot backgrounds
remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# remove axes from ggplot backgrounds
remove_axes <- theme(
  axis.title.x = element_blank(), # Remove x-axis title
  axis.title.y = element_blank(), # Remove y-axis title
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.text.y = element_blank(),  # Remove y-axis text
  axis.ticks.x = element_blank(), # Remove x-axis ticks
  axis.ticks.y = element_blank(), # Remove y-axis ticks
  axis.line.x = element_blank(),  # Remove x-axis line
  axis.line.y = element_blank()   # Remove y-axis line
)

theme_umap <- theme_classic() + 
  theme(
    text = element_text(size = 6),
    legend.position = "bottom",
    legend.text = element_text(
      size = 6,
      margin = margin(t = 0, r = 0, b = -1, l = 1)
    ),
    legend.spacing.y = unit(0, "mm"),       
    legend.spacing.x = unit(1, "mm"),       
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(0, 0, 2, 0),
    strip.text = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, size = 7.3, vjust = 0),
    legend.key.size = unit(2.5, "mm"),
    legend.background = element_blank()
  )

# remove axes from ggplot backgrounds
remove_axes <- theme(
  axis.title.x = element_blank(), # Remove x-axis title
  axis.title.y = element_blank(), # Remove y-axis title
  axis.text.x = element_blank(),  # Remove x-axis text
  axis.text.y = element_blank(),  # Remove y-axis text
  axis.ticks.x = element_blank(), # Remove x-axis ticks
  axis.ticks.y = element_blank(), # Remove y-axis ticks
  axis.line.x = element_blank(),  # Remove x-axis line
  axis.line.y = element_blank()   # Remove y-axis line
)

################
# Doublet finder
###############


paramSweep_v5 <- function(seu, PCs, assay = "RNA") {
  require(Seurat)
  require(fields)
  
  # Set pN-pK parameter ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  
  sweep.res.list <- list()
  list.ind <- 0
  
  data <- GetAssayData(seu, assay = assay, slot = "counts")
  n.cells <- ncol(data)
  
  if (n.cells > 10000) {
    real.cells <- colnames(data)[sample(n.cells, 10000)]
    data <- data[, real.cells]
  } else {
    real.cells <- colnames(data)
  }
  
  n.real.cells <- length(real.cells)
  
  for (n in seq_along(pN)) {
    message(paste("Creating artificial doublets for pN =", pN[n] * 100, "%"))
    n_doublets <- round(n.real.cells / (1 - pN[n]) - n.real.cells)
    
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste0("X", seq_len(n_doublets))
    
    data_wdoublets <- cbind(data, doublets)
    
    # Create Seurat object with doublets
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    
    # Standard preprocessing
    seu_wdoublets <- NormalizeData(seu_wdoublets)
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets)
    seu_wdoublets <- ScaleData(seu_wdoublets)
    seu_wdoublets <- RunPCA(seu_wdoublets, features = VariableFeatures(seu_wdoublets), npcs = max(PCs))
    
    # Distance matrix in PCA space
    pca.coord <- Embeddings(seu_wdoublets, reduction = "pca")
    pca.coord <- pca.coord[, PCs, drop = FALSE]
    nCells <- nrow(pca.coord)
    
    dist.mat <- fields::rdist(pca.coord)[, 1:n.real.cells]
    
    # Order distances
    for (i in 1:n.real.cells) {
      dist.mat[, i] <- order(dist.mat[, i])
    }
    
    # Trim distance matrix
    ind <- round(nCells * max(pK)) + 5
    dist.mat <- dist.mat[1:ind, ]
    
    # Compute pANN
    message("Computing pANN across all pK...")
    for (k in seq_along(pK)) {
      message(paste("pK =", pK[k]))
      pk.temp <- round(nCells * pK[k])
      pANN <- data.frame(pANN = numeric(n.real.cells))
      rownames(pANN) <- real.cells
      
      for (i in seq_len(n.real.cells)) {
        neighbors <- dist.mat[2:(pk.temp + 1), i]
        pANN$pANN[i] <- sum(neighbors > n.real.cells) / pk.temp
      }
      
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- pANN
    }
  }
  
  # Naming
  name.vec <- unlist(lapply(pN, function(pn) paste0("pN_", pn, "_pK_", pK)))
  names(sweep.res.list) <- name.vec
  
  return(sweep.res.list)
}


##################
# SPATIAL ROTATION
##################

# 1. Function with adjustable rotation
make_transcript_plot <- function(obj, img_name, gene_name, rotation_degrees = 0) {
  # Get coordinates for the specific image
  coords <- GetTissueCoordinates(obj, image = img_name)
  
  # Get the expression data
  exp_data <- FetchData(obj, vars = gene_name)
  
  # Combine into a data frame
  # We assume coords has 2 columns (x and y / imagerow and imagecol)
  df <- cbind(coords, exp_data[rownames(coords), , drop = FALSE])
  colnames(df) <- c("x_orig", "y_orig", "Expression")
  
  # --- MATHEMATICAL ROTATION ---
  # Convert degrees to radians (R functions use radians)
  # Note: To rotate CLOCKWISE, we use a negative angle in the standard formula
  angle_rad <- -rotation_degrees * (pi / 180)
  
  # Center the coordinates around (0,0) before rotating to prevent distortion
  mid_x <- mean(df$x_orig)
  mid_y <- mean(df$y_orig)
  df$x_adj <- df$x_orig - mid_x
  df$y_adj <- df$y_orig - mid_y
  
  # Standard rotation matrix formula
  df$x <- df$x_adj * cos(angle_rad) - df$y_adj * sin(angle_rad)
  df$y <- df$x_adj * sin(angle_rad) + df$y_adj * cos(angle_rad)
  # -----------------------------
  
  # Create the plot using your specific settings
  ggplot(df, aes(x = x, y = y, color = Expression)) +
    geom_point(size = 0.3, stroke = 0) + 
    scale_color_gradientn(
      colours = col_list[["spatial_slices"]],
      limits = c(global_min, global_max),
      breaks = c(global_min, global_max),
      #labels = c("Min", "Max"),
      name = "Exp."
    ) +
    coord_fixed() + 
    theme_void() + 
    ggtitle(img_name) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 6), # Matches your theme_mk_title vibe
      legend.position = "right", 
      legend.text = element_text(size = 6),
      legend.title = element_text(size =6)
    ) +
    guides(
      color = guide_colorbar(
        barwidth = 0.1, 
        barheight = 1.3
      )
    )
}

####################
# Scientific format
###################
add_sci_superscript <- function(x) {
  sapply(x, function(value) {
    numeric_value <- as.numeric(value)
    if (!is.na(numeric_value) && numeric_value != 0) {
      exponent <- floor(log10(abs(numeric_value)))
      base <- signif(numeric_value / (10^exponent), digits = 2)
      
      # Use paste0 with the Unicode multiplication sign to remove gaps
      # The * inside bquote connects elements with zero space
      as.expression(bquote(.(base) * "\u00d7" * 10^.(exponent)))
      
    } else {
      as.character(value)
    }
  }, USE.NAMES = FALSE)
}

##################################
# Violin box plot
#####################################
violin_boxplot <- function(data, feature, y_label) {
  ggplot(data, aes(x = Isolation_pub, y = .data[[feature]], fill = Isolation_pub)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.6, linewidth = 0.3) +  # Violin plot
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.1) +  # Boxplot inside
    scale_fill_manual(values = col_list[["isol_prot"]]) +  # Custom color palette
    labs(y = y_label, x = "Isolation") +
    theme_mk + remove_grid +
    theme(legend.position = "none"
    )
} 

violin_theme <- theme(
  strip.background = element_blank(),
  strip.text = element_text(size = 6),
  axis.title.x = element_blank(),
  # Your specific rotation and alignment settings
  axis.text.x = element_text(angle = 310, vjust = 0.1, hjust = 0.2), 
  axis.title.y = element_text(margin = margin(r = 5)),
  plot.title = element_blank(),
  panel.spacing = unit(1, "lines")
)

