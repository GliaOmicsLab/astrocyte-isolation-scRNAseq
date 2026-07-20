orders[["spatial_samples"]] <- c("Ctrl","1DPI", "3DPI", "7DPI")
orders[["seurat_annotated_all"]] <- c("Ast","Mg","Neu", "Opc", "Ol", "Endo", "Epen", "Peri1", "Peri2", "Vlmc", "Cp", "Oec")
orders[["samples"]] <- c("Ctrl","Tau", "MCAO")
orders[["Isolation_pub"]] <- c( "Miltenyi", "WorthGentle","WorthMech")

#-------------------------------------------------------------------------------
#Public
col_list[["column_plot"]] <- c("#eadbc7")

#Spatial
col_list[["spatial_slices"]] <- c( "#f4eff5", "#763d41")
col_list[["borders"]] <- c( "#2B1700")
col_list[["bar_plot_visium"]] <- c("Aldh1l1+" = "#f6eeef", "Gfap-high+" = "#b1686d")


#scRNA-seq
col_list[["UMAP_all"]] <- c("#fcc94f", "#e2a004", "#91C27B","#fdbcb5","#fb8072","#9E74A5","#C0A5C5","#235656","#598e8e", "#80CCCC","#666666","#B39831","#C18452")
col_list[["UMAP_astro"]] <- c("Hom" = "#9f6600", "Reac" = "#fcc94f")
col_list[["isol_prot"]] <- c("#7ecdba", "#FFCCB4","#a490bb")
col_list[["bar_plot"]] <- c("Aldh1l1+" = "#fff7e6", "Gfap+" = "#fcc94f")
col_list[["gradient_comp"]] <- c("#fef8f6","#9f6600")
col_list[["react_condition"]] <- c("#a3bccf","#6380b5","#f8f2ed")

