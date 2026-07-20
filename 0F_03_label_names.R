##########################
# Unified labels in charts
##########################
gene_labels <- list(
  #genes
  "Aldh1l1+" = bquote(italic("Aldh1l1")^"+"),
  "Gfap+" = bquote(italic("Gfap")^"+"),
  "Gfap_high+" = bquote(italic("Gfap")^"high+"),
  #isolation protocols
  "Isolation_pub" = bquote("Isol. protocol"),
  "Isolation_pub_2l" = bquote("Isol. \nprotocol"),
  #condition
  "Conditions" =  bquote("Pathology"),

  # general abbreviation
  "cell_types" = bquote("Cell type"),
  "cell_types_2l" = bquote("Cell \ntype"),
  "number_of_cells" = bquote("Number of cells"),
  "relat_prop_cells" = bquote("Relative proportion \nof cells (%)"),
  "rel_exp" = bquote("Relative \nexpression"),
  "sc_astro" = bquote("scRNA-seq data (Ast)"),
  "sc_astro_mcao" = bquote("scRNA-seq data (Ast MCAO)"),
  "our_sc_astro_mcao" = bquote("Our scRNA-seq data \n(Ast MCAO)"),
  "avrg_exp" = bquote("Avrg. \nexp."),
  "pub_data" = bquote("Public scRNA-seq data"),
  "pub_spatial" = bquote("Public  spatial RNA-seq data"^Delta),
  "qPCR"= bquote("sc-qPCR data (MCAO)"),
  "mod_score"=bquote("Module score"),
  
  # technical parameters
  "nFeature_RNA" = bquote("Number of genes \nper cell (log10)"),
  "nCount_RNA" = bquote("Number of UMIs \nper cell (log10)"),
  "percent.mt" = bquote("Mitochondrial \ngenes (%)"),
  "percent.rib" = bquote("Ribosomal \ngenes (%)"),

  # validaiton
  "Gfap_per" = bquote(italic("Actb")^"+" * italic("Aldh1l1")^"+" * italic("Gfap")^"+" ~ "cells (%)"),
  "replicate" =  bquote("Tech. \nreplicate"),
  
  #analysis
  "DEG" = bquote("DEG status")
  
)
# Helper function to get the label safely
get_label <- function(gene_name) {
  if (gene_name %in% names(gene_labels)) {
    return(gene_labels[[gene_name]])
  } else {
    # Fallback: if gene not in dictionary, just return the name as italic
    return(bquote(italic(.(gene_name))))
  }
}