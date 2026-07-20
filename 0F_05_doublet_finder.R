# DoubletFinder to identify and filter out doublets
options(future.globals.maxSize = 4 * 1024^3)

set.seed(1234)

#pK Identification (no ground-truth)
sweep.res.list <- as.list(NULL)
sweep.stats <- as.list(NULL)
bcmvn <- as.list(NULL)
bcmvn_pK_value <- as.vector(NULL)
p_bcmvn <- as.list(NULL)

seurat_SCT <- seurat_object
DefaultAssay(seurat_SCT) <- "SCT"

set.seed(1234)
suppressWarnings(sweep.res.list <- paramSweep(seurat_SCT, PCs = 1:20, sct = TRUE, num.cores = 1))

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
  
bcmvn_BC_value <- max(bcmvn$BCmetric)
bcmvn_pK_value <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric == bcmvn_BC_value,]$pK))


pdf(pdf_DF)

  p <- qplot(as.numeric(as.vector(bcmvn$pK)), bcmvn$BCmetric, geom = "line") +
    geom_vline(xintercept = bcmvn_pK_value, color = "red", linetype = 2) +
    annotate("text",
             x = if (bcmvn_pK_value < 0.2) {
               bcmvn_pK_value + 0.05
             } else {
               bcmvn_pK_value - 0.05
             },
             y = 100,
             label = paste0("Optimal pK = ",bcmvn_pK_value),
             color = "red") +
    labs(x = "pK", y = expression(BC[mvn]))
  print(p)

dev.off()

#Homotypic Doublet Proportion Estimate
annotations <- as.list(NULL)
homotypic.prop <- as.vector(NULL)
nExp_poi <- as.vector(NULL)
nExp_poi.adj <- as.vector(NULL)


annotations <- as.character(Idents(seurat_SCT))
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*length(row.names(seurat_SCT@meta.data)))  ## Assuming 5% doublet formation rate - estimated by 10x
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


rm(annotations,homotypic.prop)

set.seed(1234)
suppressWarnings(seurat_SCT <- doubletFinder_Zuz(seurat_SCT, PCs = 1:10, pN = 0.25, pK = bcmvn_pK_value, nExp = nExp_poi,
                                                       reuse.pANN = FALSE, sct = TRUE))
#--
set.seed(1234)
suppressWarnings(seurat_SCT <- doubletFinder_Zuz(seurat_SCT, PCs = 1:10, pN = 0.25, pK = bcmvn_pK_value, nExp = nExp_poi.adj,
                                                       reuse.pANN = paste0("pANN_0.25_",bcmvn_pK_value,"_",nExp_poi), sct = TRUE))


seurat_doublets_all <- data.frame(NULL)


seurat_doublets <- seurat_SCT@meta.data[,colnames(seurat_SCT@meta.data) == paste0("DF.classifications_0.25_",bcmvn_pK_value,"_",nExp_poi) |
                                                colnames(seurat_SCT@meta.data) == paste0("DF.classifications_0.25_",bcmvn_pK_value,"_",nExp_poi.adj)]
colnames(seurat_doublets) <- c("Low Confidence","High Confidence")
seurat_doublets_all <- rbind(seurat_doublets_all,seurat_doublets)


seurat_object$Doublet <- c(apply(seurat_doublets_all,1,
                                     function(x){
                                       if(x[1]==x[2])
                                         if(x[1]=="Singlet")
                                           "Singlet"
                                       else "Doublet - High Confidence"
                                       else "Doublet - Low Confidence"}))

gc()
rm(bcmvn,bcmvn_BC_value,bcmvn_pK_value,i,nExp_poi,nExp_poi.adj,p,p_bcmvn,seurat_doublets,seurat_doublets_all,seurat_SCT,sweep.res.list,sweep.stats)

