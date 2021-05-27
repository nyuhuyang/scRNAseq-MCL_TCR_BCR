########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20210525_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

nrow(df_samples)

#======1.2 load  Seurat =========================
(load(file = "data/MCL_TCR_BCR_20210525.Rda"))

meta.data = object@meta.data
for(i in 1:length(samples)){
    cells <- meta.data$orig.ident %in% samples[i]
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"project"] = df_samples$project[i]
    meta.data[cells,"date"] = df_samples$date[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()

object_list %<>% pblapply(SCTransform,method = "glmGamPoi",vars.to.regress = "percent.mt")
features <- SelectIntegrationFeatures(object.list = object_list)

options(future.globals.maxSize= object.size(object_list)*1.5)
object_list %<>% pblapply(FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
})
anchors <- FindIntegrationAnchors(object.list = object_list,reference = c(1, 2), reduction = "rpca", 
                                  dims = 1:50)
remove(object_list);GC()
# this command creates an 'integrated' data assay
object <- IntegrateData(anchorset = anchors,normalization.method = "SCT", dims = 1:50)
remove(anchors);GC()
save(object, file = "data/MCL_TCR_BCR_20210525.Rda")
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!
    
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard workflow for visualization and clustering
object %<>% NormalizeData()
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 150, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 150))
dev.off()
object %<>% RunUMAP(reduction = "pca", dims = 1:75)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:75))

object %<>% FindNeighbors(reduction = "pca", dims = 1:75)
object %<>% FindClusters(resolution = 0.8)

#=======1.9 summary =======================================
g1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = F,
                 label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
                 no.legend = T,label.size = 4, repel = T, title = "with data Integration",
                 do.print = T, do.return = T)

object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

save(object, file = "data/Lorenzo-LS6_20210408.Rda")
format(object.size(object@assays$RNA),unit = "GB")
format(object.size(object@assays$integrated),unit = "GB")

object[['RNA']] <- NULL
object[['integrated']] <- NULL
   
format(object.size(object),unit = "GB")

save(object, file = "data/MCL_TCR_BCR_SCT_20210525.Rda")

#= compare data integration==============
rm(object);GC()
load(file = "data/MCL_TCR_BCR_20210525.Rda")

# Run the standard workflow for visualization and clustering
DefaultAssay(object) <- "SCT"

object %<>% RunPCA(npcs = 150, verbose = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:150)

g2 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 0.1,label = F,
                 label.repel = T,alpha = 0.9,cols = c(Singler.colors,Singler.colors),
                 no.legend = T,label.size = 4, repel = T, title = "without Integration",
                 do.print = F, do.return = T)

jpeg(paste0(path,"UMAP_data_integration.jpeg"), units="in", width=10, height=7,res=600)
patchwork::wrap_plots(list(g2,g1), nrow = 1, ncol = 2)
dev.off()

saveRDS(object@reductions, file = "data/MCL_TCR_BCR_SCT_20210525_nointer_red.rds")


load(file = "data/MCL_TCR_BCR_20210525.Rda")
colnames(reductions[["umap"]]@cell.embeddings) %<>% paste0("old-",.)
object[["old-umap"]] <- CreateDimReducObject(embeddings = reductions[["umap"]]@cell.embeddings,
                                             key = "old-UMAP_", assay = DefaultAssay(object))
object$patient %<>% as.factor()