########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","harmony",
                   "sctransform","data.table"), function(x) {
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
df_samples <- readxl::read_excel("output/20220323/20211210_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sequence %in% "GEX")
nrow(df_samples)
df_samples$date %<>% gsub(" UTC","",.) %>% as.character()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/MCL_15_TCR_BCR_20220323.rds")
meta.data = object@meta.data
meta.data$patient = meta.data$orig.ident
for(i in 1:length(df_samples$sample.id)){
    cells <- meta.data$orig.ident %in% df_samples$patient[i]
    print(df_samples$sample.id[i])
    print(table(cells))
    meta.data[cells,"progress"] = df_samples$progress[i]
    meta.data[cells,"project"] = df_samples$project[i]
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"tissue"] = df_samples$note[i]
    meta.data[cells,"date"] = df_samples$date[i]
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$patient)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","patient",
                                "progress","project","date","tissue", "Mean.Reads.per.Cell",
                                "Number.of.Reads","Sequencing.Saturation")]
# ======= load 10X 10k PBMC data ==============================
PBMC = readRDS("../seurat_resources/10k_PBMC_5pv2_nextgem_Chromium_X/10k_PBMC_5pv2_seurat.rds")
(mito.features <- grep(pattern = "^MT-", x = rownames(PBMC), value = TRUE))
PBMC[["percent.mt"]] <- PercentageFeatureSet(object = PBMC, pattern = "^MT-")
PBMC$orig.ident = "Normal"
PBMC$patient = "Normal"
PBMC$progress = "Normal"
PBMC$tissue = "PB"
PBMC$date = "2021-08-09"
PBMC$project = "10k_PBMC_5pv2_nextgem_Chromium_X"
PBMC$Mean.Reads.per.Cell = ""
PBMC$Number.of.Reads = ""
PBMC$Sequencing.Saturation = ""
PBMC@meta.data = PBMC@meta.data[,colnames(object@meta.data)]
#======1.6 Performing SCTransform and integration =========================
set.seed(100)

object %<>% merge(PBMC)
Normal = c("SCT","RNA")[1]
if(Normal == "SCT") object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", 
                                            verbose = TRUE)

object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

saveRDS(object, file = "data/MCL_15_TCR_BCR_20220323.rds")

# Determine the ‘dimensionality’ of the dataset  =========
npcs = 100

DefaultAssay(object) <- "RNA"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:9){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 9)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"
object %<>% RunPCA(npcs = npcs, verbose = FALSE)

jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
#object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
#                                                 key = "harmonytSNE_", assay = DefaultAssay(object))

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)

resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

saveRDS(object, file = "data/MCL_15_TCR_BCR_20220323.rds")

object_SCT <- object

#=======1.9 save SCT only =======================================
format(object.size(object_SCT),unit = "GB")

format(object.size(object_SCT@assays$RNA),unit = "GB")
format(object.size(object_SCT@assays$integrated),unit = "GB")
object_SCT[['RNA']] <- NULL
object_SCT[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object_SCT),unit = "GB")
saveRDS(object_SCT, file = "data/MCL_15_SCT_TCR_BCR_20220323.rds")

