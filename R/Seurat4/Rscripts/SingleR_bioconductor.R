#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(Seurat)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

# ====== load single cell =============
object = readRDS("data/MCL_15_SCT_TCR_BCR_20220323.rds")

sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                            colData=DataFrame(object@meta.data))
rm(object);GC()

datasets = c("blue_encode_MCL","azimuth_PBMC")
print(dataset <- datasets[args])

# 2. check and prepare MCL data==============================
#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
#' 
RemoveDup <- function(mat){
    gene_id <- as.matrix(mat[,1])
    mat <- mat[,-1]
    if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
    rownames(mat) <- 1:nrow(mat)
    mat[is.na(mat)] = 0
    mat <- cbind(mat, "rowSums" = rowSums(mat))
    mat <- mat[order(mat[,"rowSums"],decreasing = T),]
    gene_id <- gene_id[as.numeric(rownames(mat))]
    remove_index <- duplicated(gene_id)
    mat <- mat[!remove_index,]
    rownames(mat) <- gene_id[!remove_index]
    return(mat[,-ncol(mat)])
}

if(dataset == "blue_encode_MCL"){
    MCL_bulk <- read.csv(file="../scRNAseq-MCL/data/RNAseq/MCL_bulk_191006.csv")
    MCL_bulk <- RemoveDup(MCL_bulk)
    B_bulk <- read.csv(file="../scRNAseq-MCL/data/RNAseq/B_bulk_191006.csv")
    B_bulk <- RemoveDup(B_bulk)
    
    B_MCL_bulk <- merge(log1p(B_bulk),log1p(MCL_bulk),
                        by="row.names",all=FALSE)
    B_MCL_bulk <- RemoveDup(B_MCL_bulk)
    
    meta.data = data.frame("label.main"= c(rep("B cells",ncol(B_bulk)),
                                           rep("MCL",ncol(MCL_bulk))),
                           "label.fine" = c(rep("B cells, PB",ncol(B_bulk)),
                                            rep("MCL",ncol(MCL_bulk))),
                           "label.ont" = c(colnames(B_bulk),colnames(MCL_bulk)))
    B_MCL_sce <- SummarizedExperiment(list(logcounts=B_MCL_bulk),
                                      colData=DataFrame(meta.data))
    
    # ====== load reference =============
    blue_encode <- BlueprintEncodeData()
    #remove = grepl("CD4|CD8|Tregs|B-cells|Monocytes",blue_encode$label.fine)
    #blue_encode = blue_encode[,!remove]
    
    #immue_exp <- DatabaseImmuneCellExpressionData()
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(B_MCL_sce),
                                     rownames(blue_encode)
    ))
    length(common)
    combine_ref = do.call("cbind", list(blue_encode[common,],
                                        B_MCL_sce[common,]))
    table(combine_ref$label.fine)
    system.time(trained <- trainSingleR(ref = combine_ref,
                                        labels=combine_ref$label.fine))
    # elapsed 4872.846 sec
}


if(dataset == "azimuth_PBMC"){
    # ======= load azimuth PBMC data ==============================
    path = "../seurat_resources/azimuth/PBMC/"
    counts <- Read10X(paste0(path, "GSE164378/GSM5008740_RNA_5P"))
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    meta.data = read.csv(paste0(path,"GSE164378/GSE164378_sc.meta.data_5P.csv"),row.names =1)
    table(rownames(meta.data) == colnames(counts))
    PBMC <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                 colData=DataFrame(meta.data))
    rm(counts);GC()
    
    # ====== conbime data =============
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(PBMC)
    ))
    length(common)
    table(PBMC$celltype.l3)
    system.time(trained <- trainSingleR(ref = PBMC[common,],
                                        labels=PBMC$celltype.l3))
}

system.time(pred <- classifySingleR(sce[common,], trained))
saveRDS(object = pred, file = paste0("output/MCL_15_SCT_TCR_BCR_20220323_",dataset,"_singleR_pred.rds"))