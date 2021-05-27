#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============
load(file = "data/MCL_TCR_BCR_20210525.Rda")

sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# 2. check and prepare MCL data==============================
#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
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
                                 #                                 rownames(immue_exp)
                                 rownames(blue_encode)
))
length(common)
combine_ref = do.call("cbind", list(blue_encode[common,],
                                    #immue_exp[common,],
                                    B_MCL_sce[common,]))
table(combine_ref$label.fine)
system.time(trained <- trainSingleR(ref = combine_ref,
                                    labels=combine_ref$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/MCL_TCR_BCR_20210525_singleR_pred.rds")