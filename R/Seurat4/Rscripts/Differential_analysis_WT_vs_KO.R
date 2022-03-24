########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files
load(file = "data/Lorenzo-LS6_20210408_SCT.Rda")
# Need 32GB
#DefaultAssay(object) = "SCT"
#Idents(object) = "Doublets"
#object <- subset(object, idents = "Singlet")
object$label.fine %<>% gsub("Macrophages activated","Macrophages",.)
Idents(object) = "label.fine"
cell.types <- c("Fibroblasts activated", "Fibroblasts","NK cells",
                "Endothelial cells","Monocytes","B cells","Granulocytes",
                "T cells","Macrophages","Cardiomyocytes")
cell.type = cell.types[args]
object <- subset(object, idents = cell.type)
Idents(object) = "conditions"
system.time(markers <- FindAllMarkers_UMI(object, 
                                       logfc.threshold = 0.01, 
                                       return.thresh = 1, 
                                            only.pos = F,latent.vars = "nFeature_SCT",
                                               test.use = "MAST"))
markers$cell.type = cell.type
if(args < 10) args = paste0("0", args)
write.csv(markers,paste0(path,args,"_FC0.01_",cell.type,".csv"))