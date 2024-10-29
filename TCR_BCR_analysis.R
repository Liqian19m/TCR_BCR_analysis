library(scRepertoire)
library(ggpubr)
library(Seurat)

####TCR analysis
samplelist <- list.dirs(path="/TCR/scTCR",full.names = T, recursive = F)
TCR_list <- list()
for (i in 1:length(samplelist)){
    sample_name <- basename(samplelist[i]) 
    filepath <- paste0("/TCR/scTCR/", sample_name)
    file <- list.files(path = filepath, pattern = "contig_annotations.csv$", full.names = TRUE)
    if (length(file) > 0) {
        TCR_list[[sample_name]] <- read.csv(file, header = TRUE, as.is = TRUE)
    } else {
        print(paste("No file found for sample:", sample_name))
    }
}

combined <- combineTCR(TCR_list,
                      ID = names(TCR_list),
                      sample=c("Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer",
                                 "Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus"),
                      cells ="T-AB")

combined_1 <- addVariable(combined, name = "CNSgroup", 
                           variables = c("NCR", "NCR", "NCR", "CR","CR","CR","NCR","NCR","NCR","NCR","NCR","NCR","CR","CR","CR","CR","CR","CR","NCR","NCR")
                           )

save(combined_1, file = "TCR_combined_1.rds")

combined_sc=combined_1
for (i in 1:length(combined_sc)) {
  barcode_i <- combined_sc[[i]]$barcode
  for (n in 1:length(barcode_i)) {

    barcode_i[n] <- sapply(barcode_i[n], function(x) paste(strsplit(x, "_")[[1]][2:3], collapse = "_"))
  }
  combined_sc[[i]]$barcode <- barcode_i
}

##combining Tcell_Seurat_object
load("DBEI_Tcell.rda")
DBEITcell_TCR <- combineExpression(combined_sc,DBEITcell , 
                            cloneCall="gene", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500)) 
save(DBEITcell_TCR, file = "DBEITcell_TCR.rds")

###################################################################################################################
###BCR analysis
samplelist <- list.dirs(path="/BCR/scBCR",full.names = T, recursive = T)
BCR_list <- list()
for (i in 1:length(samplelist)){
    sample_name <- basename(samplelist[i])  # 获取样本名
    filepath <- paste0("/BCR/scBCR/", sample_name)
    file <- list.files(path = filepath, pattern = "contig_annotations.csv$", full.names = TRUE)
    if (length(file) > 0) {
        BCR_list[[sample_name]] <- read.csv(file, header = TRUE, as.is = TRUE)
    } else {
        print(paste("No file found for sample:", sample_name))
    }
}

combined <- combineBCR(BCR_list,
                      ID = names(BCR_list),
                      sample=c("Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer",
                              "Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus","Cancer","Thrombus")
                            )

combined_1 <- addVariable(combined, name = "CNSgroup", 
                           variables = c("NCR", "NCR", "NCR", "CR","CR","CR","NCR","NCR","NCR","NCR","NCR","NCR","CR","CR","CR","CR","CR","CR","NCR","NCR")
                           )
combined_1[[1]][1:5,ncol(combined_1[[1]])] # This is showing the first 5 values of the new column added
save(combined_1, file = "BCR_combined_1.rds")

combined_sc=combined_1
for (i in 1:length(combined_sc)) {
  barcode_i <- combined_sc[[i]]$barcode
  for (n in 1:length(barcode_i)) {

    barcode_i[n] <- sapply(barcode_i[n], function(x) paste(strsplit(x, "_")[[1]][2:3], collapse = "_"))
  }
  combined_sc[[i]]$barcode <- barcode_i
}
head(combined_sc$Cancer_PNA1C)

##combining Bcell_Seurat_object

load("DBEI_Bcell.rda")
DBEIBcell <- SetIdent(DBEIBcell, value = "Bcell_type")
DBEIBcell_BCR <- combineExpression(combined_sc,DBEIBcell , 
                            cloneCall="gene+nt", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500)) 
save(DBEIBcell_BCR, file = "DBEIBcell_BCR.rds")
