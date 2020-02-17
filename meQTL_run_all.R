library(MatrixEQTL)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

for(cell in c("T","B")) {
  
  # Set Path constants
  preset.path             <- paste0("~/MatrixEQTL_Input_Data/",cell,"_Cell/")
  preset.output.path      <- paste0("~/Results_mQTL/",cell,"_Cell/")
  
  # Setup Paths for input data + run settings
  useModel                <- modelLINEAR
  SNP_file_name           <- paste0(preset.path, "Genotypes_",cell,"_Cell.csv")
  expression_file_name    <- paste0(preset.path, "Mvals_",cell,"_Cell.csv")
  covariates_file_name    <- paste0(preset.path, "Covariates",cell,"_Cell.csv")
  snps_location_file_name <- paste0(preset.path, "Genotype_Annotation_",cell,"_Cell.csv")
  gene_location_file_name <- paste0(preset.path, "Meth_Annotation_",cell,"_Cell.csv")
  output_file_name        <- tempfile()
  errorCovariance         <- numeric()
  
  # Read Genotype Data
  snps                    <- SlicedData$new()
  snps$fileDelimiter      <- ","
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows       <- 1        
  snps$fileSkipColumns    <- 1 
  snps$fileSliceSize      <- 10000 
  snps$LoadFile(SNP_file_name)
  
  
  # Read Methylation Data
  gene                    <- SlicedData$new()
  gene$fileDelimiter      <- ","
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows       <- 1 
  gene$fileSkipColumns    <- 1
  gene$fileSliceSize      <- 10000
  gene$LoadFile(expression_file_name)
  
  # Read Covariates
  cvrt                    <- SlicedData$new();
  cvrt$fileDelimiter      <- ","
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows       <- 1
  cvrt$fileSkipColumns    <- 1
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name)
  }
  
  # Read Annotation
  snpspos                 <- read_csv(snps_location_file_name) %>% as.data.frame %>% dplyr::select(-X1)
  genepos                 <- read_csv(gene_location_file_name) %>% as.data.frame %>% dplyr::select(-X1)
  
  # Setup output temp files
  output_file_name_cis    <- tempfile()
  output_file_name_tra    <- tempfile()
  
  # Cutoff thresholds
  pvOutputThreshold_cis   <- 1e-2
  pvOutputThreshold_tra   <- 1e-4
  cisDist                 <- 1e6
  
  # Run Engine
  me      <- Matrix_eQTL_main(snps                  = snps,
                              gene                  = gene,
                              cvrt                  = cvrt,
                              output_file_name      = output_file_name_tra,
                              pvOutputThreshold     = pvOutputThreshold_tra,
                              useModel              = useModel, 
                              errorCovariance       = errorCovariance, 
                              verbose               = T,
                              output_file_name.cis  = output_file_name_cis,
                              pvOutputThreshold.cis = pvOutputThreshold_cis,
                              snpspos               = snpspos, 
                              genepos               = genepos,
                              cisDist               = cisDist,
                              pvalue.hist           = "qqplot",
                              min.pv.by.genesnp     = F,
                              noFDRsaveMemory       = F)
  
  # Disconnect from temporary data structures
  unlink(output_file_name_tra)
  unlink(output_file_name_cis)
  
  # Results
  results_cis   <- me$cis$eqtls   %>% as.data.frame %>% arrange(FDR)
  results_trans <- me$trans$eqtls %>% as.data.frame %>% arrange(FDR)
  
  save(me, file = paste0(preset.output.path, "MatrixEQTL_Results.RData"), compress = T)
  write_csv(results_cis,   path = paste0(preset.output.path,"Cis_Results.csv"))
  write_csv(results_trans, path = paste0(preset.output.path,"Trans_Results.csv"))
  
}