# Lymphocyte_meQTL
Mapping lymphocyte methylation quantitative trait loci in early arthritis

- DNA_Methylation_Preprocessing.Rmd
Preprocessing of DNA methylation data generated using the Illunina Infinium MethylationEPIC kit. IDAT files are read in and samples failing quality control checks removed. CpGs that do not pass the detection p-value threshold (indistinguishable from background signal) are removed, and the data normalised using the functional normalisation method described by Fortin et al. (PMID 25599564). Following normalisation, CpGs harbouring SNPs with a minor allale frequency >0.05, those mapping to the X/Y chromosomes, and those described as cross-reactive/ multi-mapping in previous studies were removed.
Surrogate variable analysis is performed to identify latent sources of variation, and these combined into a table with disease diagnosis to use as covariates in meQTL mapping. Methylation data are adjusted for covariates for plotting purposes. 
The outputs are DNA methylation M-values (as well as data adjusted for covariates), an annotation file detailing genomic mappings for all CpGs passing QC, a table of covariates, and the sample annotation file.

- meQTL_run_all.R
This script runs the main engine for mapping meQTLs in T and B cells seperately using MatrixEQTL (Shabalin PMID 22492648). This is submitted to a cluster using the shell script mQTL_cluster.sh.

- Get_Residuals.Rmd
Calculation of DNAm and expression residuals that were inputted to the causal inference test (CIT) and the eQTM analysis. 

- CIT.Rmd
Prepares input data and performs the CIT analysis for both cell types and all traits. 
