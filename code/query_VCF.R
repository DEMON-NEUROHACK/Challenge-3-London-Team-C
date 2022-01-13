#### NEUROHACK 2022 Team: Bioinformagicians ####
## Preprocessing script for VCFs

#### Compute requirements ####
## Minimum compute requirements are based on a Macbook Pro 
## on which the 20 dummy VCFs were preprocessed (though more is better!): 
## 
## Cores: >=12
## Memory: >=32GB
## Storage: >=3GB


#### Install packages ####
## CRAN 
if(!require(remotes)) install.packages("remotes")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyr)) install.packages("tidyr")
if(!require(tibble)) install.packages("tibble")
if(!require(here)) install.packages("here")
if(!require(stringr)) install.packages("stringr")
if(!require(parallel)) install.packages("parallel")
if(!require(data.table)) install.packages("data.table")
if(!require(scales)) install.packages("scales")
if(!require(stats)) install.packages("stats")
if(!require(Matrix)) install.packages("Matrix")
## Bioconductor 
if(!require(BiocManager)) install.packages("BiocManager")
if(!require(VariantAnnotation)) BiocManager::install("VariantAnnotation")
if(!require(StructuralVariantAnnotation)) BiocManager::install("StructuralVariantAnnotation")

#### Set working directory ####
setwd(here::here())
source(here::here("code/gene_length.R"))
source(here::here("code/utils.R"))

#### Set up save directory ####
## This is the folder where all output files will be written to. 
save_dir <- here::here("data")

#### List VCFs ####
root <- "MND_ALS_Data_Central" ## Define where data is located
snp_vcfs <- list.files(file.path(root,"SNP_VCFs"),".vcf$", 
                       full.names = TRUE, recursive = TRUE)
indel_vcfs <- list.files(file.path(root,"Indel_VCFs"),".vcf$", 
                         full.names = TRUE, recursive = TRUE)
sv_vcfs <- list.files(file.path(root,"SV_VCFs"),".vcf$", 
                      full.names = TRUE, recursive = TRUE)
 
#### Source custom functions from R scripts ####
source(here::here("code/gene_length.R"))
source(here::here("code/utils.R"))


#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### SNPs ####
## Estimated time: ~20min per 20 VCFs.
mat_var <- genotype_matrix(vcfs = snp_vcfs,
                           save_prefix=file.path(save_dir,
                                                 "SNP_VCFs/by_variant/SNP_VCFs"))

grl <- gather_vcfs(vcfs = snp_vcfs)
gr_dt <- merge_vcfs(grl = grl, 
                    save_prefix = file.path(save_dir,
                                            "SNP_VCFs/by_variant/SNP_VCFs"))
gene_counts <- get_gene_counts(gr_dt = gr_dt, 
                               gene_col = "CSQ", 
                               extract_string = "ENSG")
gene_counts <- add_gene_length(gene_counts = gene_counts)
gene_counts <- add_gene_score(gene_counts = gene_counts,
                              formula_str = "GT_int ~ count + TXLENGTH + QUAL")
### Save matrices 
data_type <- "gene_scores" # gene_scores"
mat_gene <- create_gene_matrix(gene_counts = gene_counts, 
                               data_type = data_type, 
                               save_prefix = file.path(save_dir,
                                                       "SNP_VCFs/by_gene/SNP_VCFs"))



#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### SVs #### 
# Process Structural Variants (SVs) at the variant- and gene- levels. 
## Estimated time: ~30min per 20 VCFs.
#### Create variant-level matrices ####
genotypes_list <- genotype_matrix_grouped(vcfs = sv_vcfs, 
                                          grouping_var = "SVTYPE",
                                          save_prefix = file.path(save_dir,
                                                                  "SV_VCFs/by_variant/SV_VCFs"))
#### Create gene-level matrices ####
grl <- gather_vcfs(vcfs = sv_vcfs,
                   cols=c("id","SVTYPE","SVLEN","ensembl_gene_id","QUAL","GT"))
gr_dt <- merge_vcfs(grl = grl,
                    save_prefix = file.path(save_dir,
                                            "SV_VCFs/by_variant/SV_VCFs"))
gene_counts <- get_gene_counts(gr_dt = gr_dt)
gene_counts <- add_gene_length(gene_counts = gene_counts)
gene_counts <- add_gene_score(gene_counts = gene_counts)
### Save matrices 
data_type <- "gene_scores" # gene_scores"
mat_list <- create_gene_matrix(gene_counts = gene_counts,
                               grouping_var = "SVTYPE",
                               data_type = data_type, 
                               save_prefix = file.path(save_dir,
                                                       "SV_VCFs/by_gene/SV_VCFs"))


#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### ------------------------------------------ ####
#### Indels ####
## Create variant-level matrices.
## Estimated time: ~2-3 hours per 20 VCFs (after the genotype_matrix, which only takes ~10min).
mat <- genotype_matrix(vcfs = indel_vcfs, 
                       structural_variants = TRUE,
                       save_prefix = file.path(save_dir,
                                               "Indel_VCFs/by_variant/Indel_VCFs"))

#### Create gene-level matrices ####
grl <- gather_vcfs(vcfs = indel_vcfs)
gr_dt <- merge_vcfs(grl = grl, 
                    save_prefix =file.path(save_dir,
                                           "Indel_VCFs/by_variant/Indel_VCFs"))
gene_counts <- get_gene_counts(gr_dt = gr_dt, 
                               gene_col = "CSQ", 
                               extract_string = "ENSG")
gene_counts <- add_gene_length(gene_counts = gene_counts)
gene_counts <- add_gene_score(gene_counts = gene_counts)
### Save matrices 
data_type <- "gene_scores" # gene_scores"
mat_list <- create_gene_matrix(gene_counts = gene_counts,
                               grouping_var = "SVTYPE",
                               data_type = data_type, 
                               save_prefix = file.path(save_dir,
                                                       "Indel_VCFs/by_gene/Indel_VCFs"))

