NEUROHACK 2022
================
<h4>
README updated: <i>Jan-14-2022</i>
</h4>

## Team

![](figures/logo/team_logo.png)

### Members

``` r
members <- readxl::read_excel("presentations/team_members.xlsx")
data.frame(members)
```

    ##               Name                    Affiliation
    ## 1   Brian Schilder        Imperial College London
    ## 2        Yizhou Yu        University of Cambridge
    ## 3  Hanz Tantiangco       University of Sheffield 
    ## 4      Areda Elezi                
    ## 5 Siddharth Grover Indian Institute of Technology
    ## 6    Davide Spalla              Donders Institute
    ## 7        Guan Wang         University of Brighton
    ## 8   Renata Kabiljo          King’s College London
    ##                                                                                                                                                                                                           Expertise
    ## 1 Bioinformatics\r\nEvolutionary biology\r\nNeuroscience\r\nMachine learning\r\nGWAS summary statistics\r\nFine-mapping\r\neQTL colocalization\r\nMatrix decomposition\r\nSingle-cell RNA-seq\r\nR, Python, web dev
    ## 2                                                                                             Bioinfo\r\nWet lab (cell culture)\r\nML \r\nSNP data\r\nMolecular modelling / virtual drug screening / structural bio
    ## 3                                                                                        Computational chemistry\r\nDrug discovery\r\nDeep learning (Pytorch)\r\nNeuroscience\r\n(some) bioinformatics\r\nPython, R
    ## 4                                                                                                        Bioinformatics\r\nNexflow pipelines (RNAseq etc)\r\nWeb and software dev\r\nWet lab\r\nBasic ML \r\nPython
    ## 5                                                                                                                                                                    Machine Learning\r\nData Mining\r\nPython, C++
    ## 6                                                                                                                                         Machine learning\r\nNeuroscience\r\nPython\r\nData analysis/visualization
    ## 7                                                                                                                        GWAS analysis\r\nBulk RNA-seq data analysis (short-reads)\r\nPopulation genetics; genomics
    ## 8                                                                                                                                                          NGS Data Analysis\r\nPython, R\r\nGeneral Bioinformatics
    ##                                                                                                                             Roles
    ## 1                 Team lead. \r\nSV/SNP/Indel data preprocessing (variant and gene levels). \r\nPresentation creation/presenting.
    ## 2         Preprocessing of functional impact predictions from deep learning models (DeepSEA, Basenj). \r\nVirtual drug screening.
    ## 3 Assist in classifier model design.\r\nExplore dimensionality reduction/feature prioritization pre-step with PCA and autencoder.
    ## 4                                                      Search for additional data modalities.\r\nPrepare Expansion Hunter script.
    ## 5                                                                                              Assist in classifier model design.
    ## 6                                                                                              Design and train classifier model.
    ## 7                                                                              Preprocess combined SNP-to-gene model (cS2G) data.
    ## 8                                                                       Generate and preprocess retroviral insertion predictions.

### Project title

Predicting ALS drug targets using integrative multi-modal machine
learning

### Goals

1.  Identify molecular targets (genes and variants) of ALS (status and
    survival time).  
2.  Identify potential therapeutics for ALS.

### Workflow

1.  Preprocess input datasets
    -   MND\_ALS VCFs  
    -   External datasets  
2.  Filter input features  
3.  Train predictive model  
4.  Extract key genes
5.  Identify therapeutics

## Abstract

\[ALS background info\] \[our approach\] \[main findings\] [Future
directions](#future-directions)

## Introduction

### ALS

\[ALS background info: @Siddharth Grover @Areda Elezi @Guan Wang \]

### Classifier model

The classifier model is an instance of a deep neural network for binary
classification of survival status (SHORT vs LONG) from any number of
features coming from different data modalities.

The main idea is to have a model able to handle data from different data
streams, and integrate them to make prediction about phenotypic
outcomes. The architecture is design for extreme flexibility, making it
straightforward to feed any number of different data streams to the
model. Model complexity is kept to a minimum to ensure interpretability
and good functioning in the low sampling regime, but can easily be
extended to enhance the representational power of the model.

## Materials

The following datatypes were prepared as numeric participant x feature
matrices that could be fed into our machine learning model.

### Phenotypes

Participant x phenotype data for each of the 20 ALS participants were
encoded as numeric vectors and fed as input to the model. This currently
only include sex, but can easily be expanded to other categorical or
continuous traits as they become available.

### [Retroviral insertions](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/HERVK_Insertions)

Retroviral insertions were identified using the pipeline described
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/data/HERVK_Insertions/readme.md).

**Table 1**: Genotype encodings.

``` r
data.table::fread("data/genotype_encodings.csv")
```

    ##      key value
    ##  1:          0
    ##  2: NULL     0
    ##  3:          0
    ##  4:  ./.     1
    ##  5:  0/0     1
    ##  6:  1/.     2
    ##  7:  1/0     2
    ##  8:  ./1     2
    ##  9:  0/1     2
    ## 10:  1/1     3
    ## 11:    1     2
    ## 12:  1/2     4
    ## 13:  2/1     4
    ## 14:  2/2     5

### [SNPs: variant-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/SNP_VCFs/by_variant)

For each individual, genotypes from both alleles were encoded as integer
values at every SNP position (**Table 1**), following standards set by
the [`snpMatrix` data class in
R](https://rdrr.io/bioc/VariantAnnotation/man/genotypeToSnpMatrix-methods.html).
These encodings were then merged across individuals and cast into a
participant x variant matrix.

All code for preparing variant- and gene-level matrices (for SNPs, cS2G,
SVs , and Indel datatypes) can be found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/query_VCF.Rmd).

### [SNPs: gene-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/SNP_VCFs/by_gene)

Data generated in the prior variant-level step were then aggregated to
generate a gene-level dataset After numerically encoding all genotypes,
SNPs were mapped to genes using overlap/proximity-based annotations
already included in the VCF files. Genotype encoding were then averaged
to produce gene level encodings. Finally, we generated gene-level scores
for each participant by estimating the encoded genotype residuals,
applying z-score normalisation, and then rescaling all values from 0-1:,
Residuals were computed using the following formula:

-   *GT\_int*: Gene-level average of numerically encoded genotypes.
-   *count*: The number of time each gene is counted in the VCF
    annotations.
-   *SVLEN*: The length of the SV.
-   *TXLENGTH*: Gene length, computed by taking the mean length of all
    of each gene’s trancripts.
-   *QUAL*: Genotype quality score, produced during imputation.

`GT_int ~ count + SVLEN + TXLENGTH + QUAL`

### [cS2G: gene-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/cS2G/by_gene)

Unlike the aforementioned gene-level scores derived from SNPs, these
SNP-to-gene mappings were made using scores generated by the combined
SNPs-to-genes (cS2G) model [(Gazal et al., bioRxiv,
2021)](https://doi.org/10.1101/2021.08.02.21261488). cS2G integrates
data from various sources (e.g. QTLs, chromatin interactions) to create
more accurate SNP-to-gene mappings. We specifically used the predictions
generated for all SNP positions in the UK Biobank (found
[here](https://alkesgroup.broadinstitute.org/cS2G/cS2G_UKBB/)).

We incorporated cS2G predictions to generate gene-level scores as
follows:

*cS2G*: The cS2G-predicted probability of a SNP-to-gene interactions,
averaged by gene.

`GT_int ~ cS2G`

As before, the residuals were normalised and rescaled to generate
gene-level scores for each participant.

### [Indels: variant-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/Indel_VCFs/by_variant)

Similar to the SNP variants, insertions and deletions (indels) were
numerically encoded and cast into a participant x variant matrix.

### [SVs: variant-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/SV_VCFs/by_variant)

Structural variant (SV) genotypes were numerically encoded and then
split into four SV types: deletions (DEL), insertions (INS),
duplications/tandem repeats (DUP:TANDEM), and inversions (INV). Each of
these SV types were cast into their own separate participant x variant
matrices.

### [SVs: gene-level](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/SV_VCFs/by_gene)

Using the same strategy as the SNP gene-level matrices, gene-level
scores were computed for each SV type, producing a series of participant
x gene matrices.

## Methods

### Dimensionality reduction model

\[PCA/autoencoder description here by @Hanz\]

Prior to training the classifier model, dimensionality reduction and feature selection were performed on the training datasets. The purpose of this initial step is to optimise model training by only selecting the top features from the dataset. In this project, we compared PCA and autoencoders for dimensionality reduction and feature selection.

The code used to run PCA can be found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/PCA.ipynb).  
The code used to create and train the autoencoder can be found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/autoencoder_hanz.ipynb).

**Figure 1.**: Dimensionality reduction model architecture.

![](figures/models/autoencoder.png)

### Classifier model

Each input modality (datatype) was fed into an supervised
fully-connected artificial neural network (ANN), such that each
participant is a sample and each gene/variant/annotation is a feature.
The inputs of the model are initially partitioned into separate
“channels” of the model input, and are reduced in dimensionality by the
subsequent layers. Next, the reduced representation from each modality
are concatenated into a single vector (layer 3 in **Fig. 2** below).
Finally, the data is further compressed to predict which phenotypic
category each participant belongs to.

The model is currently designed to provide categorical predictions for
each sample: short-survival vs. long-survival, or ALS vs. control
(depending on the data available). However, it can also easily be
adapted to continuous phenotypic data (e.g. survival years, GWAS-derived
polygenic risk score (PRS)).

Once fully trained, the model can be interrogated to extract the most
relevant features per modality. This allows us to generate ranked lists
of genes/variants/annotations which can be used in the candidate
therapeutics prediction step.

All code used to create, train and evaluate the classifier model can be
found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/multimodal_classifier.ipynb).

**Figure 2**: Classifier model architecture.

![](figures/models/classifer.png)

### Therapeutics prediction

Once the relative importance of each gene for predicting ALS survival
have been identified, three complementary approaches will be used to
identify candidate therapeutics for ALS: 1. virtual screening, 2.
perturbation database queries, 3. literature mining.

#### Virtual screening

\[description by @Yizhou here\]

#### Perturbation database queries

#### Literature mining

## Results

### Retroviral Insertions HERV-K



### PCA and autoencoder

![](figures/PCA-autoencoder.png)

## Conclusions

1.  
2.  
3.  

## Future directions

1.  Estimate the size of repeats within a genome using Expansion Hunter by searching through a BAM/CRAM file for reads that span, flank, and are fully contained in each repeat.
2.  
3.  

<hr>

## References

> van Rheenen, W., van der Spek, R.A.A., Bakker, M.K. et al. Common and
> rare variant association analyses in amyotrophic lateral sclerosis
> identify 15 risk loci with distinct genetic architectures and
> neuron-specific biology. Nat Genet 53, 1636–1648 (2021).
> <https://doi.org/10.1038/s41588-021-00973-1>

<hr>

## Session info

<details>

``` r
utils::sessionInfo()
```

    ## R version 4.1.0 (2021-05-18)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7        fansi_1.0.0       crayon_1.4.2      utf8_1.2.2       
    ##  [5] digest_0.6.29     cellranger_1.1.0  lifecycle_1.0.1   magrittr_2.0.1   
    ##  [9] evaluate_0.14     pillar_1.6.4      rlang_0.4.12      stringi_1.7.6    
    ## [13] readxl_1.3.1      data.table_1.14.2 vctrs_0.3.8       ellipsis_0.3.2   
    ## [17] rmarkdown_2.11    tools_4.1.0       stringr_1.4.0     xfun_0.29        
    ## [21] yaml_2.2.1        fastmap_1.1.0     compiler_4.1.0    pkgconfig_2.0.3  
    ## [25] htmltools_0.5.2   knitr_1.37        tibble_3.1.6

</details>

<br>
