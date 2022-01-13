README
================
<h4>
README updated: <i>Jan-13-2022</i>
</h4>

# NEUROHACK 2022

## Team

![](figures/logo/team_logo.png)

### Members

``` r
members <- readxl::read_excel("presentations/team_members.xlsx")
# knitr::kable(members)
data.table::data.table(members)
```

    ##                Name                    Affiliation
    ## 1:   Brian Schilder        Imperial College London
    ## 2:        Yizhou Yu        University of Cambridge
    ## 3:  Hanz Tantiangco       University of Sheffield 
    ## 4:      Areda Elezi                Crick Institute
    ## 5: Siddharth Grover Indian Institute of Technology
    ## 6:    Davide Spalla              Donders Institute
    ## 7:        Guan Wang         University of Brighton
    ## 8:   Renata Kabiljo          King’s College London
    ##                                                                                                                                                                                                            Expertise
    ## 1: Bioinformatics\r\nEvolutionary biology\r\nNeuroscience\r\nMachine learning\r\nGWAS summary statistics\r\nFine-mapping\r\neQTL colocalization\r\nMatrix decomposition\r\nSingle-cell RNA-seq\r\nR, Python, web dev
    ## 2:                                                                                             Bioinfo\r\nWet lab (cell culture)\r\nML \r\nSNP data\r\nMolecular modelling / virtual drug screening / structural bio
    ## 3:                                                                                        Computational chemistry\r\nDrug discovery\r\nDeep learning (Pytorch)\r\nNeuroscience\r\n(some) bioinformatics\r\nPython, R
    ## 4:                                                                                                        Bioinformatics\r\nNexflow pipelines (RNAseq etc)\r\nWeb and software dev\r\nWet lab\r\nBasic ML \r\nPython
    ## 5:                                                                                                                                                                    Machine Learning\r\nData Mining\r\nPython, C++
    ## 6:                                                                                                                                         Machine learning\r\nNeuroscience\r\nPython\r\nData analysis/visualization
    ## 7:                                                                                                                        GWAS analysis\r\nBulk RNA-seq data analysis (short-reads)\r\nPopulation genetics; genomics
    ## 8:                                                                                                                                                          NGS Data Analysis\r\nPython, R\r\nGeneral Bioinformatics
    ##                                                                                                                              Roles
    ## 1:                 Team lead. \r\nSV/SNP/Indel data preprocessing (variant and gene levels). \r\nPresentation creation/presenting.
    ## 2:         Preprocessing of functional impact predictions from deep learning models (DeepSEA, Basenj). \r\nVirtual drug screening.
    ## 3: Assist in classifier model design.\r\nExplore dimensionality reduction/feature prioritization pre-step with PCA and autencoder.
    ## 4:                                                      Search for additional data modalities.\r\nPrepare Expansion Hunter script.
    ## 5:                                                                                              Assist in classifier model design.
    ## 6:                                                                                              Design and train classifier model.
    ## 7:                                                                              Preprocess combined SNP-to-gene model (cS2G) data.
    ## 8:                                                                       Generate and preprocess retroviral insertion predictions.

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

\[ALS background info\] \[our approach\] \[main findings\] [future
directions](#future-directions)

## Materials

The following datatypes were prepared as numeric participant x feature
matrices that could be fed into our machine learning model.

### Phenotypes

Participant x phenotype data for each of the 20 ALS participants were
encoded as numeric vectors and fed as input to the model. This currently
only include sex, but can easily be expanded to other categorical or
continuous traits as they become available.

### Retroviral insertions

Retroviral insertions were identified using the pipeline described
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/data/HERVK_Insertions/readme.md).

### SNPs: variant-level

### SNPs: gene-level

### Indels: variant-level

### SVs: variant-level

### SVs: gene-level

## Methods

### Dimensionality reduction model

\[description here\]

**Figure 1.**:

![](figures/models/autoencoder.png)

### Classifier moel

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

**Figure 2**:

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

## Conclusions

## Future directions

<hr>

# File descriptions

## [data/SNP\_VCFs/](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/SNP_VCFs)

### by\_variant/

### SNP\_VCFs.all.tsv.gz

All 20 MND ALS dummy VCFs from *SNP\_VCFs/* merged into one participant
x variant matrix.

The matrix is filled with genotypes numerically encoded as detailed
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/data/genotype_encodings.csv).

### SNP\_VCFs.merged.tsv.gz

All 20 MND ALS dummy VCFs from *SNP\_VCFs/* merged into a `data.table`
in long-format.

This is the data from which *SNP\_VCFs.all.tsv.gz* is derived.

## [data/HERVK\_Insertions/](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/tree/main/data/HERVK_Insertions)

### HERV\_K\_Insertions.txt

A dummy file with HERV\_K insertions for all 20 MND ALS dummy subjects.

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
