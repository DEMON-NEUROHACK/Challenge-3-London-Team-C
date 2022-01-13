README
================
<h4>
README updated: <i>Jan-13-2022</i>
</h4>

# NEUROHACK 2022

## Team

![](figures/logo/team_logo.pdf)

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

#### pheno\_data.tsv

Participant x phenotype data for each of the 20 ALS participants.

Column descriptions:  
*iid*: Participant ID. *PheFemale\_Sex*: Sex.

**Creator**: Brian Schilder

### Retroviral insertions

#### HERV\_K\_Insertions.txt

A participant x region dummy file with HERV\_K insertions for all 20 MND
ALS dummy subjects.

[Additional
info.](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/data/HERVK_Insertions/readme.md)

**Creator**: Renata Kabiljo

### SNPs: variant-level

#### SNP\_VCFs.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SNP\_VCFs/* merged into one participant
x variant matrix.

**Creator**: Brian Schilder

### SNPs: gene-level

#### SNP\_VCFs.gene\_scores\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SNP\_VCFs/* merged into one participant
x gene matrix.

**Creator**: Brian Schilder

### Indels: variant-level

#### Indel\_VCFs.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *Indel\_VCFs/* merged into one
participant x variant matrix.

**Creator**: Brian Schilder

### SVs: variant-level

#### SV\_VCFs.DEL.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
variant matrix. Deletion SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.DUP-TANDEM.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
variant matrix. Duplications/tandem repeat SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.INS.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
variant matrix. Insertion SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.INV.genotype\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
variant matrix. Inversion SVs only.

**Creator**: Brian Schilder

### SVs: gene-level

#### SV\_VCFs.DEL.gene\_scores\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
gene matrix. Deletion SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.DUP-TANDEM.gene\_scores\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
gene matrix. Duplication/tandem repeat SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.INS.gene\_scores\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
gene matrix. Insertion SVs only.

**Creator**: Brian Schilder

#### SV\_VCFs.INV.gene\_scores\_matrix.tsv.gz

All 20 MND ALS dummy VCFs from *SV\_VCFs/* merged into one participant x
gene matrix. Inversion SVs only.

**Creator**: Brian Schilder

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
    ##  [1] compiler_4.1.0  magrittr_2.0.1  fastmap_1.1.0   tools_4.1.0    
    ##  [5] htmltools_0.5.2 yaml_2.2.1      stringi_1.7.6   rmarkdown_2.11 
    ##  [9] knitr_1.37      stringr_1.4.0   xfun_0.29       digest_0.6.29  
    ## [13] rlang_0.4.12    evaluate_0.14

</details>

<br>
