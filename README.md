NEUROHACK 2022
================
<h4>
README updated: <i>Jan-14-2022</i>
</h4>

## Team

![](figures/logo/team_logo.png)

### Members

![](presentations/team_members.png)

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

Amyotrophic lateral sclerosis (ALS), also known as Motor Neuron Disease
(MND), is a progressive neurodegenerative disease. ASL onsets in
individuals between 55 and 70 years old, with a predominance in the male
population and a mean survival rate of 3 to 5 years \[1\]. ALS affects
the upper and lower motor neurons and it is clinically identified by
weakness in spinal and bulbar muscles with atrophy, spasticity, weight
loss and ultimately respiratory failure \[2\]. Although it is not
usually inherited from the parents, 30 genes have been linked to the
presence of the disease, with the GGGGCC repeat expansion of the
*C9orf72* gene being present in 40% of European-ancestry patients \[1\].
Current approaches for the treatment of ASL have relied on prescription
of drugs that target cellular pathways that are responsible for
neurodegeneration.

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

### [Functional annotations](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/data/snp_geno_matrix_h3k4me3_inputForML.tsv)

Deep learning models have shown great promise in predicting regulatory
effects from single nucleotide polymorphisms. In this part, we used a
deep learning model called Basenji, created by [Kelley et
al.](https://pubmed.ncbi.nlm.nih.gov/29588361/) to predict the
regulatory activity of the genes based on methylation status of histone
3 lysine 4 (H3K4M3). [Del et
al.](https://www.nature.com/articles/s41467-020-18515-4) have shown that
H3K4M3 status are significantly informative for diseases. We therefore
hypothesise that including this information as a modality in our machine
learning model will enhance its accuracy. The code developed for this
part is
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/functional_annotations_H3K4ME3_allTissues.Rmd)

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

Prior to training the classifier model, dimensionality reduction and
feature selection were performed on the training datasets. The purpose
of this initial step is to optimise model training by only selecting the
top features from the dataset. In this project, we compared PCA and
autoencoders for dimensionality reduction and feature selection.

The code used to run PCA can be found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/PCA.ipynb).  
The code used to create and train the autoencoder can be found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/autoencoder.ipynb).

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

**Figure 2**: Classifier model architecture.

![](figures/models/classifer.png)

The model is currently designed to provide categorical predictions for
each sample: short-survival vs. long-survival, or ALS vs. control
(depending on the data available). However, it can also easily be
adapted to continuous phenotypic data (e.g. survival years, GWAS-derived
polygenic risk score (PRS)).

Once fully trained, the model can be interrogated to extract the most
relevant features per modality. This allows us to generate ranked lists
of genes/variants/annotations which can be used in the candidate
therapeutics prediction step.

The model outputs consists of. 1 - a modality weigth, informing about
the magnitude and the direction fo the contribution of each modality to
the prediction perfromance:  
![](figures/models/modality_importance.png).

2 - a contribution score for each feature in each modality, that can be
used to rank feature within modality:  
![](figures/models/feature_ranking.png).

All code used to create, train and evaluate the classifier model can be
found
[here](https://github.com/DEMON-NEUROHACK/Challenge-3-London-Team-C/blob/main/code/multimodal_classifier.ipynb).

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

1.  
2.  
3.  

## Future directions

1.  
2.  
3.  

<hr>

## References

> 1.  Marisa Cappella, Pierre-François Pradat, Giorgia Querin, Maria
>     Biferi. Beyond the Traditional Clinical Trials for Amyotrophic
>     Lateral Sclerosis and The Future Impact of Gene Therapy. Journal
>     of Neuromuscular Diseases, IOS Press, 2021, 8 (1), pp.25 - 38.
>     ff10.3233/jnd-200531ff. ffhal-03346426  
> 2.  Miller RG, Mitchell JD, Moore DH. Riluzole for amyotrophic lateral
>     sclerosis (ALS)/motor neuron disease (MND). Cochrane Database of
>     Systematic Reviews 2012, Issue 3. Art. No.: CD001447. DOI:
>     10.1002/14651858.CD001447.pub3.  
> 3.  van Rheenen, W., van der Spek, R.A.A., Bakker, M.K. et al. Common
>     and rare variant association analyses in amyotrophic lateral
>     sclerosis identify 15 risk loci with distinct genetic
>     architectures and neuron-specific biology. Nat Genet 53, 1636–1648
>     (2021). <https://doi.org/10.1038/s41588-021-00973-1>

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
    ##  [1] compiler_4.1.0    magrittr_2.0.1    fastmap_1.1.0     tools_4.1.0      
    ##  [5] htmltools_0.5.2   yaml_2.2.1        stringi_1.7.6     rmarkdown_2.11   
    ##  [9] data.table_1.14.2 knitr_1.37        stringr_1.4.0     xfun_0.29        
    ## [13] digest_0.6.29     rlang_0.4.12      evaluate_0.14

</details>

<br>
