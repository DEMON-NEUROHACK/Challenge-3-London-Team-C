README
================
<h4>
README updated: <i>Jan-12-2022</i>
</h4>

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

Predicted HERV\_K insertions for all 20 MND ALS subjects. Predictions were obtained by running a tool Retroseq on 20 whole genome sequences of our ALS subjects

<hr>

## Citation

> van Rheenen, W., van der Spek, R.A.A., Bakker, M.K. et al. Common and
> rare variant association analyses in amyotrophic lateral sclerosis
> identify 15 risk loci with distinct genetic architectures and
> neuron-specific biology. Nat Genet 53, 1636–1648 (2021).
> <https://doi.org/10.1038/s41588-021-00973-1>

<hr>

# Session info

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
