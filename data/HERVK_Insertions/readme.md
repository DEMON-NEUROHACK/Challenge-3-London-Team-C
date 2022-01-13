To install and run Retroseq tool for prediction on HERV_K insertions:

```
conda install -c bioconda samtools=0.1.19
conda install -c hcc retroseq

retroseq.pl -discover -bam BAM.bam -output discoverBAM.bed -eref c3_London_c/Dummy_HERV_K_Insertions/hervk_eref.tab -id 80

retroseq.pl -call -bam BAM.bam -input discoverBAM.bed -ref REFERENCEHG19 -output BAM.vcf
```
