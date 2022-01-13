# To install and run Retroseq tool for prediction on HERV_K insertions:

```
conda install -c bioconda samtools=0.1.19
conda install -c hcc retroseq

retroseq.pl -discover -bam BAM.bam -output discoverBAM.bed -eref c3_London_c/Dummy_HERV_K_Insertions/hervk_eref.tab -id 80

retroseq.pl -call -bam BAM.bam -input discoverBAM.bed -ref REFERENCEHG19 -output BAM.vcf
```

This should result in 20 vcf files, one for each BAM.

# To aggregate these into a matrix

```python aggregate.py > HERVKMatrixforModel.txt```

Dependencies:
```
list.txt
```
which lists the samples

And the three pythong scripts in Challenge-3-London-Team-C/data/HERVK_Insertions
```
aggregate.py
combineSidewise.py
extractColumn.py
filter.py
```

# If one sample is missing:
Remove the missing sample name from list.txt and manually edit the resulting HERVKMatrixforModel.txt to add one row with the missing sample name and fill it with NA

#Copy the resulting HERVKMatrixforModel.txt in DNANexus folder c3_london_c/model_input
