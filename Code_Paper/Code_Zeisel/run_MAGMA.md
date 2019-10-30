---
title: "Run MAGMA"
output: 
  html_document:
    keep_md: true
---

We can now run MAGMA on the processed intelligence GWAS using the following commands:


```r
cd MAGMA/
int="../../Data/GWAS/Intelligence/sumstats/int.annotated_35kbup_10_down.genes.raw"
f="../../Data/random_gene_sets/synaptic_ARC_Seth_Grant.t.txt"
cell_type="../MAGMA/top10.txt"
magma --gene-results  $int --set-annot  $cell_type --out int
cd ..
```
