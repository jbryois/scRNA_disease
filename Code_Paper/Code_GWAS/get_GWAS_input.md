---
title: "GWAS input"
output: 
  html_document:
    toc: true
    toc_float: true
    keep_md: true
    highlight: pygments
---

This document shows how to get the necessary input files for the gene expression association analysis starting from GWAS summary statistics.

We will use the GWAS sumstats of [intelligence](https://www.nature.com/articles/s41588-018-0152-6) as an example.

# Get GWAS summary statistics


```bash
mkdir -p ../../Data/GWAS/Intelligence
cd ../../Data/GWAS/Intelligence
if [ ! -f SavageJansen_IntMeta_sumstats.zip ]
then
    wget https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip
    unzip SavageJansen_IntMeta_sumstats.zip
fi
```

# MAGMA processing

First MAGMA needs to be downloaded and installed. See the following [link](https://ctg.cncr.nl/software/magma)

### Filter MAF and INFO

We will only keep SNP with MAF >1% and with a minimum INFO score >0.6


```bash
cd sumstats
awk 'NR==1 || $13>0.6 && $7>0.01 && $7 <0.99' \
SavageJansen_2018_intelligence_metaanalysis.txt > int_filtered.txt
```

### Get input for MAGMA

We need to get two different files for MAGMA:

1) A file that contains SNP, pvalue, Number of samples. Ex:


```bash
SNP     P       N
rs12184267      0.3598  225955
rs12184277      0.5116  226215
rs116801199     0.7644  226626
```

2) A file that contains SNP, chr, position, position. Ex:


```bash
rs12184267      1       715265  715265
rs12184277      1       715367  715367
rs116801199     1       720381  720381
```

These can be simply generated using awk scripts. However, sometimes the sumstats do not contain the position but only the rsid. In such case the following script can be used (**fast_match.py**). 

**fast_match.py** filters out any SNP that is not in the 1000 genomes data and gets the position based on SNP ID from the 1000 genomes file.

The g1000_eur.bim file is part of the MAGMA installation (see link above).


```bash
python ../../../../Code_Paper/utils/fast_match.py \
-b ../../../../../../MAGMA/g1000_ceu_phase3/g1000_eur.bim \
-bcols 1,0,3 \
-g int_filtered.txt \
-gcols 0,10,11
```

### Get annotation for MAGMA

Now that the sumstats are ready, we can use MAGMA to generate an annotation file. 

The following command maps all SNPs to gene using a window starting 35kb upstream and finishing 10kb downstream of each gene.

The file NCBI37.3.gene.loc.extendedMHCexcluded is a filtered version of the MAGMA NCBI37.3.gene.loc file (removed extended MHC).


```bash
magma --annotate window=35,10 --snp-loc int_filtered.txt.bed \
--gene-loc ../../../NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out int.annotated_35kbup_10_down
```

### Get gene-level association

Now we can get the gene level associations.


```bash
magma --bfile ../../../../../../MAGMA/g1000_ceu_phase3/g1000_eur \
--pval int_filtered.txt.pval ncol=3 \
--gene-annot int.annotated_35kbup_10_down.genes.annot \
--out int.annotated_35kbup_10_down
```

The magma file is now ready for association with gene expression specificity.

# LDSC

LDSC needs to be downloaded and installed. See the following [link](https://github.com/bulik/ldsc/wiki)

The summary stats ready for association can then be obtained using the following command:


```bash
munge_sumstats.py \
--sumstats int_filtered.txt \
--merge-alleles ../../../../../../Software/ldsc/w_hm3.snplist \
--signed-sumstat Zscore,0 \
--N-col N_analyzed \
--out int
```

The LDSC sumstats file is now also ready.
