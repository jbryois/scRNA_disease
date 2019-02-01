# Integrating genetics with single cell RNA-seq

The paper is currently on bioRxiv and can be accessed [here](https://www.biorxiv.org/content/10.1101/528463v1).

# Requirements

1. R libraries

```
install.packages("tidyverse")
install.packages("AnnotationDbi")
install.packages("org.Hs.eg.db")
install.packages("GenABEL")
```

2. [MAGMA](https://ctg.cncr.nl/software/magma) and its auxiliary files

3. [Partitioned LDSC](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) and its auxiliary files

# Steps

The code to reproduce our results is located in this repository. The following links show the essential steps:

1) [Get GWAS summary statistics in the right format for MAGMA and LDSC](Code_Paper/Code_GWAS/get_GWAS_input.md)

2) [Get MAGMA and LDSC input for the Zeisel et al. data set](Code_Paper/Code_Zeisel/get_Zeisel_input.md)

3) [Get MAGMA and LDSC input for the GTEx et al. data set](Code_Paper/Code_GTEx/get_GTEx_input.md)

4) [Get MAGMA and LDSC input for the Skene et al. data set](Code_Paper/Code_Skene/get_Skene_input.md)

5) [Get MAGMA and LDSC input for the Habib et al. data set](Code_Paper/Code_Habib/get_Habib_input.md)

6) [Get MAGMA and LDSC input for the Saunders et al. data set](Code_Paper/Code_Saunders/get_Saunders_input.md)

# Run MAGMA

Once the GWAS sumstats are ready and the specificity files are ready, you can use the following [code](Code_Paper/Code_Zeisel/run_MAGMA.md) to test for associations using MAGMA (Zeisel data set in this example).

If you want to run your GWAS with our specificity files, you just need to get the 'Lvl5_spe_norm_no_filter.txt' files in the different MAGMA folders and run the code above.

# Run LDSC

Once the GWAS sumstats are ready and the specificity files are ready, you can use the following [code](Code_Paper/LDSC_pipeline/README.md) to test for associations using LDSC.

The code to look for heritability enrichment of the top10% most specific genes was made to be run in parallele on a SLURM cluster. It would need to be adapted if you want to run it locally.

# Expression Weighted Celltype Enrichment

The paper describing EWCE is accessible [here](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full).

EWCE code is accessible [here](https://github.com/NathanSkene/EWCE).


