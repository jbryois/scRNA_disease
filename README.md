# Integrating genetics with single cell RNA-seq

The paper is currently on bioRxiv and can be accessed [here](https://www.biorxiv.org/content/10.1101/528463v1).

The code to reproduce our results is located in this repository. The following links show the essential steps:

1) Get GWAS summary statistics in the right format for MAGMA and LDSC. [Code](Code_Paper/Code_GWAS/get_GWAS_input.md)

2) Get MAGMA and LDSC input for the Zeisel et al. data set. [Code](Code_Paper/Code_Zeisel/get_Zeisel_input.md)

3) Get MAGMA and LDSC input for the GTEx et al. data set. [Code](Code_Paper/Code_GTEx/get_GTEx_input.md)

## Run MAGMA

Once the GWAS sumstats are ready and the specificity files are ready, you can use the following [code](Code_Paper/Code_Zeisel/run_MAGMA.md) to test for associations using MAGMA (Zeisel data set in this example).

## Run LDSC

Once the GWAS sumstats are ready and the specificity files are ready, you can use the following [code](Code_Paper/LDSC_pipeline/README.md) to test for associations using LDSC.

The code to look for heritability enrichment of the top10% most specific genes was made to be run in parallele on a SLURM cluster. It would need to be adapted if you want to run it locally.

**Under construction**

To do: 
1) Add Saunders, Habib and Skene datasets processing files
2) Add LDSC get pvalue file
3) Clean a bit