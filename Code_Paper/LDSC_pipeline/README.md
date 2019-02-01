<<<<<<< HEAD
The following code is intended to be run on a cluster using SLURM.

You will need to upload the 'Bed_4LDSC2_exp_gt_0.1' folders on the cluster. LDSC will also need to be downloaded and installed on the server.

Annotation files for partitioned LDSC are also necessary.

Ex:

```
scp -r LDSC_GTEx/Bed_4LDSC2_exp_gt_0.1 user@server.com:/home/LDSC_GTEx
````

=======
>>>>>>> ea5a3788b0585a60b4af9a4e04ec7b8d8c3d3d06
# Load modules

```
module load ldsc/1.0.0 
module load bedtools/2.25.0
module load r/3.4.1
```

# Tissue Analysis

Create annotation file for top10% for all genes.

For each .bed file in the folder (extended coordinates of the top10% most specific gene in each cell type), 
the script will find SNPs that overlap with the coordinates, create a new annotation file indicate which SNP
belong to the most specific cell type and compute LD scores for the annotation.

```
cp get_annotation_ldscores_tissue_v2.sh Bed_4LDSC2_exp_gt_0.1
cp get_partitioned_h2_tissue_v2.sh Bed_4LDSC2_exp_gt_0.1 
cp get_tissue_pvalue.R Bed_4LDSC2_exp_gt_0.1 
cd Bed_4LDSC2_exp_gt_0.1 
sbatch -t 48:00:00 -n 1 -o log_get_annot_ld_scores_tissue --wrap="sh get_annotation_ldscores_tissue_v2.sh"
```

# Heritability enrichment analysis

Once the previous step is finished (all ld scores have been calculated). This script will look for heritability enrichment
for all annotations

```
cd Bed_4LDSC2_exp_gt_0.1
sbatch -t 1:00:00 -n 1 -o log_get_partitioned_finucane_h2clozuk --wrap="sh get_partitioned_h2_tissue_v2.sh int.sumstats.gz"
```