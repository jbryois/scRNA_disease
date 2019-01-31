#!/usr/bin/sh

#Julien Bryois 5.10.2017

module load ldsc/1.0.0 
module load bedtools/2.25.0
module load r/3.4.1

sumstats=$1

path_name="/nas/depts/007/sullilab/shared/partitioned_LDSC/"

#Downloaded on LDSC wiki
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
all_annotations="1000G_EUR_Phase3_baseline"

for f in *_tissue_dir 
do
echo $f
gwas_name=`basename $sumstats | cut -d "." -f 1`
echo $gwas_name
cd $f
sbatch -t 0:30:00 -n 1 -o log_${gwas_name}_${f} --wrap="ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ../${gwas_name}_${f}"
cd ..
done