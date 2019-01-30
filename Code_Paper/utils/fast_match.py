#!/usr/bin/env python
# coding: utf-8

import argparse
import time

#Times the execution of the script

start = time.time()

#Parse arguments

parser = argparse.ArgumentParser(description='Filter GWAS sumstats based on SNP name to only retain SNPs present in bim file. Produces two outputs, one with pvalues and one with positions (input for MAGMA)')
parser.add_argument('-b', '--bim',
                        help='Bim file',
                        required='True')
parser.add_argument('-bcols', '--bimcols',
                        help='Columns to use from bim file (SNP,CHR,POS)',
                        required='True')
parser.add_argument('-g', '--gwas',
                        help='GWAS file',
                        required='True')
parser.add_argument('-gcols', '--gwascols',
                        help='Columns to use from gwas file (SNP,PVALUE,N)',
                        required='True')
args = parser.parse_args()

##################
###     Main   ###
##################

my_dict = {}

#Get snp,chr,pos from bim columns argument as integer
snp_col,chr_col,pos_col = (int(x) for x in args.bimcols.split(','))

#Get snp,pvalue,n from gwas columns argument as integer. If N index is not present, n is not defined
gwas_cols = args.gwascols.split(',')
length_gwas_cols=len(gwas_cols)

if(length_gwas_cols>2):
	snp_colg,pval_colg,n_colg = (int(x) for x in gwas_cols)
else:
	snp_colg,pval_colg = (int(x) for x in gwas_cols)

print('')
print('Reading Bim file: {}'.format(args.bim))
print('Column numbers are 0 based')
print('Creating dictionary using SNP (col {}) key'.format(snp_col))
print('Values are SNP (col {}), CHR (col {}), POS (col {}), POS (col {})'.format(snp_col,chr_col,pos_col,pos_col))

#Store each line in a dictionary with the SNP name as key
with open(args.bim,'r') as file:
    for line in file:
        items = line.split()
        snp,value = items[snp_col], [items[snp_col],items[chr_col],items[pos_col],items[pos_col]]
        my_dict[snp] = value

print('')
print('Reading GWAS file: {}'.format(args.gwas))
print('Printing output files for SNPs (col {}) in dictionary'.format(snp_colg))
print('Output file 1 contains SNP (col {}), CHR (col {}), POS (col {}), POS (col {})'.format(snp_col,chr_col,pos_col,pos_col))
if(length_gwas_cols>2):
	print('Output file 2 contains SNP (col {}), Pvalue (col {}), N (col {})'.format(snp_colg,pval_colg,n_colg))
else:
	print('Output file 2 contains SNP (col {}), Pvalue (col {})'.format(snp_colg,pval_colg))

#Open GWAS file, split lines, extract snp, pvalue and n columns and write bed files + pvalues files if snp from 1000g
#was genotypes (is in the dictionary)

with open(args.gwas,'r') as file2, open(args.gwas+'.pval', 'w') as out1, open(args.gwas+'.bed', 'w') as out2:
   for index,line in enumerate(file2):
        items = line.split()
        if(length_gwas_cols>2):
        	snp,pvalue,n = items[snp_colg],items[pval_colg],items[n_colg]
        else:
        	snp,pvalue = items[snp_colg],items[pval_colg]
        if (index==0 and length_gwas_cols>2):
        	print("SNP","P","N",sep="\t",file=out1)
        elif(index==0 and length_gwas_cols==2):
        	print("SNP","P",sep="\t",file=out1)
        elif (snp in my_dict.keys() and length_gwas_cols>2):
        	print(snp,pvalue,n,sep="\t",file=out1)
        	print("\t".join(my_dict[snp]),file=out2)
        elif (snp in my_dict.keys() and length_gwas_cols==2):
        	print(snp,pvalue,sep="\t",file=out1)
        	print("\t".join(my_dict[snp]),file=out2)

##########################
### Print Running time ###
##########################

end = time.time()
run_time=end - start
print('Running time = {:.2f} seconds'.format(run_time))

