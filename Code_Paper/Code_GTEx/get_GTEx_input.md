---
title: "Single Cell Data Set - GTEx 2017"
output: 
  html_document:
    toc: true
    toc_float: true
    keep_md: true
---

# Load Data

### Load necessary libraries


```r
library(tidyverse)
```

### Load single cell dataset


```r
file="../../Data/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"
exp <- read_tsv(file,skip=2)
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   gene_id = col_character(),
##   Description = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

### Scale to 1M TPM


```r
exp_scaled <- apply(exp[-c(1,2)],2,function(x) x*1e6/sum(x))
exp <- cbind(exp[c(1,2)],exp_scaled) %>% as.tibble()
```

Only keep genes with a unique name


```r
exp <- exp %>% add_count(gene_id) %>% 
  filter(n==1) %>%
  mutate(gene_id=gsub("\\..+","",gene_id)) %>%
  gather(key = Lvl5,value=Expr,-gene_id,-Description) %>% 
  as.tibble()
```

### Load gene coordinates

Load gene coordinates and extend upstream and downstream by 100kb.

File downloaded from MAGMA website (https://ctg.cncr.nl/software/magma).

Filtered to remove extended MHC (chr6, 25Mb to 34Mb).


```r
gene_coordinates <- 
  read_tsv("../../Data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-100000<0,0,X3-100000),end=X4+100000) %>%
  select(X2,start,end,1) %>% 
  rename(chr="X2", ENTREZ="X1") %>% 
  mutate(chr=paste0("chr",chr))
```

# get ENTREZ names


```r
entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
```

```
## 
```

Only keep genes with a unique entrez and ensembl id


```r
entrez_ensembl_unique_genes_entrez <- entrez_ensembl %>% count(gene_id) %>% filter(n==1)
entrez_ensembl_unique_genes_ens <- entrez_ensembl %>% count(ensembl_id) %>% filter(n==1)
entrez_ensembl <- filter(entrez_ensembl,gene_id%in%entrez_ensembl_unique_genes_entrez$gene_id & ensembl_id %in% entrez_ensembl_unique_genes_ens$ensembl_id)
colnames(entrez_ensembl)[1] <- "ENTREZ"
colnames(entrez_ensembl)[2] <- "Gene"
gene_coordinates <- inner_join(entrez_ensembl,gene_coordinates) %>% as.tibble()
```

```
## Joining, by = "ENTREZ"
```

# Transform Data


```r
exp_lvl5 <- exp %>% rename(Expr_sum_mean=Expr,Gene=gene_id)
```

### Get only MAGMA genes 


```r
exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates)
```

```
## Joining, by = "Gene"
```

### Write dictonary for cell type names


```r
dic_lvl5 <- exp_lvl5 %>% select(Lvl5) %>% unique() %>% mutate(makenames=make.names(Lvl5))
write_tsv(dic_lvl5,"../../Data/GTEx//dictionary_cell_type_names.txt")
```

### Get summary stats on the dataset


```r
sumstats_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>%
  summarise(total_exp=sum(Expr_sum_mean),Ngenes=sum(Expr_sum_mean>0))
```

# QC

### Remove not expressed genes


```r
not_expressed <- exp_lvl5 %>% group_by(Gene) %>% summarise(total_sum=sum(Expr_sum_mean)) %>% filter(total_sum==0) %>% select(Gene) %>% unique() 
exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)
```

# Specificity Calculation

The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).

### Lvl5


```r
exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean),max=max(Expr_sum_mean))

exp_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>%
    mutate(spe_norm=GenABEL::rntransform(specificity))
```

### Functions


#### Get MAGMA input continuous


```r
magma_input <- function(d,Cell_type,spec, outputfile_name){
 d_spe <- d %>% select(ENTREZ,Cell_type,spec) %>% spread(key=Cell_type,value=spec)
 colnames(d_spe) <- make.names(colnames(d_spe))
 dir.create("MAGMA_GTEx/", showWarnings = FALSE)
 write_tsv(d_spe,paste0("MAGMA_GTEx/",Cell_type,"_",spec,"_",outputfile_name,".txt"))
}
```

#### Get LDSC input top 10%


```r
write_group  = function(df,Cell_type,outputfile_name) {
  df <- select(df,Lvl5,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC_GTEx/Bed_4LDSC2_",outputfile_name), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC_GTEx/Bed_4LDSC2_",outputfile_name,"/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}
```


```r
ldsc_bedfile <- function(d,Cell_type, outputfile_name){
  d_spe <- d %>% group_by_(Cell_type) %>% filter(specificity>=quantile(specificity,0.9)) 
  d_spe %>% do(write_group(.,Cell_type,outputfile_name))
}
```

### Write MAGMA/LDSC input files 

Keep all genes for continuous analysis.


```r
magma_input(exp_lvl5,"Lvl5","spe_norm","no_filter")
```


```r
exp_lvl5 %>% filter(max>1) %>% ldsc_bedfile("Lvl5","exp_gt_0.1")
```
