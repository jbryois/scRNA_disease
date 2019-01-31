---
title: "Single Cell Data Set - Saunders 2018"
output: 
  html_document:
    toc: true
    toc_float: true
    keep_md: true
---

The data was downloade from [here](http://dropviz.org)

# Load Data

### Load single cell dataset


```r
library(tidyverse)
file="../../Data/Saunders/metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv"
exp <- read_csv(file)
exp <- gather(exp,key = tissue_subcluster,value=Expr,-X1) %>% as.tibble() %>% rename(Gene="X1")
```

### Load Annotation


```r
annot <- read_csv("../../Data/Saunders/annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv") %>% mutate(Lvl5=make.names(tissue_subcluster)) %>% select(-X1)
```

### Load gene coordinates

Load gene coordinates and extend upstream and downstream coordinates by 100kb.

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

### Load mouse to human 1to1 orthologs

File downloaded from: http://www.informatics.jax.org/homology.shtml and parsed.


```r
m2h <- read_tsv("../../Data/m2h.txt",col_types = "iccccc") %>% 
  select(musName,entrez) %>%
  rename(Gene=musName) %>% 
  rename(ENTREZ=entrez)
```

# Transform Data

### Add Cell type names


```r
exp_lvl5 <- inner_join(exp,annot,by="tissue_subcluster") %>% 
  as.tibble() %>% ungroup() %>% rename(Expr_sum_mean=Expr)
```

### Write dictonary for cell type names


```r
write_tsv(annot,"../../Data/Saunders/dictionary_cell_type_names.txt")
```

### Get summary stats on the dataset


```r
sumstats_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>%
  summarise(total_umi=sum(Expr_sum_mean),
            Ngenes=sum(Expr_sum_mean>0)) %>% arrange(total_umi)
```

# QC

### Remove low quality cell types

Cell types should have at least 200k reads.


```r
filter_lvl5 <- sumstats_lvl5 %>% filter(total_umi>200000)
exp_lvl5 <- filter(exp_lvl5,Lvl5%in%filter_lvl5$Lvl5)
```

### Remove not expressed genes


```r
not_expressed <- exp_lvl5 %>% group_by(Gene) %>% summarise(total_sum=sum(Expr_sum_mean)) %>% filter(total_sum==0) %>% select(Gene) %>% unique() 
exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)
```

### Scale to 1M molecules

Each cell type is scaled to the same total number of molecules. 


```r
exp_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>% mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))
```

# Specificity Calculation

The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).

### Lvl5


```r
exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean),
         max=max(Expr_sum_mean))
```

### Get 1to1 orthologs and Normalise

Get 1to1 orthologs, get binned specificity and standard normalized specificity within cell type.


```r
  exp_lvl5 <- inner_join(exp_lvl5,m2h,by="Gene") %>% group_by(Lvl5) %>%
    mutate(spe_norm=GenABEL::rntransform(specificity))
```

### Functions

#### Get MAGMA input continuous


```r
magma_input <- function(d,Cell_type,spec, outputfile_name){
 d_spe <- d %>% select(ENTREZ,Cell_type,spec) %>% spread(key=Cell_type,value=spec)
 colnames(d_spe) <- make.names(colnames(d_spe))
 dir.create("MAGMA_Saunders", showWarnings = FALSE,recursive = TRUE)
 write_tsv(d_spe,paste0("MAGMA_Saunders/",Cell_type,"_",spec,"_",outputfile_name,".txt"))
}
```

#### Get LDSC input top 10%


```r
write_group  = function(df,Cell_type,outputfile_name) {
  df <- select(df,Lvl5,chr,start,end,ENTREZ)
  dir.create("LDSC_Saunders", showWarnings = FALSE,recursive = TRUE)
  dir.create(paste0("LDSC_Saunders/Bed_4LDSC2_",outputfile_name), showWarnings = FALSE)
  write_tsv(df[-1],paste0("LDSC_Saunders/Bed_4LDSC2_",outputfile_name,"/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}
```


```r
ldsc_bedfile <- function(d,Cell_type, outputfile_name){
d_spe <- d %>% inner_join(gene_coordinates,by="ENTREZ") %>% group_by_(Cell_type) %>% filter(specificity>=quantile(specificity,0.9)) 
d_spe %>% do(write_group(.,Cell_type,outputfile_name))
}
```

### Write MAGMA/LDSC input files 

Keep all genes for continuous analysis.


```r
magma_input(exp_lvl5,"Lvl5","spe_norm","no_filter")
```

Filter genes with expression below 1 in the top cell type (scaled to 1M molecules) and take top 10% most specific genes.


```r
exp_lvl5 %>% filter(max>1) %>% ldsc_bedfile("Lvl5","exp_gt_0.1")
```
