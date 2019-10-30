---
title: "Single Cell Data Set - Skene 2018"
output:
  html_document:
    keep_md: yes
---



# Load Data

### Load single cell dataset


```r
library(tidyverse)
f <- "../../Data/KI_lvl1/celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0.rda"
load(f)
exp_lvl5 <- as.data.frame(celltype_data[[1]]$all_sct) %>% rownames_to_column("Gene")
```

Only keep genes with a unique name


```r
exp_lvl5 <- exp_lvl5 %>% 
  gather(key = Lvl5,value=Expr_sum_mean,-Gene) %>%
  as.tibble()
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
  rename(Gene=musName) %>% rename(ENTREZ=entrez)
```

# Transform Data

### Write dictonary for cell type names


```r
dic_lvl5 <- select(exp_lvl5,Lvl5) %>% unique() %>% mutate(makenames=make.names(Lvl5))
write_tsv(dic_lvl5,"../../Data/KI_lvl1/dictionary_cell_type_names.txt")
```

# QC

### Scale to 1 million molecules

Each cell type is scaled to the same total number of molecules. 


```r
exp_lvl5 <- exp_lvl5 %>% group_by(Lvl5) %>% mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))
```

# Specificity Calculation

The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).

### Lvl5


```r
exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean))
```

### Get 1to1 orthologs 

Only keep genes with a 1to1 orthologs


```r
exp_lvl5 <- inner_join(exp_lvl5,m2h,by="Gene")
```

### Only keep MAGMA genes 

Only keep protein coding genes present in the MAGMA gene file.


```r
exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates)
```

### Get number of genes

Get number of genes that represent 10% of the dataset


```r
n_genes <- length(unique(exp_lvl5$ENTREZ))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
```

# Save expression profile for other processing


```r
save(exp_lvl5,file = "expression.ready.Rdata")
```

### Functions

#### Get MAGMA input top10%


```r
magma_top10 <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group_magma(.,Cell_type))
}
```


```r
write_group_magma  = function(df,Cell_type) {
  df <- select(df,Lvl5,ENTREZ)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
  write_tsv("MAGMA/top10.txt",append=T)
return(df)
}
```

#### Get LDSC input top 10%


```r
write_group  = function(df,Cell_type) {
  df <- select(df,Lvl5,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}
```


```r
ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group(.,Cell_type))
}
```

### Write MAGMA/LDSC input files 

Filter out genes with expression below 1.


```r
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% magma_top10("Lvl5")
```

```
## Warning: group_by_() is deprecated. 
## Please use group_by() instead
## 
## The 'programming' vignette or the tidyeval book can help you
## to program with group_by() : https://tidyeval.tidyverse.org
## This warning is displayed once per session.
```


```r
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("Lvl5")
```
