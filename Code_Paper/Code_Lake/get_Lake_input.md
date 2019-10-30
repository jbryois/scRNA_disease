Load Data
=========

### Load single cell dataset

    library(tidyverse)

Compute sum of UMIs
===================

Cerebellum ce ll types

    cerebellum <- read_tsv("../../Data/Lake2018/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")
    annot <- tibble(cell=colnames(cerebellum)[-1],cell_type=gsub("_cbm[0-9]_[ATGCN]+$","",colnames(cerebellum)[-1]))

    cerebellum <- gather(cerebellum,cell,UMI,Gran_cbm1_TTAATCAGTCGC:Mic_cbm9_ATGGTCGGAGAC,-Gene)
    cerebellum <- inner_join(cerebellum,annot)

    cerebellum <- cerebellum %>% group_by(Gene,cell_type) %>% summarise(sumUMI=sum(UMI))
    write_tsv(cerebellum,path = "../Data/Lake2018/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

Frontal cortex cell types

    fc <- read_tsv("../../Data/Lake2018/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

    fc <- gather(fc,cell,UMI,Ex1_fcx8_GGACGCCTTTAA:Mic_fcx13_GGGAGTAACTTG,-Gene)
    cell_type1M <- gsub("_fcx[0-9]{1,2}_[ATGCN]+$","",fc$cell[1:100000000])
    cell_type2M <- gsub("_fcx[0-9]{1,2}_[ATGCN]+$","",fc$cell[100000001:200000000])
    cell_type3M <- gsub("_fcx[0-9]{1,2}_[ATGCN]+$","",fc$cell[200000001:353993295])
    cell_types <- c(cell_type1M,cell_type2M,cell_type3M)
    fc$cell_type <- cell_types

    fc <- fc %>% group_by(Gene,cell_type) %>% summarise(sumUMI=sum(UMI))
    write_tsv(fc,path = "../../Data/Lake2018/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

Visaul cortex cell types

    vc <- read_tsv("../../Data/Lake2018/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

    vc <- gather(vc,cell,UMI,Ex1_occ17_AAGTGAGTGACC:Mic_occ24_ATTAATTGAAGA,-Gene)
    cell_type1M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[1:100000000])
    cell_type2M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[100000001:200000000])
    cell_type3M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[200000001:300000000])
    cell_type4M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[300000001:400000000])
    cell_type5M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[400000001:500000000])
    cell_type6M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[500000001:600000000])
    cell_type7M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[600000001:700000000])
    cell_type8M <- gsub("_occ[0-9]{1,2}_[ATGCN]+$","",vc$cell[700000001:804566088])


    cell_types <- c(cell_type1M,cell_type2M,cell_type3M,cell_type4M,cell_type5M,cell_type6M,cell_type7M,cell_type8M)
    vc$cell_type <- cell_types

    vc <- vc %>% group_by(Gene,cell_type) %>% summarise(sumUMI=sum(UMI))
    write_tsv(vc,path = "../../Data/Lake2018/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

Load computed UMI for each tissue
=================================

    cerebellum <- read_tsv("../../Data/Lake2018/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

    ## Parsed with column specification:
    ## cols(
    ##   Gene = col_character(),
    ##   cell_type = col_character(),
    ##   sumUMI = col_double()
    ## )

    cerebellum <- mutate(cerebellum,tissue="cerebellum")

    fc <- read_tsv("../../Data/Lake2018/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

    ## Parsed with column specification:
    ## cols(
    ##   Gene = col_character(),
    ##   cell_type = col_character(),
    ##   sumUMI = col_double()
    ## )

    fc <- mutate(fc,tissue="Frontal Cortex")

    vc <- read_tsv("../../Data/Lake2018/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017_sumUMI.txt")

    ## Parsed with column specification:
    ## cols(
    ##   Gene = col_character(),
    ##   cell_type = col_character(),
    ##   sumUMI = col_double()
    ## )

    vc <- mutate(vc,tissue="Visual Cortex")

Join the cell types from the three tissues

    exp <- rbind(cerebellum,fc,vc)

Get sum of UMI for each gene and each cell type with the same name

    exp_lvl5 <- exp %>% group_by(Gene,cell_type) %>% summarise(sumUMI=sum(sumUMI)) %>% ungroup()

Tidy the data

    exp_lvl5 <- exp_lvl5 %>% spread(cell_type,sumUMI)
    exp_lvl5[is.na(exp_lvl5)] <- 0
    exp_lvl5 <- gather(exp_lvl5,cell_type,sumUMI,Ast:Purk2)

### Load gene coordinates

Load gene coordinates and extend upstream and downstream by 100kb.

File downloaded from MAGMA website
(<a href="https://ctg.cncr.nl/software/magma" class="uri">https://ctg.cncr.nl/software/magma</a>).

Filtered to remove extended MHC (chr6, 25Mb to 34Mb).

    gene_coordinates <- read.table("../../Data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",header=F,stringsAsFactors = F) %>% 
      mutate(start=ifelse(V3-100000<0,0,V3-100000),end=V4+100000,V1=as.character(V1)) %>% 
      select(2,start,end,1) %>% 
      as.tibble() %>% 
      rename(chr="V2", ENTREZ="V1") %>% 
      mutate(chr=paste0("chr",chr))

    ## Warning: `as.tibble()` is deprecated, use `as_tibble()` (but mind the new semantics).
    ## This warning is displayed once per session.

### Write dictonary for cell type names

    exp_lvl5 <- exp_lvl5 %>% ungroup() %>% rename(Lvl5=cell_type,Expr_sum_mean=sumUMI)
    dic_lvl5 <- select(exp_lvl5,Lvl5) %>% unique() %>% mutate(makenames=make.names(Lvl5))
    write_tsv(dic_lvl5,"../../Data/Lake2018/dictionary_cell_type_names.txt")

Remove not expressed genes

    genes_to_remove <- exp_lvl5 %>% group_by(Gene) %>% summarise(sum=sum(Expr_sum_mean)) %>% filter(sum==0)
    exp_lvl5 <- filter(exp_lvl5,!Gene%in%genes_to_remove$Gene)

QC
==

Remove not expressed genes

    genes_to_remove <- exp_lvl5 %>% 
      group_by(Gene) %>% 
      summarise(sum=sum(Expr_sum_mean)) %>% 
      filter(sum==0)
    exp_lvl5 <- filter(exp_lvl5,!Gene%in%genes_to_remove$Gene)

### Scale to 1 million molecules

Each cell type is scaled to the same total number of molecules.

    exp_lvl5 <- exp_lvl5 %>% 
      group_by(Lvl5) %>% 
      mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))

Specificity Calculation
=======================

The specifitiy is defined as the proportion of total expression
performed by the cell type of interest (x/sum(x)).

### Lvl5

    exp_lvl5 <- exp_lvl5 %>% group_by(Gene) %>% 
      mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean))

### Keep only ENTREZ genes protein coding and MAGMA tested genes

    entrez2symbol <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG) %>% rename(Gene="symbol",ENTREZ="gene_id")

    ## 

    exp_lvl5 <- inner_join(exp_lvl5,entrez2symbol,by="Gene") 
    exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates,by="ENTREZ") 

### Get number of genes

Get number of genes that represent 10% of the dataset

    n_genes <- length(unique(exp_lvl5$ENTREZ))
    n_genes_to_keep <- (n_genes * 0.1) %>% round()

Save expression profile for other processing
============================================

    save(exp_lvl5,file = "expression.ready.Rdata")

### Functions

#### Get MAGMA input top10%

    magma_top10 <- function(d,Cell_type){
      d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
      d_spe %>% do(write_group_magma(.,Cell_type))
    }

    write_group_magma  = function(df,Cell_type) {
      df <- select(df,Lvl5,ENTREZ)
      df_name <- make.names(unique(df[1]))
      colnames(df)[2] <- df_name  
      dir.create(paste0("MAGMA/"), showWarnings = FALSE)
      select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
      write_tsv("MAGMA/top10.txt",append=T)
    return(df)
    }

#### Get LDSC input top 10%

    write_group  = function(df,Cell_type) {
      df <- select(df,Lvl5,chr,start,end,ENTREZ)
      dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
      write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
    return(df)
    }

    ldsc_bedfile <- function(d,Cell_type){
      d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
      d_spe %>% do(write_group(.,Cell_type))
    }

### Write MAGMA/LDSC input files

Filter out genes with expression below 1 TPM.

    exp_lvl5 %>% filter(Expr_sum_mean>1) %>% magma_top10("Lvl5")

    ## Warning: group_by_() is deprecated. 
    ## Please use group_by() instead
    ## 
    ## The 'programming' vignette or the tidyeval book can help you
    ## to program with group_by() : https://tidyeval.tidyverse.org
    ## This warning is displayed once per session.

    exp_lvl5 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("Lvl5")
