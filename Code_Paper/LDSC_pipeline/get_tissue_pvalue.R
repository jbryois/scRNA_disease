#bsub -Ip /bin/bash
library("tidyverse")
library("stringr")

files <- list.files(".",pattern="tissue_dir.results",full.names = TRUE)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
  		 makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest() %>% 
  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>% 
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
    
write_tsv(d,path="cell_types_GWAS_pvalues.txt")