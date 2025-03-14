#this script calculates additional values to pyseer results

#load R packages
library(dplyr)
library(readr)
library(stringr)

#load metadata
read_tsv('data/MP_metadata_1121.tsv')%>% 
  as.data.frame() %>% 
  `rownames<-`(unlist(.[,1]))->metadata

#load results. `pyseer_ec1.txt` is the significant kmers list from pyseer. `unitig_presabs.txt` is the presence absence of unitigs detected by unitig caller
gwas_ec1 <- read_tsv('data/pyseer_ec1.txt') %>%  mutate(`kmer length`=nchar(variant))
kmer_presabs_ec1 <- read_lines('data/unitig_presabs.txt') %>% as.list
kmer_presabs_ec1 %>% 
  `names<-`(kmer_presabs_ec1 %>% 
              lapply(function(x)str_extract(x,'^[ATGC]*')) %>% 
              unlist)->kmer_presabs_ec1
kmer_presabs_ec1 %>% lapply( function(x) x %>% str_remove('^[ATGC]* \\| ') %>%  
                               str_replace_all(':1 ',',') %>% str_remove (':1') %>% 
                               str_split_1(',')
) -> kmer_presabs_ec1

portion_ec1_2023 = sum(metadata$clade=='T1-2-EC1'&metadata$year==2023) / sum(metadata$clade=='T1-2-EC1')

gwas_ec1 %>%rowwise %>%  mutate(
  af_2023=mean( (metadata  %>% filter(year==2023,clade=='T1-2-EC1') )$strain %in% kmer_presabs_ec1[[ variant ]] ),
  af_non2023=mean( (metadata  %>% filter(year!=2023,clade=='T1-2-EC1') )$strain %in% kmer_presabs_ec1[[ variant ]] ),
  odds_ratio=af_2023/af_non2023,
  FST=1 - (portion_ec1_2023*af_2023*(1-af_2023) + (1-portion_ec1_2023)*af_non2023*(1-af_non2023))/(af*(1-af))
  #  portion_2023=mean(kmer_presabs_ec1[[ variant ]] %in% (v12_meta_1039 %>% filter(year==2023,clade=='T1-2 (EC1)') )$strain )
) %>% select(variant,af,af_2023,af_non2023,odds_ratio,FST,everything()) %>% 
  arrange(`lrt-pvalue`)-> gwas_ec1

write_tsv(gwas_ec1,file = 'output/ec1_gwas_2023.tsv')