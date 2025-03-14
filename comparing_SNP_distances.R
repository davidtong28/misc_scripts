#This is an R script that compares and test for significant differences between two groups of sample pairs in a population using Kruskal-Wallis test
#We recommend running this script step-by-step on RStudio

#load R packages

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

#Load metadata file

read_tsv('data/MP_metadata_1121.tsv')%>% 
  as.data.frame() %>% 
  `rownames<-`(unlist(.[,1]))->metadata

#First run SNP-dists on your alignment file. Then load results from SNP-dists using:

read_tsv('data/snp_dist_2023_ec1.tsv')->ec1_snp_dist
ec1_snp_dist[,-1] %>% 
  as.data.frame() %>% 
  `rownames<-`(unlist(ec1_snp_dist[,1]))->ec1_snp_dist

#make sure row and column names are identical:
all(rownames(ec1_snp_dist)==colnames(ec1_snp_dist) )


#Alternatively, provide a phylogenetic tree to generate SNP distance matrix:
# subtrees(read.tree('data/2023_ec1_rooted_long.tre') %>% 
#              ladderize())[[2]] -> ec1_2023
# 
# ec1_2023 %>% ape::cophenetic.phylo() ->ec1_phylo_dist


#create a second matrix that classifies sample pairs into two catagories, in this case, same/different locations


ec1_2023_geo<-matrix(NA, 
                     nrow =  nrow(ec1_snp_dist) , 
                     ncol =  ncol(ec1_snp_dist) )%>% 
  `rownames<-`(rownames(ec1_snp_dist)) %>%
  `colnames<-`(rownames(ec1_snp_dist))

for (i in 1:nrow(ec1_snp_dist)) {
  for (j in 1:nrow(ec1_snp_dist)) {
    ec1_2023_geo[i, j] <- ifelse(metadata[rownames(ec1_snp_dist),"location"][i] == metadata[rownames(ec1_snp_dist),"location"][j], "same", "different")
  }
}

#The two matrices should share identical row and column names, if not, change the order accordingly
all(rownames(ec1_snp_dist)==rownames(ec1_2023_geo) ) 
all(colnames(ec1_snp_dist)==colnames(ec1_2023_geo) ) 

#remove self-pairs
diag(ec1_2023_geo)<-'NA'

#merge the information of the 2 matrices to form a list of SNP distances
ec1_2023_dist_list<-rbind(data.frame('dist'=ec1_snp_dist[ec1_2023_geo =='same'],'location'='same'),
      data.frame('dist'=ec1_snp_dist[ec1_2023_geo =='different'],'location'='different'))

#draw violin plots demonstrating the distribution of 
ec1_2023_dist_list %>% 
  ggplot() +
  geom_violin(aes(x=location,y=dist),draw_quantiles = c(0.25,0.5,0.75),fill='grey')+
  ggtitle('Distribution of SNP distances between pairs of same and different location isolates')+
  ylab("SNP distances")+
  theme_classic()

ggsave(filename = 'output/geo_violin.pdf',width = 1500,height = 800,units='px',dpi = 150)

#conduct Kruskal-Wallis test on the two groups
kruskal.test(
  list(ec1_snp_dist[ec1_2023_geo =='same'],
       ec1_snp_dist[ec1_2023_geo =='different']) 
  )

#noticing an enrichment of short-distance SNPs in the same-location group, we further examine the existence of local dissemination by conducting a Fisher's exact test
#first split each group into =<7 and >7 SNP distance (selected based on the violin plot) and construct a 2x2 contingency table
ec1_2023_dist_list %>% mutate(group=case_when(dist<=7 ~ 'close',
                                              dist>7 ~ 'far'  ) 
                              ) %>%
  group_by(location,group) %>% 
  summarise(count=n(),.groups = 'keep') %>% 
  pivot_wider(names_from = group,values_from = count) %>% as.data.frame() %>% 
  `rownames<-`(.$location) %>% 
  select(-1) %>% 
  as.matrix()->SNP_contingency_table

#examine if short-distance sample pairs are distributed differently across location groups
fisher.test(SNP_contingency_table)

