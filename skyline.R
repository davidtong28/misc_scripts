#this script describes how to draw a skyline plot of different subpopulations in a time-based phylogeny using the skygrowth R package.
#We recommend running this script step-by-step in RStudio.

#load packages
library(phytools)
library(skygrowth)
library(ape)
library(ggplot2)


#First, obtain the time based trees. You may run BEAST or Bactdating for example. We used the BEAST tree in this example.
#internal population structure often confound the overall shape of a skyline plot, and is often recommended to split the whole population before plotting. There are a few options. 
#You may run BEAST separately on each subpopulation and plot them separately. This can be done in BEAST alone.
#In this demonstration, we chose to split a dated tree into subtrees.

#load metadata
read_tsv('data/MP_metadata_1121.tsv')%>% 
  as.data.frame() %>% 
  `rownames<-`(unlist(.[,1]))->metadata

#load the dated tree of the whole population, constructed by BEAST.
read.nexus('data/T12_treeanno.nex')->t12_dated

#determine the ancestral node for the EC1 subpopulation, and split the dated tree
mrca(t12_dated)['GCA_000733995','2023378']->node_EC1   #800
phytools::splitTree(t12_dated,list(node = node_EC1 , bp = 0))->t12_dated_split

#calculate effective population size variations over time for each subpopulation
pre1_dated_fit <- skygrowth.map( t12_dated_split[[1]], res = 100  )     #res specifies the number of time intervals
(pre1_dated_fit $time + max(metadata[metadata$clade=='T1-2-PRE1',"date"],na.rm = T) - max(pre1_dated_fit $time) ) -> pre1_dated_fit $time


ec1_dated_fit <- skygrowth.map( t12_dated_split[[2]], res = 100  )
(ec1_dated_fit $time+ max(metadata[metadata$clade=='T1-2-EC1',"date"],na.rm = T) - max(ec1_dated_fit $time) )-> ec1_dated_fit $time


#repeat the process for T2-2
read.nexus('data/T22_treeanno.nex')->t22_dated

mrca(t22_dated)['008_18M1176','420-2_17M1024']-> node_EC2  #160
phytools::splitTree(t22_dated,list(node = node_EC2 , bp = 0))->t22_dated_split

pre2_dated_fit <- skygrowth.map( t22_dated_split[[1]],res = 100  )
(pre2_dated_fit $time + max(metadata[metadata$clade=='T2-2-PRE2',"date"],na.rm = T) - max(pre2_dated_fit $time) )->pre2_dated_fit $time

ec2_dated_fit <- skygrowth.map( t22_dated_split[[2]],res = 100  )
(ec2_dated_fit $time + max(metadata[metadata$clade=='T2-2-EC2',"date"],na.rm = T) - max(ec2_dated_fit $time) )-> ec2_dated_fit $time



#merge all data together and draw the skyline plot
skyline_dated_merged<-rbind(        
  cbind('time'=pre1_dated_fit$time,pre1_dated_fit$ne_ci,label='T1-2-PRE1'),
  cbind('time'=ec1_dated_fit$time,ec1_dated_fit$ne_ci,label='T1-2-EC1'),
  cbind('time'=pre2_dated_fit$time,pre2_dated_fit$ne_ci,label='T2-2-PRE2'),        
  cbind('time'=ec2_dated_fit$time,ec2_dated_fit$ne_ci,label='T2-2-EC2')
) %>% 
  as.data.frame() %>% 
  mutate(across(1:4,as.numeric),
         label=factor(label,levels=unique(label)))

#input 95% HPD of each ancestral node, estimated by BEAST
root_HPDs<-data.frame(clade=c('T1-2-PRE1','T1-2-EC1','T2-2-PRE2','T2-2-EC2'),
                      x=c(pre1_dated_fit$time[1],ec1_dated_fit$time[1],pre2_dated_fit$time[1],ec2_dated_fit$time[1]),
                      xmin=c(1938.8069,1993.2346,1961.2220,2011.9660),
                      xmax=c(1962.0051,2001.1373,1983.6617,2015.1640),
                      y=c(pre1_dated_fit$ne[1],ec1_dated_fit$ne[1],pre2_dated_fit$ne[1],ec2_dated_fit$ne[1])
                      )

#draw skyline and HPDs in a single ggplot
skyline_dated_merged %>%
  rename('Clade'='label') %>% 
  ggplot()+
  geom_errorbar(mapping = aes(x=x,xmin=xmin,xmax=xmax,y=y),alpha=0.5,data = root_HPDs)+
  geom_line(mapping=aes(x=time,y=ne,color=Clade,linetype=Clade))+
  geom_ribbon(mapping=aes(x=time,ymin=nelb,ymax=neub,fill=Clade,alpha=Clade))+
  geom_errorbar(mapping = aes(x=x,xmin=xmin,xmax=xmax,y=y),data = root_HPDs)+
  scale_x_continuous(breaks = 5*( (1935/5):(2025/5) ))+
  scale_linetype_manual(values=c(1,5,1,5))+
  scale_colour_manual(values=c('#E793A8','#E793A8','#7a87f9','#7a87f9'))+
  scale_alpha_manual(values=c(0.2,0.4,0.2,0.4) )+
  scale_fill_manual(values=c('#E793A8','#E793A8','#7a87f9','#7a87f9'))+
  scale_y_log10()+
  ylab('Effective Population Size')+
  xlab('Year')+
  theme_classic()
ggsave(filename = 'output/skyline.pdf',width = 1250,height = 800,units='px',dpi = 150)
