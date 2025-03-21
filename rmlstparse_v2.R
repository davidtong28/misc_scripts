#!/usr/bin/env Rscript

#first run rmlst.sh on a desired set of assembled genomes in the bash terminal. Gather all JSON results to a folder used as the input (-f) of this script.
#We recommend running this script on the bash terminal using the Rscript command for efficiency. R and required R packages must be installed in the environment.

library(dplyr)
library(jsonlite)
library(stringr)
library(argparse)

#specifying command line arguments
parser<-ArgumentParser()

parser$add_argument("-f","--folderpath",type="character",help="path to the input folder that contains all rmlst json files")


parser$add_argument("-o","--outdir",type="character",help="output directory folder")



#generate `_sum.csv` file
report_list<-function(json_path,isolate_name="isolate"){
  rmlst_json=fromJSON(json_path)
  fields<-if( (rmlst_json$fields %>% as.data.frame %>% nrow)==0){c(NA,NA,NA)}else{
    c(rmlst_json$fields$rST %>% paste0(collapse = "; "),
      rmlst_json$fields$genus %>% paste0(collapse = "; "),
      rmlst_json$fields$species %>% paste0(collapse = "; ") )}
  taxon=if( (rmlst_json$taxon_prediction %>% as.data.frame %>% nrow)==0)
  {c("unknown",NA,NA,NA,NA)}else{
    c(ifelse((rmlst_json$taxon_prediction %>% as.data.frame %>% nrow)==1,"pure","admixed"),
      rmlst_json$taxon_prediction$taxon %>% paste0(collapse = "; "),
      rmlst_json$taxon_prediction$support %>% paste0(collapse = "; "),
      rmlst_json$taxon_prediction$rank %>% paste0(collapse = "; "),
      rmlst_json$taxon_prediction$taxonomy %>% paste0(collapse = "; ") ) 
  }
  
  alleles=if( (rmlst_json$exact_matches) %>% length==0){c(0,0,0)}else{
    c( (rmlst_json$exact_matches) %>% length,
       rmlst_json$exact_matches %>% lapply(function(x)nrow(x)>1) %>% unlist %>% sum,
       rmlst_json$exact_matches %>% lapply(FUN=function(x)nrow(x)) %>% unlist %>% mean )
  } 
  
  
  report=c(isolate_name,fields,taxon,alleles) %>% as.data.frame() %>% t %>%
    `colnames<-`(c("name","rST","genus","species","purity","taxon","support","rank","taxonomy","identified_loci","duplicated_loci","alleles_per_loci") ) %>% 
#    `row.names<-`(isolate_name)
  return(report)
}


#generate `_allele.csv` file
allele_table<-function(json_path){
  rmlst_json=fromJSON(json_path)
  exact_match<-rmlst_json$exact_matches 
  exact_match0<-exact_match %>% `names<-`(NA)
  exact_match2<-list()
  for (x in c(1:length(exact_match))){
    exact_match2[x]<-list( exact_match0[x] %>% as.data.frame() %>% mutate(allele=names(exact_match)[x])  )
  }
  exact_match3<-exact_match2%>% rbind_pages() %>% as.matrix() %>% as.data.frame()
  
  exact_match3$NA.linked_data<-exact_match3$NA.linked_data %>% lapply(function(x){as.data.frame(x) %>% select("value","frequency") %>% arrange(-frequency) %>% unlist %>%  tryCatch(error=function(x)NA) }) %>% paste0
  exact_match3 <-exact_match3 %>% rename_all(function(x)str_remove(x,"NA.")) %>% mutate_if(is.list,as.character)
  return(exact_match3 %>% select("allele","allele_id","linked_data","contig","length","start","end","orientation"))
}


#function for execution
single_report_gen=function(jsonpath,name,outdir=paste0( str_remove(jsonpath,"[^/]*$"),"rmlst_",name) ){
  
  ifelse(dir.exists(outdir),T,dir.create(outdir) )
  summary<-report_list(jsonpath,name)
  write.csv(summary,file = paste0(outdir,"/",name,"_sum.csv"),row.names = F )
  allele<-allele_table(jsonpath)
  write.csv(allele,file = paste0(outdir,"/",name,"_allele.csv") )
}



#get arguments
args <- parser$parse_args()

listjson<-list.files(path=args$f,pattern='.json',full.names=T)
listjsonshort<-list.files(path=args$f,pattern='.json',full.names=F)

for (x in 1:length(listjson)){
  single_report_gen(listjson[x],str_remove(listjsonshort[x],'.json'),args$outdir)
}
