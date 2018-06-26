##compare cancer genes for mutational data.


source("../../bin/WGSData.R")

##now store the file
#cfile='../../data/Census_allTue Jan 19 18-58-56 2016.csv'
#synStore(File(cfile,parentId='syn4984723'))

require(parallel)
require(tidyverse)
##now get
all.genes=unique(read.table('../../data/HugoGIDsToEntrez_DAVID.txt',sep='\t',header=T,as.is=T,quote='"')[,1])

impact=c("LOW",'MODERATE','HIGH')

allsoms <- synTableQuery("SELECT * FROM syn12555329 where parentId='syn5578958'") %>% as.data.frame()
print(paste('Selecting from',nrow(allsoms),'mutation files'))
allsoms=allsoms[which(!gsub('CT0+','',allsoms$patientId)%in%c("10","13")),]
print(paste('Removing patient 10 and 13 to get',nrow(allsoms),'mutation files'))

allsoms=allsoms[unlist(sapply(impact,grep,allsoms$name)),]
print(paste("Found",nrow(allsoms),'with',paste(impact,collapse=' or '),'impact'))
som.germ=getAllMutData(allsoms,filt=c('PASS'))

require(pbmcapply)
allstats<-pbmclapply(as.character(all.genes),function(x) try(getMutationStatsForGene(gene=x, som.germ = som.germ, filt=c('PASS'), impact = "HIGH", redo=T)), mc.cores = 8)
  
names(allstats)<-as.character(all.genes)

fulldf<-data.frame(Hugo_Symbol=unlist(sapply(allstats,function(x) as.character(x$Hugo_Symbol))),
                   Protein_Change=unlist(sapply(allstats,function(x) as.character(x$Protein_Change))),
                   Sample_ID=unlist(sapply(allstats,function(x) as.character(x$Sample_ID))),
                   Mutation_Status=unlist(sapply(allstats,function(x) as.character(x$Mutation_Status))),
                   Chromosome = unlist(sapply(allstats,function(x) as.character(x$Chromosome))),
                   Start_Position = unlist(sapply(allstats,function(x) as.character(x$Start_Position))),
                   End_Position = unlist(sapply(allstats,function(x) as.character(x$End_Position))),
                   Reference_Allele = unlist(sapply(allstats,function(x) as.character(x$Reference_Allele))),
                   Variant_Allele = unlist(sapply(allstats,function(x) as.character(x$Variant_Allele))),
                   Mutation_Type = unlist(sapply(allstats,function(x) as.character( x$Mutation_Type))),
                   ExAC_AF = unlist(sapply(allstats,function(x) as.character( x$ExAC_AF))))

 udf<-distinct(fulldf) 

 write.table(udf,file='allGeneMutationsInDermalsFiltered.tsv',sep='\t',row.names=F,quote=F)
write.table(subset(udf,Mutation_Status=="Germline"),file='germlineAllGeneMutationsInDermalsFiltered.tsv',sep='\t',row.names=F,quote=F)
write.table(subset(udf,Mutation_Status=="Somatic"),file='somaticAllGeneMutationsInDermalsFiltered.tsv',sep='\t',row.names=F,quote=F)

for(f in c("allGeneMutationsInDermalsFiltered.tsv","germlineAllGeneMutationsInDermalsFiltered.tsv","somaticAllGeneMutationsInDermalsFiltered.tsv")){
    synStore(File(f,parentId='syn5605256'),used=list(list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/analysis/2016-03-07/mutsWithoutComms.R',wasExecuted=TRUE)))
}
