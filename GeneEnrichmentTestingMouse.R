#gene overlaps with MGEnrichment database
#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#2025
##############################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(lsmeans)
library(nlme)
library(cowplot)
library(rJava)
library(xlsx)
library("GeneOverlap")
##############################################################################################################


setwd("/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/MGenrichment")

#####################################################
#read in DEG GR list
#####################################################

DEG_list <- read.csv("/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/MGenrichment/GRlist.csv")


#clean symbols:
DEG_list$symbol <- gsub(";.*$", "", DEG_list$symbol)
DEG_list$symbol <- trimws(DEG_list$symbol)

genelists <- DEG_list


unique(genelists$groups)
unique(genelists$listname)
unique(genelists$source)
genelists$listname2 <- paste(genelists$listname,genelists$source)
genelists$mgi_symbol <- genelists$symbol

List <- split(genelists$mgi_symbol, genelists$listname2)

List <- lapply(List, function(x) Filter(Negate(is.na), x))
  
#####################################################
#read in gene lists 
#####################################################

#DEGs from all samples for CTL vs ELA for both sexes
Alldata_both <- read.csv("/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/2025-10-12-PVNSeq/CTLvsELA_both\ _DEgenes.csv")
Alldata_both <- Alldata_both %>% filter(adj.P.Val < 0.05)

Alldata_both_gs <- na.omit(Alldata_both$genesymbol)

Alldata_both_gs <- as.character(Alldata_both_gs)

#####################################################
#all detected RNAseq genes for background
#####################################################

bg <- read.csv("/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/2025-10-12-PVNSeq/CTLvsELA_both\ _DEgenes.csv")

mm10genome <- length(unique(na.omit(bg$genesymbol)))


#####################################################
#overlaps
#modifiy function from MGEnrichment to take new input list
#####################################################


Overlap_fxn <- function(targetlistname,genelists,genomesize){
  out <- NULL
  target <- targetlistname
  inputlistname <- c("GRlist")
  
  for (i in 1:length(genelists)) { 
    
    #call gene overlaps
    go.obj <- newGeneOverlap(target,
                             genelists[[i]],
                             genome.size=genomesize)
    
    #perform test
    go.obj <- testGeneOverlap(go.obj) #returns onetailed pvalue
    #return odds ratio:
    OR <- getOddsRatio(go.obj)
    pvalue <- getPval(go.obj)
    
    #extract contingency table
    CT <- getContbl(go.obj)
    notAnotB <- CT[1,1]
    inAnotB <- CT[1,2]
    inBnotA <- CT[2,1]
    inBinA <- CT[2,2]
    
    CTlist <- cbind(notAnotB,inAnotB,inBnotA,inBinA)
    
    
    #two sided fisher's exact test
    #test <- fisher.test(CT,alternative='two.sided')
    
    #get gene list B:
    intersection <- go.obj@intersection
    intersection <- paste(as.character(intersection),collapse=", ",sep="")
    
    #get listname
    listname <- paste(names(genelists[i]))
    
    results <- cbind(listname,pvalue,OR, CTlist,intersection)
    
    names(results) <- c("listname","pvalue","OR","notAnotB","inAnotB","inBnotA","inBinA","geneids")
    out <- rbind(out,results) 
    
  }
  
  #remove first row as overlap with self:
  out <- as.data.frame(out)
  out2 <- out #[- grep(inputlistname, out$listname),] #not needed b/c targets not in database
  
  #add in targetlist name (assumes first list is the input)
  out2$targetlist <- paste(inputlistname)
  
  
  rownames(out2) <- NULL
  
  #return results
  return(out2)
}


#example from database:
targetlistname <- Alldata_both_gs
genelists <- List
genomesize <- mm10genome

test <- Overlap_fxn(targetlistname,List,genomesize)






######################################################################
#FDR correction
######################################################################

#filter for desired lists:
dfoverlaps <- test

dfoverlaps$pvalue <- as.numeric(paste(dfoverlaps$pvalue))

#adjust matrix pvalues by FDR
dfoverlaps$FDR <- p.adjust(dfoverlaps$pvalue, method='fdr')

#write.csv(dfoverlaps, "MGEnrichment2024.csv")


#add back in database info
df <-  DEG_list %>% 
dplyr:: select(-symbol) %>%
  distinct()

df$listname2 <- paste(df$listname,df$source)


dfoverlaps2 <- merge(dfoverlaps, df, by.x = "listname", by.y = "listname2", all.x=T)
dfoverlaps2 <- distinct(dfoverlaps2)



#A = target list
#B = database gene list

# % of target list found in the database gene list: overlap/size of A *100
#inBinA/(inAnotB + inBinA)*100

dfoverlaps2$inBinA <- as.numeric(dfoverlaps2$inBinA)
dfoverlaps2$inAnotB <- as.numeric(dfoverlaps2$inAnotB)
dfoverlaps2$inBnotA <- as.numeric(dfoverlaps2$inBnotA)

dfoverlaps2$Percent_target_in_database <- dfoverlaps2$inBinA/(dfoverlaps2$inAnotB + dfoverlaps2$inBinA)*100

#percent of target list shared with database gene list
dfoverlaps2$Percent_shared <- dfoverlaps2$inBinA/(dfoverlaps2$inAnotB + dfoverlaps2$inBinA + dfoverlaps2$inBnotA)*100


#filter for significant
Sig <- dfoverlaps2 #%>% filter(FDR < 0.05) %>% filter(groups != "Human Brain Disorders")
Sig$listname <- as.factor(Sig$listname.y)
Sig$targetlist <- as.factor(Sig$targetlist)
Sig$neg_log_pvalue <- -log(Sig$FDR)

write.csv(Sig, "GRlist_Enrichment2026.csv")


