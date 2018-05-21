#Decontamination Script for M2W; RPK 2017, last edited by CJC on 20180119, CJC on 20180312; RPK 20180425

library(here); library(plyr); library(biomformat); library(tidyverse); library(phyloseq);  library(vegan); library(reshape2); library(MASS)
bit=function(x){print(x[c(1:10),c(1:10)]); print(dim(x))}  #shows top-left corner of a data.frame/matrix; useful for large datasets
colMax <- function(data) sapply(data, max, na.rm = TRUE) #get column-wise maximum value
Col2RN<-function(df, x){row.names(df)<-df[,x]; df<-df[-x]; return(df)} #convert column to row.names
rowVar<-function(data) apply(data, 1, FUN=var, na.rm=T) #get row-wise variance of a df

FK_12S_biom <- import_biom("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/biom/FK_12S_NT_json_obs_md.biom",parseFunction = parse_taxonomy_default)
load("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/OccupancyModeling/Florida/ProbOcc12s_dataOUT14_13_19_FK.Rdata")
FK_12S <- merge_phyloseq(FK_12S_biom@otu_table, FK_12S_biom@tax_table,sample_data(read.csv("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/metadata/m2w_12S_FK_sam_data.csv", row.names=2)))
(READcount0<-sum(sum(FK_12S@otu_table)))
(OTUcount0<-nrow(FK_12S@otu_table))

FK_12S_prob <- subset(dataOUT,probOcc >0.8)
FK_12S_unique_OTUs <- as.character(unique(FK_12S_prob$OTUname))
FK_12S_pruned <- prune_taxa(FK_12S_unique_OTUs, FK_12S)

a <- prune_taxa(taxa_sums(FK_12S_pruned) > 0, FK_12S_pruned)
b <- as.data.frame(FK_12S_pruned@otu_table)
c <- as.data.frame(FK_12S@otu_table)

#pdf(file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/OccupancyModeling/Florida/hist_12S_prob.pdf")
#hist(dataOUT$probOcc)
#dev.off()

allOTUs<-as.data.frame(otu_table(FK_12S_pruned, taxa_are_rows=TRUE))
taxa<-as.data.frame(tax_table(FK_12S_pruned))
(READcount1<-sum(sum(allOTUs)))
(OTUcount1<-nrow(allOTUs))
Meta<-as.data.frame(FK_12S_pruned@sam_data)


###################################################
#interal DECONTAMINATION STEP using controls 
###################################################
samples_all <- allOTUs
CONTROL<-Meta[which(grepl("EB|CB|NTC|blank", Meta$sample_name)),]
ArtComm<-Meta[which(grepl("ArtComm", row.names(Meta))),]  #select relevant samples from the metadata

samples_con <- samples_all[,colnames(samples_all)%in%row.names(CONTROL)] #controls dataset
samples_art <- samples_all[,colnames(samples_all)%in%row.names(ArtComm)] #artComm dataset
samples_env <- samples_all[,!colnames(samples_all)%in%c(colnames(samples_con), colnames(samples_art))] #environmental samples dataset

contam_proportion_df=data.frame(samples_con) #set up df for proportion contamination calculation

contam_proportion_df=scale(contam_proportion_df, center=F, scale=colSums(contam_proportion_df)) #for each dup/OTU in controls, calc frequency
contam_proportion_df[is.na(contam_proportion_df)]<-0 
#NOTE from Collin [although I'm not sure this still applies]: the above line doesn't make sense for 12s, since there are almost no reads in the ArtComms. Better to estimate proportion of contam based on N reads / OTU in the NTCs, then scale this by the mean number of reads in the environmental samples, and use the result as a portion of likely contamination.
#contam_proportion_df=scale(contam_proportion_df, center=F, scale=rep(mean(colSums(samples_env)), times=9))

contam_proportion_vec =apply(contam_proportion_df, 1, max)  #take max of the frequencies of probable contam for each OTU.  this is the max proportional contribution of contamination to each of these OTUs/dups in the dataset. Accordingly, you will then subtract this proportion from each dup/OTU occurrence in the dataset.  A different (better?) way to do this would be to fit a likelihood distribution to the observed proportions of the dups/OTUs in the controls, and come up with the most likely proportion attributable to contamination.  But, in most cases, the most likely contribution is at or near zero (using a beta distrib), and in many cases the models don't converge anyway.


#Step 1: remove the fraction of each OTU abundace attributable to contamination
DECONTAM = samples_env
for (i in 1:ncol(DECONTAM)){
  nreads=sum(DECONTAM[,i])
  contam_proportion_vec_forSample= ceiling(contam_proportion_vec*nreads)
  DECONTAM[,i]= samples_env[,i]-contam_proportion_vec_forSample
}
DECONTAM[DECONTAM <0]<-0  #assign zero to any reads for which there are a negative (corrected) nuMBer of reads
DECONTAM = DECONTAM[rowSums(DECONTAM)>0,]  #remove reads that no longer occur in the data



#Step 2: remove control taxon OTU from environmental samples in which it occurs
#DECONTAM <- DECONTAM[!(rownames(DECONTAM) %in% control_taxon),]

#NOTE: if your control taxon or ArtComm includes species that are truly in the sampled environment,
#remove any other OTU that *predominantly* occurs in control samples 
#tail(samples_con[order(rowSums(samples_con)),], 20)  #visualize, if desired
#DECONTAM <- DECONTAM[!(rownames(DECONTAM) %in% names(which(rowSums(samples_con)>rowSums(samples_env)))),]

#auxiliary cleanup
#DECONTAM<-DECONTAM[,colSums(DECONTAM)>10000] #low abundance outlier
#DECONTAM <-DECONTAM[order(rowSums(DECONTAM), decreasing=T),]

(READcount2<-sum(sum(DECONTAM)))
(OTUcount2<-nrow(DECONTAM))

distmat
#######################################################
#Third decontamination: Outlier removal, for PCR replicates that are very different from their sister replicates
#######################################################
(replicates<-substr(colnames(DECONTAM),1,nchar(colnames(DECONTAM))-1)) #truncate names to group unique bottles of water
distmat<-vegdist(t(DECONTAM))
distList=list(NA) ; distList.tri=list(NA); index=1
for (i in unique(replicates)){
  rowMatch<-which(replicates%in%i)
  distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
  distList.tri[[index]]<-	distList[[index]][upper.tri(	distList[[index]])]
  index=index+1
}

#visualize, if desired 
#hist(unlist(distList.tri))

normparams=fitdistr(unlist(distList.tri), "normal")$estimate  #fit distribution to bray-curtis dissimilarities. lognormal, beta, normal, etc, may have less-good fits
probs=pnorm(unlist(distList.tri), normparams[1], normparams[2])
outliers =which(probs>0.95)
minOutlier<-min(unlist(distList.tri)[outliers]) #minimum outlier value
#minOutlier<-0.5  #here, if desired, instead of fitting a distribution you can a hard cutoff because of the distribution visualized

#remove outliers
distList<-distList[-which(lapply(distList, length)==1)] 	
namesOutliers=list(NA)
for (i in 1:length(distList)){
  namesOutliers[[i]]<-intersect(
    names(which(colSums(distList[[i]]>=minOutlier)>0)),
    names(which.max(rowMeans(distList[[i]])))
  )
}
DECONTAM<-DECONTAM[,-match(unlist(namesOutliers),colnames(DECONTAM))]
DECONTAM<-DECONTAM[rowSums(DECONTAM)>0,]

(READcount3<-sum(sum(DECONTAM)))
(OTUcount3<-nrow(DECONTAM))

############write out
library(magrittr)
DECONTAM_biom_FK_12S<-prune_taxa(row.names(DECONTAM), FK_12S_pruned)%>%
  prune_samples(colnames(DECONTAM), .)

sample_names(DECONTAM_biom_FK_12S)
FK_12S_m2w <- subset_samples(DECONTAM_biom_FK_12S, description == "m2w")
sample_names(FK_12S_m2w)

write.csv(FK_12S_m2w@otu_table, "/Users/anni_djurhuus/Documents/Projects/FKON/m2w/google_drive/data/MEGAN6/post_decont/FK/12S/FK_12S_otu_table.csv")  #is this the best way to save biom files?


DECONTAM<-data.frame(DECONTAM, allOTUs[match(row.names(DECONTAM), row.names(allOTUs)),64:65])
#write.csv(DECONTAM, "M2W_DECONTAM_12s_FK_20171101.csv")

SummaryTableOut_FK_12S<-t(data.frame(c(READcount0, OTUcount0),c(READcount1, OTUcount1),c(READcount2, OTUcount2),c(READcount3, OTUcount3)))
colnames(SummaryTableOut_FK_12S)<-c("Reads","OTUs")
row.names(SummaryTableOut_FK_12S)<-c("Raw", "PostOccupancy","PostSubtraction","final")

#######################################################