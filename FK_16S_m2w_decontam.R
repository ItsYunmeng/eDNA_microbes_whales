#Decontamination Script for M2W; RPK 2017, last edited by CJC on 20180119, CJC on 20180312; RPK 20180425

library(here); library(plyr); library(biomformat); library(tidyverse); library(phyloseq);  library(vegan); library(reshape2); library(MASS)
bit=function(x){print(x[c(1:10),c(1:10)]); print(dim(x))}  #shows top-left corner of a data.frame/matrix; useful for large datasets
colMax <- function(data) sapply(data, max, na.rm = TRUE) #get column-wise maximum value
Col2RN<-function(df, x){row.names(df)<-df[,x]; df<-df[-x]; return(df)} #convert column to row.names
rowVar<-function(data) apply(data, 1, FUN=var, na.rm=T) #get row-wise variance of a df

FK_16S_biom <- import_biom("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/biom/M2W_16S_FK_json_obs_md.biom",parseFunction = parse_taxonomy_default)
#sample_names(FK_16S_biom)
load("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/OccupancyModeling/Florida/ProbOcc16s_dataOUT14_24_38_FK.Rdata")
FK_16S <- merge_phyloseq(FK_16S_biom@otu_table, FK_16S_biom@tax_table,sample_data(read.csv("/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/metadata/m2w_16S_FK_sam_data.csv", row.names=2, check.names = FALSE)))
  (READcount0<-sum(sum(FK_16S@otu_table)))
  (OTUcount0<-nrow(FK_16S@otu_table))

hist(dataOUT$probOcc)

FK_16S_prob <- subset(dataOUT,probOcc >0.8)
FK_16S_unique_OTUs <- as.character(unique(FK_16S_prob$OTUname))
FK_16S_pruned <- prune_taxa(FK_16S_unique_OTUs, FK_16S)

#a <- as.data.frame(FK_16S@tax_table)
#FK_16S_unique_OTUs

allOTUs<-as.data.frame(otu_table(FK_16S_pruned, taxa_are_rows=TRUE))
(READcount1<-sum(sum(allOTUs)))
(OTUcount1<-nrow(allOTUs))
Meta<-as.data.frame(FK_16S_pruned@sam_data)


###################################################
#interal DECONTAMINATION STEP using controls 
###################################################
samples_all <- allOTUs
CONTROL<-Meta[which(grepl("EB|CB|NTC|blank", row.names(Meta))),]
ArtComm<-Meta[which(grepl("MOCK", row.names(Meta))),]  #select relevant samples from the metadata

samples_con <- samples_all[,colnames(samples_all)%in%row.names(CONTROL)] #controls dataset
samples_art <- samples_all[,colnames(samples_all)%in%row.names(ArtComm)] #artComm dataset
samples_env <- samples_all[,!colnames(samples_all)%in%c(colnames(samples_con), colnames(samples_art))] #environmental samples dataset

contam_proportion_df=data.frame(samples_con) #set up df for proportion contamination calculation

contam_proportion_df=scale(contam_proportion_df, center=F, scale=colSums(contam_proportion_df)) #for each dup/OTU in controls, calc frequency
contam_proportion_df[is.na(contam_proportion_df)]<-0 
#NOTE from Collin [although I'm not sure this still applies]: the above line doesn't make sense for 12s, since there are almost no reads in the ArtComms. Better to estimate proportion of contam based on N reads / OTU in the NTCs, then scale this by the mean nuMBer of reads in the environmental samples, and use the result as a portion of likely contamination.
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
DECONTAM <- DECONTAM[!(rownames(DECONTAM) %in% names(which(rowSums(samples_con)>rowSums(samples_env)))),]

#auxiliary cleanup
#DECONTAM<-DECONTAM[,colSums(DECONTAM)>10000] #low abundance outlier
DECONTAM <-DECONTAM[order(rowSums(DECONTAM), decreasing=T),]

(READcount2<-sum(sum(DECONTAM)))
(OTUcount2<-nrow(DECONTAM))


single_samples <- as.data.frame(cbind(c(DECONTAM$WS0116a),c(DECONTAM$WS0615a), c(DECONTAM$WS0415b),c(DECONTAM$WS0415c)))
row.names(single_samples) <- row.names(DECONTAM)
single_samples$a <- row.names(DECONTAM)

#single_samples <- cbind(c(DECONTAM$WS0116a),c(DECONTAM$WS0615a))
colnames(single_samples) <- c("WS0116a","WS0615a","WS0415b","WS0415c", "a")
#colnames(single_samples) <- c("WS0116a","WS0615a")

DECONTAM <- DECONTAM[,-c(15,20, 24,25)] #Remove single samples
#DECONTAM <- DECONTAM[,-c(15,20)] #Remove single samples

#######################################################
#Third decontamination: Outlier removal, for PCR replicates that are very different from their sister replicates
#######################################################
(replicates<-substr(colnames(DECONTAM),1,nchar(colnames(DECONTAM))-1)) #truncate names to group unique bottles of water
distmat<-vegdist(t(DECONTAM))
distList=list(NA) ; distList.tri=list(NA); index=1
for (i in unique(replicates)){
  rowMatch<-which(replicates%in%i)
  distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
  distList.tri[[index]]<-	distList[[index]][upper.tri(distList[[index]])]
  index=index+1
}

#visualize, if desired 
#hist(unlist(distList.tri))

normparams=fitdistr(unlist(distList.tri), "normal")$estimate  #fit distribution to bray-curtis dissimilarities. lognormal, beta, normal, etc, may have less-good fits
probs=pnorm(unlist(distList.tri), normparams[1], normparams[2])
outliers =which(probs>0.95)
#minOutlier<-min(unlist(distList.tri)[outliers]) #minimum outlier value
minOutlier<-0.5  #here, if desired, instead of fitting a distribution you can a hard cutoff because of the distribution visualized

#remove outliers
#distList<-distList[-which(lapply(distList, length)==1)] 	
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

#DECONTAM <- cbind(DECONTAM, single_samples)
#install.packages("compare")
#library(compare)
#comparison <- compare(DECONTAM,single_samples,allowAll=TRUE)
library(dplyr)

#join samples that did not have triplicates!
DECONTAM$a <- row.names(DECONTAM)
joined <- inner_join(DECONTAM, single_samples, by = c("a"))
row.names(joined) <- joined$a
DECONTAM <- joined[,-c(15)] #remove a

DECONTAM_biom_FK_16S<-prune_taxa(row.names(DECONTAM), FK_16S_pruned)%>%
  prune_samples(colnames(DECONTAM), .)

sample_names(DECONTAM_biom_FK_16S)
FK_16S_m2w <- subset_samples(DECONTAM_biom_FK_16S, description == "m2w")
sample_names(FK_16S_m2w)

write.csv(FK_16S_m2w@otu_table, "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/post_decont/FK/16S/FK_16S_otu_table.csv")  #is this the best way to save biom files?
write.csv(FK_16S_m2w@tax_table, "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/post_decont/FK/16S/FK_16S_taxa_table.csv")  #is this the best way to save biom files?
write.csv(FK_16S_m2w@sam_data, "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/post_decont/FK/16S/FK_16S_sam_table.csv")  #is this the best way to save biom files?


#DECONTAM<-data.frame(DECONTAM, allOTUs[match(row.names(DECONTAM), row.names(allOTUs)),64:65])
#write.csv(DECONTAM, "M2W_DECONTAM_12s_FK_20171101.csv")

SummaryTableOut_FK_16S<-t(data.frame(c(READcount0, OTUcount0),c(READcount1, OTUcount1),c(READcount2, OTUcount2),c(READcount3, OTUcount3)))
colnames(SummaryTableOut_FK_16S)<-c("Reads","OTUs")
row.names(SummaryTableOut_FK_16S)<-c("Raw", "PostOccupancy","PostSubtraction","final")

#######################################################

