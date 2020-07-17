source("Illumina_functions.R")
source("pipeline_functions.R")
library(swamp)
library(WriteXLS)
library(limma)
library(sva)

## reading in the files, because maimages reads only Cy3 and Cy5 "other" column needs to be developed for Alexa
## genepix.median is the image analysis program
## other.columns is character vector of names of other columns to be read containing spot-specific information
dat <- read.maimages(dir()[grep(".gpr", dir())], source="genepix.median", other.columns=c(B="F488 Median", Rsd="F635 SD", Gsd="F532 SD", Bsd="F488 SD"))
gal <- read.delim("gal_to_GEO.txt", sep="\t", header=T)
# last row of the gal file is empty
gal = gal[-nrow(gal),]

# phenodata - pheno is uploaded again later using the eUTOPIA function
pd <- read.table("GSE92899_metadata.txt", sep="\t", header=T)

#gather the CyDyes together
idxR = which(pd$dye == "Cy5")
idxG = which(pd$dye == "Cy3")
idxO = which(pd$dye == "a488")

dat2 <- cbind(dat$R[,pd$filenames[idxR]], dat$G[,pd$filenames[idxG]], dat$other[[1]][,pd$filenames[idxO]])
rownames(dat2) = dat$genes$ID

sdmat2 <- cbind(dat$other[[2]][,pd$filenames[idxR]],dat$other[[3]][,pd$filenames[idxG]],dat$other[[4]][,pd$filenames[idxO]])
rownames(dat2) = dat$genes$ID

colnames(dat2) <- pd$GSM
par(mar=rep(4,4))
boxplot(log2(dat2), las=2, cex=0.7)

### sample negCtrl 6h area4 has problems, rmv the sample from all the data
dat3 <- dat2[,-which(colnames(dat2) %in% "GSM2439596")]
sdmat3 <- sdmat2[,-which(colnames(dat2) %in% "GSM2439596")]
pd3 <- pd[-which(colnames(dat2) %in% "GSM2439596"),]

boxplot(log2(dat3), las=2, cex=0.7)

##take only the probe intensities to work with, leaving out the neg and pos controls
#probeID = gdata::startsWith(str = dat$genes$ID,pattern = "A_")
dat4 <- dat3[which(gal$ControlType=="probe"),]
sdmat4 <- sdmat3[which(gal$ControlType=="probe"),]

#dat4 <- dat3[which(probeID),]
#sdmat4 <- sdmat3[which(probeID),]
gal4 <- gal[which(gal$ControlType=="probe"),]
boxplot(log2(dat4), las=2, cex=0.4)

##############
# 1st filter #
##############


##build empty matrix and give a score 1 (good) if the probes intensity is in quantile of 75% *1,5 greater than negative control to execute the real intensities from the background
score.neg <- matrix(0, nrow=dim(dat4)[1], ncol=dim(dat4)[2])
for(i in 1:dim(dat4)[2]){
  thres <- NULL
  thres <- quantile(dat3[which(gal$ControlType=="neg"),i],.50)*1.5
  #thres <- quantile(dat3[which(probeID==FALSE),i],.75)*1.5
  score.neg[which(dat4[,i]>thres),i] <- 1
}

## 1st filter for background 
##sum rows horizontally to get a score with min. 0 and max. 59 (no. of samples)
neg.sumrow <- apply(score.neg, 1, sum)
perc <- 50
dat5 <- dat4[neg.sumrow>=round(dim(dat4)[2]*perc/100),]
sdmat5 <- sdmat4[neg.sumrow>=round(dim(dat4)[2]*perc/100),]
gal5 <- gal4[neg.sumrow>=round(dim(dat4)[2]*perc/100),]
boxplot(log2(dat5), las=2, cex=0.7)

##############
# 2nd filter #
##############

## build another empty matrix and give the score 1 if the sd of the samples having smaller value than threshold
##
score.sd <- matrix(0, nrow=dim(dat5)[1], ncol=dim(dat5)[2])
for(i in 1:dim(dat5)[2]){
  thres <- NULL
  thres <- quantile(sdmat5[,i],.90)
  score.sd[which(dat5[,i]<thres),i] <- 1
}

### 2nd filter for problems in probes intensities standard deviation

sd.sumrow <- apply(score.sd, 1, sum)
dat6 <- dat5[sd.sumrow>=round(dim(dat5)[2]*perc/100),]
gal6 <- gal5[sd.sumrow>=round(dim(dat5)[2]*perc/100),]
boxplot(log2(dat6), las=2, cex=0.7)



##########################
# normalizing and ComBat #
##########################


##normalizing the probedata
norm <- normalizeQuantiles(log2(dat6))
boxplot(norm, las=2, cex=0.7)
#plotMDS(norm, gene.selection="common", col=as.numeric(pd3$treatment))

# Probe annotation

annot <- read.delim("AllAnnotations_026652_D_AA_20190204.txt", sep="\t", header=T, row.names=1, stringsAsFactors = FALSE)
probes <- as.character(unique(gal6$ID))
p2eg <- matrix(nrow=length(probes), ncol=2)
colnames(p2eg) <- c("probeID", "EnsemblID")


for(i in 1:length(probes)){
  p2eg[i,1] <- probes[i]
  p2eg[i,2] <- annot$EnsemblID[which(rownames(annot)==probes[i])]
}
p2eg <- as.data.frame(p2eg)

unieg <- unique(p2eg$EnsemblID)

norm2 <- matrix(nrow=length(unieg), ncol=ncol(norm))
rownames(norm2) <- unieg
colnames(norm2) <- colnames(norm)

for(i in 1:length(unieg)){
  ind <- which(p2eg$EnsemblID==unieg[i])
  for(j in 1:dim(norm)[2]){
    norm2[i,j] <- median(norm[ind,j])
  }
}

norm = norm2


# PCA ON NORMALIZED

library(swamp)

pheno = load_pheno(inputFile = "GSE92899_metadata.txt",columnID = 2)
pheno$slide = as.factor(pheno$slide)
pheno$array = as.factor(pheno$array)
pheno$dose = as.numeric(pheno$dose)
pheno$time_point = as.numeric(pheno$time_point)
# remove the removed outlier
pheno = pheno[-which(rownames(pheno) %in% "GSM2439596"),]

# set directory for images
imageDir = ""
plotPrefix = "before_batch"

plot_dendrogram(data = norm, phFactor = pheno, plotMethod = "correlation",imageDir = imageDir, plotPrefix=plotPrefix)

# MDS
plot_mds(data = norm, phFactor = pheno,mdsColor = "group",mdsLabel = "group", imageDir=imageDir,plotPrefix = plotPrefix,plotWidth = 10, plotHeight =10)

# Prince plot
prince_plot(data=norm,phFactor=pheno, npc = 10,plotHeight=10,plotWidth=10,plotPrefix=plotPrefix,imageDir,plotMargins = c(5,7))

# confounding plot
confounding_plot(phFactor = pheno, imageDir=imageDir,plotPrefix=plotPrefix)


#sva
#sva_res = sva_f(data=norm,phFactor=pheno, varI = "group", coVar = NULL)

#if the re are surrogate variables, attach them to the phenodata table!  
#svaSV <- sva_res$svaSV
#svaSVc <- sva_res$svaSVc
#pheno <- cbind(pheno, svaSV, svaSVc)
#View(pheno)

#pheno$svaD.1 = as.factor(pheno$svaD.1)
#pheno$svaD.2 = as.factor(pheno$svaD.2)
#pheno$svaD.3 = as.factor(pheno$svaD.3)

#confounding_plot(phFactor = pheno, imageDir=imageDir,plotPrefix="after_sva")
#prince_plot(data=norm,phFactor=pheno, npc = 10,plotHeight=10,plotWidth=10,plotPrefix="after_sva",imageDir,plotMargins = c(5,7))

#combat
combat_res = combat_f(data=norm, phFactor=pheno, npc = 10, varI="group", coVar = NULL, batch=c("array", "dye"))

plot_mds(data = combat_res, phFactor = pheno, mdsColor = "group",mdsLabel = "group", imageDir=imageDir,plotPrefix = "after_combat",plotWidth = 10, plotHeight =10)

prince_plot(data=combat_res,phFactor=pheno, npc = 10, plotHeight=10,plotWidth=10,plotPrefix="after_combat",imageDir,plotMargins = c(5,7))
plot_dendrogram(data = combat_res, phFactor = pheno, plotMethod = "correlation",imageDir = imageDir, plotPrefix="after_combat")


table(pheno$group)

# limma analysis ,
limma_res = run_limma_function(combat_res, pheno, pvAdjMethod="fdr",varI="group", coVar=c("array", "dye"),  comps = c("Baytube_6-control_6", "Baytube_24-control_24", "Fullerene_6-control_6",
                                                                                                                      "Fullerene_24-control_24", "Graphite_6-control_6", "Graphite_24-control_24", 
                                                                                                                      "rCNT_6-control_6", "rCNT_24-control_24", "SES_6-control_6", "SES_24-control_24", 
                                                                                                                      "tCNT_6-control_6", "tCNT_24-control_24"))
#limma_res = run_limma_function(norm, pheno, pvAdjMethod="fdr",varI="group",  comps = c("Jurkat_GO-Jurkat_control", "Jurkat_GONH2-Jurkat_control"))
limma_filtered = list()

for (i in 1:length(limma_res)){
  fillimma = limma_res[[i]][abs(limma_res[[i]]$logFC) >= 0.58 & limma_res[[i]]$adj.P.Val < 0.05, ]
  if (nrow(fillimma) > 0){
    limma_filtered[[length(limma_filtered)+1]] <- fillimma
    names(limma_filtered)[[length(limma_filtered)]] <- names(limma_res)[[i]]   
  }
}

str(limma_filtered)
str(limma_res)

geoid = "GSE92899"
write.table(combat_res, file = paste(imageDir, "Corrected_Expression_Matrix_", geoid, ".txt", sep = ""), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(norm, file = paste(imageDir, "Normalized_Expression_Matrix_", geoid, ".txt", sep = ""), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
WriteXLS(limma_res, ExcelFileName = paste(imageDir, paste(geoid, "Unfiltered_DEG.xlsx", sep = "_"), sep = ""), row.names = FALSE, col.names = TRUE)
WriteXLS(limma_filtered, ExcelFileName = paste(imageDir, paste(geoid, "Filtered_DEG.xlsx", sep = "_"), sep = ""), row.names = FALSE, col.names = TRUE)
