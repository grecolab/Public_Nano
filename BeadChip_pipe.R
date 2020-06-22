library(lumi)
library(illuminaHumanv3.db)
library(illuminaHumanv4.db)
library(illuminaRatv1.db)
library(illuminaMousev2.db)
source("Illumina_functions.R")
source("pipeline_functions.R")
library(swamp)
library(WriteXLS)

# Data sets to preprocess
todo = c("GSE99929","GSE99929_THP1", "GSE99929_Jurkat", "GSE20692", "GSE30180", "GSE117056", "GSE30178", "GSE100500") 
x = 1

geoid = todo[x]

# Raw data path
inputFile = paste("/home/MOKA/laura/NanoSolveIT_data_preprocessing/Illumina_raw/", paste(geoid, "non-normalized.txt", sep = "_"), sep = "") #others non-normalized

# Create an output folder
dir.create(paste("/home/MOKA/laura/NanoSolveIT_data_preprocessing/Illumina_results/", paste(geoid, "results", sep = "_"), sep = ""))
outputFolder = paste("/home/MOKA/laura/NanoSolveIT_data_preprocessing/Illumina_results/", paste(geoid, "results/", sep = "_"), sep = "")

x.lumi <- lumiR(inputFile,lib.mapping = "lumiRatIDMapping",convertNuID=FALSE)  #change mouse to human or rat according to species

#Manifest file
#mani = "/home/MOKA/laura/NanoSolveIT_data_preprocessing/Illumina_raw/HumanHT-12_V4_0_R2_15002873_B.txt"

#Read idat files into an expression set
#x.lumi <- lumiR.idat(files = filenames, path = "/home/MOKA/laura/NanoSolveIT_data_preprocessing/Illumina_raw/GSE117056_RAW/", probeID = "ProbeID", QC = FALSE, manifestfile = mani, memory = "-Xmx4g")

# Remove outlier based on position in x.lumi@phenoData@data
#which(x.lumi@phenoData@data == "24H_A549_Si(0)_14")
#x.lumi <- x.lumi[,-c(3)]

#quality check plots before preprocessing
make_plots(x.lumi, flag = "before", outputFolder)

# background correction
datB <- lumiB(x.lumi)   
#log2 transformation
datT <- lumiT(datB, method="log2")
#quantile normalization
datN <- lumiN(datT, method="quantile")

#quality check plots after preprocessing
make_plots(datN, flag = "after", outputFolder)

#from probes to ensemblID
dataMatrix <- exprs(datN)
# filter probes
selDataMatrix = filtering(x.lumi, dataMatrix, th = 0)

# Change mapping based on platform
norm = probe_summarization(selDataMatrix, MappingDB = illuminaRatv1ENSEMBL)


pheno = load_pheno(inputFile = paste("/home/MOKA/laura/NanoSolveIT_data_preprocessing/", geoid, "_pheno.txt", sep = ""), columnID = "GSM")

pheno = pheno[colnames(norm),]
#pheno = pheno[,-7]

imageDir = outputFolder
plotPrefix = "before_batch"

pheno$dose = as.numeric(pheno$dose)
pheno$GSM = as.factor(pheno$GSM)
pheno$treatment = as.factor(pheno$treatment)
pheno$group = as.factor(pheno$group)
pheno$slide = as.factor(pheno$slide)
#pheno$array = as.factor(pheno$array)

# hierarlchical clustering with color bars
plot_dendrogram(data = norm, phFactor = pheno, plotMethod = "correlation",imageDir = imageDir, plotPrefix=plotPrefix)

# MDS
plot_mds(data = norm, phFactor = pheno,mdsColor = "group",mdsLabel = "group", imageDir=imageDir,plotPrefix = plotPrefix,plotWidth = 10, plotHeight =10)

# Prince plot
prince_plot(data=norm,phFactor=pheno, npc = 10,plotHeight=10,plotWidth=10,plotPrefix=plotPrefix,imageDir,plotMargins = c(5,7))

# confounding plot
confounding_plot(phFactor = pheno, imageDir=imageDir,plotPrefix=plotPrefix)

#sva
sva_res = sva_f(data=norm,phFactor=pheno, varI = "group", coVar = NULL)

#if the re are surrogate variables, attach them to the phenodata table!  
svaSV <- sva_res$svaSV
svaSVc <- sva_res$svaSVc
pheno <- cbind(pheno, svaSV, svaSVc)
View(pheno)

pheno$svaD.1 = as.factor(pheno$svaD.1)
#pheno$svaD.2 = as.factor(pheno$svaD.2)
#pheno$svaD.3 = as.factor(pheno$svaD.3)

confounding_plot(phFactor = pheno, imageDir=imageDir,plotPrefix="after_sva")
prince_plot(data=norm,phFactor=pheno, npc = 10,plotHeight=10,plotWidth=10,plotPrefix="after_sva",imageDir,plotMargins = c(5,7))

#combat
combat_res = combat_f(data=norm,phFactor=pheno, npc = 10,varI="group", coVar = NULL, batch=c("slide"))

plot_mds(data = combat_res, phFactor = pheno, mdsColor = "group",mdsLabel = "group", imageDir=imageDir,plotPrefix = "after_combat",plotWidth = 10, plotHeight =10)

prince_plot(data=combat_res,phFactor=pheno, npc = 10, plotHeight=10,plotWidth=10,plotPrefix="after_combat",imageDir,plotMargins = c(5,7))

table(pheno$group)

# limma analysis, corrected batches as covariates
limma_res = run_limma_function(combat_res, pheno, pvAdjMethod="fdr",varI="group", coVar=c("slide"),  comps = c("Ag10_Citrate-Control_Citrate", "Ag75_Citrate-Control_Citrate", 
                                                                                                               "Ag10_PVP-Control_HBSS", "Ag75_PVP-Control_HBSS"))
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

write.table(combat_res, file = paste(outputFolder, "Corrected_Expression_Matrix_", geoid, ".txt", sep = ""), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(norm, file = paste(outputFolder, "Normalized_Expression_Matrix_", geoid, ".txt", sep = ""), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
WriteXLS(limma_res, ExcelFileName = paste(outputFolder, paste(geoid, "Unfiltered_DEG.xlsx", sep = "_"), sep = ""), row.names = FALSE, col.names = TRUE)
WriteXLS(limma_filtered, ExcelFileName = paste(outputFolder, paste(geoid, "Filtered_DEG.xlsx", sep = "_"), sep = ""), row.names = FALSE, col.names = TRUE)

x = x+1
