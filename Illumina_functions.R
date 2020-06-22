#' this function import and clean the phenodata file
#' @param inputFile path to phenodata file
#' @param columnID colnames in the phenodata file containing the sample IDs
#' @return a dataframe containing the phenodata
#' 
load_pheno = function(inputFile, columnID){
  phTable <- read.delim(inputFile)
  rownames(phTable) = phTable[,columnID]
  
  coltypes <- unlist(lapply(phTable, class))
  print("coltypes")
  print(coltypes)
  print(table(coltypes))
  
  #Check specifically for column with logical data type values
  coltypes.logical.idx <- which(coltypes=="logical")
  if(length(coltypes.logical.idx)>0){
    print("Phenotype contains logical columns!")
    print(coltypes.logical.idx)
    for(idx in as.vector(coltypes.logical.idx)){
      print("Updating logical to character!")
      print(colnames(phTable)[idx])
      phTable[,idx] <- as.character(phTable[,idx])
    }
    coltypes <- unlist(lapply(phTable, class))
  }
  
  #Identify the indices for character columns an non-character columns
  coltypes.charOnly.idx <- which(coltypes=="character")
  coltypes.nonChar.idx <- which(!coltypes=="character")
  coltypes.charOnly.len <- length(coltypes.charOnly.idx)
  coltypes.nonChar.len <- length(coltypes.nonChar.idx)
  
  #Check mixed type columns and remove them
  remInfo <- 0
  remStr <- ""
  if(coltypes.charOnly.len>0){
    #Subset and get character type columns
    phTable.charOnly <- phTable[, coltypes.charOnly.idx, drop=F]
    print("dim(phTable.charOnly)")
    print(dim(phTable.charOnly))
    
    #Evaluate character column to contain numeric data type values
    numCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.numeric(col))))){1}else{0}}))
    #Evaluate character column to contain integer data type values
    intCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.integer(col))))){1}else{0}}))
    #Evaluate character column to contain double data type values
    doubleCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.double(col))))){1}else{0}}))
    allCheck <- numCheck+intCheck+doubleCheck
    
    #Check if character columns contain other data type values
    if(all(allCheck==0)){
      checkFailed <- names(allCheck[allCheck==0])
      remInfo <- 1
      remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
      if(coltypes.nonChar.len>0){
        phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
        phTable.comb <- phTable.nonChar
      }else{
        remStr <- paste0(remStr, "No column survived filtering!!! Please define phenotype data columns with singular data type.")
        
        taskManager$addError('error2', remStr)
        taskManager$end()
        
        return(NULL)
      }
    }else{
      if(any(allCheck==0)){
        checkFailed <- names(allCheck[allCheck==0])
        phTable.charOnly <- phTable.charOnly[,-which(colnames(phTable.charOnly) %in% checkFailed)]
        remInfo <- 1
        remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
      }
      
      print("str(phTable) -- before:")
      print(str(phTable))
      
      #Fix spaces and '-' special character
      phTable.charOnly <- as.data.frame(apply(phTable.charOnly, 2, function(x){res<-trimws(x); res<-gsub(" +", " ", res, perl=T); res<-gsub("[ -]", "_", res); return(res)}), stringsAsFactors=FALSE)
      if(coltypes.nonChar.len>0){
        phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
        print("dim(phTable.nonChar)")
        print(dim(phTable.nonChar))
        phTable.comb <- data.frame(phTable.charOnly, phTable.nonChar, stringsAsFactors=FALSE)
      }else{
        phTable.comb <- phTable.charOnly
      }
    }
    
    #Identify original order of the columns in the user provided file
    colOrgIdx <- sapply(colnames(phTable.comb), function(x){which(colnames(phTable) %in% x)})
    #Reorder the cleaned table to match the original column order
    phTable <- phTable.comb[,names(colOrgIdx[order(colOrgIdx)]), drop=F]
  }
  
  print("str(phTable) -- check:")
  print(str(phTable))
  
  #Identify columns with single level data
  nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
  nrlevels.singular <- which(nrlevels==1)
  
  #Remove columns with single level data
  if(length(nrlevels.singular)>0){
    remInfo <- 1
    remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
    if(length(nrlevels.singular)==ncol(phTable)){
      remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
      
      taskManager$addError('error2', remStr)
      taskManager$end()
      
      return(NULL)
    }
    col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
    phTable <- phTable[,-col2rem, drop=F]
  }
  
  
  print("str(phTable) -- after:")
  print(str(phTable))
  
  toRem = c()
  for(i in 1:ncol(phTable)){
    if(sum(is.na(phTable[,i])) == nrow(phTable)){
      toRem = c(toRem, i)
    }
  }
  
  if(length(toRem)>1){
    phTable = phTable[,-toRem]
  }
  
  # #Create phTable variable class table
  # phClassTable <- phTable %>% dplyr::summarise_all(class) %>% tidyr::gather(Variable, Class) %>% dplyr::mutate(Type=ifelse(Class=="character", "factor", "vector"))
  return(phTable)
}


#' this function takes in input the filtered matrix with probe as rownames and the mapping database 
#' @param selDataMatrix filtered matrix with probe as rownames
#' @param MappingDB mapping database
#' @return a a matrix with gene ensemble rownames that are the result of a mean probe summarization 
#' 

probe_summarization = function(selDataMatrix, MappingDB = illuminaHumanv3ENSEMBL){
  probeList  = rownames(selDataMatrix)  
  
  probes2ENSEMBL<-data.frame(Gene=unlist(mget(x = probeList,envir = MappingDB)))
  probes2ENSEMBL = data.frame("Probes" = rownames(probes2ENSEMBL), "Genes" = probes2ENSEMBL$Gene)
  toRem = which(is.na(probes2ENSEMBL[,"Genes"]))
  
  # remove probes that do not map to ensemble gene id
  if(length(toRem)>0) probes2ENSEMBL=probes2ENSEMBL[-toRem,]
  toRem2 = which((probes2ENSEMBL$Probes %in% rownames(selDataMatrix))==FALSE)
  if(length(toRem2)>0) probes2ENSEMBL=probes2ENSEMBL[-toRem2,]
  probes2ENSEMBL = data.frame(Probes = as.character(probes2ENSEMBL$Probes), Genes = as.character(probes2ENSEMBL$Genes))
  
  uniEG = unique(probes2ENSEMBL$Genes)  
  norm <- matrix(nrow=length(uniEG), ncol=ncol(selDataMatrix))
  colnames(norm) <- colnames(selDataMatrix)
  rownames(norm) <- uniEG
  
  pb = txtProgressBar(min = 1,max = length(uniEG),style = 3)
  for(i in 1:length(uniEG)){
    ind <- which(probes2ENSEMBL$Genes==uniEG[i])
    probi = probes2ENSEMBL[ind,1]
    if(length(probi)>1) {
      norm[i,] = colMeans(selDataMatrix[probi,])
    }else{
      norm[i,] = selDataMatrix[probi,]
      
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  return(norm)
}

#' this function take in input the lumi object, the expression matrix and the filtering threshold on the probe detection value
#' it removes from the expression matrix probes with a threshold lower than the one specified by the user
#' @param lumiObj lumi object 
#' @param dataMatrix expression matrix with probes as rownames
#' @param th threshold, default value = 0
#' @return a a matrix with gene ensemble rownames that are the result of a mean probe summarization 
#'
filtering = function(lumiObj,dataMatrix,  th = 0){
  presentCount <- detectionCall(lumiObj)
  selDataMatrix <- dataMatrix[presentCount > 0,] 
  return(selDataMatrix)
}


#' this function produces density, MA, sample relation, outliers, mds and boxplot starting from the lumi object.
#' this plot can be used for quality controls
#' it takes in input the lumi object, a flag that can be before or after normalization to attach to the file names, the folder where output will be stored
#' @param lumiObj lumi object 
#' @param flag string that will be attached to the name of the file
#' @param outputFolder path of a folder where to store the pdf of the plots
#'
make_plots = function(lumiObj, flag = "before", outputFolder){
  
  pdf(paste(outputFolder,"density_",flag,".pdf",sep=""))
  plot(lumiObj, what='density')
  dev.off()
  
  #pdf(paste(outputFolder,"MA_",flag,".pdf",sep=""))
  #MAplot(lumiObj, smoothScatter=TRUE)
  #dev.off()
  
  pdf(paste(outputFolder,"sampleRelation_",flag,".pdf",sep=""))
  plot(lumiObj, what='sampleRelation')
  dev.off()
  
  pdf(paste(outputFolder,"outlier_",flag,".pdf",sep=""))
  plot(lumiObj, what='outlier')
  dev.off()
  
  pdf(paste(outputFolder,"mds_",flag,".pdf",sep=""))
  plotSampleRelation(lumiObj, method='mds')
  dev.off()
  
  pdf(paste(outputFolder,"boxplot_",flag,".pdf",sep=""))
  plot(lumiObj, what="boxplot")
  dev.off()
}

#' this function produce a dendrogram of the sample and annotate it with color bars coming from the phenodata
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param plotMethod: the distance method for the clustering. default="correlation". hcluster from the package amap is used and method must be one of "euclidean", "maximum", "manhattan", "canberra" "binary" "pearson", "correlation", "spearman" or "kendall".
#' @param imageDir path to the directory where to store the file
#' @param plotPrefix string attached as prefix to the file name
#' @param plotWidth width of the plot default = 10
#' @param plotHeight height of the plot default = 5
#'

plot_dendrogram = function(data, phFactor, plotMethod = "correlation",imageDir, plotPrefix, plotWidth = 10, plotHeight =5){
  
  hcPlotFile <- file.path(imageDir, paste0(plotPrefix, "_hcPlot.pdf"))
  pdf(file=hcPlotFile, height=plotHeight, width=plotWidth)
  swamp::hca.plot(g=data, o=phFactor, method=plotMethod)
  dev.off()
}


#' this function produce a mds of the sample and uses phenodata to annotate it
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param topNum is an integer number specifying the number of components
#' @param mdsColor colnames of the phenodata used to color the sample labels on the plot
#' @param mdsLabel colnames of the phenodata used to plot labels on the plot
#' @param plotTitle title of the mds plot
#' @param imageDir path to the directory where to store the file
#' @param plotPrefix string attached as prefix to the file name
#' @param plotWidth width of the plot default = 10
#' @param plotHeight height of the plot default = 10
#'
plot_mds = function(data, phFactor,topNum = 10, mdsColor, mdsLabel, plotTitle = "",imageDir,plotPrefix,plotWidth = 10, plotHeight =10){
  if((mdsColor %in% colnames(phFactor)) == FALSE){
    error("mdsColor has to be one of the colnames of the phenodata ")
  }
  colVec <- get_color_palette(iVec=phFactor[,mdsColor], asFactor=TRUE)
  
  mdsPlotFile <- file.path(imageDir, paste0(plotPrefix, "_mdsPlot.pdf"))
  pdf(file=mdsPlotFile, height=plotHeight, width=plotWidth)
  limma:::plotMDS(data, top=topNum, labels=as.character(phFactor[,mdsLabel]), col=colVec, gene.selection="common", main=plotTitle)
  dev.off()
}


#' this function produce a princeplot of the expression matrix
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param npc is an integer number specifying the number of components
#' @param imageDir path to the directory where to store the file
#' @param plotPrefix string attached as prefix to the file name
#' @param plotWidth width of the plot default = 10
#' @param plotHeight height of the plot default = 10
#' @param plotMargins vector specifying the margins of the plot. Default = c(5,7). see parameter margins in function prince.plot{swamp}
#' @param note boolean to plot or not pvalues in the heatmap. default = TRUE. see parameter note in function prince.plot{swamp}
#' @param notecex font cex . default = 1. see parameter notecex in function prince.plot{swamp}
#'

prince_plot = function(data,phFactor, npc = 10,plotHeight=10,plotWidth=10,plotPrefix,imageDir,plotMargins = c(5,7), note = TRUE, notecex = 1){
  if(npc>ncol(data)){
    npc <- ncol(data)
  }
  
  pr <- swamp::prince(data, phFactor, top=npc)
  
  # generate the prince plot
  princePlotFile <- file.path(imageDir, paste0(plotPrefix, "_princePlot.pdf"))
  pdf(file=princePlotFile, height=plotHeight, width=plotWidth)
  swamp::prince.plot(prince=pr, margins=plotMargins, note=note, notecex=notecex)
  dev.off()
}

#' this function produce a princeplot of the expression matrix
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param npc is an integer number specifying the number of components
#' @param imageDir path to the directory where to store the file
#' @param plotPrefix string attached as prefix to the file name
#' @param plotWidth width of the plot default = 10
#' @param plotHeight height of the plot default = 10
#' @param plotMargins vector specifying the margins of the plot. Default = c(5,7). see parameter margins in function prince.plot{swamp}
#' @param note boolean to plot or not pvalues in the heatmap. default = TRUE. see parameter note in function prince.plot{swamp}
#' @param notecex font cex . default = 1. see parameter notecex in function prince.plot{swamp}
#'
confounding_plot = function(phFactor, imageDir,plotPrefix,plotHeight=10,plotWidth=10,plotMargins = c(5,7), note = TRUE, notecex = 1){
  confPlotFile <- paste(imageDir, plotPrefix,"_confoundingPlot.pdf",sep = "")
  pdf(file=confPlotFile, height=plotHeight, width=plotWidth)
  swamp::confounding(phFactor, margins=plotMargins, note=note, notecex=notecex)
  dev.off()
}

#' this function perform sva analysis to identify hidden variables
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param varI colnames of the phenodata specifying the variable of interests
#' @param coVar vector of strings specifying the covariates
#' @return a list containing the surrogate variables in discrete and continuous format
#'
sva_f = function(data,phFactor,varI, coVar = NULL){
  
  if((varI %in% colnames(phFactor))==FALSE){
    print("unknown variable of interest")
    return(NULL)
  }
  
  if(is.null(coVar)){
    coVar <- NULL
  }else{
    if(sum(coVar %in% colnames(phFactor))>1){
      print("one or more covariates not preset in phenodata table")
      return(NULL)
    }
    coVar <- as.list(coVar)
  }
  
  batchCorVar <- list(var.int=varI, covariates=coVar, batches=NULL)
  print(str(batchCorVar))
  batches.sva <- get.sva.batch.effects(comb.data=data, pd=phFactor, vars=batchCorVar, cmd.ui=FALSE)
  assoc.cutoff <- 0.05
  sv.filt.logic <- apply(batches.sva$pd, 2, function(x) x[batchCorVar$var.int]<assoc.cutoff)
  sv.filt.names <- names(sv.filt.logic[sv.filt.logic==F])
  
  print(sv.filt.names)
  if(length(sv.filt.names)>0){
    if(length(sv.filt.names)>1){
      #Check for confounded sva variables
      sva.assoc.mat <- assoc.var.int(batches.sva$sv[sv.filt.names], batches.sva$sv[sv.filt.names])
      colnames(sva.assoc.mat) <- rownames(sva.assoc.mat)
      sva.assoc.mat[lower.tri(sva.assoc.mat, diag=T)] <- NA
      sva.assoc.DF <- as.data.frame(as.table(sva.assoc.mat), stringsAsFactors=FALSE)
      sva.assoc.DF <- sva.assoc.DF[-which(is.na(sva.assoc.DF$Freq)),]
      print(dim(sva.assoc.DF))
      print(head(sva.assoc.DF))
      rowSel <- which(sva.assoc.DF$Freq<0.01)
      if(length(rowSel)>0){
        sva.assoc.DF.conf <- sva.assoc.DF[rowSel,]
        sva.rm.names <- unique(sva.assoc.DF.conf[,2])
        sv.filt.names <- sv.filt.names[-which(sv.filt.names %in% sva.rm.names)]
      }
    }
    #Add selected sva sv and perform combat
    sva.sv.filt <- batches.sva$sv[,sv.filt.names]
    sva.svc.filt <- batches.sva$svc[,sv.filt.names]
    print(class(sva.sv.filt))
    if(!is.data.frame(sva.sv.filt) && !is.matrix(sva.sv.filt)){
      sva.sv.filt <- data.frame(sva.sv.filt)
      sva.svc.filt <- data.frame(sva.svc.filt)
    }
    colnames(sva.sv.filt) <- paste("svaD",c(1:ncol(sva.sv.filt)),sep=".")
    colnames(sva.svc.filt) <- paste("svaC",c(1:ncol(sva.svc.filt)),sep=".")
    print(class(sva.sv.filt))
    print(dim(sva.sv.filt))
    print(head(sva.sv.filt))
    
    svaSV <- sva.sv.filt
    svaSVc <- sva.svc.filt
    svaStep <- TRUE
    return(list(svaSV=svaSV, svaSVc=svaSVc))
  }else{
    svaStep <- FALSE
    
    errStr <- paste0("All surrogate variables are confounded with variable of interest! Proceed to ComBat with known variables OR rerun SVA with different model!")
    print(errStr)
    
    return(NULL)
  }
}

#' this function perform combat analysis to remove known batches
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phFactor phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param npc is an integer number specifying the number of components
#' @param varI colnames of the phenodata specifying the variable of interests
#' @param coVar vector of strings specifying the covariates
#' @param batch colnamesf of the phenodata specifying the batch variable
#' @return a list containing the surrogate variables in discrete and continuous format
#'
combat_f = function(data,phFactor, npc = 10,varI, coVar = NULL, batch){
  if((varI %in% colnames(phFactor))==FALSE){
    print("unknown variable of interest")
    return(NULL)
  }
  
  if((batch %in% colnames(phFactor))==FALSE){
    print("unknown batch variable")
    return(NULL)
  }
  
  if(npc>ncol(data)){
    npc <- ncol(data)
  }
  
  if(is.null(coVar)){
    coVar <- NULL
  }else{
    if(sum(coVar %in% colnames(phFactor))>1){
      print("one or more covariates not preset in phenodata table")
      return(NULL)
    }
    coVar <- as.list(coVar)
  }
  
  batchCorVar <- list(var.int=varI, covariates=coVar, batches=batch)
  
  
  comb.data <- remove.batch.effects(comb.data =data,pd = phFactor,num.pc = npc, vars = batchCorVar, method="Combat", plot=FALSE, verbose=TRUE)
  return(comb.data)
}

#' this function run limma analysis to identify differentially expressed genes
#' @param data expression matrix with genes on the rows and samples on the columns
#' @param phTable phenodata dataframe with samples on the rows and variables on the columns. N.B. rownames of phenodata has to be the same of colnames of expression data
#' @param pvAdhMethod character string specifying p-value adjustment method. Possible values are "none", "BH", "fdr" (equivalent to "BH"), "BY" and "holm". See p.adjust for details.
#' @param varI colnames of the phenodata specifying the variable of interests
#' @param coVar vector of strings specifying the covariates
#' @param comps vector containint the contrasts
#' @return a list with toptable result for each contrast
#'
run_limma_function = function(data, phTable, pvAdjMethod,varI, coVar=NULL,  comps = c("crystalline_silica_800-control_0")){
  
  if(is.null(coVar)){
    coVar <- NULL
  }else{
    coVar <- as.list(coVar)
  }
  
  if((varI %in% colnames(phTable))=="FALSE"){
    error("variable of interest not in phenodata")
    return(NULL)
  }
  
  if(sum((coVar %in% colnames(phTable))==FALSE)>1){
    error("covariate/s not in phenodata")
    return(NULL)
  }
  
  
  batchCorVar <- list(var.int=varI, covariates=coVar) ## No batch for limma model
  print(str(batchCorVar))
  
  des <- build.model.matrix(pd = phTable, intercept=-1, var.int = batchCorVar$var.int, covariates = batchCorVar$covariates, verbose=T)
  ne <- limma::nonEstimable(des)
  if(!is.null(ne)){
    errStr <- paste0("Coefficients not estimable : ",paste(ne,collapse=", "),"\n\nPlease check your limma model definition!!! Possibly contains confounders.")
    print(errStr)
    return(NULL)
  }
  
  deg.list <- diff.gene.expr(data, des, contrasts=comps, pvalue=1, fcvalue=0, p.adjust.method="fdr", annot=NULL, plot=FALSE, verbose=TRUE)
  return(deg.list)
}