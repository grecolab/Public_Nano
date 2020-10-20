### Get and filter counts.

get_raw_counts <- function(bam.files, annotation, paired.end, ncores=2) {
  
  if(is.null(bam.files)){stop("Error: please provide one or more .bam files!")}
  if(is.null(annotation)){stop("Error: please provide an annotation file in .gtf or .gff format!")}
  if(!class(paired.end)=="logical"){stop("Error: indicate whether the reads are paired end")}
  if(!class(ncores)=="numeric"){stop("Error: please input a numeric value!")}
  
  raw.counts1 <- Rsubread::featureCounts(bam.files, annot.ext=annotation, isPairedEnd=paired.end, isGTFAnnotationFile = TRUE, nthreads = ncores)
  
  total.counts.matrix <- raw.counts1[[1]]
  
  genelengths <- raw.counts1$annotation$Length
  names(genelengths) <- raw.counts1$annotation$GeneID
  
  
  return(list(total.counts.matrix, genelengths))
}




filter_low_counts <- function(counts.matrix, conditions, method = "cpm", normalized=FALSE, depth=NULL, cpm=1, p.adj = "fdr"){
  
  if(is.null(counts.matrix)){stop("Error: please provide a numeric count matrix!")}
  if(is.null(conditions)){stop("Error: please provide a factor or a vector indicating the conditions!")}
  if(!method %in% c("cpm", "wilcoxon", "proportion")) {stop("Error: Please type in one of valid methods!")}
  
  
  if (method=="cpm"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 1, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
    
  }else if(method=="wilcoxon"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 2, cv.cutoff = 100, p.adj = p.adj)
    
  }else if(method=="proportion"){
    
    if(is.null(depth)){stop("Error: indicate a numeric vector indicating per sample library depths")}
    if(!class(depth)=="numeric"){stop("Error: please provide the depth argument with a numeric vector!")}
    ### Compute librarary depth
    
    
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, depth = depth, method = 3, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
    
  }
  return(filtered.counts)
}
