# Usage
# fName <- "C:/Users/heba/Downloads/dataset_158_9.txt"
# median <- read.csv(fName, sep="", na.strings="", stringsAsFactors=FALSE)
# dnaList <- median$x
# k <- ?
# t1 <- Sys.time(); x <- medianString(dnaList, k); t2 <- Sys.time(); (t2-t1)
# x


source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/neighbors.R"))

# *********************************
# *********** LIBRARIES ***********
# *********************************
library(data.table)


# *********************************
# *********** CONSTANTS ***********
# *********************************


# *********************************
# *********** VARIABLES ***********
# *********************************
MS_VERBOSE <- TRUE


# *********************************
# *********** FUNCTIONS ***********
# *********************************
motifs <- function(dnaList, k) {
  
  # get all k-length patterns
  kmers <- attributes(frequencyArray(k)[[PATTERN_IDX]])[[1]]
  
  # get the max d "mismatch" value, with respect to k
  maxd <- (k-1)
  
  out <- c()
  while (length(kmers) > 0) {
    
    pattern <- kmers[1]
    
    for (d in 1:maxd) {
      
      if (d == 1) pNeighbors <- neighbors(pattern, d)
      else pNeighbors <- neighborsExactMatch(pattern, d)
      
      for (n in pNeighbors) {
        
        # invariant: (|idx| > 0) --> (isUnique(idx) == TRUE)
        idx <- grep(n, dnaList, fixed = TRUE)
        
        if (length(idx) > 0) {
          
          out <- c(out, n)
          dnaList <- dnaList[-idx]
          
          if (length(dnaList) == 0)
            return(out)
          
        }  # END if
        
      }  # END for
      
      kmers <- setdiff(kmers, pNeighbors)
      
    }  # END for
    
  }  # END while
  
  return(out)
  
}  # END motifs

d <- function(pattern, dnaList) {
  
  maxd <- (nchar(pattern) - 1)
  res  <- 0
  
  for (d in 1:maxd) {
    
    if (d == 1) pNeighbors <- neighbors(pattern, d)
    else pNeighbors <- neighborsExactMatch(pattern, d)
    
    for (n in pNeighbors) {
      
      # invariant: (|idx| > 0) --> (isUnique(idx) == TRUE)
      idx <- grep(n, dnaList, fixed = TRUE)
      
      if (length(idx) > 0) {
        
        if (n != pattern)
          res <- (res + (d * length(idx)))
        
        dnaList <- dnaList[-idx]
        
        if (length(dnaList) == 0)
          return(res)
        
      }  # END if
      
    }  # END for
  
  }  # END for
  
  return(res)
  
}  # END d

medianString <- function(dnaList, k) {
  
  # get all k-length patterns
  kmers <- attributes(frequencyArray(k)[[PATTERN_IDX]])[[1]]

  out <- data.table(pattern = character(), d = numeric())
  for (pattern in kmers) {
    
    if (MS_VERBOSE) print(pattern)
    out <- rbind(out, list(pattern, d(pattern, dnaList)))
    
  }  # END for
  
  return(out[out$d == min(out$d),]$pattern)
  
}  # END medianString
