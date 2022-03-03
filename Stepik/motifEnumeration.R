source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/neighbors.R"))


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************


# *********************************
# *********** VARIABLES ***********
# *********************************


# *********************************
# *********** FUNCTIONS ***********
# *********************************
enumerateMotifs <- function(dna, k, d) {
  
  numDna   <- length(dna)
  lenDna   <- nchar(dna[1])
  numKmers <- (lenDna-k+1)
  allNeighbors <- list()
  
  for (i in 1:numDna) {
    
    nxtDna <- dna[i]
    dnaNeighbors <- c()
    
    for (j in 1:numKmers) {
      
      kmer <- substr(nxtDna, j, (j+k-1))
      dnaNeighbors <- c(dnaNeighbors, neighbors(kmer, d))
      
    }  # END for
    
    dnaNeighbors <- dnaNeighbors[!duplicated(dnaNeighbors)]
    allNeighbors[[i]] <- dnaNeighbors
    
  }  # END for
  
  res <- Reduce(intersect, allNeighbors)
  
  return(res)
  
}  # END enumerateMotifs
