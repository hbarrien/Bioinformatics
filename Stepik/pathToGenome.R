# Usage:
# path_in <- read.csv("C:/Users/heba/Downloads/dataset_198_3.txt", sep="", na.strings="", stringsAsFactors=FALSE)
# path <- path_in$kmer
# my_genome <- pathToGenome(path)

pathToGenome <- function(path) {
  
  assembledGenome <- ""

  if (is.null(path) || (length(path) == 0)) return(assembledGenome)
  if (length(path) == 1) return(path[1])
  
  numKmers   <- length(path)
  numSymbols <- nchar(path[1])
  
  # At this point,  (|path| > 1)
  # Loop invariant: (|path[idx]| == numSymbols)
  for (idx in 1:(numKmers-1)) {
    
    # Check invariant
    if (nchar(path[idx]) != numSymbols)
      stop(paste0("Error: kmer ", path[idx], " has different length."))
           
    if (!(substr(path[idx], 2, numSymbols) == substr(path[idx+1], 1, (numSymbols-1))))
      stop(paste0("Error: no match between kmer: ", idx, " and kmer: ", (idx+1)))
    
    assembledGenome <- paste0(assembledGenome, substr(path[idx], 1, 1))
    if (idx == (numKmers-1)) assembledGenome <- paste0(assembledGenome, path[idx+1])
    
  }  # END for
  
  return(assembledGenome)
  
}  # END pathToGenome