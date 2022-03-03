# Usage:
# k <- 9
# x <- kUniversalCircularString(k)


source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/stringReconstruction.R"))


# *********************************
# *********** LIBRARIES ***********
# *********************************
library(R.utils)


# *********************************
# *********** FUNCTIONS ***********
# *********************************
binaryStrings <- function(k) {
  
  return(intToBin(c(0:(2^k-1))))
  
}  # END binaryStrings


kUniversalCircularString <- function(k) {
  
  # Get the string resulting from the Euler circuit
  patterns <- binaryStrings(k)
  ePath    <- stringReconstruction(patterns)
  
  # Convert the string resulting from the Euler circuit to a circular string
  lastKmer   <- substr(ePath, nchar(ePath)-k+1, nchar(ePath))
  lastSymbol <- substr(lastKmer, nchar(lastKmer), nchar(lastKmer))
  circularString <- paste0(lastSymbol, substr(ePath, 1, nchar(ePath)-k))
  
  return(circularString)
  
}  # END kUniversalCircularString
