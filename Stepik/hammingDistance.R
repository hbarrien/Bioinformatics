# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
EMPTY_STRING <- ""


# *********************************
# *********** VARIABLES ***********
# *********************************


# *********************************
# *********** FUNCTIONS ***********
# *********************************
hammingDistance <- function(p, q) {
  
  if (is.null(p) || (length(p) == 0) ||
      is.null(q) || (length(q) == 0))
    return(NULL)
  
  if (length(p) != length(q))
    return(NULL)
  
  sP <- strsplit(p, EMPTY_STRING)[[1]]
  sQ <- strsplit(q, EMPTY_STRING)[[1]]
  
  hD <- (sP == sQ)
  
  return(sum(hD == FALSE))
  
}  # END hammingDistance

# Usage:
# fName <- paste0(getwd(), "/aodna.txt"); dna <- readChar(fName, file.size(fName)); fName <- paste0(getwd(), "/aodna_out.txt"); f <- file(fName); writeChar(as.character(aO), f); close(f)
approximateOccurrences <- function(pattern, dna, d, zeroBasedIdx = TRUE) {
  
  aO <- c()
  lenPattern  <- nchar(pattern)
  lenDnaParse <- (nchar(dna)-lenPattern+1)
  
  for (i in 1:lenDnaParse) {
    
    subPattern <- substr(dna, i, (i+lenPattern-1))
    hD <- hammingDistance(pattern, subPattern)
    if (is.null(hD)) return(NULL)
    
    if (hD <= d) aO <- c(aO, i)
    
  }  # END for
  
  if (zeroBasedIdx)
    aO <- (aO - 1)
  
  return(aO)
  
}  # END approximateOccurrences

countMismatches <- function(pattern, dna, d, zeroBasedIdx = TRUE) {
  
  return(length(approximateOccurrences(pattern, dna, d, zeroBasedIdx)))
  
}  # END countMismatches