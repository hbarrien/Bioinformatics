# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
A <- 65
C <- 67
G <- 71
T <- 84

SKEW_TABLE <- integer(T)
SKEW_TABLE[C] <- -1
SKEW_TABLE[G] <- 1


# *********************************
# *********** VARIABLES ***********
# *********************************
SK_VERBOSE <- 1

# *********************************
# *********** FUNCTIONS ***********
# *********************************
skew <- function(dna) {
  
  splitDna  <- strsplit(dna, "")[[1]]
  skewCount <- 0
  
  ret <- c(skewCount)
  ret <- c(ret, as.integer(lapply(splitDna, function(x) {skewCount <<- (skewCount + SKEW_TABLE[utf8ToInt(x)]); return(skewCount)})))
  
  return(ret)
  
}  # END skew

# Usage:
# fName <- paste0(getwd(), "/_skew.txt"); dna <- readChar(fName, file.size(fName)); getMinSkewPos(dna)
getMinSkewPos <- function(dna, zeroBasedIdx = TRUE) {
  
  s <- skew(dna)
  
  if (zeroBasedIdx) minSkewPos <- (which(s == min(s)) - 1)
  else minSkewPos <- which(s == min(s))
    
  return(minSkewPos)
  
}  # END getMinSkewPos
