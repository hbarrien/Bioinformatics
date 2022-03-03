# Usage:
# gappedPatterns <- read.csv("C:/Users/heba/Downloads/dataset_6206_4.txt", sep="", na.strings="", stringsAsFactors=FALSE)
# gappedPatterns <- gappedPatterns$rP
# k <- 50
# d <- 200
# x <- stringSpelledByGappedPatterns(gappedPatterns, k, d)


source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/pathToGenome.R"))


# *********************************
# *********** CONSTANTS ***********
# *********************************
READ_PAIR_SEPARATOR_IN  <- "\\|"
READ_PAIR_SEPARATOR_OUT <- "|"

FIRST_PATTERNS  <- "firstPatterns"
SECOND_PATTERNS <- "secondPatterns"

SPACE <- " "
EMPTY <- ""


# *********************************
# *********** FUNCTIONS ***********
# *********************************
separateReadPair <- function(readPair) {
  
  readPair  <- gsub(SPACE, EMPTY, readPair)
  readPairs <- strsplit(readPair, READ_PAIR_SEPARATOR_IN)[[1]]
  
  return(readPairs)
  
}  # END separateReadPair


separateReadPairs <- function(readPairs) {
  
  out <- list(firstPatterns=NULL, secondPatterns=NULL)
  
  firstPatterns  <- c()
  secondPatterns <- c()
  
  for (readPair in readPairs) {
    
    sepReadPairs <- separateReadPair(readPair)
    
    firstPatterns  <- c(firstPatterns, sepReadPairs[1])
    secondPatterns <- c(secondPatterns, sepReadPairs[2])
    
  }  # END for
  
  out[[FIRST_PATTERNS]]  <- firstPatterns
  out[[SECOND_PATTERNS]] <- secondPatterns
  
  return(out)
  
}  # END separateReadPairs

 
prefixPaired <- function(readPair) {
  
  readPairs <- separateReadPair(readPair)
  lenRead   <- nchar(readPairs[1])
  pfix      <- paste0(substr(readPairs[1], 1, (lenRead-1)), READ_PAIR_SEPARATOR_OUT, substr(readPairs[2], 1, (lenRead-1)))
  
  return(pfix)
  
}  # END prefixPaired


suffixPaired <- function(readPair) {
  
  readPairs <- separateReadPair(readPair)
  lenRead   <- nchar(readPairs[1])
  sfix      <- paste0(substr(readPairs[1], 2, lenRead), READ_PAIR_SEPARATOR_OUT, substr(readPairs[2], 2, lenRead))
  
  return(sfix)
  
}  # END suffixPaired


stringSpelledByGappedPatterns <- function(gappedPatterns, k, d) {
  
  sepReadPairs <- separateReadPairs(gappedPatterns)
  prefixString <- pathToGenome(sepReadPairs[[FIRST_PATTERNS]])
  suffixString <- pathToGenome(sepReadPairs[[SECOND_PATTERNS]])
  outString    <- prefixString
  
  splitPrefix <- strsplit(prefixString, EMPTY)[[1]]
  splitSuffix <- strsplit(suffixString, EMPTY)[[1]]
  
  for (i in (k+d+1):length(splitPrefix)) {
    
    if (splitPrefix[i] != splitSuffix[(i-k-d)]) 
      return(NULL)
    
  }  # END for
  
  outString <- paste0(outString, substr(suffixString, (nchar(suffixString)-k-d+1), nchar(suffixString)))
  
  return(outString)
  
}  # END stringSpelledByGappedPatterns
