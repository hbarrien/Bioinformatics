source(paste0(getwd(), "/Bioinformatics/Strings/0003 complementDNAStrand/complementDNAStrand.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/PatternToNumber.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/FrequencyArray_v2.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/hammingDistance.R"))

# Usage: Grand Challenge
# k <- 7
# L <- 500
# t <- 3
# d <- 1
# fName <- paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/Salmonella_enterica.txt")
# S_enterica <- inputFromFasta(fName)
# dna <- getDnaFromFasta(S_enterica)
# Optional: s <- skew(dna); plot(s)
# minSkewPos <- getMinSkewPos(dna, FALSE)[1]
# wL <- substr(dna, minSkewPos, (minSkewPos + (2*L) - 1))
# wC <- substr(dna, (minSkewPos - trunc((2*L)/2) + 1), (minSkewPos + trunc((2*L)/2) - 1))
# wR <- substr(dna, (minSkewPos - (2*L) + 1), minSkewPos)
# For all: wL, wC, wR : do findClump(), frequentWordsWithMismatches(), frequentWordsWithMismatchesAndComplements.
# Then compare results and extract the repeating terms as soutions


# *********************************
# *********** LIBRARIES ***********
# *********************************
library(seqinr)


# *********************************
# *********** CONSTANTS ***********
# *********************************
FWM_VERBOSE  <- TRUE
EMPTY_STRING <- ""


# *********************************
# *********** VARIABLES ***********
# *********************************


# *********************************
# *********** FUNCTIONS ***********
# *********************************

# ************** I/O **************
inputFromFasta <- function(fastaFName) {
  
  input <- read.fasta(file = fastaFName, as.string = TRUE)
  return(input)
  
}  # END inputFromFasta

getDnaFromFasta <- function(fastaObj) {
  
  return(toupper(fastaObj[[1]][1]))
  
}  # END getDnaFromFasta


# ********* DNA PROCESSING *********
frequentWordsWithMismatches <- function(dna, k, d) {
  
  if (FWM_VERBOSE) print("patterns: creating frequency array...")
  fa <- frequencyArray(k)
  
  lastWordStartPos <- (nchar(dna)-k+1)
  
  for (i in 1:lastWordStartPos) {
    
    if (FWM_VERBOSE) print(paste0("patterns: ", i, " of ", lastWordStartPos))
    
    startPos <- i
    stopPos  <- (i+k-1)
    pattern  <- substr(dna, startPos, stopPos)
    code     <- (patternToNumber(pattern) + 1)
    
    numMismatches <- countMismatches(pattern, dna, d)
    
    if (fa[[FREQ_IDX]][code] < numMismatches)
      fa[[FREQ_IDX]][code] <- numMismatches
    
  }  # END for
  
  idx <- which(fa[[FREQ_IDX]] == max(fa[[FREQ_IDX]]))
  res <- attributes(fa[[PATTERN_IDX]][idx])
  
  return(res)
  
}  # END frequentWordsWithMismatches

frequentWordsWithMismatchesAndComplements <- function(dna, k, d) {
  
  if (FWM_VERBOSE) print("patterns: creating frequency array...")
  fa <- frequencyArray(k)
  
  lastWordStartPos <- (nchar(dna)-k+1)
  
  for (i in 1:lastWordStartPos) {
    
    if (FWM_VERBOSE) print(paste0("patterns: ", i, " of ", lastWordStartPos))
    
    startPos <- i
    stopPos  <- (i+k-1)
    pattern  <- substr(dna, startPos, stopPos)
    code     <- (patternToNumber(pattern) + 1)
    
    # Calculate number of mismatches for the current kmer (i.e., pattern)
    numMismatchesKmer <- countMismatches(pattern, dna, d, FALSE)
    
    # Calculate number of mismatches for the reverse complement
    reverseComplement <- complementDNAString(pattern)
    numMismatchesRevComplement <- countMismatches(reverseComplement, dna, d, FALSE)
    
    if ((numMismatchesKmer > 0) && (numMismatchesRevComplement > 0))
      fa[[FREQ_IDX]][code] <- (numMismatchesKmer + numMismatchesRevComplement)
      
  }  # END for
  
  idx <- which(fa[[FREQ_IDX]] == max(fa[[FREQ_IDX]]))
  
  res <- c()
  lapply(idx, function(x) { 
                kmer <- attributes(fa[[PATTERN_IDX]][x])[[1]]
                revComplement <- complementDNAString(kmer)
                res <<- c(res, kmer, revComplement)
              })

  return(res)
  
}  # END frequentWordsWithMismatchesAndComplements