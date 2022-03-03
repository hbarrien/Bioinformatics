source(paste0(getwd(), "/Bioinformatics/Strings/0009 enumerateK-mersLexico/enumerateK-mersLexicographically_v3.R"))


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
FA_VERBOSE <- "TRUE"

b4letters <- c("A", "C", "G", "T")
names(b4letters) <- c(0, 1, 2, 3)

PATTERN_IDX <- "kmers"
FREQ_IDX    <- "freq"
WINDOW_IDX  <- "windows"
START_IDX   <- "start"
STOP_IDX    <- "stop"


# *********************************
# *********** VARIABLES ***********
# *********************************
lastk   <- -1
lenP    <- 0
kmerIdx <- NULL


# *********************************
# *********** FUNCTIONS ***********
# *********************************
frequencyArray <- function(k) {
  
  if (is.null(k) || !is.numeric(k) || (k < 1))
    return(NULL)
  
  if (k != lastk) {
    
    # p is lexicographically ordered. Indexing is base 1
    p <- lexicoPermutations(b4letters, k)
    
    lastk   <<- k
    lenP    <<- nrow(p)
    kmerIdx <<- seq(lenP)
    kmers   <- c()
    
    # Heuristic: lenP> 10000: write intermediate output to disk. Faster than adding to a vector in memory
    if (lenP > 10000) {
      
      fName <- paste0(getwd(),"/kmers.txt")
      f <- file(fName,"wa")
      writeLines("kmer",con=f, sep="\n")
      
    }  # END if
    
    for (i in 1:lenP) {
      
      print(paste0("frequencyArray: ", i, " of ", lenP))
      kmer <- paste0(p[i,], collapse = "")
      
      if (lenP > 10000) writeLines(kmer,con=f, sep="\n")
      else kmers <- c(kmers, kmer)
      
    }  # END for
    
    if (lenP > 10000) {
      
      close(f)
      x <- read.csv(fName)
      file.remove(fName)
      kmers <- as.character(x$kmer)
      
    }  # END if
    
    names(kmerIdx) <<- kmers
    
  }  # END if
  
  # Always create fresh counter sequences
  freq         <- integer(lenP)
  freqWindows  <- integer(lenP)
  startPos     <- integer(lenP)
  stopPos      <- integer(lenP)
  
  freqArray <- list(kmers=kmerIdx, freq=freq, windows=freqWindows, start=startPos, stop=stopPos)
  
  return(freqArray)
  
}  # END frequencyArray

computingFrequencies <- function(text, k) {
  
  if (FA_VERBOSE) print("computingFrequencies: creating frequency array...")
  freqArray <- frequencyArray(k)
  
  lenText <- (nchar(text)-k+1)
  for (i in 1:lenText) {
    
    if (FA_VERBOSE) print(paste0(i, " of ", lenText))
    pattern <- substr(text, i, (i+k-1))
    
    j <- (patternToNumber(pattern) + 1)
    freqArray[[2]][j] <- (freqArray[[2]][j] + 1) 
    
  }  # END for
  
  if (FA_VERBOSE) print("Done.")
  return(freqArray)
  
}  # END computingFrequencies
