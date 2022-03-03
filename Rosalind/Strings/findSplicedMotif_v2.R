# Usage (quick test):
# dna <- "ACGTACGTGACG"
# motif <- "GTA"
# x <- createMotifIndices(dna, motif)
# y <- getSplicedMotifIndexSequences(dna, motif)
#
# Timed test:
# dna <- "CAATGGCTCGGCCCTCTAGGCTGCGAGTATCGAACGTGCTCTGATATAATCGCCTTGACTGATGCGCGAGGCGTGGGTCCATCCACATTGGATCCTGACTAACCGTTTGGCCGTTTAGTCGTCTAATAAGGCATCACGTCTAGCCTTCCTCTGAGTGTCCCTATGTGCCAAGGTCCACAGGTAATTTTTAATGGTATCTGTTGTAGGCGGCGTTCCTTCCGGGCTGATTTATTGTACACAGACGTCCGAGAGGAACTATAACCGTGGGTGCCGGATTGAGCATTTCTCCCACAGCGTCTCACCCATAGTATAATACGTTTATCCAAAAGAAACAGCTCCAACTCCGCATGTATGTCTTCGAAGGGAAGAATTTTATGTTGTCCAATGACCTCTTACGCCTTAATCGACTCTAGGACCTTGAATCGCGTCCACCATCTAGCAGCCGCCCTGGACCTCGTGGAGGGGGCCGATCGTTTAGTGTCCCGCACAACATCAGCACATCGTTAAGGTTACTTTTAGAGGGTGGAACCTGATGGCACCTGGTAACAGTCCTCTTGCTGAGGTTCAACTAGTGAGCAGTGATTTTCATAACGACCTTGTGTATACGGTGCATACTAACGCATGCATATTCCGCATGGCGCCCGACAGCGACCATATCTAATGATGGATGCGTAAACCGTAAGGAGAAGGAGACTTCCCCGCTGCTCATCCATGAGAGACCCATATCCGGTTCGATGCGGCCACACGTGCATCGCAATTCCGGACATCTCGGGGAGTTGGAGATTTCTGGAGTAGGGCTTCCCCGTACTCGTAAAGCACGATGAGGCATGCAAAGGTTAGCACTTATGTAAATGTATATGGTCGATAGCACTTCTGATGCGAACTTTAGAGGAGGATTGCGTTTAGTCCGCTAGCTCTACTTCGGGGTACACAGTATTACCTCGGAACCTCTGATCTAAGAGGCGGTACTACGATCGCTTGAAT"
# motif <- "TATTCCTCCAGGGGTCAACCATGCCTACAGTCAATCACTGGATTTGTCTAATGATGCGTCGTTAAGGAAACCA"
# x <- createMotifIndices(dna, motif)
# t1 <- Sys.time(); y <- getSplicedMotifIndexSequences(dna, motif); t2 <- Sys.time(); (t2-t1)


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
SM_VERBOSE   <- TRUE
EMPTY_STRING <- ""


# *********************************
# *********** VARIABLES ***********
# *********************************
out <- NULL

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

getMotifFromFasta <- function(fastaObj) {
  
  return(toupper(fastaObj[[2]][1]))
  
}  # END getMotifFromFasta

readInputFromFasta <- function(fastaFName) {
  
  input <- inputFromFasta(fastaFName)
  dna   <- getDnaFromFasta(input)
  motif <- getMotifFromFasta(input)
  
  return(list(dna, motif))
  
} # END readInputFromFasta


# ************** spliced motif **************
createMotifIndices <- function(dna, motif) {
  
  splitDna    <- strsplit(dna, EMPTY_STRING)[[1]]
  splitMotif   <- strsplit(motif, EMPTY_STRING)[[1]]
  motifIndices <-c()
    
  for (chMotif in splitMotif) {
    
    mIdx <- grep(chMotif, splitDna, fixed = TRUE)
    if (length(mIdx) == 0) return(NULL)
    
    motifIndices <- c(motifIndices, list(mIdx))
    
  }  # END for
  
  return(motifIndices)
  
}  # END createMotifIndices

# PRECONDITION
# For every sublist in motifIndices: (|sublist| >= 1) AND (sublist sorted ASC)
computeSplicedMotifIndexSequences <- function(motifIndices, listIdx = 0, idxSeq = NULL) {
  
  if (listIdx == 0) {
    
    listIdx <- (listIdx+1)
    for (i in motifIndices[[listIdx]]) {
      
      idxSeq <- c(i)
      out <- computeSplicedMotifIndexSequences(motifIndices, (listIdx+1), idxSeq)
      
    }  # END for
    
    # POSTCONDITION
    # ((|out| MOD |motifIndices|) == 0)
    return(out)
    
  }  # END if
  
  lastIdx <- tail(idxSeq, 1)
  mIdx    <- motifIndices[[listIdx]][motifIndices[[listIdx]] > lastIdx]
  lenMotifIndices <- length(motifIndices)
  
  for (j in mIdx) {
    
    nxtIdxSeq <- idxSeq
    nxtIdxSeq[length(nxtIdxSeq)+1] <- j
    
    if (length(nxtIdxSeq) == lenMotifIndices) {
      
      out[[length(out)+1]] <<- nxtIdxSeq
      if (SM_VERBOSE) print(nxtIdxSeq)
      
    } else {
      
      computeSplicedMotifIndexSequences(motifIndices, (listIdx+1), nxtIdxSeq)
      
    }  # END if
    
  }  # END for
  
  return(out)
  
}  # END computeSplicedMotifIndexSequences

getSplicedMotifIndexSequences <- function(dna, motif) {
  
  out <<- list()
  motifIndices <- createMotifIndices(dna, motif)
  if (is.null(motifIndices)) return(NULL)
  
  x <- computeSplicedMotifIndexSequences(motifIndices)
  return(x)
  
}  # END getSplicedMotifIndexSequences
