# Problem
# A common substring of a collection of strings is a substring of every member of the collection. 
# We say that a common substring is a longest common substring if there does not exist a longer 
# common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but 
# it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" 
# and "AACCGTATA".
#
# Note that the longest common substring is not necessarily unique; for a simple example, "AA" and 
# "CC" are both longest common substrings of "AACC" and "CCAA".
# 
# Given: A collection of k (k <= 100) DNA strings of length at most 1 kbp each in FASTA format.
#  
# Return: A longest common substring of the collection. (If multiple solutions exist, you may return 
# any single solution.)


# ###################
# ## ASSUMPTIONS ####
# ###################
# 1. all DNA strings in k are of the same length (improvement: if not, choose the shortest string)


# ###################
# #### LIBRARIES ####
# ###################
library(Biostrings)


# ###################
# #### CONSTANTS ####
# ###################
DEFAULT_FILE_PATH <- "/Bioinformatics/Strings/0008 findSharedMotif/"


# ###################
# ##### TESTING #####
# ###################
TESTING <- FALSE
motifs  <- NULL
p <- NULL
q <- NULL


# ###################
# #### FUNCTIONS ####
# ###################

# #### I/O FUNCTIONS ####
# readInput
#   Access: private
readInput <- function(inputSource) {
  
  if (is.null(inputSource))
    return(NULL)
  
  wd <- getwd()
  input_filePath <- paste0(wd, DEFAULT_FILE_PATH, inputSource)  
  
  dna <- readDNAStringSet(input_filePath)
  return(dna)
  
}  # END readInput

# #### MOTIF FUNCTIONS ####
# computeSharedMotifs
# Access: private
computeSharedMotifs <- function(dnaStrings) {
  
  # get the first dna string, to be used as the dna substring pattern creator
  leadDnaString <- dnaStrings[1]
  
  # set parse indices
  fIdx <- 1  # first index
  lIdx <- 2  # last index
  
  numDnaStrings    <- length(dnaStrings)
  lenLeadDnaString <- nchar(leadDnaString)
  outSharedMotifs  <- c()
  longestMotif     <- 0
  
  while (lIdx <= lenLeadDnaString) {
  
    # get the next substring
    nxtDnaPattern <- substr(leadDnaString, fIdx, lIdx)
    
    if (TESTING) print(nxtDnaPattern)
    
    found <- (sum(grepl(nxtDnaPattern, dnaStrings, fixed=TRUE) == TRUE) == numDnaStrings)
    
    if (!found) {
      
      # reset first and last indices to start a new substring
      if ((lIdx - fIdx) == 1) {
        
        fIdx <- lIdx
        lIdx <- (lIdx + 1)
        
      } else {
        
        fIdx <- (fIdx + 1)
        
      }  # END if
        
    } else {
      
      lenNxtDnaPattern <- nchar(nxtDnaPattern)
      
      if (lenNxtDnaPattern >= longestMotif) {
      
        # add the shared dna substring (i.e., motif) to the output
        outSharedMotifs <- c(outSharedMotifs, nxtDnaPattern)
        
        longestMotif <- lenNxtDnaPattern
        if (TESTING) print(longestMotif)
        
      }  # END if
      
      # increment the last index, in order to obtain a new substring
      lIdx <- (lIdx + 1)
      
    } # END if
    
  }  # END while
  
  finalResult <- NULL
  if (length(outSharedMotifs) >= 1) {
    
    if (length(outSharedMotifs) == 1) {
      
      finalResult <- outSharedMotifs[1]
      
    } else {
    
      # eliminate redundancy
      outSharedMotifs <- unique(outSharedMotifs)
      
      # get the lengths of all resulting dna substrings
      p <<- nchar(outSharedMotifs)
      
      # get the positions of the longest dna substrings
      q <<- which(p == max(p))
      
      # get the longest strings having equal length
      finalResult <- outSharedMotifs[q]
      
      motifs <<- outSharedMotifs
      
    }  # END if
    
  }  # END if
  
  return(finalResult)
    
}  # END computeSharedMotifs

# getSharedMotifs
# Access: public
getSharedMotifs <- function(inputSource) {
  
  if (is.null(inputSource))
    return(NULL)
  
  dnaStrings <- readInput(inputSource)
  dnaStrings <- as.character(as.vector(dnaStrings))
  
  if (is.null(dnaStrings) || (length(dnaStrings) < 2)) 
    return(NULL)
  
  return(computeSharedMotifs(dnaStrings))
  
}  # END getSharedMotifs

# testGetSharedMotifs
# Access: public
testGetSharedMotifs <- function(inputSource) {
  
  t0 <- Sys.time()
  result <- getSharedMotifs(inputSource)
  t1 <- Sys.time()
  
  print(t1-t0)
  return(result)
  
}  # END testGetSharedMotifs