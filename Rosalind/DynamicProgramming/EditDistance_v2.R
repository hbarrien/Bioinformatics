## path <- "C:/Users/heba/Downloads/rosalind_edit.txt"


# ################# LIBRARIES #################
library(seqinr)


# ################# CONSTANTS #################
EMPTY <- ""


# ################# VARIABLES #################
p <- NULL
q <- NULL


# ################# FUNCTIONS #################
readFastaInput <- function(path) {
  
  return(read.fasta(file = path, as.string = TRUE))
  
}  # END readFastaInput


ed <- function(p, q) {

  ps <- strsplit(p, EMPTY)[[1]]
  qs <- strsplit(q, EMPTY)[[1]]
  
  pLen <- length(ps)
  qLen <- length(qs)
  
  m <- matrix(data = 0, nrow = (pLen+1), ncol = (qLen+1))
  
  for (i in 1:pLen+1) m[i, 1] <- i-1
  for (j in 1:qLen+1) m[1, j] <- j-1
  
  for (j in 2:ncol(m)) {
    
    for (i in 2:nrow(m)) {
      
      substitutionCost <- ifelse(ps[i-1] == qs[j-1], 0, 1)
      
      m[i, j] <- min(m[i-1, j] + 1,
                     m[i, j-1] + 1,
                     m[i-1, j-1] + substitutionCost)
      
    }  # END for
    
  }  # END for
  
  return(m[nrow(m), ncol(m)])
    
}  # END ed


editDistance <- function(path) {
  
  input <- readFastaInput(path)
  p <<- input[[1]][1]
  q <<- input[[2]][1]
  
  return(ed(p, q))
  
}  # END editDistance
