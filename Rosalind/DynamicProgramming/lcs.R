## path <- "C:/Users/heba/Downloads/rosalind_lcsq.txt"

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


lcs_len <- function(p, q) {
  
  pLen <- nchar(p)
  qLen <- nchar(q)
  
  m <- matrix(data = 0, nrow = (pLen+1), ncol = (qLen+1))
  
  for (i in pLen:1) {
    
    for (j in qLen:1) {
      
      if ((i > pLen) || (j > qLen)) m[i, j] <- 0
      else if (substr(p, i, i) == substr(q, j, j)) m[i, j] <- (1 + m[i+1, j+1])
      else m[i, j] <- max(m[i+1, j], m[i, j+1])
      
    }  # END for
    
  }  # END for
  
  return(m)
  
}  # END lcs_len


lcs_getSeq <- function(p, q, m) {
  
  s <- EMPTY
  i <- 1
  j <- 1
  
  pLen <- nchar(p)
  qLen <- nchar(q)
  
  while ((i <= pLen) && (j <= qLen)) {
    
    x <- substr(p, i, i)
    y <- substr(q, j, j)
  
    if (x == y) {
      
      s <- paste0(s, substr(p, i, i))
      i <- (i+1)
      j <- (j+1)
      
    } else if (m[i+1,j] >= m[i,j+1]) {
      
      i <- (i+1)
      
    } else {
      
      j <- (j+1)
      
    }  # END if
    
  }  # END while
  
  return(s)
  
}  # END lcs_getSeq


lcs <- function(path) {
  
  input <- readFastaInput(path)
  p <<- input[[1]][1]
  q <<- input[[2]][1]
  
  m <- lcs_len(p, q)

  return(lcs_getSeq(p, q, m))
  
}  # END lcs
