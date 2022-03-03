
L_PAREN <- "("
R_PAREN <- ")"
PIPE_SYMBOL <- "|"

pairedComposition <- function(k, d, txt) {
  
  lenTxt <- nchar(txt)
  lenLongKmer <- (2*k+d)
  posLastLongKmer <- (lenTxt-lenLongKmer+1)
  
  longKmers <- c()
  
  for (i in 1: posLastLongKmer) {
    
    longKmer  <- substr(txt, i, (i+lenLongKmer-1))
    longKmers <- c(longKmers, longKmer)
    
  }  # END for
  
  longKmers   <- sort(longKmers)
  pairedKmers <- ""
  
  for (longKmer in longKmers) {
    
    kmer1 <- substr(longKmer, 1, k)
    kmer2 <- substr(longKmer, (k+d+1), nchar(longKmer))
    
    pairedKmers <- paste0(pairedKmers, L_PAREN, kmer1, PIPE_SYMBOL, kmer2, R_PAREN)
    
  }  # END for
  
  return(pairedKmers)
  
}  # END pairedComposition
