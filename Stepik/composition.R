# Usage:
# txt <- gsub("[^ACGT]", "", readChar("C:/Users/heba/Downloads/dataset_197_3.txt", file.size("C:/Users/heba/Downloads/dataset_197_3.txt")))
# k <- 100
# x <- composition(txt, k)
# writeChar(x, "C:/Users/heba/Downloads/dataset_197_3_OUT.txt", nchars = nchar(x))

source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/patternToNumber.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/numberToPattern.R"))

composition <- function(txt, k, sortLexico = FALSE) {
  
  lenTxt   <- nchar(txt)
  numKmers <- (lenTxt-k+1)
  out      <- c()
  
  # k == 25 is an estimated maximum value that numberToPattern() can handle to produce
  # correct results. Values of k > 25 will result in unpredictable output. Hence, the
  # warning message
  doSort <- sortLexico
  if (sortLexico && (k > 25)) {
    print(paste0("k: ", k, " is too large. The output will be unsorted."))
    doSort <- FALSE
  }
  
  for (idx in 1:numKmers) {
    
    kmer <- substr(txt, idx, (idx+k-1))
    out  <- c(out, ifelse(doSort, patternToNumber(kmer), kmer))
    
  }  # END for
  
  if (!doSort) return(out)
  
  out <- sort(out)
  out <- unlist(lapply(out, function(x) { return((numberToPattern(x, k))) } ))
  
  return(out)
  
}  # END composition