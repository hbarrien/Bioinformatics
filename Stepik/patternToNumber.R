
BASE <- 4

b4codes <- c(0, 1, 2, 3)
names(b4codes) <- c("A", "C", "G", "T")


patternToNumber <- function(pattern) {
  
  l <- strsplit(pattern, "")[[1]]
  codes    <- as.integer(lapply(l, function(x) { return(b4codes[x][[1]]) } ))
  lenCodes <- length(codes)
  
  nCode <- 0
  for (i in 1:(lenCodes-1))
    nCode <- ((nCode + codes[i]) * BASE)
  
  nCode <- (nCode + codes[lenCodes])
  
  return(nCode)
  
}  # END patternToNumber