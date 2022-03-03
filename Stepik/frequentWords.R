# FrequentWords - pp. 9
library(data.table)

# Version 1
# fPath <- paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/PatternCount.R")
# source(fPath)
#
# frequentWords <- function(txt, k) {
#   
#   fw <- data.table(w = character(), f = integer())
#   
#   lenTxt <- nchar(txt)
#   for (i in 1:(lenTxt-k+1)) {
#     
#     pattern <- substr(txt, i, (i+k-1))
#     fw <- rbind(fw, list(pattern, patternCount(txt, pattern)))
#                        
#   }  # END for
#   
#   fw <- unique(fw)
#   fw <- fw[order(f, decreasing = TRUE)]
#   
#   fValue <- fw[1,]$f
#   fw <- fw[fw$f == fValue,]
#   
#   return(fw$w)
#   
# }  # END frequentWords


frequentWords <- function(txt, k, t = 0) {
  
  fw <- data.table(w = character(), f = integer())
  
  lenTxt <- nchar(txt)
  for (i in 1:(lenTxt-k+1)) {
    
    pattern <- substr(txt, i, (i+k-1))
    fw <- rbind(fw, list(pattern, 0))
    
  }  # END for
  
  fw <- fw[order(w, decreasing = FALSE),]
  
  lenFW <- nrow(fw)
  currWord   <- fw[1,]$w
  countWords <- 1
  
  for (i in 2:lenFW) {
    
    if (fw[i]$w == currWord) {
      countWords <- (countWords+1)
      next()
    }
    
    fw[(i-1),]$f <- countWords
    currWord   <- fw[i,]$w
    countWords <- 1
    
  }  # END for
  
  maxWC <- t
  if (t == 0) maxWC <- max(fw$f)
  
  return(fw[fw$f == maxWC,]$w)
  
}  # END frequentWords
