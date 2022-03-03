# PatternCount - pp. 8 
patternCount <- function(txt, pattern) {
  
  count <- 0
  lenTxt <- nchar(txt)
  lenPattern <- nchar(pattern)
  
  for (i in 1:(lenTxt-lenPattern+1)) {
    
    if (substr(txt, i, (i+lenPattern-1)) == pattern)
      count <- (count + 1)
    
  }  # END for
  
  return(count)
  
}  # END patternCount