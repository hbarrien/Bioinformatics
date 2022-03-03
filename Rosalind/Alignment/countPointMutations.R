# Evolution as a Sequence of Mistakes
# 
# Problem
# Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), 
# is the number of corresponding symbols that differ in s and t.
# 
# Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
# 
# Return: The Hamming distance dH(s,t).
# 
# Sample Dataset
# GAGCCTACTAACGGGAT
# CATCGTAATGACGGCCT
# 
# Sample Output
# 7

FALSE_POS <- 1
TRUE_POS  <- 2

hammingDistance <- function(s, t) {
  
  if (is.null(s) || is.null(t) || (length(s) != length(t)))
    return(NULL)
  
  s_nbr <- utf8ToInt(as.character(s))
  t_nbr <- utf8ToInt(as.character(t))
  
  ncDiff <- (s_nbr == t_nbr)
  return(table(ncDiff)[[FALSE_POS]])
  
}  # END hammingDistance