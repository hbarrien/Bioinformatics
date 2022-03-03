# Problem
# A permutation of length n is an ordering of the positive integers {1, 2, ..., n}. 
# For example, p = (5,3,2,1,4) is a permutation of length 5.
# 
# Given: A positive integer n <= 7.
# 
# Return: The total number of permutations of length n, followed by a list of all 
# such permutations (in any order).


# ###################
# #### CONSTANTS ####
# ###################
OUTPUT_PATH <- "/Bioinformatics/Combinatorics/0002 enumerateGeneOrders/perms.txt"


# ###################
# ### GLOBAL VARS ###
# ###################
permsFile <- NULL


# ###################
# #### FUNCTIONS ####
# ###################

# #### I/O ####
createPermsOutput <- function() {
  
  wd <- getwd()
  output_filePath <-  paste0(wd, OUTPUT_PATH)
  permsFile       <<- file(output_filePath, open = "wt")
  
  return(permsFile)
  
}  # END createPermsOutput

writePermsOutput <- function(permsMatrix) {
  
  write.csv(permsMatrix, permsFile)
  
}  # END writePermsOutput

destroyPermsOutput <- function() {
  
  close(permsFile)
  
}  # END destroyPermsOutput


# #### PERMUTATIONS ####
createAlphabet <- function(n) {
  
  return(c(1:n))
  
}  # END createAlphabet

numberPermutations <- function(n) {
  
  r <- n
  return(factorial(n) / factorial(n - r))
  
}  # numberPermutations

# courtesy of: https://www.r-bloggers.com/learning-r-permutations-and-combinations-with-base-r/
computePermutations <- function(v) {
  
  n <- length(v)
  
  if (n == 1) {
    
    return(v)
    
  } else {
    
    X <- NULL
    
    for (i in 1:n) 
      X <- rbind(X, cbind(v[i], computePermutations(v[-i])))
    
    return(X)
    
  }  # END if
  
}  # END computePermutations

permutations <- function(n) {
  
  if (is.null(n) || !is.numeric(n) || (n < 1)) 
    return(NULL)
  
  alphabet <- createAlphabet(n)
  perms    <- computePermutations(alphabet)
  
  # write to a file, in order to be able to copy output to Rosalind
  createPermsOutput()
  writePermsOutput(perms)
  destroyPermsOutput()
  
  return(perms)
  
}  # permutations
