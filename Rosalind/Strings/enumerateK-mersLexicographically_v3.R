# Problem
# Assume that an alphabet A has a predetermined order; that is, we write the alphabet 
# as a permutation A=(a1,a2,...,ak), where a1<a2<...<ak. For instance, the English 
# alphabet is organized as (A,B,...,Z).
#
# Given two strings s and t having the same length n, we say that s precedes t in the 
# lexicographic order (and write s<Lext) if the first symbol s[j] that doesn't match 
# t[j] satisfies sj<tj in A.
#
# Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive 
# integer n (n <= 10).
#
# Return: All strings of length n that can be formed from the alphabet, ordered lexicographically 
# (use the standard order of symbols in the English alphabet).

# ###################
# #### LIBRARIES ####
# ###################
library(data.table)
library(gtools)


# ###################
# #### CONSTANTS ####
# ###################
VERBOSE_1   <- "TRUE"
OUTPUT_PATH <- "/Bioinformatics//Strings/0009 enumerateK-mersLexico/k-mers.txt"


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
  
  write.csv(permsMatrix, permsFile, row.names = FALSE)
  
}  # END writePermsOutput

destroyPermsOutput <- function() {
  
  close(permsFile)
  
}  # END destroyPermsOutput


# #### PERMUTATIONS ####
createAlphabet <- function(alphabet) {

  indices <- c(1:length(alphabet))
  names(alphabet) <- indices
  
  return(alphabet)
  
}  # END createAlphabet

computePermutations <- function(alphabet, r) {
  
  v <- createAlphabet(alphabet)
  p <- permutations(length(alphabet), r, repeats.allowed = TRUE)
  
  nCols <- ncol(p)
  out   <- data.table()
  
  for (i in 1:nCols) {
    
    if (VERBOSE_1) print(paste0("computePermutations: level ", i))
    c <- as.character(lapply(p[,i], function(x) { return(v[x][[1]]) } ))
    out <- cbind(out, c)
    
  }  # END for
  
  if (VERBOSE_1) print("computePermutations: done")
  return(out)
  
}  # END computePermutations

lexicoPermutations <- function(alphabet, r, writeFile = FALSE) {
  
  if (is.null(alphabet) || is.null(r) || !is.numeric(r)) 
    return(NULL)

  if (VERBOSE_1) print("lexicoPermutations: calling computePermutations...")
  perms <- computePermutations(alphabet, r)
  
  # write to a file, in order to be able to copy output to Rosalind
  if (writeFile) {
    
    createPermsOutput()
    writePermsOutput(perms)
    destroyPermsOutput()
    
  }  # END if
  
  if (VERBOSE_1) print("lexicoPermutations: done.")
  return(perms)
  
}  # END lexicoPermutations
