# Problem
# Say that we have strings s=s1s2,...,sm and t=t1t2,...,tn with m < n. 
# Consider the substring t'=t[1:m]. We have two cases:
#
# If s=t', then we set s < Lext because s is shorter than t (e.g., APPLE < APPLET).
# Otherwise, s != t'. We define s < Lext if s < Lext' and define s> Lext if s>Lext' 
# (e.g., APPLET<LexARTS because APPL < LexARTS).
#
# Given: A permutation of at most 12 symbols defining an ordered alphabet A and a 
# positive integer n (n <= 4).
#
# Return: All strings of length at most n formed from A, ordered lexicographically. 
# (Note: As in "Enumerating k-mers Lexicographically", alphabet order is based on 
# the order in which the symbols are given.)


# ###################
# #### LIBRARIES ####
# ###################
library(data.table)


# ###################
# #### CONSTANTS ####
# ###################
A <- 65
C <- 67
G <- 71
T <- 84

LOWERCASE_A <- 97
OUTPUT_PATH <- "/Bioinformatics/Strings/0010 orderStringsLexico/ordered_str.txt"


# ###################
# ### GLOBAL VARS ###
# ###################
outFile <- NULL
alphabetCoding <- NULL


# ###################
# #### FUNCTIONS ####
# ###################

# #### I/O ####
createOutput <- function() {
  
  wd <- getwd()
  output_filePath <-  paste0(wd, OUTPUT_PATH)
  outFile       <<- file(output_filePath, open = "wt")
  
  return(outFile)
  
}  # END createOutput

writeOutput <- function(permsMatrix) {
  
  write.csv(permsMatrix, outFile)
  
}  # END writeOutput

destroyOutput <- function() {
  
  close(outFile)
  
}  # END destroyOutput

# #### PERMUTATIONS ####
numberPermutations <- function(n, r) {
  
  nPerms <- 0
  
  for (i in 1:r)
    nPerms <- (nPerms + n^i)
  
  return(nPerms)
  
}  # numberPermutations

# createAlphabet: using the alphabet argument as base, creates a new alphabet
# with each entry replicated r times
createAlphabet <- function(alphabet, r) {
  
  return(rep(alphabet,each=r))
  
}  # END createAlphabet

# createAlphabetCoding: assigns an ascending numeric code to each alphabet entry,
# starting from the first one
createAlphabetCoding <- function(a) {
  
  alphabetCoding <<- NULL
  
  # use lowercase A to avoid repeating code conflicts in both columns
  code <- LOWERCASE_A  
  
  for (i in a) {
    
    alphabetCoding <<- rbind(alphabetCoding, cbind(utf8ToInt(i), code))
    code <- (code + 1)
    
  }  # END for
  
  colnames(alphabetCoding) <<- c("letter.int.value", "code")
  
}  # END createAlphabetCoding

computePermutations <- function(v, r) {
  
  n <- length(v)
  
  if (r == 1) {
    
    return(v)
    
  } else {
    
    X <- NULL
    
    for (i in 1:n)
      X <- rbind(X, cbind(v[i], computePermutations(v[-i], (r-1))))
    
    return(unique(X))
    
  }  # END if
  
}  # END computePermutations

computeVariableLengthPermutations <- function(alphabet, r) {
  
  variableLengthPerms <- data.table(word=character(), code=character())
  
  # for r == 1, the return structure is a vector. thus, treat independently
  perms <- computePermutations(alphabet, 1)
  
  for (p in perms) {
    code <- alphabetCoding[which(alphabetCoding == utf8ToInt(p)), "code"]
    variableLengthPerms <- rbind(variableLengthPerms, list(p, intToUtf8(code)))
  }
  
  # for r > 1, the return structure is a matrix
  i <- 2
  while (i <= r) {
    
    v <- createAlphabet(alphabet, r)
    perms <- computePermutations(v, i)
    
    for (j in 1:nrow(perms)) {
      
      word <- ""
      code <- ""
      lapply(perms[j,], 
             function(x) { word <<- paste0(word, x)
                           code <<- paste0(code, intToUtf8(alphabetCoding[which(alphabetCoding == utf8ToInt(x)), "code"])) })
      
      variableLengthPerms <- rbind(variableLengthPerms, list(word, code))
      
    }  # END for
    
    i <- (i + 1)
    
  }  # END while
  
  return(variableLengthPerms[order(code), "word"])

}  # END computeVariableLengthPermutations

# permutations: main function
permutations <- function(alphabet, r) {
  
  if (is.null(alphabet) || (length(alphabet) == 0) || is.null(r) || !is.numeric(r) || !((r >= 1) && (r <= length(alphabet))))
    return(NULL)
  
  createAlphabetCoding(alphabet)
  perms <- computeVariableLengthPermutations(alphabet, r)
  
  # write to a file, in order to be able to copy output to Rosalind
  createOutput()
  writeOutput(perms)
  destroyOutput()
  
  return(perms)
  
}  # permutations
