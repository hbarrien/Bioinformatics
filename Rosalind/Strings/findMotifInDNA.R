# Problem
# Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols 
# in s (as a result, t must be no longer than s).
#
# The position of a symbol in a string is the total number of symbols found to its left, including itself 
# (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). 
# The symbol at position i of s is denoted by s[i].
#
# A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions 
# of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".
#
# The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations 
# in s if it occurs more than once as a substring of s (see the Sample below).
#
# Given: Two DNA strings s and t (each of length at most 1 kbp).
# Return: All locations of t as a substring of s.

library(stringr)

BEGIN_POS <- 1
END_POS   <- 1
motifPositions <- c()

findMotifStartPositions <- function(motif, dna, numCharsToLeft) {
  
  if (nchar(motif) > nchar(dna))
    return(motifPositions)
  
  locationPositions <- str_locate(dna, motif)
  
  if (is.na(locationPositions[BEGIN_POS])) {
  
    return(motifPositions)
    
  } else {
    
    numCharsToLeft <- (numCharsToLeft + locationPositions[BEGIN_POS])
    motifPositions <<- c(motifPositions, numCharsToLeft)
    dna <- substr(dna, (locationPositions[BEGIN_POS] + 1), nchar(dna))
  }
  
  return(findMotifStartPositions(motif, dna, numCharsToLeft))
  
}  # END findMotifStartPositions

# PRECONDITION
# (dna is not  null) && (dna is not empty) && (motif is not null) && (motif is not empty) &&
# (dna is a valid string) && (motif is a valid string)
#
# RETURNS
# empty vector    - precondition not met OR motif not initially found in dna
# nonempty vector - position numbers of the beginning of motif in dna
findMotifInDNA <- function(motif, dna) {

  motifPositions <<- c()    
  numCharsToLeft <- 0
  
  # check precondition
  if (is.null(motif) || (nchar(motif) == 0) ||
      is.null(dna)   || (nchar(dna) == 0))
    return(motifPositions)
  
  # check precondition
  if (grepl("[^ACGT]", dna) || grepl("[^ACGT]", motif))
    return(motifPositions)

  return(findMotifStartPositions(motif, dna, numCharsToLeft))
  
}  # END findMotifInDNA
