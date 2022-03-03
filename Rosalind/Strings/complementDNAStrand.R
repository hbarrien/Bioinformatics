# Problem
# In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.
# The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, 
# then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").
# 
# Given: A DNA string s of length at most 1000 bp.
# Return: The reverse complement sc of s.
# Sample Dataset
# AAAACCCGGT
# Sample Output
# ACCGGGTTTT
#
# Usage with large input strings (Book: pp. 12):
# dna <- readChar("C:/Users/heba/Downloads/dataset.txt", file.size("C:/Users/heba/Downloads/dataset.txt")-2); complementDNAString(dna)

# A: 65, C: 67, G:71, T: 84
# create reverse values for the numeric representation of each base
reverse <- vector(length = 84)
reverse[65] <-  19
reverse[84] <- -19
reverse[67] <-  4
reverse[71] <- -4

complementDNAString <- function(dna) {
  
  # convert dna string to integer vector; reverse the integer vector
  revDNA <- rev(utf8ToInt(dna))
  
  # complement reverse dna
  compDNA <- sapply(revDNA, function(x) {x + reverse[x]})
  
  # convert integer vector back to string
  compDNA <- intToUtf8(compDNA)
  
  return(compDNA)
  
}  # complementDNAString
