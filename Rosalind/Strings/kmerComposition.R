# Usage
# library(seqinr)
# dna <- read.fasta("C:/Users/heba/Downloads/rosalind_kmer.txt", as.string = TRUE, forceDNAtolower = FALSE)[[1]][1]
# k <- XX
# x <- kmerComposition(dna, k)
# writeChar(as.character(x), "C:/Users/heba/Downloads/rosalind_kmer_out.txt")

# PROBLEM
# For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.
# 
# Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times 
# that the mth k-mer (with respect to the lexicographic order) appears in s.
#
# Given: A DNA string s in FASTA format (having length at most 100 kbp).
#
# Return: The 4-mer composition of s.


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# ************* SOURCE ************
# *********************************
source(paste0(getwd(),"/Bioinformatics/Coursera-Stepik/Code challenges/composition.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/FrequencyArray_v2.R"))


# *********************************
# *********** FUNCTIONS ***********
# *********************************
kmerComposition <- function(dna, k) {
  
  fArray <- frequencyArray(k)
  kmers  <- composition(dna, k)
  
  lapply(kmers, function(kmer) {
    
    fArray$freq[[fArray$kmers[[kmer]]]] <<- (fArray$freq[[fArray$kmers[[kmer]]]] + 1)
    
  })
  
  return(fArray$freq)
  
}  # END kmerComposition