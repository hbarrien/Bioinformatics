# Exercise break, pp. 68 of Book

# Nucleotide alphabet
A <- c("A", "C", "G", "T")

# Probability of one nucleotide appearing in a given sequence
pNucleotide <- 1/length(A)

# Number of nucleotides in each proposed string
N <- 1000

# In the exercise, kmer ::== 9mer
k <- 9

# Number of kmers in a string with N nucleotides
nKmers <- (N-k+1)

# Probability of one kmer in a given string
pKmer  <- (pNucleotide**k)

# Number of occurences of kmer within the set of nKmers
nOccKmer <- (pKmer * nKmers)

# Number of occurrences of kmer in the proposed 500 sets of nKmers
nOccKmer * 500
