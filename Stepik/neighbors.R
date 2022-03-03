source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/FrequencyArray_v2.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/hammingDistance.R"))

# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************


# *********************************
# *********** VARIABLES ***********
# *********************************


# *********************************
# *********** FUNCTIONS ***********
# *********************************
neighbors <- function(pattern, d) {
  
  k  <- nchar(pattern)
  fa <- frequencyArray(k)
  
  kmers <- attributes(fa[[PATTERN_IDX]])[[1]]
  n     <- as.logical(lapply(kmers, function(x) { return(hammingDistance(pattern, x) <= d) } ))
  
  return(kmers[n == TRUE])
  
}  # END neighbors

neighborsExactMatch <- function(pattern, d) {
  
  k  <- nchar(pattern)
  fa <- frequencyArray(k)
  
  kmers <- attributes(fa[[PATTERN_IDX]])[[1]]
  n     <- as.logical(lapply(kmers, function(x) { return(hammingDistance(pattern, x) == d) } ))
  
  return(kmers[n == TRUE])
  
}  # END neighborsExactMatch
