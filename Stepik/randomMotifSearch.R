# Usage:
# dnaStrands <- c("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA")
# k <- 8; t <- 5; n <- 3; x <- randomMotifSearchDriver(dnaStrands, k, t, n)

source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Coursera-Stepik/Code challenges/greedyMotifSearch.R")


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
SEED              <- 123
EMPTY_STRING      <- ""
DEFAULT_FILE_PATH <- "C:/Users/heba/Downloads/"

# bX : b stands for "base". So, bA = base A
bA <- "A"
bC <- "C"
bG <- "G"
bT <- "T"


# *********************************
# *********** VARIABLES ***********
# *********************************
RMS_VERBOSE <- FALSE
RUN_VERBOSE <- TRUE


# *********************************
# *********** FUNCTIONS ***********
# *********************************
randomKmer <- function(dna, k) {
  
  if (is.null(dna) || (nchar(dna) < 2) || is.null(k) || !is.numeric(k) || (k < 1))
    return(NULL)
  
  numKmers <- (nchar(dna)-k+1)
  i        <- sample(1:numKmers, 1)
  kmer     <- substr(dna, i, (i+k-1))
  
  return(kmer)
  
}  # END randomKmer

randomMotifs <- function(dnaStrands, k, t) {
  
  rMotifs <- c()
  
  for (i in 1:t) {
    
    nxtDna <- dnaStrands[i]
    rKmer  <- randomKmer(nxtDna, k)
    rMotifs <- c(rMotifs, rKmer)
    
  }  # END for
  
  return(rMotifs)
  
}  # END randomMotifs

randomMotifSearch <- function(dnaStrands, k, t, setSeed = FALSE) {
  
  if (setSeed) set.seed(SEED)
  
  motifs          <- randomMotifs(dnaStrands, k, t)
  bestMotifs      <- motifs
  scoreBestMotifs <- score(bestMotifs)
  
  if (RMS_VERBOSE) {
    print("initial bestMotifs:")
    print(bestMotifs)
    print(paste0("score initial bestMotifs: ", scoreBestMotifs))
  }
  
  iterations <- 1
  
  while (TRUE) {
    
    motifMatrix   <- createMotifMatrix(motifs)
    motifCount    <- createMotifCount(motifMatrix, TRUE)
    profileMatrix <- createMotifProfile(motifCount, t, TRUE)
    
    motifs <- c()
    lapply(dnaStrands, function(x) {
      probableKmer <- profileMostProbableKmer(x, k, profileMatrix)
      motifs       <<- c(motifs, probableKmer)
    })
    
    if (RMS_VERBOSE) {
      print(paste0("motifs[", iterations, "]:"))
      print(motifs)
      
      print(paste0("score(motifs): ", score(motifs)))
      print(paste0("score(bestMotifs): ", score(bestMotifs)))
      print("========")
    }
    
    scoreMotifs <- score(motifs)
    
    if (scoreMotifs >= score(bestMotifs))
      break()
    
    bestMotifs      <- motifs
    scoreBestMotifs <- scoreMotifs
    iterations      <- (iterations + 1)
    
  }  # END while
  
  return(list(bestMotifs=bestMotifs, score=scoreBestMotifs))
  
}  # END randomMotifSearch

randomMotifSearchDriver <- function(dnaStrands, k, t, n, setSeed = FALSE) {
  
  bestRandom <- NULL
  
  for (i in 1:n) {
    
    if(RUN_VERBOSE) print(paste0("run: ", i))
    
    res <- randomMotifSearch(dnaStrands, k, t)
    
    if (is.null(bestRandom)) {
      bestRandom <- res
      next()
    }
    
    if (res$score < bestRandom$score)
      bestRandom <- res
    
  }  # END for
  
  return(bestRandom$bestMotifs)
  
}  # END randomMotifSearchDriver
