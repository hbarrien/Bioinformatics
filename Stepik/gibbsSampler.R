# Usage:
# k <- 8
# t <- 5
# N <- 100
# dnaStrands <- c("CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA")
# t1 <- Sys.time(); gibbsSampler(dnaStrands, k, t, N); t2 <- Sys.time(); (t2-t1)
#
# OPTIMIZED USAGE: 
# dnaStrands <- (readInput("dataset_163_4.txt"))$x
# k <- 15
# t <- 20
# N <- 2000
# t1 <- Sys.time()
# initMotifs <- randomMotifSearchDriver(dnaStrands, k, t, 200)
# gibbsSampler(dnaStrands, k, t, N, initMotifs)
# t2 <- Sys.time()
# (t2-t1)

source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Coursera-Stepik/Code challenges/randomMotifSearch.R")


# *********************************
# *********** LIBRARIES ***********
# *********************************
library(R.utils)


# *********************************
# *********** CONSTANTS ***********
# *********************************


# *********************************
# *********** VARIABLES ***********
# *********************************
GIBBS_VERBOSE <- FALSE
RUN_GIBBS_VERBOSE <- TRUE


# *********************************
# *********** FUNCTIONS ***********
# *********************************
gibbsSampler <- function(dnaStrands, k, t, N, initMotifs = NULL,  setSeed = FALSE) {
  
  if (setSeed) set.seed(SEED)
  
  if (is.null(initMotifs)) {
    motifs <- randomMotifs(dnaStrands, k, t)
  } else {
    motifs <- initMotifs
  }
  
  bestMotifs      <- motifs
  motifPrs        <- numeric(t)
  bestMotifPrs    <- motifPrs
  scoreBestMotifs <- score(bestMotifs)
  
  if (GIBBS_VERBOSE) {
    print("initial bestMotifs:")
    print(bestMotifs)
    print(paste0("score initial bestMotifs: ", scoreBestMotifs))
  }
  
  for (j in 1:N) {
    
    if (RUN_GIBBS_VERBOSE) print(paste0("run: ", j))
    
    currMotifs <- motifs
    
    i        <- sample(1:t, 1)
    delMotif <- motifs[i]
    motifs   <- motifs[-i]
    
    if (GIBBS_VERBOSE) {
      print(paste0("removed i:", i, ". altered motifs:"))
      print(motifs)
      print("dnaStrand to recalculate:")
      print(dnaStrands[i])
    }
    
    motifMatrix   <- createMotifMatrix(motifs)
    motifCount    <- createMotifCount(motifMatrix, TRUE)
    profileMatrix <- createMotifProfile(motifCount, t, TRUE)
    
    probableMotif <- profileMostProbableKmer(dnaStrands[i], k, profileMatrix, includePr = TRUE)
    
    if (probableMotif[["kmer"]] == delMotif) {
      motifs <- insert(motifs, i, delMotif)
      next()
    }
    
    if (probableMotif[["pr"]] < motifPrs[i]) {
      motifs <- insert(motifs, i, delMotif)
      next()
    }
    
    # At this point, ((probableMotif[["kmer"]] != delMotif) && (probableMotif[["pr]] >= delMotif))
    if (score(currMotifs) <= score(motifs))
      next()
    
    # At this point, (score(motifs) < score(currMotifs))
    motifs <- insert(motifs, i, probableMotif[["kmer"]])
    motifPrs[i] <- probableMotif[["pr"]]
    
    if (GIBBS_VERBOSE) {
      print(paste0("probableMotif:", probableMotif, ". restored motifs:"))
      print(motifs)
      
      print(paste0("score(motifs): ", score(motifs)))
      print(paste0("score(bestMotifs): ", score(bestMotifs)))
      print("========")
    }
    
    scoreMotifs <- score(motifs)
    
    if (scoreMotifs < scoreBestMotifs) {
      
      bestMotifs      <- motifs
      bestMotifPrs    <- motifPrs
      scoreBestMotifs <- scoreMotifs
      
    }  # END if
    
  }  # END for
  
  return(list(bestMotifs=bestMotifs, score=scoreBestMotifs))
  
}  # END gibbsSampler
