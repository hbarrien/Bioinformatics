# Usage:
# dnaStrands <- c("GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG")
# k <- 3
# t <- 5
# x <- greedyMotifSearch(dnaStrands, k, t)


source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Coursera-Stepik/Code challenges/profileMostProbableKmer.R")
source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Coursera-Stepik/Code challenges/hammingDistance.R")

# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
EMPTY_STRING      <- ""
DEFAULT_FILE_PATH <- "C:/Users/heba/Downloads/"

# bX : b stands for "base". So, bA = base A
bA <- "A"
bC <- "C"
bG <- "G"
bT <- "T"

A_IDX <- 1
C_IDX <- 2
G_IDX <- 3
T_IDX <- 4


# *********************************
# *********** VARIABLES ***********
# *********************************
GMS_VERBOSE <- FALSE


# *********************************
# *********** FUNCTIONS ***********
# *********************************

createMotifMatrix <- function(dnaStrands) {
  
  lenDnaStrand <- nchar(dnaStrands[1])
  motifMatrix   <- data.frame(matrix(ncol = lenDnaStrand))
  
  for (dnaStrand in dnaStrands) {
    
    dnaStrandSplit <- as.list(strsplit(dnaStrand, EMPTY_STRING)[[1]])
    motifMatrix    <- rbind(motifMatrix, dnaStrandSplit)
    
  }  # END for
  
  # delete the first row, which by default was created with NAs
  motifMatrix   <- motifMatrix[-1,]
  
  return(motifMatrix)
  
}  # END createMotifMatrix

createMotifCount <- function(motifMatrix, pseudocount = FALSE) {
  
  nCols        <- ncol(motifMatrix)
  motifCount   <- data.frame(matrix(nrow = 4, ncol = nCols))
  motifCount[] <- lapply(motifCount, function(x) {as.numeric(x)})
  
  # in order to avoid absolute 0.0 and 1.0 probabilities, add 1
  # to every count. See: Laplace's Rule, p. 89
  pCount <- ifelse(pseudocount, 1, 0)
  
  for (j in 1:nCols) {
    
    t <- table(motifMatrix[,j])
    a <- attributes(t)[["dimnames"]][[1]]
    
    motifCount[A_IDX,j] <- (ifelse((length(a[a==bA]) > 0), t[[bA]], 0) + pCount)
    motifCount[C_IDX,j] <- (ifelse((length(a[a==bC]) > 0), t[[bC]], 0) + pCount)
    motifCount[G_IDX,j] <- (ifelse((length(a[a==bG]) > 0), t[[bG]], 0) + pCount)
    motifCount[T_IDX,j] <- (ifelse((length(a[a==bT]) > 0), t[[bT]], 0) + pCount)
    
  }  # END for
  
  return(motifCount)
  
}  # END createMotifCount

createMotifProfile <- function(motifCount, t, pseudocount = FALSE) {
  
  nCols          <- ncol(motifCount)
  motifProfile   <- data.frame(matrix(nrow = 4, ncol = nCols))
  motifProfile[] <- lapply(motifProfile, function(x) {as.numeric(x)})
  
  ct <- sum(motifCount[,1])
  denominator <- ifelse(pseudocount, ct, t)
  
  for (j in 1:nCols) {
    
    motifProfile[A_IDX,j] <- (motifCount[A_IDX,j][[1]] / denominator)
    motifProfile[C_IDX,j] <- (motifCount[C_IDX,j][[1]] / denominator)
    motifProfile[G_IDX,j] <- (motifCount[G_IDX,j][[1]] / denominator)
    motifProfile[T_IDX,j] <- (motifCount[T_IDX,j][[1]] / denominator)
    
  }  # END for
  
  return(motifProfile)
  
}  # END createMotifProfile

getConsensusString <- function(motifMatrix) {
  
  numCols   <- ncol(motifMatrix)
  consensus <- ""
  
  for (j in 1:numCols) {
    
    nxtCol    <- motifMatrix[,j]
    baseCount <- table(nxtCol)
    
    if (length(baseCount[baseCount == max(baseCount)]) > 1) {
      nxtLetter <- attributes(baseCount[baseCount == max(baseCount)])$dimnames$nxtCol[1]
    } else { 
      nxtLetter <- attributes(baseCount[baseCount == max(baseCount)])
    }
    
    consensus <- paste0(consensus, nxtLetter)
    
  }  # END for
  
  return(consensus)
  
}  # END getConsensusString

score <- function(motifs) {
  
  if (is.null(motifs) || (length(motifs) == 0))
    return(-1)
  
  motifMatrix <- createMotifMatrix(motifs)
  numMotifs   <- nrow(motifMatrix)
  consensus   <- getConsensusString(motifMatrix)
  sc <- 0
  
  for (i in 1:numMotifs) {
    
    nxtRow   <- motifMatrix[i,]
    nxtMotif <- ""
    
    lapply(nxtRow, function(x) { nxtMotif <<- paste0(nxtMotif, x) })
    sc <- (sc + hammingDistance(consensus, nxtMotif))
        
  }  # END for
  
  return(sc)
  
}  # END score

greedyMotifSearch <- function(dnaStrands, k, t, pseudocount = FALSE) {
  
  bestMotifs  <- as.character(lapply(dnaStrands, function(x) {return(substr(x, 1, k))}))
  firstDna    <- dnaStrands[1]
  lenFirstDna <- nchar(firstDna)
  nKmers      <- (lenFirstDna-k+1)
  
  for (i in 1:nKmers) {
    
    motifs <- c(substr(firstDna, i, (i+k-1)))
    
    for (j in 2:t) {
      
      motifMatrix  <- createMotifMatrix(motifs)
      motifCount   <- createMotifCount(motifMatrix, pseudocount)
      motifProfile <- createMotifProfile(motifCount, t, pseudocount)
      
      nxtDna <- dnaStrands[j]
      probableKmer <- profileMostProbableKmer(nxtDna, k, motifProfile)
      
      if (is.null(probableKmer))
        stop("greedyMotifSearch: an error occurred.")
      
      if (length(probableKmer) == 0)
        next()
      
      motifs <- c(motifs, probableKmer)
      
    }  # END for
    
    if (GMS_VERBOSE) {
     
      print("bestMotifs:")
      print(bestMotifs)
      print(paste0("score(bestMotifs): ", score(bestMotifs)))
      print("motifs:")
      print(motifs)
      print(paste0("score(motifs): ", score(motifs)))
      print("=============================")
      
    }  # END if
    
    if (score(motifs) < score(bestMotifs))
      bestMotifs <- motifs
    
  }  # END for
  
  return(bestMotifs)
  
}  # END greedyMotifSearch
