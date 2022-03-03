# Usage:
  # t1 <- Sys.time()
  # patterns <- read.csv("C:/Users/heba/Downloads/dataset_200_8.txt", sep="", na.strings="", stringsAsFactors=FALSE)$kmer
  # x <- deBruijinGraphAdjacencyListFromPatterns(patterns)
  # writeOverlap(x)
  # t2 <- Sys.time()
  # (t2-t1)


source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/deBruijinGraph.R"))

createNodeListFromPatterns <- function(patterns) {
  
  nodeList <- c()
  
  for (pattern in patterns) {
    
    pfix <- prefix(pattern)
    sfix <- suffix(pattern)
    
    nodeList <- c(nodeList, pfix, sfix)
    
  }  # END for
  
  return(sort(unique(nodeList)))
  
}  # END createNodeListFromPatterns

computeAdjancencyFromPatterns <- function(patterns, dbMatrix) {
  
  for (pattern in patterns) {
    
    pfix <- prefix(pattern)
    sfix <- suffix(pattern)
    
    dbMatrix[pfix, sfix] <- (dbMatrix[pfix, sfix]+1)
    
  }  # END for
  
  return(dbMatrix)
  
}  # END computeAdjancencyFromPatterns

deBruijinGraphAdjacencyListFromPatterns <- function(patterns) {
  
  print("Creating node list...")
  nodeList <- createNodeListFromPatterns(patterns)
  
  print("Creating matrix...")
  dbMatrix <- deBruijinMatrix(nodeList)
  
  print("Computing adjacency...")
  dbMatrix <- computeAdjancencyFromPatterns(patterns, dbMatrix)
  
  print("Creating deBruijin graph...")
  dbGraph  <- createDeBruijinGraph(dbMatrix)
  
  return(dbGraph)
  
}  # END deBruijinGraphAdjacencyListFromPatterns
