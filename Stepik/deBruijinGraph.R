# Usage:
# txt <- "ACGGAGCGTACGAAGAGGTTTTACACAGCCGGGCGAGATCAAGGTCATCCCGATGCTGTGTAGCTAGATCTGTATTGCCTTCGTTGATTCCCAACGTTATTATGCGTATCGATTTCGACTACCGGTTGTAGGGTTTCGGACGACTACCCCTCCCGGGCATACTCTTAGTCACTTCTACCCGGCATATCGCACACTCACCCAATAGGGAATCAGTTGGGCAAAGTCAGAGACAAACGCCCGAATGGGATAGACTGCTACTATATAGGGCAACGGCATGAGCCCGTTTGCGGAATCTACCATTCGGTAGATTAACGCAATAGGGACTTGTAGTCGGAATACCCAACTTGGGCCTGCAATTACGGATCGCGAGGAACCCGGAAATTCACCAACTAGCCCTTAGTGTTGAGATCGGATGCTACGAAAATCACGGTGAATCGGTAGGTACAGTCAGTATGCCATACACCAATTAAGAGCGCCGAGATGTACGCCACTCACTCGCTTTGCTACTGTGAGACTCGCGCCGGATCCGCTGTACCTTCATATCCTAGCAAACCCTGTTGTATGGTCGGCCCGGTGCGTGCGCCATGTCATCAGACGTAGAAGGCCTGCCGGCATGACTGACGCACTAGCCGAGAGAATACAACTATGCTACCAGACGCCAAACTATTTACGCGGTCAACATGTTTAACCAGTCTCCATTATTATCATCTTCGCTAATCAGCGTGTCTTTTCCCTAGCACCGTCGTGTCGGTCGATAGATTATGATGAAGAGGTTGATGAAAATGCTAGGCGACACCCCTGTACGTGGTGGGAGGATCGGGTTGGTCCCCAATGAAGTTTCTCCGTTGCTACATTGTTGTGAGCCTCGACCCAAATCCTCTCTATCAATTGGTTAATTATTTGTCAATTACCCCGATTGCTCTGTCAAATGGTCTTCTGCATTAGCTAGTTGTACTATCCCGAGCGCTACACCCGCACTTACGAGCAATGACTCCCTTTGCGGGCTGCAGCGCTGAGTCTAGTGAACGGCTGGGTGGACAAACGAATGTACGCCAATAGCTCTCCTCTTCTTAGTTCGTGTCGTATTCAAGGGTGGCCTTCACGATCACAGGCAAGGTAGCGTCCTTATGTAGTCACTACTTCGAGTTTGAATAGCTACTGTTAAATACCCATTGATCAATGCGTAATTCGACTACCGTTGTATTTACAATGTGAAACTGCCGGGAGGTGATAGGCCTAAATGGTTCATGTTAGGCTCGGTCGTAACTTTCCTCGCTAAAGATATGACTAAGTAAAGGACTATTCGCTTAGGGTGTGGTACACGCACATTAACTGCGCGTCGAGGGTAGCAAGAGTTGGATTACAAACGGTTCCTGGAGACGCTTTTGGGAACAAGCCGGCTAGAGCAACCTGGTAACAGAATCGGTGGTTTCCATAAAGAGCGAAGACGCCATTATTGTTGCTAAAGGCTGCGGCTCAAACTTGGATCCGTATGTATGAACTAGTAGCTCGCCGTAGATTGCCAGGTCGTGGGGCAAGGTGTCGATCCTACATTTTGGTAAGTCTTCGAGAGAAGGGGCATAATTATCGTGGCGTCTAGTTGCATTAGCACGATGCGGAGGCCGTACAGACGACGGAGAAGGGAATAAGGTTTCGTGGCTACACACCAGTCGGAAGGTAGTCGTTTCCCATGGTCAAGTTTGCCCTGCGTGATCAGTATTCCCCTATGCGAAAAGTGGAAACGTAGGTTCGGAAACCAGGGAGGCGGAGGTTGAAGACGGCAAGCACTGTATTTAAGAAAGACGACTCCGATAGATCCGCACCCAGTAAAACGTTATGCGCTCAGATGGCAGACTGGGTAGTGGTTGAGTGACATAGATCTATTCGAATTTGCGTCATCCTGAACCCTTAGACACCGATGCATGAAAACACGTATTACGGGCTATAACATGGAGACCCTTAGCAATCGGATCAAATAGCATTTGACTCGCTCTCTTATG"
# k <- 12
# t1 <- Sys.time()
# x <- deBruijinGraphAdjacencyList(txt, k)
# writeOverlap(x)
# t2 <- Sys.time()
# (t2-t1)

source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/overlapGraph.R"))

createNodeList <- function(txt, k) {
  
  k <- (k-1)
  nodeList <- c()
  numNodes <- (nchar(txt)-k+1)
  
  for (i in 1:numNodes) {
    
    nxtNode  <- substr(txt, i, (i+k-1))
    nodeList <- c(nodeList, nxtNode)
    
  }  # END for
  
  return(nodeList)
  
}  # END createNodeList

deBruijinMatrix <- function(nodeList) {
  
  uniqueNodeList  <- unique(nodeList)
  deBruijinMatrix <- matrix(0L, nrow = length(uniqueNodeList), ncol = length(uniqueNodeList), dimnames = list(uniqueNodeList, uniqueNodeList))
  
  return(deBruijinMatrix)
  
  
}  # END deBruijinMatrix

computeAdjancency <- function(nodeList, deBruijinMatrix) {
  
  lenNodeList <- (length(nodeList) - 1)
  
  for (idx in 1:lenNodeList) {
    
    r  <- nodeList[idx]
    cl <- nodeList[idx+1]
    
    deBruijinMatrix[r, cl] <- (deBruijinMatrix[r, cl]+1)
    
  }  # END for
  
  return(deBruijinMatrix)
  
}  # END computeAdjancency

# Precondition: dim(deBruijinMatrix) == (n x n)
createDeBruijinGraph <- function(deBruijinMatrix, sortList = FALSE) {
  
  kmers <- row.names(deBruijinMatrix)
  if (sortList) kmers <- sort(kmers)
  
  deBruijinGraph <- list()
  graphIdx <- 1
  nxtPath  <- NULL
  
  for (rKmer in kmers) {
    
    nxtPath <- new("Node")
    nxtPath@kmer <- rKmer
    
    j <- 1
    for (cKmer in kmers) {
      
      if (deBruijinMatrix[rKmer, cKmer] == 0)
        next()
      
      for (i in 1:deBruijinMatrix[rKmer, cKmer]) {
        
        nxtPath@adjacencyList[[j]] <- cKmer
        j <- (j+1)
        
      }  # END for
      
    }  # END for
    
    if (length(nxtPath@adjacencyList) > 0) {
      
      deBruijinGraph[[graphIdx]] <- nxtPath
      graphIdx <- (graphIdx+1)
      
    }  # END if
    
  }  # END for
  
  return(deBruijinGraph)
  
}  # END createDeBruijinGraph

deBruijinGraphAdjacencyList <- function(txt, k, sortList = FALSE) {
  
  nodeList <- createNodeList(txt, k)
  dbMatrix <- deBruijinMatrix(nodeList)
  dbMatrix <- computeAdjancency(nodeList, dbMatrix)
  dbGraph  <- createDeBruijinGraph(dbMatrix, sortList)
  
  return(dbGraph)
  
}  # END deBruijinGraphAdjacencyList
