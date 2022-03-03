
setClass("Node", slots=list(kmer="character", adjacencyList="list"))

COMMA <- ","
SPACE <- " "
ADJACENT_TO <- "->"

# Precondition: ((pattern is not null) && (|pattern| > 1))
prefix <- function(pattern) {
  
  return(substring(pattern, 1, (nchar(pattern)-1)))
  
}  # END prefix

# Precondition: ((pattern is not null) && (|pattern| > 1))
suffix <- function(pattern) {
  
  return(substring(pattern, 2, nchar(pattern)))
  
}  # END suffix

printOverlap <- function(adjacencyGraph) {
  
  numNodes <- length(adjacencyGraph)
  
  for (i in 1:numNodes) {
    
    node <- adjacencyGraph[[i]]
    lenAdjacenyList <- length(node@adjacencyList)
    
    adjList <- ""
    for (j in 1:lenAdjacenyList) {
      
      adjList <- paste0(adjList, node@adjacencyList[[j]])
      if (j < lenAdjacenyList) adjList <- paste0(adjList, COMMA, SPACE)
      
    }  # END for
    
    print(paste0(node@kmer, SPACE, ADJACENT_TO, SPACE, adjList))
    
  }  # END for
  
}  # END printOverlap

writeOverlap <- function(adjacencyGraph) {
  
  f <- file("C:/Users/heba/Downloads/overlap_graph_1_out.txt", "a")
  
  numNodes <- length(adjacencyGraph)
  
  for (i in 1:numNodes) {
    
    node <- adjacencyGraph[[i]]
    lenAdjacenyList <- length(node@adjacencyList)
    
    adjList <- ""
    for (j in 1:lenAdjacenyList) {
      
      adjList <- paste0(adjList, node@adjacencyList[[j]])
      if (j < lenAdjacenyList) adjList <- paste0(adjList, COMMA, SPACE)
      
    }  # END for
    
    writeLines(paste0(node@kmer, SPACE, ADJACENT_TO, SPACE, adjList), f)
    
  }  # END for
  
  close(f)
  
}  # END writeOverlap

# Precondition: (patterns is not null) && (|patterns| > 0)
overlap <- function(patterns) {
  
  numPatterns    <- length(patterns)
  adjacencyGraph <- c()
  numSymbols     <- nchar(patterns[1])
  
  # Loop invariant: (|patterns[idx]| == numSymbols)
  for (i in 1:numPatterns) {
    
    print(i)
    
    # Check invariant
    if (nchar(patterns[i]) != numSymbols)
      stop(paste0("Error: kmer ", patterns[i], " has different length."))
    
    node <- new("Node")
    node@kmer  <- patterns[i]
    nodeSuffix <- suffix(node@kmer)
    
    j <- 1
    for (k in 1:numPatterns) {
      
      evalKmer <- patterns[k]
      
      if (nodeSuffix == prefix(evalKmer)) {
        node@adjacencyList[[j]] <- evalKmer
        j <- (j+1)
      }
      
    }  # END for
    
    if (length(node@adjacencyList) > 0)
      adjacencyGraph <- c(adjacencyGraph, node)
    
  }  # END for
  
  return(adjacencyGraph)
  
}  # END overlap