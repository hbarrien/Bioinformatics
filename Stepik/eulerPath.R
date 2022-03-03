# Usage:
# t1 <- Sys.time()
# fPath <- "C:/Users/heba/Downloads/[file_name].txt"
# x <- eulerPath(fPath, formatOutput=TRUE)
# t2 <- Sys.time()
# (t2-t1)
#
# Author        : Herbert Barrientos
# Author contact: hpbarr@gmail.com
# Creation date : 2020-06-20
# Description   : Implementation of the Hierholzer algorithm for directed graphs
# Source        : https://www-m9.ma.tum.de/graph-algorithms/hierholzer/index_en.html
#
# Improvements (current state of this program is to satisfy a course's programming exercise):
#   a. include preconditions
#   b. include postconditions
#   c. include loop invariants
#   d. code optimization (for loops, tour construction, validation checks)
#   e. enhance documentation (program headers with author, version, scope, etc.; argument descriptions, inline comments)


# *********************************
# *********** LIBRARIES ***********
# *********************************
library(dplyr)


# *********************************
# *********** CONSTANTS ***********
# *********************************
ADJACENCY_LIST_POINTER <- " -> "     # Symbol used to read and process input
ADJACENCY_POINTER_OUT  <- "->"       # Symbol used to create output

NODE       <- "NODE"                 # Column name
IN_DEGREE  <- "IN"                   # Column name
OUT_DEGREE <- "OUT"                  # Column name
VISITED    <- "-1"                   # Indicator for visited edges

# Constants used to check for Euler Path and Circuit
TOO_MANY_IN_OR_OUT_EDGES     <- -40  # Graph has too many incoming or outgoing edges
MORE_THAN_ONE_START_NODE     <- -30  # Graph has more than one start node
MORE_THAN_ONE_END_NODE       <- -20  # Graph has more than one end node
GRAPH_NOT_STRONGLY_CONNECTED <- -10  # Graph has at least one disconncted node, or OUT/IN values are invalid
NO_EULER_PATH                <- 0    # Graph does NOT satisfy conditions for Euler Path (unknown reasons)
GRAPH_HAS_EULER_CIRCUIT      <- 1    # Graph has an Euler path which is also an Euler circuit (any node can be a start and an end node)
GRAPH_HAS_EULER_PATH         <- 2    # Graph has an Euler path (it has a start node and an end node)

COMMA <- ","


# *********************************
# *********** VARIABLES ***********
# *********************************
adjacencyList <- NULL                # Table representing a graph as an adjacency list
nodeDegrees   <- NULL                # Table containing the in-degrees and out-degrees of each node in adjacencyList


# *********************************
# *********** FUNCTIONS ***********
# *********************************

# Read input from text file: a graph represented as an adjacency list
# Input format per line: STRING -> STRING[, STRING,...,STRING]
#
# Returns:
#  Table with n columns - first column contains a node name; from the second column
#  on are the adjacent nodes. Columns not filled with node names contain NA
readGraphAdjacencyList <- function(fPath) {
  
  f <- file(fPath)
  adjList <- readLines(f)
  close(f)
  
  adjList <- gsub(ADJACENCY_LIST_POINTER, COMMA, adjList)
  adjList <- gsub(" +", "", adjList)
  adjList <- strsplit(adjList, COMMA)
  
  nodes <- c()
  nRowsAdjList   <- length(adjList)
  maxColsAdjList <- max(lengths(adjList))
  cNames <- c(NODE)
  
  # Ensure that all nodes are included in the NODE column
  mainNodes <- c()
  for (i in 1:nRowsAdjList)
    mainNodes <- c(mainNodes, adjList[[i]][1])
  
  missingFromMain <- c()
  for (i in 1:nRowsAdjList) {
    
    p <- which(!(adjList[[i]][-1] %in% mainNodes))
    if (length(p) > 0) missingFromMain <- c(missingFromMain, adjList[[i]][p+1])
    
  }  # END for
  
  if (length(missingFromMain) > 0) {
    
    lapply(missingFromMain, function(x) {adjList <<- c(adjList, x)})
    nRowsAdjList <- length(adjList)
    
  }  # END if
  
  for (i in 1:(maxColsAdjList-1))
    cNames <- c(cNames, paste0("V", i))
  
  adjListOut <- matrix(NA, nrow = nRowsAdjList, ncol = maxColsAdjList)
  colnames(adjListOut) <- cNames
  
  for (i in 1:nRowsAdjList) {
    
    adjListOut[i, NODE] <- adjList[[i]][1]
    
    for (j in 2:length(adjList[[i]]))
      adjListOut[i, cNames[j]] <- adjList[[i]][j]
    
  }  # END for
  
  adjListOut <- as.data.frame(adjListOut)
  adjListOut <- (adjListOut %>% mutate_all(as.character))
  
  return(adjListOut)
  
}  # END readGraphAdjacencyList


# Calculates the in-degree and out-degree of every node in global variable adjacencyList
#
# Returns:
#   Table of the form:
#     [row name] - a node name
#     IN         -  column having the number of incoming edges for node name (in-degree)
#     OUT        - column having the number of outgoing edges for node name (out-degree)
calculateNodeDegrees <- function() {
  
  nodes <- adjacencyList$NODE
  nodeDegrees <- as.data.frame(matrix(0L, nrow = nrow(adjacencyList), ncol = 2, dimnames = list(nodes, c(IN_DEGREE, OUT_DEGREE))))
  
  nRowsAdjacencyList <- nrow(adjacencyList)
  nColsAdjacencyList <- ncol(adjacencyList)
  
  for (i in 1:nRowsAdjacencyList) {
    
    outNode <- adjacencyList[i,]$NODE
    
    for (j in 1:(nColsAdjacencyList-1)) {
      
      cName  <- paste0("V", j)
      inNode <- adjacencyList[i, cName]
      
      if(is.na(inNode)) break()
      
      nodeDegrees[outNode, OUT_DEGREE] <- (nodeDegrees[outNode, OUT_DEGREE]+1)
      nodeDegrees[inNode, IN_DEGREE]   <- (nodeDegrees[inNode, IN_DEGREE]+1)
      
    }  # END for
    
  }  # END for
  
  return(nodeDegrees)
  
}  # END calculateNodeDegrees


# Determines if the graph in global variable nodeDegrees is balanced or not
#
# Returns:
#   TRUE  - the graph is balanced
#   FALSE -  otherwise
isGraphBalanced <- function(nodeDegrees) {
  
  return(all((nodeDegrees$IN > 0) & (nodeDegrees$OUT > 0) & (nodeDegrees$IN == nodeDegrees$OUT)))
  
}  # END isGraphBalanced


# Determines whether or not the graph in global variable nodeDegrees contains an Euler path or an Euler circuit
#
# Returns (see value descriptions in the CONSTANTS section):
#   GRAPH_NOT_STRONGLY_CONNECTED - graph does not have an Euler path or circuit
#   TOO_MANY_IN_OR_OUT_EDGES     - graph does not have an Euler path or circuit
#   MORE_THAN_ONE_START_NODE     - graph does not have an Euler path or circuit
#   MORE_THAN_ONE_END_NODE       - graph does not have an Euler path or circuit
#   NO_EULER_PATH                - graph does not have an Euler path or circuit
#   GRAPH_HAS_EULER_CIRCUIT      - graph has an Euler circuit
#   GRAPH_HAS_EULER_PATH         - graph has an Euler circuit
graphHasEulerPath <- function() {
  
  startNode <- 0
  endNode   <- 0
 
  nRowsgraphDegree <- nrow(nodeDegrees)
  
  for (i in 1:nRowsgraphDegree) {
    
    if ((nodeDegrees[i,]$OUT < 1) && (nodeDegrees[i,]$IN < 1)) {
      
      return(GRAPH_NOT_STRONGLY_CONNECTED)
      
    } else if (((nodeDegrees[i,]$OUT - nodeDegrees[i,]$IN) > 1) ||
               ((nodeDegrees[i,]$IN - nodeDegrees[i,]$OUT) > 1)) {
      
      return((TOO_MANY_IN_OR_OUT_EDGES))
      
    } else if ((nodeDegrees[i,]$OUT - nodeDegrees[i,]$IN) == 1) {
      
      startNode <- (startNode+1)
      if (startNode > 1) return(MORE_THAN_ONE_START_NODE) 
      
    } else if ((nodeDegrees[i,]$IN - nodeDegrees[i,]$OUT) == 1) { 
      
      endNode <- (endNode+1)
      if (endNode > 1) return(MORE_THAN_ONE_END_NODE)
    } 
    
  }  # END for
  
  if ((startNode == 0) && (endNode == 0)) return(GRAPH_HAS_EULER_CIRCUIT)
  if ((startNode == 1) && (endNode == 1)) return(GRAPH_HAS_EULER_PATH)
  
  # Unknown reasons...
  return(NO_EULER_PATH)
  
}  # END graphHasEulerPath


# Finds the starting node for graph traversal
# Precondition: nodeDegrees represents an Euler Path or an Euler Circuit
#
# Returns:
#   A start node name
findStartNode <- function() {
  
  if (graphHasEulerPath() == GRAPH_HAS_EULER_PATH)
    return(row.names(nodeDegrees[(nodeDegrees$OUT - nodeDegrees$IN) == 1,]))
  
  # At this point, graph MUST be an Euler Circuit. Choose the first node
  # with one entry edge and one exit edge
  startNode <- row.names(nodeDegrees[((nodeDegrees$OUT == 1) & (nodeDegrees$IN == 1)),])[1]
  
  # If not found, select a node with the minimum number of edges
  if (is.na(startNode) || is.null(startNode))
    startNode <- row.names(nodeDegrees[min(nodeDegrees$OUT) == min(nodeDegrees$IN),])[1]
  
  return(startNode)
  
}  # END findStartNode


# Takes the first node, in argument tour, having at least one unvisited edge
#
# Returns: 
#   A node with at least one unvisited edge
#   NULL if all nodes in tour have no visited edges
nodeWihUnvisitedEdge <- function(tour) {
  
  for (node in tour) {
    
    nodes <- as.list(adjacencyList[adjacencyList$NODE == node,])
    
    # First element is node. From the second on are adjacent nodes
    for (i in 2:length(nodes)) {
      
      if (is.na(nodes[[i]]))
        break()
      
      if (nodes[[i]] != VISITED)
        return(nodes[[NODE]])
      
    }  # END for
    
  }  # END for
  
  return(NULL)
  
}  # END nodeWihUnvisitedEdge


# Selects an unvisited edge, represented by currNode->u, where u is a node that
# has not been marked as visited, and updates global variables adjacencyList
# and nodeDegrees
#
# Returns: 
#   u    - the node on the other side of the taken edge
#   NULL - no unvisited edge was found
takeUnvisitedEdge <- function(currNode) {
  
  r <- adjacencyList[adjacencyList$NODE == currNode,]
  
  for (i in 2:ncol(r)) {
    
    u <- r[[i]]
    
    if (!is.na(u) && (u != VISITED)) {
      
      r[[i]] <- VISITED
      
      adjacencyList[adjacencyList$NODE == currNode,] <<- r
      nodeDegrees[currNode,]$OUT <<- (nodeDegrees[currNode,]$OUT - 1)
      
      return(u)
      
    }  # END if
    
  }  # END for
  
  return(NULL)
  
}  # END takeUnvisitedEdge


# Updates argument tour by replacing the first occurrence of matchNode
# with subTour
#
# Returns:
#   Updated version of argument tour
integrateTours <- function(stuck, matchNode, subTour, tour)  {
  
  newTour <- c()
  oneTime <- FALSE
  
  if (stuck) {
    
    if (subTour[1] == tour[length(tour)]) return(c(tour, subTour[2:length(subTour)]))
    else stop(paste0("Error. No connecting element found for matchNode: ", matchNode))
    
  }  # END if
  
  for (node in tour) {
    
    if ((node != matchNode) || oneTime) {
      
      newTour <- c(newTour, node)
      next()
      
    }  # END if
    
    newTour <- c(newTour, subTour)
    oneTime <- TRUE
    
  }  # END for
  
  return(newTour)
  
}  # END integrateTours


# Determines whether or not an Euler path or circuit has already been found
#
# Returns:
#   TRUE  - an Euler path or circuit has been found, i.e., the process is complete
#   FALSE - otherwise
isEulerPathComplete <- function() {
  
  return(all(nodeDegrees$OUT == 0))
  
}  # END isEulerPathComplete


# Calculates an Euler path or an Euler circuit
#
# Precondition: global variables adjacenyList and nodeDegrees contain the representation
# of a graph. Furthermore, both global variables represent an Euler path or an Euler circuit. 
# The process stops execution if an error is found
#
# Returns
#   A vector containing the graph's traversal ordering for an Euler path or an Euler circuit
computeEulerPath <- function() {
  
  startNode <- findStartNode()
  tour      <- c(startNode)
  subTour   <- c()
  
  repeat {
    
    currNode  <- nodeWihUnvisitedEdge(tour)
    startNode <- currNode
    
    if (is.null(currNode))
      stop("Error: No node with unvisited edge.")
    
    subTour <- c(subTour, startNode)
    
    stuck <- FALSE
    repeat{
      
      u <- takeUnvisitedEdge(currNode)
      
      # Path traversal got stuck
      if (is.null(u)) {
        
        stuck <- TRUE
        break()
        
      }  # END if
      
      subTour  <- c(subTour, u)
      currNode <- u
      
      if (startNode == currNode)
        break()
      
    }  # END repeat
    
    tour    <- integrateTours(stuck, startNode, subTour, tour)
    subTour <- c()
    
    if (isEulerPathComplete())
      break()
    
  }  # END repeat
  
  return(tour)
  
}  # END computeEulerPath


# Formats output
#
# Returns:
#   A string of the form: "NODE->NODE->...->NODE"
formatEulerPath <- function(ePath) {
  
  ePathFormatted <- ""
  
  for (node in ePath)
    ePathFormatted <- ifelse((ePathFormatted == ""), 
                             paste0(ePathFormatted, node), 
                             paste0(ePathFormatted, ADJACENCY_POINTER_OUT, node))
  
  return(ePathFormatted)
  
}  # END formatEulerPath


# Main function
#
# Returns
#   Value of variable g, if the graph in memory does not represent an Euler path or circuit
#   The unformatted output returned by computeEulerPath()
#   Formatted output
eulerPath <- function(fPath, formatOutput = FALSE) {
  
  adjacencyList <<- readGraphAdjacencyList(fPath)
  nodeDegrees   <<- calculateNodeDegrees()
  
  g <- graphHasEulerPath()
  
  if ((g != GRAPH_HAS_EULER_CIRCUIT) && (g != GRAPH_HAS_EULER_PATH))
    return(g)
  
  ePath <- computeEulerPath()
  
  if (formatOutput)
    ePath <- formatEulerPath(ePath)
  
  return(ePath)
  
}  # END eulerPath
