source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/PatternToNumber.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/FrequencyArray_v2.R"))

# Usage: t1 <- Sys.time(); fw <- patterns(genome, k, L, T); t2 <- Sys.time(); (t2-t1)
#        t1 <- Sys.time(); x <- findClump(genome, k, L, t); t2 <- Sys.time(); (t2 - t1)


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
FC_VERBOSE <- TRUE


# *********************************
# *********** VARIABLES ***********
# *********************************
rPath <- paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/dataset.txt")


# *********************************
# *********** FUNCTIONS ***********
# *********************************
patterns <- function(txt, k, L, t) {

  if (FC_VERBOSE) print("patterns: creating frequency array...")
  fa <- frequencyArray(k)
  
  lastWordStartPos <- (nchar(txt)-k+1)

  for (i in 1:lastWordStartPos) {
    
    if (FC_VERBOSE) print(paste0("patterns: ", i, " of ", lastWordStartPos))
    
    startPos <- i
    stopPos  <- (i+k-1)
    pattern  <- substr(txt, startPos, stopPos)
    code     <- (patternToNumber(pattern) + 1)
    
    if (fa[[FREQ_IDX]][code] == 0) {
      
      fa[[FREQ_IDX]][code]  <- 1
      fa[[START_IDX]][code] <- startPos
      fa[[STOP_IDX]][code]  <- stopPos
      next()
      
    }  # END if
    
    # At this point: (fa[[FREQ_IDX]][code] > 0)
    distance <- (stopPos - fa[[START_IDX]][code])
    if (distance <= L) {
      
      fa[[FREQ_IDX]][code] <- (fa[[FREQ_IDX]][code] + 1)
      
      if (fa[[FREQ_IDX]][code] == t)
        fa[[WINDOW_IDX]][code] <- (fa[[WINDOW_IDX]][code] + 1)
        
      next()
      
    }  # END if
    
    # At this point, (distance > L)
    # Start with this kmer as a new point of reference
    fa[[FREQ_IDX]][code]  <- 1
    fa[[START_IDX]][code] <- startPos
    fa[[STOP_IDX]][code]  <- stopPos
    
  }  # END for

  return(fa)
  
}  # END patterns

findClump <- function(genome, k, L, t) {
  
  # genome == NULL: when input is too large to be passed as a string
  if (is.null(genome))
    genome <- readChar(rPath, nchars = file.size(rPath))
  
  fa <- patterns(genome, k, L, t)
  
  result <- attributes(fa[[PATTERN_IDX]][fa[[WINDOW_IDX]] > 0])
  return(result[[1]])
  
}  # END findClump