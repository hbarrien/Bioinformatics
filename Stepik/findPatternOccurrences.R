# PP. 13
# Usage: p <- "CTTGATCAT"; input <- readChar("C:/Users/heba/Downloads/dataset.txt", nchars = file.size("C:/Users/heba/Downloads/dataset.txt")-2 ); start <- Sys.time(); x <- findPatternOccurrences(p, input, TRUE); stop <- Sys.time(); (stop-start)

library(data.table)

computePatternOccurrences <- function(pattern, genome) {

  lenPattern <- nchar(pattern)

  firstPatternNucleotide <- utf8ToInt(pattern)[1]
  numericGenome <- utf8ToInt(genome)

  t <- data.table(nucleotide = integer(), idx = integer())
  t <- rbind(t, list(numericGenome, seq(length(numericGenome))))
  t <- t[t$nucleotide == firstPatternNucleotide,]
  t$s <- list(lapply(t$idx, function(x) { return(substr(genome, x, ((x+lenPattern-1))))} ))

  return(t[t$s == pattern,]$idx)

}  # END computePatternOccurrences

findPatternOccurrences <- function(pattern, genome, zeroBasedIndexing) {
  
  patternOccurrences <- computePatternOccurrences(pattern, genome)
  if (zeroBasedIndexing) return(as.integer(lapply(patternOccurrences, function(x) { return(x-1)})))
  
  return(patternOccurrences)
  
}  # END findPatternOccurrences
