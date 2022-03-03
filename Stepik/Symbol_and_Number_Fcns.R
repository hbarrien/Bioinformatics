SYMBOLS             <- c("A", "C", "G", "T")
SYMBOL_CODES        <- c(0, 1, 2, 3)

SYMBOL_TO_NUMBER <- SYMBOL_CODES
names(SYMBOL_TO_NUMBER) <- SYMBOLS

NUMBER_TO_SYMBOL <- SYMBOLS
names(NUMBER_TO_SYMBOL) <- SYMBOL_CODES


patternPrefix <- function(pattern) {
  
  p <- substr(pattern, 1, (nchar(pattern) - 1))
  
  return(p)
  
}  # END patternPrefix

lastSymbol <- function(pattern) {
  
  lenPattern <- nchar(pattern)
  symbol <- substr(pattern, lenPattern, lenPattern)
  
  return(symbol)
  
}  # END lastSymbol

symbolToNumber <- function(symbol) {
  
  return(SYMBOL_TO_NUMBER[symbol][[1]])
  
}  # END symbolToNumber

numberToSymbol <- function(index) {
  
  return(NUMBER_TO_SYMBOL[index+1][[1]])
  
}  # END numberToSymbol

patternToNumber <- function(pattern) {
  
  if (nchar(pattern) == 0) return(0)
  
  symbol <- lastSymbol(pattern)
  prefix <- patternPrefix(pattern)
  return(4 * patternToNumber(prefix) + symbolToNumber(symbol))
  
}  # END patternToNumber

numberToPattern <- function(index, k) {

  if (k == 1) return(numberToSymbol(index))
  
  prefixIndex <- trunc(index / 4)
  r <- (index %% 4)
  
  symbol <- numberToSymbol(r)
  prefixPattern <- numberToPattern(prefixIndex, (k - 1))
  return(paste0(prefixPattern, symbol))
  
}  # END numberToPattern

