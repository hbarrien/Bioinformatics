
# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
VERBOSE <- "TRUE"
BASE    <- 4

b4letters <- c("A", "C", "G", "T")
names(b4letters) <- c(0, 1, 2, 3)


# *********************************
# *********** FUNCTIONS ***********
# *********************************

# IMPORTANT NOTE (N1): the %% operator will return incorrect values 
# for large values of 'index'. The gmp library will help transforming 
# a value in E notation to decimal. For instance: x <- as.bigz(y), where 
# y is a large int. To date (2020-06-06), I have not found a suitable
# replacement for %%, when the numerator is large...
# library(gmp)
numberToPattern <- function(index, k) {
  
  # N1: index <- as.bigz(index)
  intPart <- 0
  res <- c()
  
  for (i in 1:k) {
    
    # N1: print(class(index))
    intPart <- trunc(index/BASE)
    res     <- c(res, (index %% BASE))
    index   <- intPart
    
  }  # END for
  
  # N1: res <- as.integer(res)
  # N1: print(res)
  
  res <- rev(res)
  outStr <- as.character(lapply(res, function(x) { return(b4letters[x+1][[1]]) } ))
  
  return(paste0(outStr, collapse = ""))
  
}  # END numberToPattern
