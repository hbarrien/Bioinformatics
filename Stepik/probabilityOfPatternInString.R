mChoosek <- function(m, k) {
  
  return(lfactorial(m) / (lfactorial(k) * lfactorial(m - k)))
  
}  # END mChoosek

Pr <- function(N, A, k, t) {
  
  return(mChoosek(((N-t) * (k-1)), t) / '^'(A, (t*k)))
  
}  # END Pr
