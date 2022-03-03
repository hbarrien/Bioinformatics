# Problem
# A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences 
# can be finite or infinite. Two examples are the finite sequence (pi,sqrt(2), 0, pi) and the infinite 
# sequence of odd numbers (1, 3, 5, 7, 9, ...). We use the notation an to represent the nth term of a sequence.
#
# A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. 
# In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive 
# the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal 
# to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs 
# alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence 
# relation Fn=F(n-1)+F(n-2) (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it 
# was known to Indian mathematicians over two millennia ago.
# 
# When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation 
# to generate terms for progressively larger values of n. This problem introduces us to the computational technique of 
# dynamic programming, which successively builds up solutions by using the answers to smaller cases.
# 
# Given: Positive integers n <= 40 and k <= 5.
# Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each 
# generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
# Sample Dataset
# 5 3
# Sample Output
# 19

# memoi: cache for previous calculations
memoi <- c(1, 1)

# last_k: keep track of the last k value input by the user. The purpose of this value is that, if the new k entered 
# is different from the last computed k, then memoi has to be reset to its initial values
last_k <- -1

# calcRabbits
#
# Scope: private
calcRabbits <- function(n, d, k) {
  
  if (n <= length(memoi))
    return(memoi[n])
  
  if (d == (length(memoi) + 1)) {
  
    newRabbits <- (memoi[d-1] + (k * memoi[d-2]))
    memoi <<- c(memoi, newRabbits)
    
    if (n == d)
      return(memoi[n])
    else
      return(calcRabbits(n, n, k))
  }
  
  return(calcRabbits(n, (d-1), k))
  
}  # END calcRabbits

# numRabbits
# 
# Scope: public
#
# Precondition: (n is int) & (k is int) & (n > 0) & (k > 0)
# 
# Returns: -1: error
#          >0: computed number of rabbits
numRabbits <- function(n, k) {
  
  # check precondition 1: must use is.double() to check for parameter type, as 1, 1.0, -1.0 are treated as doubles
  if (!(is.double(n) && is.double(k)))
    return(-1)
  
  # check precondition 2: ensure n and k are integers
  if (((n - trunc(n)) > 0) || ((k - trunc(k)) > 0))
    return(-1)
  
  # check precondition 3: ensure (n > 0) && (k > 0)
  if ((n < 1) || (k < 1))
    return(-1)
  
  # reset memoi if parameter k is different from the last computed k, set new value for last_k
  if (k != last_k) {
    print(paste0("changing... ", "last_k: ", last_k, ". k: ", k))
    memoi  <<- c(1, 1)
    last_k <<- k
  }
  
  # d: number that starts equal to n, and decreases by 1 with each recursive 
  # call when: (d != (length(memoi) + 1)). Used as index for memoi
  d <- n
  return(calcRabbits(n, d, k))
  
}  # END numRabbits
