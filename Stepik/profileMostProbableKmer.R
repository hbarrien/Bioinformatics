# Usage:
# dnaStrand <- "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
# k <- 5
# 0.2 0.2 0.3 0.2 0.3
# 0.4 0.3 0.1 0.5 0.1
# 0.3 0.3 0.5 0.2 0.4
# 0.1 0.2 0.1 0.1 0.2
# x <- profileMostProbableKmer(dnaStrand, k, NULL, "dataset.txt")
# x
#
# Sample Output:
# CCGAG


# *********************************
# *********** LIBRARIES ***********
# *********************************
library(Rcpp)


# *********************************
# *********** CONSTANTS ***********
# *********************************
DEFAULT_FILE_PATH <- "C:/Users/heba/Downloads/"
EMPTY_STRING      <- ""


# *********************************
# *********** VARIABLES ***********
# *********************************
PMFK_VERBOSE <- FALSE


# *********************************
# *********** FUNCTIONS ***********
# *********************************

# ************** I/O **************
readInput <- function(inputSource) {
  
  if (is.null(inputSource))
    return(NULL)
  
  input_filePath <- paste0(DEFAULT_FILE_PATH, inputSource)  
  
  input <- read.csv(input_filePath, sep=",", na.strings="", stringsAsFactors=FALSE)
  return(input)
  
}  # END readInput


# ******** DATA PROCESSING ********
# Returns: kmer
cppFunction('Rcpp::CharacterVector computeKmerProbabilitiesCpp(Rcpp::CharacterVector dnaStrand, int k, NumericMatrix profileMatrix) {

            int nKmers  = (dnaStrand.size()-k+1);
            int rIdx;
            double p;
            double bestp = 0.0;
            Rcpp::CharacterVector out(1);
            
            for (int i = 0; i < nKmers; i++) {
            
              int m = 0;
              Rcpp::CharacterVector nxtKmer(k);

              for (int j = i; j <= (i+k-1); j++) {

                nxtKmer[m] = dnaStrand[j];
                m++;

              }  // END for

              p = 1.0;
            
              for (int j = 0; j < k; j++) {
            
                Rcpp::String nxtLetter(nxtKmer[j]);

                if (nxtLetter == "A") {
                  rIdx = 0;
                } else {
                  if (nxtLetter == "C") {
                    rIdx = 1;
                  } else {
                    if (nxtLetter == "G") {
                      rIdx = 2;
                    } else {
                      rIdx = 3;
                    }
                  }
                }
                
                double nxtP = profileMatrix(rIdx, j);
                p = (p * nxtP);

              }  // END for
            
              if (p > bestp) {

                out = nxtKmer;
                bestp = p;

              }  // END if

            }  // END for
            
            return out;
            
        }')

# Returns: <kmer, pr>
cppFunction('Rcpp::List computeKmerProbabilitiesKmerPrCpp(Rcpp::CharacterVector dnaStrand, int k, NumericMatrix profileMatrix) {

            int nKmers  = (dnaStrand.size()-k+1);
            int rIdx;
            double p;
            double bestp = 0.0;
            Rcpp::CharacterVector out(1);
            
            for (int i = 0; i < nKmers; i++) {
            
              int m = 0;
              Rcpp::CharacterVector nxtKmer(k);
            
              for (int j = i; j <= (i+k-1); j++) {
            
                nxtKmer[m] = dnaStrand[j];
                m++;
            
              }  // END for
            
              p = 1.0;
            
              for (int j = 0; j < k; j++) {
            
                Rcpp::String nxtLetter(nxtKmer[j]);
            
                if (nxtLetter == "A") {
                  rIdx = 0;
                } else {
                  if (nxtLetter == "C") {
                    rIdx = 1;
                  } else {
                    if (nxtLetter == "G") {
                      rIdx = 2;
                    } else {
                      rIdx = 3;
                    }
                  }  
                }
            
                double nxtP = profileMatrix(rIdx, j);
                p = (p * nxtP);
            
              }  // END for
            
              if (p > bestp) {
            
                out = nxtKmer;
                bestp = p;
            
              }  // END if
            
            }  // END for
            
            return Rcpp::List::create(Rcpp::Named("kmer") = out, Rcpp::Named("pr") = bestp);
            
            }')

profileMostProbableKmer <- function(dnaStrand, k, profileMatrix, profileMatrixFName = NULL, includePr = FALSE) {
  
  if (is.null(dnaStrand) || (nchar(dnaStrand) == 0) ||
      is.null(k) || !is.numeric(k) || (k < 1))
    return(NULL)
  
  if (is.null(profileMatrix)) {
    
    if (is.null(profileMatrixFName) || (nchar(profileMatrixFName) == 0)) return(NULL)
    profileMatrix <- readInput(profileMatrixFName)
    
  }  # END if
  
  if (!includePr) {

    pr <- computeKmerProbabilitiesCpp(strsplit(dnaStrand, EMPTY_STRING)[[1]], k, as.matrix(profileMatrix))
    
    if (length(pr) == 0) return(NULL)
    
    res <- ""
    lapply(pr, function(x) { res <<- paste0(res, x) })
    
  } else {
    
    pr  <- computeKmerProbabilitiesKmerPrCpp(strsplit(dnaStrand, EMPTY_STRING)[[1]], k, as.matrix(profileMatrix))
    
    res <- ""
    lapply(pr[[1]], function(x) { res <<- paste0(res, x) })
    
    res <- list(kmer=res, pr=pr[[2]])
    
  }  # END if
  
  return(res)
  
}  # END profileMostProbableKmer


# ==================================================================
# ORIGINAL
#
# ### LIBRARIES ###
# library(data.table)
# library(compiler)
#
# ### CONSTANTS ####
# A <- 65
# C <- 67
# G <- 71
# T <- 84
# 
# ROW_IDX    <- c(integer(84))
# ROW_IDX[A] <- 1
# ROW_IDX[C] <- 2
# ROW_IDX[G] <- 3
# ROW_IDX[T] <- 4
#
# ### VARIABLES ###
# computeKmerProbabilitiesCmp <- NULL
#
# ### FUNCTIONS ###
# getRowIdx <- function(letter) {
#   
#   return(ROW_IDX[utf8ToInt(letter)])
#   
# }  # END getRowIdx
#
# computeKmerProbabilities <- function(dnaStrand, k, profileMatrix) {
# 
#   if(PMFK_VERBOSE) print("computeKmerProbabilities")
# 
#   lenText <- nchar(dnaStrand)
#   nKmers  <- (lenText-k+1)
#   out     <- data.table(kmer=character(nKmers), p=numeric(nKmers))
# 
#   for (i in 1:nKmers) {
# 
#     nxtKmer <- substr(dnaStrand, i, (i+k-1))
#     nxtKmerSplit <- strsplit(nxtKmer, EMPTY_STRING)[[1]]
#     p <- 1
# 
#     for (j in 1:k) {
# 
#       nxtLetter <- nxtKmerSplit[j]
#       nxtP      <- profileMatrix[getRowIdx(nxtLetter), j]
#       p         <- (p * nxtP)
# 
#       if (p == 0.0) break()
# 
#     }  # END for
# 
#     if(PMFK_VERBOSE) print(paste0("kmer: ", nxtKmer, ", p = ", p))
#     out[i,]$kmer <- nxtKmer
#     out[i,]$p <- p
# 
#   }  # END for
# 
#   return(out)
# 
# }  # END computeKmerProbabilities
#
# profileMostProbableKmer <- function(dnaStrand, k, profileMatrix, profileMatrixFName = NULL) {
#   
#   if (is.null(dnaStrand) || (nchar(dnaStrand) == 0) ||
#       is.null(k) || !is.numeric(k) || (k < 1))
#     return(NULL)
#   
#   if (is.null(profileMatrix)) {
#     
#     if (is.null(profileMatrixFName) || (nchar(profileMatrixFName) == 0)) return(NULL)
#     profileMatrix <- readInput(profileMatrixFName)
#     
#   }  # END if
#   
#   if (is.null(computeKmerProbabilitiesCmp))
#     computeKmerProbabilitiesCmp <<- cmpfun(computeKmerProbabilities)
#   
#   pr <- computeKmerProbabilitiesCmp(dnaStrand, k, profileMatrix)
#   if (nrow(pr) == 0) return(NULL)
#   
#   return(list(pr[pr$p == max(pr$p),]))
#   
# }  # END profileMostProbableKmer
# ==================================================================
