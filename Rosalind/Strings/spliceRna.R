# Usage (for testing purtposes):
# fastaFName <- paste0(getwd(), "/Bioinformatics/Strings/0012 spliceRna/fasta.txt")
# fastaInput <- readInputFromFasta(fastaFName)
# mRna <- getmRna(fastaInput[[1]], fastaInput[[2]])
# mRna

source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Strings/0002 transcribeDNAintoRNA/transcribeDNAintoRNA.R")
source("C:/Users/heba/Desktop/HBC/Work/Projects/R/R/R-3.6.2/wd/Bioinformatics/Strings/0011 translateRnaIntoProtein/translateRnaIntoProtein.R")


# *********************************
# *********** LIBRARIES ***********
# *********************************


# *********************************
# *********** CONSTANTS ***********
# *********************************
SR_VERBOSE   <- TRUE
EMPTY_STRING <- ""


# *********************************
# *********** VARIABLES ***********
# *********************************


# *********************************
# *********** FUNCTIONS ***********
# *********************************

# ************** I/O **************
inputFromFasta <- function(fastaFName) {
  
  input <- read.fasta(file = fastaFName, as.string = TRUE)
  return(input)
  
}  # END inputFromFasta

getDnaFromFasta <- function(fastaObj) {
  
  return(toupper(fastaObj[[1]][1]))
  
}  # END getDnaFromFasta

getIntronsFromFasta <- function(fastaObj) {
  
  introns <- c()
  for (idx in 2:length(fastaObj)) {
    
    introns <- toupper(c(introns, fastaObj[[idx]][1]))
    
  }  # END for
  
  return(introns)
  
}  # END getIntronsFromFasta

readInputFromFasta <- function(fastaFName) {
  
  input <- inputFromFasta(fastaFName)
  dna   <- getDnaFromFasta(input)
  introns <- getIntronsFromFasta(input)
  
  return(list(dna, introns))
  
} # END readInputFromFasta


# ************** mRna **************
spliceDna <- function(dna, introns) {
  
  exons <- dna
  
  for (intron in introns)
    exons <- gsub(intron, EMPTY_STRING, exons, fixed = TRUE)
  
  return(exons)
  
}  # END spliceDna

getmRna <- function(dna, introns) {
  
  exons  <- spliceDna(dna, introns)
  rna  <- transcribeDnaToRna(exons)
  mRna <- computeRnaIntoProteinTranslation(rna)
  
  return(mRna)
  
}  # END getmRna
