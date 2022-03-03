# Problem
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English 
# alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from 
# these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along 
# with DNA strings and RNA strings.
#
# The RNA codon table dictates the details regarding the encoding of specific codons into the 
# amino acid alphabet.
# 
# Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
# 
# Return: The protein string encoded by s.


library(data.table)


# ################################
# ########## CONSTANTS ###########
# ################################
RNA_CODON_LENGTH     <- 3
TRANSLATION_STOP_UAA <- "UAA"
TRANSLATION_STOP_UAG <- "UAG"
TRANSLATION_STOP_UGA <- "UGA"
TRANSLATION_STOP     <- paste0("(", TRANSLATION_STOP_UAA, "|", TRANSLATION_STOP_UAG, "|",TRANSLATION_STOP_UGA, ")")


# ################################
# ########## VARIABLES ###########
# ################################
translation_table <- NULL



# ################################
# #### TRANSLATION TABLE FILE ####
# ################################
fPath <- paste0(getwd(), "/Bioinformatics/Strings/0011 translateRnaIntoProtein/")
fName <- paste0(fPath, "translation_table.txt")
translation_table <- as.data.table(read.csv(fName, stringsAsFactors=FALSE))


# ################################
# ########## FUNCTIONS ###########
# ################################
computeRnaIntoProteinTranslation <- function(rna) {

  protein <- ""
  
  for(i in seq(from=1, to=nchar(rna), by=RNA_CODON_LENGTH)) {
    
    codon <- substr(rna, i, (i+RNA_CODON_LENGTH-1))
    
    if (!length(translation_table[translation_table$rna.codon == codon,]$amino.acid))
      stop(paste0("Unrecognized rna codon: ", codon))
    
    if (grepl(TRANSLATION_STOP, codon))
      break()
  
    protein <- paste0(protein, translation_table[translation_table$rna.codon == codon,]$amino.acid)
  
  }  # END for
  
  return(protein)
  
}  # END computeRnaIntoProteinTranslation

translateRnaIntoProtein <- function(rnaFileName) {
  
  if (is.null(rnaFileName) || (rnaFileName == ""))
    stop("Bad input.")

  fName <- paste0(fPath, rnaFileName)
  rna <- readChar(fName, file.info(fName)$size)
  
  return(computeRnaIntoProteinTranslation(rna))
  
}  # END translateRnaIntoProtein
