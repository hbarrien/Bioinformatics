# Problem
# The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. 
# For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string 
# has the same GC-content.
# 
# DNA strings must be labeled when they are consolidated into a database. A commonly used method of string 
# labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', 
# followed by some labeling information. Subsequent lines contain the string itself; the first line to begin 
# with '>' indicates the label of the next string.
#
# In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where 
# "xxxx" denotes a four-digit code between 0000 and 9999.
# 
# Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
# Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. 
# Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see 
# the note on absolute error below.
#
# Note on Absolute Errorclick to collapse
# We say that a number x is within an absolute error of y to a correct solution if x is within y of the 
# correct solution. For example, if an exact solution is 6.157892, then for x to be within an absolute error 
# of 0.001, we must have that |x−6.157892|<0.001, or 6.156892<x<6.158892.
#
# Error bounding is a vital practical tool because of the inherent round-off error in representing decimals 
# in a computer, where only a finite number of decimal places are allotted to any number. After being compounded 
# over a number of operations, this round-off error can become evident. As a result, rather than testing whether 
# two numbers are equal with x=z, you may wish to simply verify that |x−z| is very small.
#
# The mathematical field of numerical analysis is devoted to rigorously studying the nature of computational approximation.

library(data.table)

# file path to the data input file
fpath <- "C:/Users/heba/Desktop/Bioinformatics/01 Strings/0005 identifyDNA_GC-Content/rosalind_gc.txt"

# numeric values for each letter (each representing a nucleobase)
A <- 65
C <- 67
G <- 71
T <- 84

# variables holding extracted data from the input file
data_label <- ""
dna_data   <- ""

# data table containing processed dna strings
rosalind_data <- data.table(data_label=character(), dna_data=character(), count_A=integer(), count_C=integer(), 
                            count_G=integer(), count_T=integer(), GC_content=double())

# input data file
text_data <- read.csv(fpath, header = FALSE)

for (i in 1:nrow(text_data)) {
  
  input_row <- as.character(text_data[i,1])
  
  # Follow FASTA format to obtain complete input data for one output row
  if (grepl(">", input_row)) {
    data_label <- gsub(">", "", input_row)
    next()
  } else {
    dna_data<- input_row
  }
  
  # convert dna character string to a numeric vector
  dna_data_numeric <- utf8ToInt(as.character(dna_data))

  # count number of instances of each letter and store results in countDNA  
  countDNA <- integer(84)
  invisible(sapply(dna_data_numeric, function(x) {countDNA[x] <<- (countDNA[x] + 1)}))
  
  # calculate percentage of GC-content
  gc_content <- ((countDNA[G] + countDNA[C]) / length(dna_data_numeric)) * 100
  
  # create a new row with all the results and add it to the output table
  rosalind_data <- rbind(rosalind_data, list(data_label, dna_data, countDNA[A], countDNA[C], countDNA[G], countDNA[T], gc_content))
  
}  # END for

# obtain the dna strand with the highes GC-content
highest_gc_content <- 1
for (i in 2: nrow(rosalind_data)) {
  
  if (rosalind_data[i,]$GC_content > rosalind_data[highest_gc_content,]$GC_content)
    highest_gc_content <- i
}

# Display the final result
cat(paste(rosalind_data[highest_gc_content,]$data_label, round(rosalind_data[highest_gc_content,]$GC_content, digits=6), sep="\n"))
