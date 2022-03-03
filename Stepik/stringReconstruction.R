# Usage:
# t1 <- Sys.time()
# patterns <- read.csv("C:/Users/heba/Downloads/dataset_203_7.txt", sep="", na.strings="", stringsAsFactors=FALSE)$kmer
# x <- stringReconstruction(patterns)
# t2 <- Sys.time()
# (t2-t1)


# ===============
# GENOME ASSEMBLY
# ===============

source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/deBruijinGraphFromKmerComposition.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/eulerPath.R"))
source(paste0(getwd(), "/Bioinformatics/Coursera-Stepik/Code challenges/pathToGenome.R"))

stringReconstruction <- function(patterns) {
  
  print("de Bruijn...")
  db <- deBruijinGraphAdjacencyListFromPatterns(patterns)
  writeOverlap(db)
  
  print("Creating Euler path...")
  fpath <- "C:/Users/heba/Downloads/overlap_graph_1_out.txt"
  path  <- eulerPath(fpath)
  file.remove(fpath)
  
  print("Converting Euler path to genome string...")
  txt <- pathToGenome(path)
  
  return(txt)
  
}  # END stringReconstruction
