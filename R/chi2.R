chi2fishercpp <- function(inputfilepath, outputfilepath, th){ 
  library(gdsfmt)
  library(Rcpp)
  #th <-  7
  #sourceCpp("/home/dan/Workspace/R/chiSquare_test/chiSquare_test/run_test.cpp")
  #inputfilepath <- "/home/dan/Workspace/R/1000G_gds/1000_genomes_FINGBR.gds"
  #outputfilepath <- "/home/dan/Workspace/R/KekteevaAngiraChi2/data"
  f <- openfn.gds(filename = inputfilepath, readonly = T)
  phenotype <-as.factor(as.character(read.gdsn(index.gdsn(f, "phenotype")) ) )
  phen_levels_number <- nlevels(phenotype)
  phenotype_length <- length(phenotype)
  
  phenotype <- as.numeric(as.factor(as.numeric(phenotype)))
  genotype <- read.gdsn(index.gdsn(f,"genotype"), c(1,1), c(-1,-1))
  genotypes_num_cols <- as.integer(NCOL(genotype))
  
  closefn.gds(f)
  start.time <- Sys.time()
  
  df <- run_chisquare_test(phenotype_length =as.integer( phenotype_length), phen_levels_number = as.integer(phen_levels_number), genotype_num_cols = as.integer(genotypes_num_cols), phenotype_vector_r = phenotype, genotype_matrix_r = genotype, threads = th)
  
  end.time <- Sys.time()
  filename <- "rcppOutput.csv"
  dir.create(outputfilepath, showWarnings = F)
  outcsv <- file.path(outputfilepath, filename)
  
  write.csv(df, file = outcsv)
  
  time.taken <- end.time - start.time
  time.taken
}
