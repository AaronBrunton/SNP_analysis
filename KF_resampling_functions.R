

library(dartR)
library(tidyverse)
library(reshape2)
library(cowplot)

# dat should be a character string to an rds object, nreps is the number of replicates to run.. 

SNP_optimization <- function(dat, nreps){
# Read in data
  data <- readRDS(dat)
  if(methods::is(data,"genind")){
    # Convert to genlight
    dat <- data
    print("Data is in genind format, moving onto resampling procedure")
  } else if(methods::is(data, "genlight")){
    print("Data is in genlight format, converting to genind before resampling procedure")
    dat <- dartR::gl2gi(data)
  }

# Record the number of SNPs, round down to nearest 100 so we don't throw an error in our loop 
nSNPs <- ncol(dat@tab)
nSNPs <- floor(nSNPs/100)*100

# Convert back to genlight for empirical calculations
Emp_gl <- dartR::gi2gl(dat)

### Calculate empirical statistics 
Emp_het <- dartR::gl.report.heterozygosity(Emp_gl, plot.out = F)

# Extract columsn we are intersted in by name
cols <- c("Ho", "He", "uHe", "FIS")
Emp_het <- Emp_het[,which(colnames(Emp_het) %in% cols)]

# How many boostraps do we want, technically its not bootstrapping because we are not creating a dataset the same size as our original, so technically its just resampling with replacement
boots <- nreps
samp_size <- seq(from = 100, to= nSNPs, by = 100)
Sampled_loc <- list() # Make a list to store the locus names in each boostrap sample, maybe also do some kind of seed??
Boot_res <- list()


# Just test against different SNP thresholds, includes all individuals 
for (j in 1:length(samp_size)) {
 # Set progress bar 
pb <- txtProgressBar(min = 0, max = (length(samp_size)-1), initial = 0, style = 3) 
  Tmp_res <- list()
  Tmp_loc <- list()
for(i in 1:boots){
  sampled_data <-dat[,sample(nLoc(dat), size = samp_size[j], replace = T)]
  Tmp_loc[[i]] <- colnames(sampled_data@tab) # Get locus names considered in the sampling replicate
  names(Tmp_loc)[[i]] <- paste("Sample", i, "loci", samp_size[j], sep = '_')
  # Convert for summary stat estimation
  gl_mod <- dartR::gi2gl(sampled_data, verbose = 0)
  # Estimate stats
  Sub_het <- dartR::gl.report.heterozygosity(gl_mod, plot.out = F)
  # Extract columsn we are intersted in by name
  cols <- c("Ho", "He", "uHe", "FIS")
  Sub_het <- Sub_het[,which(colnames(Sub_het) %in% cols)]
  rownames(Sub_het) <- paste("Samp", i, "est", samp_size[j], sep = '_')
  
  # Calculate error (actual stat - subsample estimate)
  Error <- Sub_het - Emp_het
  rownames(Error) <- paste("Samp", i, "error", samp_size[j], sep = '_')
  MSE <- Error^2
  rownames(MSE) <- paste("Samp", i, "MSE", samp_size[j], sep = '_')
  Sub_est <- rbind(Sub_het, Error, MSE)
  Sub_est$Stat <- c("Est", "Error", "MSE")
  Sub_est$nSNP <- samp_size[j]
  
  # Put into our list of results
  Tmp_res[[i]] <- Sub_est
  # Get rid of anything just to be safe 
  remove(cols, Sub_het, Error, MSE, Sub_est, gl_mod, sampled_data)
  }
  Boot_res[[j]] <- Tmp_res
  names(Boot_res)[j] <- samp_size[j]
  Sampled_loc[[j]] <- Tmp_loc
  names(Sampled_loc)[j] <- samp_size[j]
  remove(Tmp_res, Tmp_loc)
  setTxtProgressBar(pb,j)
}
Final_res_ind <- list(Boot_res, Sampled_loc, Emp_het)
names(Final_res_ind) <- c("SNP_Resampling_Results", "Sampled_loci", "Empirical_Estimates")
return(Final_res_ind)
}



# Test different SNP and individual combination
# Make a list to store the locus names in each boostrap sample, maybe also do some kind of seed??

Ind_optimization <- function(dat, nreps){

Sampled_Ind <- list()
Boot_res_ind <- list()
Ind_SNP_res <- list()

# Read in data
data <- readRDS(dat)
if(methods::is(data,"genind")){
  # Convert to genlight
  dat <- data
  print("Data is in genind format, moving onto resampling procedure")
} else if(methods::is(data, "genlight")){
  print("Data is in genlight format, converting to genind before resampling procedure")
  dat <- dartR::gl2gi(data)
}

# Convert back to genlight for empirical calculations
Emp_gl <- dartR::gi2gl(dat)


### Calculate empirical statistics 
Emp_het <- dartR::gl.report.heterozygosity(Emp_gl, plot.out = F)

# Extract columsn we are intersted in by name
cols <- c("Ho", "He", "uHe", "FIS")
Emp_het <- Emp_het[,which(colnames(Emp_het) %in% cols)]

# Record the number of individuals, round down to nearest 5 so we don't throw an error in our loop 
nInds <- nrow(dat@tab)
nInds <- floor(nInds/5)*5

Ind_samp_size <- seq(from = 5, to= nInds, by = 5)

boots <- nreps

# Just test against different individual thresholds, includes all SNPs
for (j in 1:length(Ind_samp_size)) {
  Tmp_res <- list()
  Tmp_ind <- list()
  pb <- txtProgressBar(min = 0, max = (length(Ind_samp_size)-1), initial = 0, style = 3) 
  for(i in 1:boots){
    sampled_data <-dat[sample(nInd(dat), size = Ind_samp_size[j], replace = T),]
    Tmp_ind[[i]] <- rownames(sampled_data@tab) # Get locus names considered in the sampling replicate
    names(Tmp_ind)[[i]] <- paste("Sample", i, "ind", Ind_samp_size[j], sep = '_')
    # Convert for summary stat estimation
    gl_mod <- dartR::gi2gl(sampled_data, verbose = 0)
    # Estimate stats
    Sub_het <- dartR::gl.report.heterozygosity(gl_mod, plot.out = F)
    # Extract columsn we are intersted in by name
    cols <- c("Ho", "He", "uHe", "FIS")
    Sub_het <- Sub_het[,which(colnames(Sub_het) %in% cols)]
    rownames(Sub_het) <- paste("Samp", i, "est", Ind_samp_size[j], sep = '_')
    
    # Calculate error (actual stat - subsample estimate)
    Error <- Sub_het - Emp_het
    rownames(Error) <- paste("Samp", i, "error", Ind_samp_size[j], sep = '_')
    MSE <- Error^2
    rownames(MSE) <- paste("Samp", i, "MSE", Ind_samp_size[j], sep = '_')
    Sub_est <- rbind(Sub_het, Error, MSE)
    Sub_est$Stat <- c("Est", "Error", "MSE")
    Sub_est$nInd <- Ind_samp_size[j]
    
    # Put into our list of results
    Tmp_res[[i]] <- Sub_est
    # Get rid of anything just to be safe 
    remove(cols, Sub_het, Error, MSE, Sub_est, gl_mod, sampled_data)
  }
  Boot_res_ind[[j]] <- Tmp_res
  names(Boot_res_ind)[j] <- Ind_samp_size[j]
  Sampled_Ind[[j]] <- Tmp_ind
  names(Sampled_Ind)[j] <- Ind_samp_size[j]
  remove(Tmp_res, Tmp_ind)
  setTxtProgressBar(pb,j)
}
Final_res_ind <- list(Boot_res_ind, Sampled_Ind, Emp_het)
names(Final_res_ind) <- c("Individual_Resampling_Results", "Sampled_Individuals", "Empirical_Estimates")
return(Final_res_ind)
}
