

library(dartR)
library(tidyverse)
library(reshape2)
library(cowplot)

# Read in data
genind_data <- readRDS("bfly_data.rds")

# Record the number of SNPs, round down to nearest 100 so we don't throw an error in our loop 
nSNPs <- ncol(genind_data@tab)
nSNPs <- floor(nSNPs/100)*100
# Record the number of individuals, round down to nearest 5 so we don't throw an error in our loop 
nInds <- nrow(genind_data@tab)
nInds <- floor(nInds/5)*5

# Convert to genlight
gl <- gi2gl(genind_data)

### Calculate empirical statistics 
Emp_het <- gl.report.heterozygosity(gl, plot.out = F)

# Extract columsn we are intersted in by name
cols <- c("Ho", "He", "uHe", "FIS")
Emp_het <- Emp_het[,which(colnames(Emp_het) %in% cols)]

# How many boostraps do we want, technically its not bootstrapping because we are not creating a dataset the same size as our original, so technically its just resampling with replacement
boots <- 100
samp_size <- seq(from = 100, to= nSNPs, by = 100)
Ind_samp_size <- seq(from = 5, to= nInds, by = 5)
Sampled_loc <- list() # Make a list to store the locus names in each boostrap sample, maybe also do some kind of seed??
Boot_res <- list()


# Just test against different SNP thresholds, includes all individuals 
for (j in 1:length(samp_size)) {
  Tmp_res <- list()
  Tmp_loc <- list()
for(i in 1:boots){
  sampled_data <-genind_data[,sample(nLoc(genind_data), size = samp_size[j], replace = T)]
  Tmp_loc[[i]] <- colnames(sampled_data@tab) # Get locus names considered in the sampling replicate
  names(Tmp_loc)[[i]] <- paste("Sample", i, "loci", samp_size[j], sep = '_')
  # Convert for summary stat estimation
  gl_mod <- gi2gl(sampled_data, verbose = 0)
  # Estimate stats
  Sub_het <- gl.report.heterozygosity(gl_mod, plot.out = F)
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
}

# Look for any duplicate samples 
Dup_SNP_check <- lapply(Sampled_loc, as.data.frame)
for(i in 1:length(Dup_SNP_check)){
  tmp <- table(duplicated(Dup_SNP_check[[i]]))
  if(names(tmp) == TRUE) {
    print(paste('Duplicates found in ', names(Dup_SNP_check)[i], sep = ""))
  }
}

remove(Dup_SNP_check)

### Investigate results 
# Bind list into a single data frame
nSNP_res <- bind_rows(Boot_res)

# Look at some stats, we group the data by the statistic (error, MSE, estiamte) and by the number of SNPs in the run
Stats <- nSNP_res %>% group_by(Stat, nSNP) %>% summarise(across(c(Ho, He, uHe, FIS),mean)) %>%
  ungroup()

nSNP_res$Stats <- as.factor(nSNP_res$Stat)
nSNP_res$nSNP <- as.numeric(nSNP_res$nSNP)

Ho_stat <- nSNP_res %>% select(Stat,nSNP, Ho)
He_stat <- nSNP_res %>% select(Stat,nSNP, He)
uHe_stat <- nSNP_res %>% select(Stat,nSNP, uHe)
FIS_stat <- nSNP_res %>% select(Stat,nSNP, FIS)

# Set new names for our plots 
Ho_stat$Stat[Ho_stat$Stat == "Est"] <- 'Observed Heterozygosity'
Ho_stat$Stat[Ho_stat$Stat == "MSE"] <- 'Squared Error'
Ho_stat$Stat <- as.factor(Ho_stat$Stat)
# Reorder levels for plotting purposes
Ho_stat$Stat <- factor(Ho_stat$Stat, levels = c('Observed Heterozygosity', 'Error', 'Squared Error'))
# Check to make sure it worked
levels(Ho_stat$Stat)

## Let's plot it for observed heterozygosity, we can add best fit line using geom_smooth(color = 'black')
Ho_SNPs <- ggplot(data = Ho_stat, aes(x = nSNP, y = Ho, color = Stat)) + geom_point(alpha = 0.5) +
  scale_color_manual(values = c('Squared Error' = '#f46d43', 'Error' = '#a50026', 
                                'Observed Heterozygosity' = '#4575b4')) + 
  theme_classic() + xlab('Number of SNPS') + ylab(NULL) + facet_wrap( ~ Stat, scales = "fixed") + 
  theme(legend.position="none")
                              


# Test different SNP and individual combination
# Make a list to store the locus names in each boostrap sample, maybe also do some kind of seed??
Sampled_Ind <- list()
Boot_res_ind <- list()
Ind_SNP_res <- list()

# Just test against different individual thresholds, includes all SNPs
for (j in 1:length(Ind_samp_size)) {
  Tmp_res <- list()
  Tmp_ind <- list()
  for(i in 1:boots){
    sampled_data <-genind_data[sample(nInd(genind_data), size = Ind_samp_size[k], replace = T),]
    Tmp_ind[[i]] <- rownames(sampled_data@tab) # Get locus names considered in the sampling replicate
    names(Tmp_ind)[[i]] <- paste("Sample", i, "ind", Ind_samp_size[j], sep = '_')
    # Convert for summary stat estimation
    gl_mod <- gi2gl(sampled_data, verbose = 0)
    # Estimate stats
    Sub_het <- gl.report.heterozygosity(gl_mod, plot.out = F)
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
}

# Look for any duplicate samples 
Dup_Ind_check <- lapply(Sampled_Ind, as.data.frame)
for(i in 1:length(Dup_Ind_check)){
  tmp <- table(duplicated(Dup_Ind_check[[i]]))
  if(names(tmp) == TRUE) {
    print(paste('Duplicates found in ', names(Dup_Ind_check)[i], sep = ""))
  }
}

remove(Dup_Ind_check)

### Investigate results 
# Bind list into a single data frame
nInd_res <- bind_rows(Boot_res_ind)

# Look at some stats, we group the data by the statistic (error, MSE, estiamte) and by the number of SNPs in the run
Stats_Ind <- nInd_res %>% group_by(Stat, nInd) %>% summarise(across(c(Ho, He, uHe, FIS),mean)) %>%
  ungroup()

nInd_res$nInd <- as.numeric(nInd_res$nInd)

Ho_stat_ind <- nInd_res %>% select(Stat,nInd, Ho)
He_stat_ind <- nInd_res %>% select(Stat,nInd, He)
uHe_stat_ind <- nInd_res %>% select(Stat,nInd, uHe)
FIS_stat_ind <- nInd_res %>% select(Stat,nInd, FIS)

# Set new names for our plots 
Ho_stat_ind$Stat[Ho_stat_ind$Stat == "Est"] <- 'Observed Heterozygosity'
Ho_stat_ind$Stat[Ho_stat_ind$Stat == "MSE"] <- 'Squared Error'
Ho_stat_ind$Stat <- as.factor(Ho_stat_ind$Stat)
# Reorder levels for plotting purposes
Ho_stat_ind$Stat <- factor(Ho_stat_ind$Stat, levels = c('Observed Heterozygosity', 'Error', 'Squared Error'))
# Check to make sure it worked
levels(Ho_stat$Stat)

## Let's plot it for observed heterozygosity, we can add best fit line using geom_smooth(color = 'black')
Ho_Inds <- ggplot(data = Ho_stat_ind, aes(x = nInd, y = Ho, color = Stat)) + geom_point(alpha = 0.5) +
  scale_color_manual(values = c('Squared Error' = '#f46d43', 'Error' = '#a50026', 
                                'Observed Heterozygosity' = '#4575b4')) + 
  theme_classic() + xlab('Number of Individuals') + ylab(NULL) + facet_wrap( ~ Stat, scales = "fixed") + 
  theme(legend.position="none")


## Make a multi panel of the nSNP and nIND plots

plot_grid(Ho_SNPs, Ho_Inds, nrow = 2)
                          