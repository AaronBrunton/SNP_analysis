
# List all of our rds files 
Gen_data <- list.files(path = '.', pattern = "*.rds")

# Load functions
source("KF_resampling_functions.R")

# Make lists to store results 
SNP_res <- list()
Ind_res <- list()
# Run for each rds object, this will take a while since we have a lot of species and since we should probably do quite a rew replicates
for (i in 1:length(Gen_data)) {
  SNP_res[[i]] <- SNP_optimization(dat = Gen_data[i], nreps = 1000)
  Ind_res[[i]] <- Ind_optimization(dat = Gen_data[i], nreps = 1000)
}
