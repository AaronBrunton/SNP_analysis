# Compile list of genind objects and randomly sample for diversity analysis

# setwd("/mnt/Rcache/ajb050/R/SNP_analysis/snp_data_files")

# Read in genind datasets for each species
file_names <- c('ArWolf_data.rds', 'bfly_data.rds', 'bluejay.rds', 'centro_data.rds','coffs_genind.rds','DH_Lizard_data.rds','egk_genind.rds',
                'Cherry_data.rds', 'jawfish.rds', 'kpen_data.rds') # feel free to add more genind object s if you want

# Create an empty list to store the new datasets
new_datasets <- list()

# Loop through each file
for (i in 1:length(file_names)) {
  # Read in the genind dataset from the RDS file
  genind_data <- readRDS(file_names[i])
  
  # Randomly sample 30 individuals
  n_ind_subsample <- 30
  sampled_data <- genind_data[sample(nInd(genind_data), size = n_ind_subsample), ]
  
  # Replace NAs with 0
  sampled_data@tab[is.na(sampled_data@tab)] <- 0
  
  # Convert 'pop' to factor class with a new name (e.g., PopA, PopB, etc.)
  new_pop_name <- paste0("Pop", letters[i])  # This will create PopA, PopB, PopC, ...
  sampled_data$pop <- factor(rep(new_pop_name, nInd(sampled_data)))
  
  # Store the sampled data with the new name in the list
  new_datasets[[new_pop_name]] <- sampled_data
}

# You can access the new datasets like this:
# new_datasets$PopA, new_datasets$PopB, ...
wolf_dat <- new_datasets$Popa
butterfly_dat <- new_datasets$Popb
bluejay_dat <- new_datasets$Popc
centro_dat <- new_datasets$Popd
coffs_dat <- new_datasets$Pope
DHLiz_data <- new_datasets$Popf
egk_dat <- new_datasets$Popg
cherry_dat <- new_datasets$Poph
jawfish_dat <- new_datasets$Popi
kpeng_dat <- new_datasets$Popj


