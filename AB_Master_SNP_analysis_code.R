# Multi-species SNP analysis - AB Aug 2023

# List of data for each species (e.g., popA_species1, popA_species2, ...)
# rename object in list so we know what they are in the large results table

# Define species names
species_names <- c("wolf", "butterfly", "bluejay", "centro", "coffs", "DHLiz", "egk", "cherry", "jawfish", "kpeng")
# Create a list of genind objects for each species
genind_species_list <- list(wolf_dat, butterfly_dat, bluejay_dat, centro_dat, coffs_dat, DHLiz_data, egk_dat, cherry_dat,
                            jawfish_dat, kpeng_dat)

# Load required packages
library(adegenet)
library(knitr)

set.seed(123)

# Function to calculate genetic diversity and FIS
calculateGeneticDiversity <- function(pop) {
  n_ind <- nInd(pop)  # Number of individuals
  
  # Calculate observed heterozygosity (Ho)
  Ho <- colSums(tab(pop) == 1) / n_ind
  
  # Calculate allele frequencies (p and q)
  p <- colSums(tab(pop) == 0) / (n_ind * 2)
  q <- 1 - p
  
  # Calculate expected heterozygosity (He)
  He <- 1 - (p^2 + q^2)
  
  # Calculate unbiased expected heterozygosity (uHe)
  uHe <- He * (2 * n_ind / (2 * n_ind - 1))
  
  # Calculate inbreeding coefficient (FIS)
  FIS <- (mean(He) - mean(Ho)) / mean(He)
  
  # Return the results
  return(data.frame(Observed_Heterozygosity = Ho,
                    Expected_Heterozygosity = He,
                    Unbiased_Expected_Heterozygosity = uHe,
                    FIS = FIS))
}

# Define parameters
replicates <- c(100, 200, 400, 500, 750)  # Number of resampling replicates
sample_sizes <- c(2, 4, 6, 8, 10, 15, 20, 25, 30)  # Number of individuals per population
snp_counts <- c(50, 100, 200, 500, 1000, 1500, 2000)  # Number of SNPs

# Create a data frame to store the multiresults
multiresults <- data.frame(Species = character(),
                           Replicates = integer(),
                           Sample_Size = integer(),
                           SNP_Count = integer(),
                           Observed_Heterozygosity = numeric(),
                           Expected_Heterozygosity = numeric(),
                           Unbiased_Expected_Heterozygosity = numeric(),
                           FIS = numeric())

# Perform resampling and calculate genetic diversity for different parameter combinations
for (i in seq_along(genind_species_list)) {
  current_genind <- genind_species_list[[i]]
  current_species <- species_names[i]
  
  for (replicate_count in replicates) {
    for (n_ind in sample_sizes) {
      for (n_snp in snp_counts) {
        # Create resampled population
        pop_resample <- current_genind[sample(nInd(current_genind), size = n_ind), ]
        
        # Randomly select SNPs
        pop_resample <- pop_resample[, sample(nLoc(pop_resample), size = n_snp)]
        
        # Calculate genetic diversity and FIS
        diversity <- calculateGeneticDiversity(pop_resample)
        
        # Add the multiresults to the data frame
        multiresults <- rbind(multiresults, data.frame(Species = current_species,
                                                       Replicates = replicate_count,
                                                       Sample_Size = n_ind,
                                                       SNP_Count = n_snp,
                                                       Observed_Heterozygosity = mean(diversity$Observed_Heterozygosity, na.rm = TRUE),
                                                       Expected_Heterozygosity = mean(diversity$Expected_Heterozygosity, na.rm = TRUE),
                                                       Unbiased_Expected_Heterozygosity = mean(diversity$Unbiased_Expected_Heterozygosity, na.rm = TRUE),
                                                       FIS = mean(diversity$FIS, na.rm = TRUE)))
      }
    }
  }
}

# Print the multiresults as a table
kable(multiresults, caption = "Genetic Diversity and FIS Estimates for Multiple Species")
head(multiresults)
