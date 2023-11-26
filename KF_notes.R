# Purpose: Explore the calculations of diversity statistics 
# Date: October 9th, 2023
# Author(s): Keaka Farleigh
# Email: farleik@miamioh.edu

### Questions/Notes
# 1. Why is the horned lizard data different than others? Maybe just an upload error/difference?
# 2. I modified line 6 of the subsample code to simply list the files based on their extension
# 3. On line 21 of the subsample code you replace all NAs with 0, but our calculation of allele frequencies (p specifically) considers 0's right? Will converting NAs to 0 bias our estimates?
# 4. Are 0's always the dominant allele in the genind objects? If not, should we calculate if 0 or 2 have the higher relative frequency? 
# 5. The current version of the allele frequency calculations only consider homozygotes, we need to also include heterozygotes, see the code below. 


### Load in data 
genind_data <- readRDS("bfly_data.rds")

# Isolate genetic data 
sampled_data <- genind_data[sample(nInd(genind_data), size = 30), ]

# Get the original data (NAs are NA)
Gen_org <- sampled_data@tab

# Modify and extract modified data
sampled_data@tab[is.na(sampled_data@tab)] <- 0
Gen_mod <- sampled_data@tab

# Calculate p and q the original way with and without coding NAs
# KF code to account for NAs
p1 <- colSums(Gen_org[,1:4] == 0, na.rm = T)/((nrow(Gen_org) - colSums(is.na(Gen_org[,1:4])))* 2)

# Original code to caldulate p
p2 <- colSums(Gen_mod[,1:4] == 0) / (nrow(Gen_mod) * 2)

# Compare them 
p1
p2

## Identify whether 0 or 2 have the higher frequencies

dom_allele <- list()
for (i in 1:4){
  freq0 <- ((sum(Gen_org[,i] == 0, na.rm = T)*2) + sum(Gen_org[,i] == 1, na.rm = T))/((nrow(Gen_org) - sum(is.na(Gen_org[,i])))* 2)
  freq2 <- ((sum(Gen_org[,i] == 2, na.rm = T)*2) + sum(Gen_org[,i] == 1, na.rm = T))/((nrow(Gen_org) - sum(is.na(Gen_org[,i])))* 2)
  if(freq0 >= freq2){
    dom_allele[[i]] <- 0
  } else {
    dom_allele[[i]] <- 2
  }
}

# Then you could calculate p and q
p <- list()
q <- list()
for(i in 1:4){
  # Count the number of homozygote individuals with the dominant allele, p_hom = frequency of the dominant allele in the homozygotes, 
  p_hom <- (sum(Gen_org[,i] == dom_allele[[i]], na.rm = T)*2)/((nrow(Gen_org) - sum(is.na(Gen_org[,i])))*2)
  p_het <- sum(Gen_org[,i] == 1, na.rm = T)/((nrow(Gen_org) - sum(is.na(Gen_org[,i]))) * 2)
  p[[i]] <- p_hom + p_het
  q[[i]] <- 1-p[[i]]
  remove(p_hom)
  remove(p_het)
}

