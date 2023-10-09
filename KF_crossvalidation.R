# Purpose: Perform cross validation
# Date: October 9th, 2023
# Author(s): Keaka Farleigh
# Email: farleik@miamioh.edu

### Questions/Notes
# 1. Why is the horned lizard data different than others?
# 2. I modified line 6 of the subsample code to simply list the files based on their extension
# 3. It seems like there are two columns per SNP in the tab tables, right? Its nitpicky but shouldn't we just consider one column to get the number of het per SNP?

### Load you pacakges
library(tidyverse)
library(caret)

### Load in data 
Dat <- readRDS("bfly_data.rds")

# Isolate genetic data 
Gen <- Dat[sample(nInd(Dat), size = 30), ]

