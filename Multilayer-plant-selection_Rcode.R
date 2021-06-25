# Methods for the selection of the most important plants in the Norwood farm network #
# Code collated and adapted by Fred M. Windsor (fredric.windsor@newcastle.ac.uk)     #
# The origins of all code are referenced in line and in the main text of the study   #



#### SET-UP #### 

rm(list=ls())
setwd("../")



#### LIBRARIES #### 

library(magrittr); library(bipartite); library(tidyverse); library(igraph); library(gridExtra); 
library(lhs); library(plyr); library(centiserve)




#### FUNCTIONS ####

source("Functions/network_functions.R") # Adapted from Pilosof et al. (2017)
source("Functions/ga.R") # From M'Gonigle et al. (2017)
source("Functions/control.R") # From M'Gonigle et al. (2017)
source("Functions/objective_functions.R") # Adapted from M'Gonigle et al. (2017)
source("Functions/stability_functions.R") # From Sauve et al. (2016)




#### DATA INPUT ####

# Read in the data for the norwood farm network
dframe1 <- read.csv("nore2_aggregated.csv", header = T)

# Sort out the shorthand guild names for the following manipulations
dframe1$lower.names <- substr(dframe1$lower, 1, 4)
dframe1$upper.names <- substr(dframe1$upper, 1, 4)

# Subset the data into different groups (plant-pollinators and plant-herbivore-parasitoids)
plant.poll_edgelist <- dframe1 %>% filter(upper.names == "02FV" | upper.names == "12BF")
plant.herb_edgelist <- dframe1 %>% filter(upper.names == "03AP")
plant.lm.para_edgelist <- dframe1 %>% filter(upper.names == "06MI")

# Remove non plant-polinator interactions for ease 
plant.poll_edgelist <- subset(plant.poll_edgelist, lower.names == "01PL")
plant.herb_edgelist <- subset(plant.herb_edgelist, lower.names == "01PL")
plant.lm.para_edgelist <- subset(plant.lm.para_edgelist, lower.names == "01PL")
multilayer_edgelist <- rbind(plant.poll_edgelist, plant.herb_edgelist, plant.lm.para_edgelist)

# Convert the edgelists into matrices
plant.poll_mat <- edgelist2Matrix(plant.poll_edgelist)
plant.herb_mat <- edgelist2Matrix(plant.herb_edgelist)
plant.lm.para_mat <- edgelist2Matrix(plant.lm.para_edgelist)
multilayer_mat <- edgelist2Matrix(multilayer_edgelist)

# Transpose matrices (required that plants are columns for analyses)
plant.poll_mat <- t(plant.poll_mat)
plant.herb_mat <- t(plant.herb_mat)
plant.lm.para_mat <- t(plant.lm.para_mat)
multilayer_mat <- t(multilayer_mat)

# Convert the matrix to binary for initial assessment
plant.poll_mat_bin <- binarise(plant.poll_mat)
plant.lm.para_mat_bin <- binarise(plant.lm.para_mat)
plant.herb_mat_bin <- binarise(plant.herb_mat)
multilayer_mat_bin <- binarise(multilayer_mat)




#### GENETIC ALGORITHM ####

### PLANT-POLLINATOR NETWORK ### 

## WEIGHTED ##

# Store the matrix as "v.mat"
v.mat <- plant.poll_mat

# Optimising "abundance" for 5 species mixes
mix.a <- find.mix(f=abundance, k=5)

# Optimising "richness" for 5 species mixes
mix.r <- find.mix(f=richness, k=5)

# Optimising both "abundance" and "richness" for 5 species
mix.ar <- find.mix(f=abundance.richness, k=5)

# Compare "abundance" model scores for the above mixes (higher is better)
score(c(abundance), mix.a); score(c(abundance), mix.r); score(c(abundance), mix.ar)

# Compare "richness" model scores for the above mixes (higher is better)
score(c(richness), mix.a); score(c(richness), mix.r); score(c(richness), mix.ar)

# Compare "abundance" and "richness" model scores for the above mixes (higher is better)
score(c(abundance.richness), mix.a); score(c(abundance.richness), mix.r); score(c(abundance.richness), mix.ar)

# Store the data 
mixes_abundance_weighted <- cbind(poll_abun = mix.a)
mixes_richness_weighted <- cbind(poll_rich = mix.r)
mixes_both_weighted <- cbind(poll_both = mix.ar)



### PLANT-LEAF MINER PARASITOID NETWORK ### 

# Store the matrix as "v.mat"
v.mat <- plant.lm.para_mat

# Optimising "abundance" for 5 species mixes
mix.a <- find.mix(f=abundance, k=5)

# Optimising "richness" for 5 species mixes
mix.r <- find.mix(f=richness, k=5)

# Optimising both "abundance" and "richness" for 5 species
mix.ar <- find.mix(f=abundance.richness, k=5)

# Compare "abundance" model scores for the above mixes (higher is better)
score(c(abundance), mix.a); score(c(abundance), mix.r); score(c(abundance), mix.ar)

# Compare "richness" model scores for the above mixes (higher is better)
score(c(richness), mix.a); score(c(richness), mix.r); score(c(richness), mix.ar)

# Compare "abundance" and "richness" model scores for the above mixes (higher is better)
score(c(abundance.richness), mix.a); score(c(abundance.richness), mix.r); score(c(abundance.richness), mix.ar)

# Store the data 
mixes_abundance_weighted <- cbind(mixes_abundance_weighted, lmpara_abun = mix.a)
mixes_richness_weighted <- cbind(mixes_richness_weighted, lmpara_rich = mix.r)
mixes_both_weighted <- cbind(mixes_both_weighted, lmpara_both = mix.ar)



### PLANT-HERBIVORE NETWORK ### 

# Store the matrix as "v.mat"
v.mat <- plant.herb_mat

# Optimising "abundance" for 5 species mixes
mix.a <- find.mix(f=abundance, k=5)

# Optimising "richness" for 5 species mixes
mix.r <- find.mix(f=richness, k=5)

# Optimising both "abundance" and "richness" for 5 species
mix.ar <- find.mix(f=abundance.richness, k=5)

# Compare "abundance" model scores for the above mixes (higher is better)
score(c(abundance), mix.a); score(c(abundance), mix.r); score(c(abundance), mix.ar)

# Compare "richness" model scores for the above mixes (higher is better)
score(c(richness), mix.a); score(c(richness), mix.r); score(c(richness), mix.ar)

# Compare "abundance" and "richness" model scores for the above mixes (higher is better)
score(c(abundance.richness), mix.a); score(c(abundance.richness), mix.r); score(c(abundance.richness), mix.ar)

# Store the data 
mixes_abundance_weighted <- cbind(mixes_abundance_weighted, herb_abun = mix.a)
mixes_richness_weighted <- cbind(mixes_richness_weighted, herb_rich = mix.r)
mixes_both_weighted <- cbind(mixes_both_weighted, herb_both = mix.ar)



### MULTILAYER NETWORK ### 

# Set as multilayer matrix
v.mat <- multilayer_mat_bin

# Store the matrices as "m1", "m2" and "m3"
m1 <- plant.poll_mat
m2 <- plant.lm.para_mat
m3 <- plant.herb_mat

# Optimising "richness" for 5 species mixes (additive method)
mix.r.multilayer.add5 <- find.mix(f=multilayer.additive.richness, k=5)

# Optimising "abundance" for 5 species mixes (additive method)
mix.a.multilayer.add5 <- find.mix(f=multilayer.additive.abundance, k=5)

# Optimising "abundance" for 5 species mixes (additive method)
mix.ar.multilayer.add5 <- find.mix(f=multilayer.additive.both, k=5)




#### SENSITIVITY ANALYSIS ####

### SET UP SAMPLES FOR SENSITIVITY ANALYSIS ###

# Latin hypercube sampling (n = 1000) for five parameters (s, N, p.mutate, p.recomb, p.sex)
samples.LHS <- randomLHS(1000, 5)

# Convert the hypercube samples to a sensible distrbution for variables (we want the algorithms to converge in reasonable time!)
samples.sensitivity <- matrix(nrow=nrow(samples.LHS), ncol=ncol(samples.LHS))
colnames(samples.sensitivity) <- c("s","N", "p.mutate","p.rec","p.sex")
samples.sensitivity[,1] <- qunif(samples.LHS[,1], min = 1, max = 10)
samples.sensitivity[,2] <- qunif(samples.LHS[,2], min = 10, max = 1000)
samples.sensitivity[,3] <- qunif(samples.LHS[,3], min = 0, max = 0.1)
samples.sensitivity[,4] <- qunif(samples.LHS[,4], min = 0.25, max = 0.75)
samples.sensitivity[,5] <- qunif(samples.LHS[,5], min = 0, max = 0.5)

# Round the values to sensible numbers for the parameters 
samples.sensitivity[,1] <- round_any(samples.sensitivity[,1], 0.5)
samples.sensitivity[,2] <- round_any(samples.sensitivity[,2], 1)
samples.sensitivity[,3] <- round_any(samples.sensitivity[,3], 0.01)
samples.sensitivity[,4] <- round_any(samples.sensitivity[,4], 0.01)
samples.sensitivity[,5] <- round_any(samples.sensitivity[,5], 0.01)

# These samples change each time - therefore we also provide the data frame for the LHS we present in the study
# samples.sensitivity <- read.csv("GA_LHS_samples.csv)

# A word of warning - this section was run using a HPC suite (contact fredric.windsor@ncl.ac.uk for SLURM script and R directory)
# The code is not parallelised so a serial run for the bipartite network sensitivity analyses take approx. 6 hours per network
# The multilayer genetic algorithm takes approx. 5.5 days 
# If anyone wants to parallelise the code to increase speed please get in touch with: fredric.windsor@ncl.ac.uk



### PLANT-POLLINATOR NETWORK ###

# Set up two matrices to store the data (1. identities of plants species, 2. Species richness of different taxa)
mix.sensitivity.plant.poll_identity <- data.frame(Plant1 = rep(NA, 1000), Plant2 = rep(NA, 1000), Plant3 = rep(NA, 1000), Plant4 = rep(NA, 1000), Plant5 = rep(NA, 1000))
mix.sensitivity.plant.poll_n <- data.frame(Pollinator_n = rep(NA, 1000), Parasitoid_n = rep(NA, 1000), Herbivore_n = rep(NA, 1000))

# Set data source
v.mat <- plant.poll_mat

# Use the sampled parameters in the GA and store the results in a matrix
for (i in 1:nrow(samples.sensitivity)){
  mix <- find.mix(f=richness, 
                  k = 5, 
                  s = samples.sensitivity[i,"s"], 
                  N=samples.sensitivity[i,"N"], 
                  p.mutate = samples.sensitivity[i,"p.mutate"], 
                  p.rec = samples.sensitivity[i,"p.rec"], 
                  p.sex = samples.sensitivity[i,"p.sex"])
  
  # Species identities
  mix.sensitivity.plant.poll_identity$Plant1[i] <- mix[1]; mix.sensitivity.plant.poll_identity$Plant2[i] <- mix[2]; 
  mix.sensitivity.plant.poll_identity$Plant3[i] <- mix[3]; mix.sensitivity.plant.poll_identity$Plant4[i] <- mix[4]
  mix.sensitivity.plant.poll_identity$Plant5[i] <- mix[5]
  
  # Species richness
  mix.sensitivity.plant.poll_n$Pollinator_n[i] <- table(rowSums(multilayer_mat[1:257,mix])>0)["TRUE"]
  mix.sensitivity.plant.poll_n$Parasitoid_n[i] <- table(rowSums(multilayer_mat[286:381,mix])>0)["TRUE"]
  mix.sensitivity.plant.poll_n$Herbivore_n[i] <- table(rowSums(multilayer_mat[258:380,mix])>0)["TRUE"]
  
  print(i)
  
}

# Merge the sample data and the outputs from the sensitivity analysis 
mix.sensitivity.plant.poll_combined <- cbind(samples.sensitivity, mix.sensitivity.plant.poll_identity, mix.sensitivity.plant.poll_n)



### PLANT-PARASITOID NETWORK ###

# Set up two matrices to store the data (1. identities of plants species, 2. Species richness of different taxa)
mix.sensitivity.plant.lm.para_identity <- data.frame(Plant1 = rep(NA, 1000), Plant2 = rep(NA, 1000), Plant3 = rep(NA, 1000), Plant4 = rep(NA, 1000), Plant5 = rep(NA, 1000))
mix.sensitivity.plant.lm.para_n <- data.frame(Pollinator_n = rep(NA, 1000), Parasitoid_n = rep(NA, 1000), Herbivore_n = rep(NA, 1000))

# Set data source
v.mat <- plant.lm.para_mat

# Use the sampled parameters in the GA and store the results in a matrix
for (i in 1:nrow(samples.sensitivity)){
  mix <- find.mix(f=richness, 
                  k = 5, 
                  s = samples.sensitivity[i,"s"], 
                  N=samples.sensitivity[i,"N"], 
                  p.mutate = samples.sensitivity[i,"p.mutate"], 
                  p.rec = samples.sensitivity[i,"p.rec"], 
                  p.sex = samples.sensitivity[i,"p.sex"])
  
  # Species identities
  mix.sensitivity.plant.lm.para_identity$Plant1[i] <- mix[1]; mix.sensitivity.plant.lm.para_identity$Plant2[i] <- mix[2]; 
  mix.sensitivity.plant.lm.para_identity$Plant3[i] <- mix[3]; mix.sensitivity.plant.lm.para_identity$Plant4[i] <- mix[4]
  mix.sensitivity.plant.lm.para_identity$Plant5[i] <- mix[5]
  
  # Species richness
  mix.sensitivity.plant.lm.para_n$Pollinator_n[i] <- table(rowSums(multilayer_mat[1:257,mix])>0)["TRUE"]
  mix.sensitivity.plant.lm.para_n$Parasitoid_n[i] <- table(rowSums(multilayer_mat[286:381,mix])>0)["TRUE"]
  mix.sensitivity.plant.lm.para_n$Herbivore_n[i] <- table(rowSums(multilayer_mat[258:380,mix])>0)["TRUE"]
  
  print(i)
  
}

# Merge the sample data and the outputs from the sensitivity analysis 
mix.sensitivity.plant.lm.para_combined <- cbind(samples.sensitivity, mix.sensitivity.plant.lm.para_identity, mix.sensitivity.plant.lm.para_n)



### PLANT-HERBIVORE NETWORK ###

# Set up two matrices to store the data (1. identities of plants species, 2. Species richness of different taxa)
mix.sensitivity.plant.herb_identity <- data.frame(Plant1 = rep(NA, 1000), Plant2 = rep(NA, 1000), Plant3 = rep(NA, 1000), Plant4 = rep(NA, 1000), Plant5 = rep(NA, 1000))
mix.sensitivity.plant.herb_n <- data.frame(Pollinator_n = rep(NA, 1000), Parasitoid_n = rep(NA, 1000), Herbivore_n = rep(NA, 1000))

# Set data source
v.mat <- plant.herb_mat

# Use the sampled parameters in the GA and store the results in a matrix
for (i in 1:nrow(samples.sensitivity)){
  mix <- find.mix(f=richness, 
                  k = 5, 
                  s = samples.sensitivity[i,"s"], 
                  N=samples.sensitivity[i,"N"], 
                  p.mutate = samples.sensitivity[i,"p.mutate"], 
                  p.rec = samples.sensitivity[i,"p.rec"], 
                  p.sex = samples.sensitivity[i,"p.sex"])
  
  # Species identities
  mix.sensitivity.plant.herb_identity$Plant1[i] <- mix[1]; mix.sensitivity.plant.herb_identity$Plant2[i] <- mix[2]; 
  mix.sensitivity.plant.herb_identity$Plant3[i] <- mix[3]; mix.sensitivity.plant.herb_identity$Plant4[i] <- mix[4]
  mix.sensitivity.plant.herb_identity$Plant5[i] <- mix[5]
  
  # Species richness
  mix.sensitivity.plant.herb_n$Pollinator_n[i] <- table(rowSums(multilayer_mat[1:257,mix])>0)["TRUE"]
  mix.sensitivity.plant.herb_n$Parasitoid_n[i] <- table(rowSums(multilayer_mat[286:381,mix])>0)["TRUE"]
  mix.sensitivity.plant.herb_n$Herbivore_n[i] <- table(rowSums(multilayer_mat[258:380,mix])>0)["TRUE"]
  
  print(i)
  
}

# Merge the sample data and the outputs from the sensitivity analysis 
mix.sensitivity.plant.herb_combined <- cbind(samples.sensitivity, mix.sensitivity.plant.herb_identity, mix.sensitivity.plant.herb_n)



### MULTILAYER NETWORK ###

# Set up two matrices to store the data (1. identities of plants species, 2. Species richness of different taxa)
mix.sensitivity.plant.multi_identity <- data.frame(Plant1 = rep(NA, 1000), Plant2 = rep(NA, 1000), Plant3 = rep(NA, 1000), Plant4 = rep(NA, 1000), Plant5 = rep(NA, 1000))
mix.sensitivity.plant.multi_n <- data.frame(Pollinator_n = rep(NA, 1000), Parasitoid_n = rep(NA, 1000), Herbivore_n = rep(NA, 1000))

# Set aDATA SOURCE
v.mat <- multilayer_mat

# Store the matrices as "m1", "m2" and "m3"
m1 <- plant.poll_mat
m2 <- plant.lm.para_mat
m3 <- plant.herb_mat

# Use the sampled parameters in the GA and store the results in a matrix
for (i in 1:nrow(samples.sensitivity)){
  mix <- find.mix(f=multilayer.additive.richness, 
                  k = 5, 
                  s = samples.sensitivity[i,"s"], 
                  N=samples.sensitivity[i,"N"], 
                  p.mutate = samples.sensitivity[i,"p.mutate"], 
                  p.rec = samples.sensitivity[i,"p.rec"], 
                  p.sex = samples.sensitivity[i,"p.sex"])
  
  # Species identities
  mix.sensitivity.plant.multi_identity$Plant1[i] <- mix[1]; mix.sensitivity.plant.multi_identity$Plant2[i] <- mix[2]; 
  mix.sensitivity.plant.multi_identity$Plant3[i] <- mix[3]; mix.sensitivity.plant.multi_identity$Plant4[i] <- mix[4]
  mix.sensitivity.plant.multi_identity$Plant5[i] <- mix[5]
  
  # Species richness
  mix.sensitivity.plant.multi_n$Pollinator_n[i] <- table(rowSums(multilayer_mat[1:257,mix])>0)["TRUE"]
  mix.sensitivity.plant.multi_n$Parasitoid_n[i] <- table(rowSums(multilayer_mat[286:381,mix])>0)["TRUE"]
  mix.sensitivity.plant.multi_n$Herbivore_n[i] <- table(rowSums(multilayer_mat[258:285,mix])>0)["TRUE"]
  
  print(i)
  
}

# Merge the sample data and the outputs from the sensitivity analysis 
mix.sensitivity.plant.multi_combined <- cbind(samples.sensitivity, mix.sensitivity.plant.multi_identity, mix.sensitivity.plant.multi_n)



### SPECIES RICHNESS FROM DIFFERENT ALGORITHMS ###

## We want this data stored for plotting (Line 1241)

# Add columns to each table to indicate the mix that the data was derived from
mix.sensitivity.plant.poll_combined$mix <- rep("Pollinator", 1000)
mix.sensitivity.plant.herb_combined$mix <- rep("Herbivore", 1000)
mix.sensitivity.plant.lm.para_combined$mix <- rep("Parasitoid", 1000)
mix.sensitivity.plant.multi_combined$mix <- rep("Multilayer", 1000)

# Join dataframes together for plotting
mix.richness <- rbind(mix.sensitivity.plant.poll_combined,mix.sensitivity.plant.herb_combined,mix.sensitivity.plant.lm.para_combined, mix.sensitivity.plant.multi_combined)

# Order the factors of the mix
mix.richness$mix <- ordered(mix.richness$mix, levels = c("Pollinator", "Parasitoid", "Herbivore", "Multilayer"))




#### NETWORK PROPERTIES OF PLANT MIXTURES ####

# Set up a dataframe to store the results (to compare bipartite and compound)
mix.properties <- data.frame(Mix = c("Pollinator", "Parasitoid", "Herbivore", "Multilayer"), Pollinator_n = NA, Parasitoid_n = NA, Herbivore_n = NA)

### SPECIES RICHNESS ###

## Pollinator plant mix

# Pollinator species richness
mix.properties$Pollinator_n[1] <- table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,1]])>0)["TRUE"]

# Parasitoid species richness
mix.properties$Parasitoid_n[1] <- table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,1]])>0)["TRUE"]

# Herbivore species richness
mix.properties$Herbivore_n[1] <- table(rowSums(multilayer_mat[258:285,mixes_both_weighted[,1]])>0)["TRUE"]


## Parasitoid plant mix

# Pollinator species richness
mix.properties$Pollinator_n[2] <- table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,2]])>0)["TRUE"]

# Parasitoid species richness
mix.properties$Parasitoid_n[2] <- table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,2]])>0)["TRUE"]

# Herbivore species richness
mix.properties$Herbivore_n[2] <- table(rowSums(multilayer_mat[258:285,mixes_both_weighted[,2]])>0)["TRUE"]


## Herbivore plant mix

# Pollinator species richness
mix.properties$Pollinator_n[3] <- table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,3]])>0)["TRUE"]

# Parasitoid species richness
mix.properties$Parasitoid_n[3] <- table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,3]])>0)["TRUE"]

# Herbivore species richness
mix.properties$Herbivore_n[3] <- table(rowSums(multilayer_mat[258:285,mixes_both_weighted[,3]])>0)["TRUE"]


## Multilayer plant mix

# Pollinator species richness
mix.properties$Pollinator_n[4] <- table(rowSums(multilayer_mat[1:257,mix.r.multilayer.add5])>0)["TRUE"]

# Parasitoid species richness
mix.properties$Parasitoid_n[4] <- table(rowSums(multilayer_mat[286:381,mix.r.multilayer.add5])>0)["TRUE"]

# Herbivore species richness
mix.properties$Herbivore_n[4] <- table(rowSums(multilayer_mat[258:285,mix.r.multilayer.add5])>0)["TRUE"]


## Difference from optimal (calculate differentce from optimal as dictated by the bipartite network algortihms)
mix.properties$Pollinator_dfo <- mix.properties$Pollinator_n - mix.properties$Pollinator_n[1]
mix.properties$Parasitoid_dfo <- mix.properties$Parasitoid_n - mix.properties$Parasitoid_n[2]
mix.properties$Herbivore_dfo <- mix.properties$Herbivore_n - mix.properties$Herbivore_n[3]

# Reorder mixes for plotting
mix.properties$Mix <- ordered(mix.properties$Mix, levels = c("Pollinator", "Parasitoid", "Herbivore", "Multilayer"))



### DEGREE ### 

## Pollinator plant mix 
mean(degree((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,1]])))))

## Parasitoid plant mix
mean(degree((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,2]])))))

## Herbivore plant mix
mean(degree((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,3]])))))

## Multilayer plant mix
mean(degree((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mix.r.multilayer.add5])))))



### KATZ CENTRALITY ###

## Pollinator plant mix 
mean(katzcent((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,1]]))), alpha = 0.05))

## Parasitoid plant mix
mean(katzcent((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,2]]))), alpha = 0.05))

## Herbivore plant mix
mean(katzcent((graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,3]]))), alpha = 0.05))

## Multilayer plant mix
mean(katzcent(graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mix.r.multilayer.add5])), alpha = 0.05))




### CONNECTANCE ### 

## Pollinator plant mix 
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]])))

## Parasitoid plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,2]])))

## Herbivore plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,3]])))

## Multilayer plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.add5])))




### ROBUSTNESS ### (from Bane et al. 2018)

## Generate networks to work from

rob.plant.poll_mat <- multilayer_mat[1:257,]
rob.plant.herb_mat <- multilayer_mat[258:285,]
rob.plant.lm.para_mat <- multilayer_mat[286:381,]


## Pollinator plant mix 

# Remove plants that are not shared between layers
commonplants.poll <- mixes_both_weighted[,1]
plant.poll.poll_mat <- rob.plant.poll_mat[,commonplants.poll]
plant.herb.poll_mat <- rob.plant.herb_mat[,commonplants.poll]
plant.lm.para.poll_mat <- rob.plant.lm.para_mat[,commonplants.poll]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_Rvalues <- NULL
herb_Rvalues <- NULL
lmpara_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll.poll_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll.poll_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll.poll_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll.poll_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb.poll_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para.poll_mat # make copy of original plant-leaf miner parasitoid matrix
  
  PLANT<-colSums(pl.p_mat)    # save plant degrees
  POLL<-rowSums(pl.p_mat)     # save polinator degrees
  HERB<-rowSums(pl.h_mat)     # save herbivore degrees
  LMPARA<-rowSums(pl.lmp_mat) # save leaf miner parasitoid degrees
  
  poll_survivors<-NULL   # create survivors to save pollinator counts
  herb_survivors<-NULL   # create survivors to save herbivore counts
  lmpara_survivors<-NULL # create survivors to save herbivore counts
  
  plantdeaths<-NULL              # create survivors to save pollinator counts
  triggerhat<-colnames(pl.p_mat) # create hat with plant names to pick from
  pastplantdeaths<-NULL          # set plant deaths to zero
  
  ntriggerhat<-NULL # set the number in the hat to zero
  i=0               # used to iterate in a for loop later
  j=j+1             # used to iterate in a for loop later
  ttally=0          # set ttally to 0
  
  while(sum(colSums(pl.p_mat))>0){ # while there are still plant nodes left in the networks
    
    if(length(ntriggerhat)>0){              # if trigger is greater than 0 (i.e. >1 plant has been made extinct) then...
      ntriggerhat<-sample(ntriggerhat)      # sample the vector of trigger values              
      cascadelength[c]<-length(ntriggerhat) # record cascade value in column c, the number of trigger values
      c=c+1                                 # set c to the next value (e.g. current value + 1)
      
      for (m in 1:length(ntriggerhat)){ # then for each value recorded in the trigger hat
        ntrigger<-ntriggerhat[m]        # set the column of ntrigger to the correct number
        pl.p_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        pl.h_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        
        pastplantdeaths<-c(pastplantdeaths,ntrigger) # record the plant deaths as the value of plants and triggers lost
        i=i+1                                        # set i to the next value (e.g. current value + 1)
        
        poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators
        pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct
        poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # set survivors as pollinators with more than 1 observation
        
        herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check herbivores
        pl.h_mat[herb.ext,]<-0                                      # make herbivores extinct
        herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # set survivors as herbivores with more than 1 observation
        
        lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check leaf miner parasitoids
        pl.lmp_mat[lmpara.ext,]<-0                                          # make leaf miner parasitoids extinct
        lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # set survivors as leaf miner parasitoids with more than 1 observation
        
        trigmat[j,i]<-1                                           # record the run in the trigmat (shows that this section was used)
        plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # set number of plant survivors 
        exttimesC[k,i]<-ntrigger                                  # record the number of primary extinctions that have occurred
        
      }
    }
    else{
      rtrigger<-sample((names(colSums(pl.p_mat)[(colSums(pl.p_mat))>0])),1) # sample a trigger from present plant taxa (those which are not extinct)
      
      pl.p_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.h_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.lmp_mat[,rtrigger]<-0 # set the trigger to zero (e.g. make it extinct)
      
      pastplantdeaths<-c(pastplantdeaths,rtrigger) # names of the past plant deaths plus the new plant death (rtrigger)
      i=i+1                                        # set i to the next value (e.g. current value + 1)
      
      poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # save number of survivors
      
      herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.h_mat[herb.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # save number of survivors
      
      lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.lmp_mat[lmpara.ext,]<-0                                          # make pollinators extinct (those which have lost over 50% of individuals)
      lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # save number of survivors
      
      ttally<-ttally+1                                          # set the tally to the next value
      trigmat[j,i]<-0                                           # identify that nothing was removed
      plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # record plant survivors
      exttimesC[k,i]<-rtrigger                                  # record the extinct node (a matrix of the order of species removal)
      
    }
    
    p.ext<-which(((PLANT-colSums(pl.p_mat))/PLANT)>=threshold)                                  # record the plants that have gone extinct (e.g. lost 50% of abundance)
    ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext))) # remove the plants that                                            
    
  }
  
  triggertally[k]<-ttally                                 # record the number of triggers used
  poll_survivorsx<-c(nrow(pl.p_mat),poll_survivors)       # record the number of pollinator survivor and the total number of starting taxa
  herb_survivorsx<-c(nrow(pl.h_mat),herb_survivors)       # record the number of herbivore survivor and the total number of starting taxa
  lmpara_survivorsx<-c(nrow(pl.lmp_mat),lmpara_survivors) # record the number of leaf miner parasitoid survivor and the total number of starting taxa
  
  poll_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_Rvalues)
median(herb_Rvalues)
median(lmpara_Rvalues)


## Parasitoid plant mix

# Remove plants that are not shared between layers
commonplants.para <- mixes_both_weighted[,2]
plant.poll.para_mat <- rob.plant.poll_mat[,commonplants.para]
plant.herb.para_mat <- rob.plant.herb_mat[,commonplants.para]
plant.lm.para.para_mat <- rob.plant.lm.para_mat[,commonplants.para]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_Rvalues <- NULL
herb_Rvalues <- NULL
lmpara_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll.para_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll.para_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll.para_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll.para_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb.para_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para.para_mat # make copy of original plant-leaf miner parasitoid matrix
  
  PLANT<-colSums(pl.p_mat)    # save plant degrees
  POLL<-rowSums(pl.p_mat)     # save polinator degrees
  HERB<-rowSums(pl.h_mat)     # save herbivore degrees
  LMPARA<-rowSums(pl.lmp_mat) # save leaf miner parasitoid degrees
  
  poll_survivors<-NULL   # create survivors to save pollinator counts
  herb_survivors<-NULL   # create survivors to save herbivore counts
  lmpara_survivors<-NULL # create survivors to save herbivore counts
  
  plantdeaths<-NULL              # create survivors to save pollinator counts
  triggerhat<-colnames(pl.p_mat) # create hat with plant names to pick from
  pastplantdeaths<-NULL          # set plant deaths to zero
  
  ntriggerhat<-NULL # set the number in the hat to zero
  i=0               # used to iterate in a for loop later
  j=j+1             # used to iterate in a for loop later
  ttally=0          # set ttally to 0
  
  while(sum(colSums(pl.p_mat))>0){ # while there are still plant nodes left in the networks
    
    if(length(ntriggerhat)>0){              # if trigger is greater than 0 (i.e. >1 plant has been made extinct) then...
      ntriggerhat<-sample(ntriggerhat)      # sample the vector of trigger values              
      cascadelength[c]<-length(ntriggerhat) # record cascade value in column c, the number of trigger values
      c=c+1                                 # set c to the next value (e.g. current value + 1)
      
      for (m in 1:length(ntriggerhat)){ # then for each value recorded in the trigger hat
        ntrigger<-ntriggerhat[m]        # set the column of ntrigger to the correct number
        pl.p_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        pl.h_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        
        pastplantdeaths<-c(pastplantdeaths,ntrigger) # record the plant deaths as the value of plants and triggers lost
        i=i+1                                        # set i to the next value (e.g. current value + 1)
        
        poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators
        pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct
        poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # set survivors as pollinators with more than 1 observation
        
        herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check herbivores
        pl.h_mat[herb.ext,]<-0                                      # make herbivores extinct
        herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # set survivors as herbivores with more than 1 observation
        
        lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check leaf miner parasitoids
        pl.lmp_mat[lmpara.ext,]<-0                                          # make leaf miner parasitoids extinct
        lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # set survivors as leaf miner parasitoids with more than 1 observation
        
        trigmat[j,i]<-1                                           # record the run in the trigmat (shows that this section was used)
        plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # set number of plant survivors 
        exttimesC[k,i]<-ntrigger                                  # record the number of primary extinctions that have occurred
        
      }
    }
    else{
      rtrigger<-sample((names(colSums(pl.p_mat)[(colSums(pl.p_mat))>0])),1) # sample a trigger from present plant taxa (those which are not extinct)
      
      pl.p_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.h_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.lmp_mat[,rtrigger]<-0 # set the trigger to zero (e.g. make it extinct)
      
      pastplantdeaths<-c(pastplantdeaths,rtrigger) # names of the past plant deaths plus the new plant death (rtrigger)
      i=i+1                                        # set i to the next value (e.g. current value + 1)
      
      poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # save number of survivors
      
      herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.h_mat[herb.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # save number of survivors
      
      lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.lmp_mat[lmpara.ext,]<-0                                          # make pollinators extinct (those which have lost over 50% of individuals)
      lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # save number of survivors
      
      ttally<-ttally+1                                          # set the tally to the next value
      trigmat[j,i]<-0                                           # identify that nothing was removed
      plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # record plant survivors
      exttimesC[k,i]<-rtrigger                                  # record the extinct node (a matrix of the order of species removal)
      
    }
    
    p.ext<-which(((PLANT-colSums(pl.p_mat))/PLANT)>=threshold)                                  # record the plants that have gone extinct (e.g. lost 50% of abundance)
    ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext))) # remove the plants that                                            
    
  }
  
  triggertally[k]<-ttally                                 # record the number of triggers used
  poll_survivorsx<-c(nrow(pl.p_mat),poll_survivors)       # record the number of pollinator survivor and the total number of starting taxa
  herb_survivorsx<-c(nrow(pl.h_mat),herb_survivors)       # record the number of herbivore survivor and the total number of starting taxa
  lmpara_survivorsx<-c(nrow(pl.lmp_mat),lmpara_survivors) # record the number of leaf miner parasitoid survivor and the total number of starting taxa
  
  poll_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_Rvalues)
median(herb_Rvalues)
median(lmpara_Rvalues)


## Herbivore plant mix

# Remove plants that are not shared between layers
commonplants.herb <- mixes_both_weighted[,3]
plant.poll.herb_mat <- rob.plant.poll_mat[,commonplants.herb]
plant.herb.herb_mat <- rob.plant.herb_mat[,commonplants.herb]
plant.lm.para.herb_mat <- rob.plant.lm.para_mat[,commonplants.herb]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_Rvalues <- NULL
herb_Rvalues <- NULL
lmpara_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll.herb_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll.herb_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll.herb_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll.herb_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb.herb_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para.herb_mat # make copy of original plant-leaf miner parasitoid matrix
  
  PLANT<-colSums(pl.p_mat)    # save plant degrees
  POLL<-rowSums(pl.p_mat)     # save polinator degrees
  HERB<-rowSums(pl.h_mat)     # save herbivore degrees
  LMPARA<-rowSums(pl.lmp_mat) # save leaf miner parasitoid degrees
  
  poll_survivors<-NULL   # create survivors to save pollinator counts
  herb_survivors<-NULL   # create survivors to save herbivore counts
  lmpara_survivors<-NULL # create survivors to save herbivore counts
  
  plantdeaths<-NULL              # create survivors to save pollinator counts
  triggerhat<-colnames(pl.p_mat) # create hat with plant names to pick from
  pastplantdeaths<-NULL          # set plant deaths to zero
  
  ntriggerhat<-NULL # set the number in the hat to zero
  i=0               # used to iterate in a for loop later
  j=j+1             # used to iterate in a for loop later
  ttally=0          # set ttally to 0
  
  while(sum(colSums(pl.p_mat))>0){ # while there are still plant nodes left in the networks
    
    if(length(ntriggerhat)>0){              # if trigger is greater than 0 (i.e. >1 plant has been made extinct) then...
      ntriggerhat<-sample(ntriggerhat)      # sample the vector of trigger values              
      cascadelength[c]<-length(ntriggerhat) # record cascade value in column c, the number of trigger values
      c=c+1                                 # set c to the next value (e.g. current value + 1)
      
      for (m in 1:length(ntriggerhat)){ # then for each value recorded in the trigger hat
        ntrigger<-ntriggerhat[m]        # set the column of ntrigger to the correct number
        pl.p_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        pl.h_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        
        pastplantdeaths<-c(pastplantdeaths,ntrigger) # record the plant deaths as the value of plants and triggers lost
        i=i+1                                        # set i to the next value (e.g. current value + 1)
        
        poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators
        pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct
        poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # set survivors as pollinators with more than 1 observation
        
        herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check herbivores
        pl.h_mat[herb.ext,]<-0                                      # make herbivores extinct
        herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # set survivors as herbivores with more than 1 observation
        
        lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check leaf miner parasitoids
        pl.lmp_mat[lmpara.ext,]<-0                                          # make leaf miner parasitoids extinct
        lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # set survivors as leaf miner parasitoids with more than 1 observation
        
        trigmat[j,i]<-1                                           # record the run in the trigmat (shows that this section was used)
        plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # set number of plant survivors 
        exttimesC[k,i]<-ntrigger                                  # record the number of primary extinctions that have occurred
        
      }
    }
    else{
      rtrigger<-sample((names(colSums(pl.p_mat)[(colSums(pl.p_mat))>0])),1) # sample a trigger from present plant taxa (those which are not extinct)
      
      pl.p_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.h_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.lmp_mat[,rtrigger]<-0 # set the trigger to zero (e.g. make it extinct)
      
      pastplantdeaths<-c(pastplantdeaths,rtrigger) # names of the past plant deaths plus the new plant death (rtrigger)
      i=i+1                                        # set i to the next value (e.g. current value + 1)
      
      poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # save number of survivors
      
      herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.h_mat[herb.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # save number of survivors
      
      lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.lmp_mat[lmpara.ext,]<-0                                          # make pollinators extinct (those which have lost over 50% of individuals)
      lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # save number of survivors
      
      ttally<-ttally+1                                          # set the tally to the next value
      trigmat[j,i]<-0                                           # identify that nothing was removed
      plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # record plant survivors
      exttimesC[k,i]<-rtrigger                                  # record the extinct node (a matrix of the order of species removal)
      
    }
    
    p.ext<-which(((PLANT-colSums(pl.p_mat))/PLANT)>=threshold)                                  # record the plants that have gone extinct (e.g. lost 50% of abundance)
    ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext))) # remove the plants that                                            
    
  }
  
  triggertally[k]<-ttally                                 # record the number of triggers used
  poll_survivorsx<-c(nrow(pl.p_mat),poll_survivors)       # record the number of pollinator survivor and the total number of starting taxa
  herb_survivorsx<-c(nrow(pl.h_mat),herb_survivors)       # record the number of herbivore survivor and the total number of starting taxa
  lmpara_survivorsx<-c(nrow(pl.lmp_mat),lmpara_survivors) # record the number of leaf miner parasitoid survivor and the total number of starting taxa
  
  poll_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_Rvalues)
median(herb_Rvalues)
median(lmpara_Rvalues)


## Multilayer (additive) plant mix

# Remove plants that are not shared between layers
commonplants.add <- mix.r.multilayer.add5
plant.poll.add_mat <- rob.plant.poll_mat[,commonplants.add]
plant.herb.add_mat <- rob.plant.herb_mat[,commonplants.add]
plant.lm.para.add_mat <- rob.plant.lm.para_mat[,commonplants.add]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_Rvalues <- NULL
herb_Rvalues <- NULL
lmpara_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll.add_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll.add_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll.add_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll.add_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb.add_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para.add_mat # make copy of original plant-leaf miner parasitoid matrix
  
  PLANT<-colSums(pl.p_mat)    # save plant degrees
  POLL<-rowSums(pl.p_mat)     # save polinator degrees
  HERB<-rowSums(pl.h_mat)     # save herbivore degrees
  LMPARA<-rowSums(pl.lmp_mat) # save leaf miner parasitoid degrees
  
  poll_survivors<-NULL   # create survivors to save pollinator counts
  herb_survivors<-NULL   # create survivors to save herbivore counts
  lmpara_survivors<-NULL # create survivors to save herbivore counts
  
  plantdeaths<-NULL              # create survivors to save pollinator counts
  triggerhat<-colnames(pl.p_mat) # create hat with plant names to pick from
  pastplantdeaths<-NULL          # set plant deaths to zero
  
  ntriggerhat<-NULL # set the number in the hat to zero
  i=0               # used to iterate in a for loop later
  j=j+1             # used to iterate in a for loop later
  ttally=0          # set ttally to 0
  
  while(sum(colSums(pl.p_mat))>0){ # while there are still plant nodes left in the networks
    
    if(length(ntriggerhat)>0){              # if trigger is greater than 0 (i.e. >1 plant has been made extinct) then...
      ntriggerhat<-sample(ntriggerhat)      # sample the vector of trigger values              
      cascadelength[c]<-length(ntriggerhat) # record cascade value in column c, the number of trigger values
      c=c+1                                 # set c to the next value (e.g. current value + 1)
      
      for (m in 1:length(ntriggerhat)){ # then for each value recorded in the trigger hat
        ntrigger<-ntriggerhat[m]        # set the column of ntrigger to the correct number
        pl.p_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        pl.h_mat[,ntrigger]<-0          # set the triggers to zero (e.g. make them extinct)
        
        pastplantdeaths<-c(pastplantdeaths,ntrigger) # record the plant deaths as the value of plants and triggers lost
        i=i+1                                        # set i to the next value (e.g. current value + 1)
        
        poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators
        pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct
        poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # set survivors as pollinators with more than 1 observation
        
        herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check herbivores
        pl.h_mat[herb.ext,]<-0                                      # make herbivores extinct
        herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # set survivors as herbivores with more than 1 observation
        
        lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check leaf miner parasitoids
        pl.lmp_mat[lmpara.ext,]<-0                                          # make leaf miner parasitoids extinct
        lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # set survivors as leaf miner parasitoids with more than 1 observation
        
        trigmat[j,i]<-1                                           # record the run in the trigmat (shows that this section was used)
        plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # set number of plant survivors 
        exttimesC[k,i]<-ntrigger                                  # record the number of primary extinctions that have occurred
        
      }
    }
    else{
      rtrigger<-sample((names(colSums(pl.p_mat)[(colSums(pl.p_mat))>0])),1) # sample a trigger from present plant taxa (those which are not extinct)
      
      pl.p_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.h_mat[,rtrigger]<-0   # set the trigger to zero (e.g. make it extinct)
      pl.lmp_mat[,rtrigger]<-0 # set the trigger to zero (e.g. make it extinct)
      
      pastplantdeaths<-c(pastplantdeaths,rtrigger) # names of the past plant deaths plus the new plant death (rtrigger)
      i=i+1                                        # set i to the next value (e.g. current value + 1)
      
      poll.ext<-which(((POLL-rowSums(pl.p_mat))/POLL)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.p_mat[poll.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      poll_survivors[i]<-length(which(rowSums(pl.p_mat)>0))       # save number of survivors
      
      herb.ext<-which(((HERB-rowSums(pl.h_mat))/HERB)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.h_mat[herb.ext,]<-0                                      # make pollinators extinct (those which have lost over 50% of individuals)
      herb_survivors[i]<-length(which(rowSums(pl.h_mat)>0))       # save number of survivors
      
      lmpara.ext<-which(((LMPARA-rowSums(pl.lmp_mat))/LMPARA)>=threshold) # check pollinators (which have lost more than 50% of individuals)
      pl.lmp_mat[lmpara.ext,]<-0                                          # make pollinators extinct (those which have lost over 50% of individuals)
      lmpara_survivors[i]<-length(which(rowSums(pl.lmp_mat)>0))           # save number of survivors
      
      ttally<-ttally+1                                          # set the tally to the next value
      trigmat[j,i]<-0                                           # identify that nothing was removed
      plantsurvivorsC[j,i]<-length(which(colSums(pl.p_mat)==0)) # record plant survivors
      exttimesC[k,i]<-rtrigger                                  # record the extinct node (a matrix of the order of species removal)
      
    }
    
    p.ext<-which(((PLANT-colSums(pl.p_mat))/PLANT)>=threshold)                                  # record the plants that have gone extinct (e.g. lost 50% of abundance)
    ntriggerhat<-c(setdiff(names(p.ext),pastplantdeaths),setdiff(pastplantdeaths,names(p.ext))) # remove the plants that                                            
    
  }
  
  triggertally[k]<-ttally                                 # record the number of triggers used
  poll_survivorsx<-c(nrow(pl.p_mat),poll_survivors)       # record the number of pollinator survivor and the total number of starting taxa
  herb_survivorsx<-c(nrow(pl.h_mat),herb_survivors)       # record the number of herbivore survivor and the total number of starting taxa
  lmpara_survivorsx<-c(nrow(pl.lmp_mat),lmpara_survivors) # record the number of leaf miner parasitoid survivor and the total number of starting taxa
  
  poll_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_Rvalues)
median(herb_Rvalues)
median(lmpara_Rvalues)




### STABILITY ### 

## Pollinator plant mix

# Create adjacency matrix
multilayer.poll_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,1]]))))

# Create a Jacobian matrix
multilayer.poll_jac.mat <- jacobian_binary(multilayer.poll_adj.mat)

# Calculate the stability
stability(multilayer.poll_jac.mat, s2 = 13)


## Parasitoid plant mix

# Create adjacency matrix
multilayer.para_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,2]]))))

# Create a Jacobian matrix
multilayer.para_jac.mat <- jacobian_binary(multilayer.para_adj.mat)

# Calculate the stability
stability(multilayer.para_jac.mat, s2 = 12)


## Herbivore plant mix

# Create adjacency matrix
multilayer.herb_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mixes_both_weighted[,3]]))))

# Create a Jacobian matrix
multilayer.herb_jac.mat <- jacobian_binary(multilayer.herb_adj.mat)

# Calculate the stability
stability(multilayer.herb_jac.mat, s2 = 10)


## Multilayer plant mix

# Create adjacency matrix
multilayer.add_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(bipartite::empty(multilayer_mat[,mix.r.multilayer.add5]))))

# Create a Jacobian matrix
multilayer.add_jac.mat <- jacobian_binary(multilayer.add_adj.mat)

# Calculate the stability
stability(multilayer.add_jac.mat, s2 = 13)




#### FIGURE 2 ####

# Plot species richness of pollinators for each mix
plot2a <- ggplot(aes(x=mix, y = Pollinator_n, fill = mix), data = mix.richness) + 
  geom_violin() + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), 
        legend.position = "NA", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("gold","darkblue","yellowgreen","purple3")) +
  ylab("Species richness (n)") + 
  coord_cartesian(ylim = c(0,200)) + 
  ggtitle("a")
plot2a

# Plot species richness of parasitoids for each mix
plot2b <- ggplot(aes(x=mix, y = Parasitoid_n, fill = mix), data = mix.richness) + 
  geom_violin() + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), 
        legend.position = "NA", axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("gold","darkblue","yellowgreen","purple3")) +
  coord_cartesian(ylim = c(0,80)) +
  ggtitle("b")
plot2b

# Plot species richness of herbivores for each mix
plot2c <- ggplot(aes(x=mix, y = Herbivore_n, fill = mix), data = mix.richness) + 
  geom_violin() + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), 
        legend.position = "NA", axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("gold","darkblue","yellowgreen","purple3")) +
  coord_cartesian(ylim = c(0,15)) +
  scale_y_continuous(breaks = c(0,5,10,15)) + 
  ggtitle("c")
plot2c

# Multiplot
grid.arrange(plot2a, plot2b, plot2c, nrow = 1, bottom = "Plant species mix")




#### FIGURE 3 ####

# Calculating the number of plant species across the sensitivity analyses
plant.S_poll <- unique(unlist(mix.sensitivity.plant.poll_combined[,6:10]))
plant.S_herb <- unique(unlist(mix.sensitivity.plant.herb_combined[,6:10]))
plant.S_para <- unique(unlist(mix.sensitivity.plant.lm.para_combined[,6:10]))

mean(c(length(plant.S_poll), length(plant.S_herb), length(plant.S_para)))

plant.S_multi <- unique(unlist(mix.sensitivity.plant.multi_combined[,6:10]))
length(plant.S_multi)

# Get summaries of the counts for each species of plant in each mix
plant.poll.identity_summary <- data.frame(table(unlist(mix.sensitivity.plant.poll_combined[,6:10])))
plant.herb.identity_summary <- data.frame(table(unlist(mix.sensitivity.plant.herb_combined[,6:10])))
plant.lm.para.identity_summary <- data.frame(table(unlist(mix.sensitivity.plant.lm.para_combined[,6:10])))
plant.multi.identity_summary <- data.frame(table(unlist(mix.sensitivity.plant.multi_combined[,6:10])))

# Alter plant names for plotting 
plant.poll.identity_summary$plant <- substring(plant.poll.identity_summary$Var1, first = 9)
plant.herb.identity_summary$plant <- substring(plant.herb.identity_summary$Var1, first = 9)
plant.lm.para.identity_summary$plant <- substring(plant.lm.para.identity_summary$Var1, first = 9)
plant.multi.identity_summary$plant <- substring(plant.multi.identity_summary$Var1, first = 9)

# Sort the dataframes by frequency
plant.poll.identity_summary <- transform(plant.poll.identity_summary, plant = reorder(plant, -Freq))
plant.herb.identity_summary <- transform(plant.herb.identity_summary, plant = reorder(plant, -Freq))
plant.lm.para.identity_summary <- transform(plant.lm.para.identity_summary, plant = reorder(plant, -Freq))
plant.multi.identity_summary <- transform(plant.multi.identity_summary, plant = reorder(plant, -Freq))

# Subset to the top 15 
plant.poll.identity_summary_sub <- plant.poll.identity_summary[with(plant.poll.identity_summary, order(-Freq)),][1:12,]
plant.herb.identity_summary_sub <- plant.herb.identity_summary[with(plant.herb.identity_summary, order(-Freq)),][1:12,]
plant.lm.para.identity_summary_sub <- plant.lm.para.identity_summary[with(plant.lm.para.identity_summary, order(-Freq)),][1:12,]
plant.multi.identity_summary_sub <- plant.multi.identity_summary[with(plant.multi.identity_summary, order(-Freq)),][1:12,]

# Plotting the frequency that species occur in the optimal mixes
plot3a <- ggplot(aes(x=plant, y=Freq), data = data.frame(plant.poll.identity_summary_sub)) + 
  stat_summary(geom = "bar", fun.y = "identity", fill = "gold", colour = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"), plot.margin = unit(c(5,5,5,10), "mm")) +
  xlab("") + 
  ylab("") +
  ggtitle("a")
plot3a

plot3b <- ggplot(aes(x=plant, y=Freq), data = data.frame(plant.herb.identity_summary_sub)) + 
  stat_summary(geom = "bar", fun.y = "identity", fill = "yellowgreen", colour = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"), plot.margin = unit(c(5,5,5,10), "mm")) +
  xlab("") + 
  ylab("") + 
  ggtitle("b")
plot3b

plot3c <- ggplot(aes(x=plant, y=Freq), data = data.frame(plant.lm.para.identity_summary_sub)) + 
  stat_summary(geom = "bar", fun.y = "identity", fill = "darkblue", colour = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"), plot.margin = unit(c(5,5,5,10), "mm")) +
  xlab("") + 
  ylab("") + 
  ggtitle("c")
plot3c

plot3d <- ggplot(aes(x=plant, y=Freq), data = data.frame(plant.multi.identity_summary_sub)) + 
  stat_summary(geom = "bar", fun.y = "identity", fill = "purple3", colour = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1, angle = 45), axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"), plot.margin = unit(c(5,5,5,10), "mm")) +
  xlab("") + 
  ylab("") + 
  ggtitle("d")
plot3d

grid.arrange(plot3a, plot3b, plot3c, plot3d, nrow = 2, bottom = "Plant Species", left = "Frequency in optimal mix (n)")




#### FIGURE 4 ####

## BIPARTITE SUBNETWORK PLANT MIX NETWORK

# Pollinator plant mix
plant.poll_graph <- graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]]))

# Set the node properties
V(plant.poll_graph)$code <- substr(names(V(plant.poll_graph)), 1, 4)
V(plant.poll_graph)$color <- V(plant.poll_graph)$code

# Set the colours to the appropriate values
V(plant.poll_graph)$color <- gsub("01PL", "forestgreen", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("03AP", "yellowgreen", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("06MI", "darkblue", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("02FV", "gold", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("12BF", "gold", V(plant.poll_graph)$color)

# Plot the network
plot(plant.poll_graph, vertex.label = NA, vertex.size = 6, vertex.color = V(plant.poll_graph)$color)
legend("topleft", legend = "a", bty = "n", cex = 2, inset = c(0.005,0.005))
legend("topleft", legend = c("Pollinators", "Herbivores", "Parasitoids", "Plants"), 
       col = unique(V(plant.poll_graph)$color), bty = "n", pch = 20, pt.cex = 3, cex = 1, horiz = F, inset = c(0.06,0.15))


## MULTILAYER NETWORK PLANT MIX NETWORK 

# Multilayer plant mix 
plant.multilayer_graph <- graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.add5]))

# Set the node properties
V(plant.multilayer_graph)$code <- substr(names(V(plant.multilayer_graph)), 1, 4)
V(plant.multilayer_graph)$color <- V(plant.multilayer_graph)$code

# Set the colours to the appropriate values
V(plant.multilayer_graph)$color <- gsub("01PL", "forestgreen", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("03AP", "yellowgreen", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("06MI", "darkblue", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("02FV", "gold", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("12BF", "gold", V(plant.multilayer_graph)$color)

# Set layout
par(mfrow=c(1,2))

# Plot the network
plot(plant.poll_graph, vertex.label = NA, vertex.size = 6, vertex.color = V(plant.poll_graph)$color)
legend("topleft", legend = "a", bty = "n", cex = 2, inset = c(-0.23,0.005))
legend("topleft", legend = c("Pollinators", "Herbivores", "Parasitoids", "Plants"), 
       col = unique(V(plant.poll_graph)$color), bty = "n", pch = 20, pt.cex = 3, cex = 1, horiz = F, inset = c(-0.1,0.15))

# Plot the network
plot(plant.multilayer_graph, vertex.label = NA, vertex.size = 6, vertex.color = V(plant.multilayer_graph)$color)
legend("topleft", legend = "b", bty = "n", cex = 2, inset = c(-0.23,0.005))


## SUPPLEMENTARY FIGURE 1 ##

plota <- ggplot() + 
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = s), size = 0.8, width = 0.1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = s), size = 0.8, width = 0.1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = s), size = 0.8, width = 0.1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = s), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = s), size = 0.8, width = 0.1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(50,200), xlim = c(1,10)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") +
  ggtitle("a")
plota

plotb <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = s), size = 0.8, width = 0.1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = s), size = 0.8, width = 0.1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = s), size = 0.8, width = 0.1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = s), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = s), size = 0.8, width = 0.1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(40,80), xlim = c(1,10)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("s (selection strength)") + 
  ylab("") + 
  ggtitle("b")
plotb

plotc <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = s), size = 0.8, width = 0.1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = s), size = 0.8, width = 0.1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = s), size = 0.8, width = 0.1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = s), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = s), size = 0.8, width = 0.1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(0,15), xlim = c(1,10)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") + 
  ggtitle("c")
plotc

plotd <- ggplot() + 
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = N), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = N), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = N), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = N), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(30,200), xlim = c(0,1000)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") +
  ggtitle("d")
plotd

plote <- ggplot() +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = N), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = N), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = N), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = N), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(25,80), xlim = c(0,1000)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("N (initial population size)") + 
  ylab("") + 
  ggtitle("e")
plote

plotf <- ggplot() +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = N), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = N), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = N), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = N), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(0,15), xlim = c(0,1000)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") + 
  ggtitle("f")
plotf

plotg <- ggplot() + 
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.mutate), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.mutate), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.mutate), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.mutate), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.mutate), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(75,200), xlim = c(0,0.1)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("Species richness (n)") +
  ggtitle("g")
plotg

ploth <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.mutate), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.mutate), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.mutate), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.mutate), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.mutate), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(35,80), xlim = c(0,0.1)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("p.mutate (proportion of mutation)") + 
  ylab("") + 
  ggtitle("h")
ploth

ploti <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.mutate), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.mutate), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.mutate), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.mutate), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.mutate), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(0,15), xlim = c(0,0.1)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") + 
  ggtitle("i")
ploti

plotj <- ggplot() + 
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.rec), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.rec), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.rec), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.rec), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.rec), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(70,200), xlim = c(0.25,0.75)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") +
  ggtitle("j")
plotj

plotk <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.rec), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.rec), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.rec), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.rec), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.rec), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(40,80), xlim = c(0.25,0.75)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("p.rec (proportion of recombination)") + 
  ylab("") + 
  ggtitle("k")
plotk

plotl <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.rec), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.rec), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.rec), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.rec), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.rec), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(0,15), xlim = c(0.25,0.75)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") + 
  ggtitle("l")
plotl

plotm <- ggplot() + 
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.sex), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.sex), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.sex), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Pollinator_n, x = p.sex), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Pollinator_n, x = p.sex), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(70,200), xlim = c(0,0.5)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") +
  ggtitle("m")
plotm

plotn <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.sex), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.sex), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.sex), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Parasitoid_n, x = p.sex), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Parasitoid_n, x = p.sex), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(40,80), xlim = c(0,0.5)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("p.sex (proportion of sexual reproduction)") + 
  ylab("") + 
  ggtitle("n")
plotn

ploto <- ggplot() +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1.5, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.sex), size = 0.8, width = 0.001, colour = "gold", data = mix.sensitivity.plant.poll_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1.5, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.sex), size = 0.8, width = 0.001, colour = "darkblue", data = mix.sensitivity.plant.lm.para_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1.5, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.sex), size = 0.8, width = 0.001, colour = "yellowgreen", data = mix.sensitivity.plant.herb_combined) +
  stat_summary(geom = "line", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "point", fun.y = mean, aes(y = Herbivore_n, x = p.sex), size = 1.5, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, aes(y = Herbivore_n, x = p.sex), size = 0.8, width = 0.001, colour = "purple3", data = mix.sensitivity.plant.multi_combined) +
  coord_cartesian(ylim = c(0,20), xlim = c(0,0.5)) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 12), legend.position = "NA") + 
  xlab("") + 
  ylab("") + 
  ggtitle("o")
ploto

grid.arrange(plota, plotb, plotc, plotd, plote, plotf, plotg, ploth, ploti, plotj, plotk, plotl, plotm, plotn, ploto, nrow = 5)
