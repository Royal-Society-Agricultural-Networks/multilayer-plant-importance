# Methods for the selection of the most important plants in the Norwood farm network #
# Code collated and adapted by Fred M. Windsor (fredric.windsor@newcastle.ac.uk)     #
# The origins of all code are referenced in line and in the main text of the study   #



#### SET-UP #### 

rm(list=ls())
setwd("")



#### LIBRARIES #### 

library(magrittr); library(bipartite); library(tidyverse); library(igraph); library(gridExtra)



#### FUNCTIONS ####

source("Functions/network_functions.R") # Adapted from Pilosof et al. (2017)
source("Functions/ga.R") # From M'Gonigle et al. (2017)
source("Functions/control.R") # From M'Gonigle et al. (2017)
source("Functions/objective_functions.R") # From M'Gonigle et al. (2017)
source("Functions/stability_functions.R") # From Sauve et al. (2016)




#### DATA INPUT ####

# Read in the data for the norwood farm network
dframe1 <- read.csv("Data/nore2_aggregated.csv", header = T)

# Sort out the shorthand guild names for the following manipulations
dframe1$lower.names <- substr(dframe1$lower, 1, 4)
dframe1$upper.names <- substr(dframe1$upper, 1, 4)

# Subset the data into different groups (plant-pollinators and plant-herbivore-parasitoids)
plant.poll_edgelist <- dframe1 %>% filter(upper.names == "02FV" | upper.names == "12BF")
plant.herb_edgelist <- dframe1 %>% filter(upper.names == "03AP")
plant.lm.para_edgelist <- dframe1 %>% filter(upper.names == "06MI")
#herb.para_edgelist <- dframe1 %>% filter(upper.names == "04PR" | upper.names == "05SE")

# Remove non plant-polinator interactions for ease 
plant.poll_edgelist <- subset(plant.poll_edgelist, lower.names == "01PL")
plant.herb_edgelist <- subset(plant.herb_edgelist, lower.names == "01PL")
plant.lm.para_edgelist <- subset(plant.lm.para_edgelist, lower.names == "01PL")
#herb.para_edgelist <- subset(plant.lm.para_edgelist, lower.names != "01PL")
multilayer_edgelist <- rbind(plant.poll_edgelist, plant.herb_edgelist, plant.lm.para_edgelist)

# Convert the edgelists into matrices
plant.poll_mat <- edgelist2Matrix(plant.poll_edgelist)
plant.herb_mat <- edgelist2Matrix(plant.herb_edgelist)
plant.lm.para_mat <- edgelist2Matrix(plant.lm.para_edgelist)
#herb.para_mat <- edgelist2Matrix(herb.para_edgelist)
multilayer_mat <- edgelist2Matrix(multilayer_edgelist)

# Transpose matrices (required that plants are columns for analyses)
plant.poll_mat <- t(plant.poll_mat)
plant.herb_mat <- t(plant.herb_mat)
plant.lm.para_mat <- t(plant.lm.para_mat)
#herb.para_mat <- t(herb.para_mat)
multilayer_mat <- t(multilayer_mat)




#### GENETIC ALGORITHM ####

### PLANT-POLLINATOR NETWORK ### 

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

# Write the data into .csv files
write.csv(mixes_abundance_weighted, "Tables/Plant-mixes_abundance_k5.csv")
write.csv(mixes_richness_weighted, "Tables/Plant-mixes_richness_k5.csv")
write.csv(mixes_both_weighted, "Tables/Plant-mixes_both_k5.csv")



### MULTILAYER NETWORK ### 

## BINARY ##

# Store the matrices as "m1", "m2" and "m3"
m1 <- plant.poll_mat
m2 <- plant.lm.para_mat
m3 <- plant.herb_mat

# Optimising "richness" for 2 species mixes (additive method)
mix.r.multilayer.add2 <- find.mix(f=multilayer.additive.richness, k=2)

# Optimising "richness" for 2 species mixes (multiplicative method)
mix.r.multilayer.multi2 <- find.mix(f=multilayer.multiplicative.richness, k=2)

# Optimising "richness" for 5 species mixes (additive method)
mix.r.multilayer.add5 <- find.mix(f=multilayer.additive.richness, k=5)

# Optimising "richness" for 5 species mixes (multiplicative method)
mix.r.multilayer.multi5 <- find.mix(f=multilayer.multiplicative.richness, k=5)

# Optimising "richness" for 10 species mixes (additive method)
mix.r.multilayer.add10 <- find.mix(f=multilayer.additive.richness, k=10)

# Optimising "richness" for 10 species mixes (multiplicative method)
mix.r.multilayer.multi10 <- find.mix(f=multilayer.multiplicative.richness, k=10)

# Optimising "abundance" for 2 species mixes (additive method)
mix.a.multilayer.add2 <- find.mix(f=multilayer.additive.abundance, k=2)

# Optimising "abundance" for 2 species mixes (multiplicative method)
mix.a.multilayer.multi2 <- find.mix(f=multilayer.multiplicative.abundance, k=2)

# Optimising "abundance" for 5 species mixes (additive method)
mix.a.multilayer.add5 <- find.mix(f=multilayer.additive.abundance, k=5)

# Optimising "abundance" for 5 species mixes (multiplicative method)
mix.a.multilayer.multi5 <- find.mix(f=multilayer.multiplicative.abundance, k=5)

# Optimising "abundance" for 10 species mixes (additive method)
mix.a.multilayer.add10 <- find.mix(f=multilayer.additive.abundance, k=10)

# Optimising "abundance" for 10 species mixes (multiplicative method)
mix.a.multilayer.multi10 <- find.mix(f=multilayer.multiplicative.abundance, k=10)

# Optimising "richness.abundance" for 2 species mixes (additive method)
mix.ar.multilayer.add2 <- find.mix(f=multilayer.additive.both, k=2)

# Optimising "richness.abundance" for 2 species mixes (multiplicative method)
mix.ar.multilayer.multi2 <- find.mix(f=multilayer.multiplicative.both, k=2)

# Optimising "richness.abundance" for 5 species mixes (additive method)
mix.ar.multilayer.add5 <- find.mix(f=multilayer.additive.both, k=5)

# Optimising "richness.abundance" for 5 species mixes (multiplicative method)
mix.ar.multilayer.multi5 <- find.mix(f=multilayer.multiplicative.both, k=5)

# Optimising "richness.abundance" for 10 species mixes (additive method)
mix.ar.multilayer.add10 <- find.mix(f=multilayer.additive.both, k=10)

# Optimising "abundance" for 10 species mixes (multiplicative method)
mix.ar.multilayer.multi10 <- find.mix(f=multilayer.multiplicative.both, k=10)




#### NETWORK PROPERTIES OF PLANT MIXTURES ####

### SPECIES RICHNESS ###

## Pollinator plant mix

# Pollinator species richness
table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,1]])>0)["TRUE"]

# Parasitoid species richness
table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,1]])>0)["TRUE"]

# Herbivore species richness
table(rowSums(multilayer_mat[258:380,mixes_both_weighted[,1]])>0)["TRUE"]


## Parasitoid plant mix

# Pollinator species richness
table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,2]])>0)["TRUE"]

# Parasitoid species richness
table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,2]])>0)["TRUE"]

# Herbivore species richness
table(rowSums(multilayer_mat[258:380,mixes_both_weighted[,2]])>0)["TRUE"]


## Herbivore plant mix

# Pollinator species richness
table(rowSums(multilayer_mat[1:257,mixes_both_weighted[,3]])>0)["TRUE"]

# Parasitoid species richness
table(rowSums(multilayer_mat[286:381,mixes_both_weighted[,3]])>0)["TRUE"]

# Herbivore species richness
table(rowSums(multilayer_mat[258:380,mixes_both_weighted[,3]])>0)["TRUE"]


## Multilayer (additive) plant mix

# Pollinator species richness
table(rowSums(multilayer_mat[1:257,mix.r.multilayer.add])>0)["TRUE"]

# Parasitoid species richness
table(rowSums(multilayer_mat[286:381,mix.r.multilayer.add])>0)["TRUE"]

# Herbivore species richness
table(rowSums(multilayer_mat[258:380,mix.r.multilayer.add])>0)["TRUE"]


## Multilayer (multiplicative) plant mix

# Pollinator species richness
table(rowSums(multilayer_mat[1:257,mix.r.multilayer.multi])>0)["TRUE"]

# Parasitoid species richness
table(rowSums(multilayer_mat[286:381,mix.r.multilayer.multi])>0)["TRUE"]

# Herbivore species richness
table(rowSums(multilayer_mat[258:380,mix.r.multilayer.multi])>0)["TRUE"]



### DEGREE ### 

## Pollinator plant mix 
mean(degree((graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]])))))

## Parasitoid plant mix
mean(degree((graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,2]])))))

## Herbivore plant mix
mean(degree((graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,3]])))))

## Multilayer (additive) plant mix
mean(degree((graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.add5])))))

## Multilayer (multiplicative) plant mix
mean(degree((graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.multi5])))))



### CONNECTANCE ### 

## Pollinator plant mix 
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]])))

## Parasitoid plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,2]])))

## Herbivore plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,3]])))

## Multilayer (additive) plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.add])))

## Multilayer (multiplicative) plant mix
edge_density(graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.multi])))



### ROBUSTNESS ###

## Generate networks to work from

# Subsets of the multilayer network
rob.plant.poll_mat <- multilayer_mat[1:257,]
rob.plant.herb_mat <- multilayer_mat[258:380,]
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


## Multilayer (multiplicative) plant mix

# Remove plants that are not shared between layers
commonplants.multi <- mix.r.multilayer.multi5
plant.poll.multi_mat <- rob.plant.poll_mat[,commonplants.multi]
plant.herb.multi_mat <- rob.plant.herb_mat[,commonplants.multi]
plant.lm.para.multi_mat <- rob.plant.lm.para_mat[,commonplants.multi]

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
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll.multi_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll.multi_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll.multi_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll.multi_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb.multi_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para.multi_mat # make copy of original plant-leaf miner parasitoid matrix
  
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
multilayer.poll_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]]))))

# Create a Jacobian matrix
multilayer.poll_jac.mat <- jacobian_binary(multilayer.poll_adj.mat)

# Calculate the stability
stability(multilayer.poll_jac.mat, s2 = 13)


## Parasitoid plant mix

# Create adjacency matrix
multilayer.para_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,2]]))))

# Create a Jacobian matrix
multilayer.para_jac.mat <- jacobian_binary(multilayer.para_adj.mat)

# Calculate the stability
stability(multilayer.para_jac.mat, s2 = 12)


## Herbivore plant mix

# Create adjacency matrix
multilayer.herb_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,3]]))))

# Create a Jacobian matrix
multilayer.herb_jac.mat <- jacobian_binary(multilayer.herb_adj.mat)

# Calculate the stability
stability(multilayer.herb_jac.mat, s2 = 10)


## Multilayer (additive) plant mix

# Create adjacency matrix
multilayer.add_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.add5]))))

# Create a Jacobian matrix
multilayer.add_jac.mat <- jacobian_binary(multilayer.add_adj.mat)

# Calculate the stability
stability(multilayer.add_jac.mat, s2 = 8)


## Multilayer (multiplicative) plant mix

# Create adjacency matrix
multilayer.multi_adj.mat <- as.matrix(as_adjacency_matrix(graph_from_incidence_matrix(empty(multilayer_mat[,mix.r.multilayer.multi5]))))

# Create a Jacobian matrix
multilayer.multi_jac.mat <- jacobian_binary(multilayer.multi_adj.mat)

# Calculate the stability
stability(multilayer.multi_jac.mat, s2 = 5)


### PLOTTING A STABLE AND UNSTABLE NETWORK (FIGURE 1)

## uNSTABLE NETWORK

# Pollinator plant mix
plant.poll_graph <- graph_from_incidence_matrix(empty(multilayer_mat[,mixes_both_weighted[,1]]))

# Set the node properties
V(plant.poll_graph)$code <- substr(names(V(plant.poll_graph)), 1, 4)
V(plant.poll_graph)$color <- V(plant.poll_graph)$code

# Set the colours to the appropriate values
V(plant.poll_graph)$color <- gsub("01PL", "darkgreen", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("03AP", "yellowgreen", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("06MI", "darkblue", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("02FV", "gold", V(plant.poll_graph)$color)
V(plant.poll_graph)$color <- gsub("12BF", "gold", V(plant.poll_graph)$color)

# Plot the network
plot(plant.poll_graph, vertex.label = NA, vertex.size = 6, vertex.color = V(plant.poll_graph)$color)
legend("topleft", legend = c("Pollinators", "Herbivores", "Parasitoids", "Plants"), 
       col = unique(V(plant.poll_graph)$color), bty = "n", pch = 20, pt.cex = 3, cex = 1, horiz = F, inset = c(0.005,0.005))


## STABLE NETWORK 

# Multilayer plant mix 
plant.multilayer_graph <- graph_from_incidence_matrix(empty(multilayer_mat[,mix.ar.multilayer.multi5]))

# Set the node properties
V(plant.multilayer_graph)$code <- substr(names(V(plant.multilayer_graph)), 1, 4)
V(plant.multilayer_graph)$color <- V(plant.multilayer_graph)$code

# Set the colours to the appropriate values
V(plant.multilayer_graph)$color <- gsub("01PL", "darkgreen", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("03AP", "yellowgreen", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("06MI", "darkblue", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("02FV", "gold", V(plant.multilayer_graph)$color)
V(plant.multilayer_graph)$color <- gsub("12BF", "gold", V(plant.multilayer_graph)$color)

# Plot the network
plot(plant.multilayer_graph, vertex.label = NA, vertex.size = 6, vertex.color = V(plant.multilayer_graph)$color)
legend("topleft", legend = c("Pollinators", "Herbivores", "Parasitoids", "Plants"), 
       col = unique(V(plant.multilayer_graph)$color), bty = "n", pch = 20, pt.cex = 3, cex = 1, horiz = F, inset = c(0.005,0.005))




#### SIMULATING INCREASES AND DECREASES IN IMPORTANT PLANT SPECIES ####

### SCENARIO 0 - NO CHANGE TO THE ABUNDANCE OF PLANTS IMPORTANT FOR ANY PROCESS

## Simulate the robustness of the network to random species removals

# Remove plants that are not shared between layers
commonplants <- Reduce(intersect, list(colnames(plant.poll_mat),colnames(plant.herb_mat),colnames(plant.lm.para_mat)))
plant.poll_mat <- plant.poll_mat[,commonplants]
plant.herb_mat <- plant.herb_mat[,commonplants]
plant.lm.para_mat <- plant.lm.para_mat[,commonplants]

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
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll_mat      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb_mat      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para_mat # make copy of original plant-leaf miner parasitoid matrix
 
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

# Histograms of R values for groups
plot1a <- ggplot(aes(x = poll_Rvalues), data = data.frame(poll_Rvalues)) + geom_histogram(color = "black") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(poll_Rvalues), color = "red") + ggtitle("Pollinators")
plot1b <- ggplot(aes(x = herb_Rvalues), data = data.frame(herb_Rvalues)) + geom_histogram(color = "black") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(herb_Rvalues), color = "red") + ggtitle("Herbivores")
plot1c <- ggplot(aes(x = lmpara_Rvalues), data = data.frame(lmpara_Rvalues)) + geom_histogram(color = "black") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(lmpara_Rvalues), color = "red") +  ggtitle("Parasitoids")

# Arrange the plots
grid.arrange(plot1a, plot1b, plot1c, nrow = 1)



### SCENARIO 1 - INCREASE THE ABUNDANCE OF PLANTS IMPORTANT FOR POLLINATORS

## Increase the relative abundance of important plant species for pollinators

# Calculate the proprortion of important plant species within the data
plant_total.abun <- sum(colSums(multilayer_mat))
plant.poll_imp.abun <- sum(colSums(multilayer_mat[,imp.plant.poll_names]))
round(plant.poll_imp.abun/plant_total.abun,2) # Important plants account for 24% of the total plant abundance

# Create a matrix to create an altered proportion of important plants
multilayer_mat.poll.inc <- multilayer_mat

# Increase the abundance of important plants for pollinators by +10%
multilayer.poll_inc <- multilayer_mat[,imp.plant.poll_names] * 1.7 # Generates a 10% proportional increase in the important plants for pollinators
multilayer_mat.poll.inc[,imp.plant.poll_names] <- multilayer.poll_inc

# Check that a 10% increase in abundance has been achieved
multilayer_total.abun.inc <- sum(colSums(multilayer_mat.poll.inc))
multilayer_imp.abun.inc <- sum(colSums(multilayer_mat.poll.inc[,imp.plant.poll_names]))
round(multilayer_imp.abun.inc/multilayer_total.abun.inc,2) # Important plants account for 34% of the total plant abundance


## Simulate the robustness of the network to random species removals

# Split out the individual layers from the multilayer network
plant.poll_mat.poll.inc <- multilayer_mat.poll.inc[1:257,]
plant.herb_mat.poll.inc <- multilayer_mat.poll.inc[258:285,]
plant.lm.para_mat.poll.inc <- multilayer_mat.poll.inc[286:381,]

# Remove plants that are not shared between layers
plant.poll_mat.poll.inc <- plant.poll_mat.poll.inc[,commonplants]
plant.herb_mat.poll.inc <- plant.herb_mat.poll.inc[,commonplants]
plant.lm.para_mat.poll.inc <- plant.lm.para_mat.poll.inc[,commonplants]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_poll.inc_Rvalues <- NULL
herb_poll.inc_Rvalues <- NULL
lmpara_poll.inc_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll_mat.poll.inc      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb_mat.poll.inc      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para_mat.poll.inc # make copy of original plant-leaf miner parasitoid matrix
  
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
  
  poll_poll.inc_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_poll.inc_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_poll.inc_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_poll.inc_Rvalues)
median(herb_poll.inc_Rvalues)
median(lmpara_poll.inc_Rvalues)

# Percentage change in the robustness (R values)
1-(median(poll_Rvalues)/median(poll_poll.inc_Rvalues))
1-(median(herb_Rvalues)/median(herb_poll.inc_Rvalues))
1-(median(lmpara_Rvalues)/median(lmpara_poll.inc_Rvalues))

# Histograms of R values for groups
plot2a <- ggplot(aes(x = poll_poll.inc_Rvalues), data = data.frame(poll_poll.inc_Rvalues)) + geom_histogram(color = "black", fill = "gold") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(poll_poll.inc_Rvalues), color = "red") + ggtitle("")
plot2b <- ggplot(aes(x = herb_poll.inc_Rvalues), data = data.frame(herb_poll.inc_Rvalues)) + geom_histogram(color = "black", fill = "gold") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(herb_poll.inc_Rvalues), color = "red") + ggtitle("")
plot2c <- ggplot(aes(x = lmpara_poll.inc_Rvalues), data = data.frame(lmpara_poll.inc_Rvalues)) + geom_histogram(color = "black", fill = "gold") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(lmpara_poll.inc_Rvalues), color = "red") +  ggtitle("")

# Arrange the plots
grid.arrange(plot2a, plot2b, plot2c, nrow = 1)



### SCENARIO 2 - INCREASE THE ABUNDANCE OF PLANTS IMPORTANT FOR PARASITOIDS

## Increase the relative abundance of important plant species for parasitoids

# Calculate the proprortion of important plant species within the data
plant.lmpara_imp.abun <- sum(colSums(multilayer_mat[,imp.plant.lmpara_names]))
round(plant.lmpara_imp.abun/plant_total.abun,2) # Important plants account for 37% of the total plant abundance

# Create a matrix to create an altered proportion of important plants
multilayer_mat.lmpara.inc <- multilayer_mat

# Increase the abundance of important plants for pollinators by +10%
multilayer.lmpara_inc <- multilayer_mat[,imp.plant.lmpara_names] * 1.5 # Generates a 10% proportional increase in the important plants for pollinators
multilayer_mat.lmpara.inc[,imp.plant.lmpara_names] <- multilayer.lmpara_inc

# Check that a 10% increase in abundance has been achieved
multilayer_total.abun.inc <- sum(colSums(multilayer_mat.lmpara.inc))
multilayer_imp.abun.inc <- sum(colSums(multilayer_mat.lmpara.inc[,imp.plant.lmpara_names]))
round(multilayer_imp.abun.inc/multilayer_total.abun.inc,2) # Important plants account for 47% of the total plant abundance


## Simulate the robustness of the network to random species removals

# Split out the individual layers from the multilayer network
plant.poll_mat.lmpara.inc <- multilayer_mat.lmpara.inc[1:257,]
plant.herb_mat.lmpara.inc <- multilayer_mat.lmpara.inc[258:285,]
plant.lm.para_mat.lmpara.inc <- multilayer_mat.lmpara.inc[286:381,]

# Remove plants that are not shared between layers
plant.poll_mat.lmpara.inc <- plant.poll_mat.lmpara.inc[,commonplants]
plant.herb_mat.lmpara.inc <- plant.herb_mat.lmpara.inc[,commonplants]
plant.lm.para_mat.lmpara.inc <- plant.lm.para_mat.lmpara.inc[,commonplants]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_lmpara.inc_Rvalues <- NULL
herb_lmpara.inc_Rvalues <- NULL
lmpara_lmpara.inc_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll_mat.lmpara.inc      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb_mat.lmpara.inc      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para_mat.lmpara.inc # make copy of original plant-leaf miner parasitoid matrix
  
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
  
  poll_lmpara.inc_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_lmpara.inc_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_lmpara.inc_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_lmpara.inc_Rvalues)
median(herb_lmpara.inc_Rvalues)
median(lmpara_lmpara.inc_Rvalues)

# Percentage change in the robustness (R values)
1-(median(poll_Rvalues)/median(poll_lmpara.inc_Rvalues))
1-(median(herb_Rvalues)/median(herb_lmpara.inc_Rvalues))
1-(median(lmpara_Rvalues)/median(lmpara_lmpara.inc_Rvalues))

# Histograms of R values for groups
plot3a <- ggplot(aes(x = poll_lmpara.inc_Rvalues), data = data.frame(poll_lmpara.inc_Rvalues)) + geom_histogram(color = "black", fill = "orange") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(poll_lmpara.inc_Rvalues), color = "red") + ggtitle("")
plot3b <- ggplot(aes(x = herb_lmpara.inc_Rvalues), data = data.frame(herb_lmpara.inc_Rvalues)) + geom_histogram(color = "black", fill = "orange") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(herb_lmpara.inc_Rvalues), color = "red") + ggtitle("")
plot3c <- ggplot(aes(x = lmpara_lmpara.inc_Rvalues), data = data.frame(lmpara_lmpara.inc_Rvalues)) + geom_histogram(color = "black", fill = "orange") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(lmpara_lmpara.inc_Rvalues), color = "red") +  ggtitle("")

# Arrange the plots
grid.arrange(plot3a, plot3b, plot3c, nrow = 1)



### SCENARIO 3 - DECREASE THE ABUNDANCE OF PLANTS IMPORTANT FOR HERBIVORES

## Increase the relative abundance of important plant species for herbivores

# Calculate the proprortion of important plant species within the data
plant.herb_imp.abun <- sum(colSums(multilayer_mat[,imp.plant.herb_names]))
round(plant.herb_imp.abun/plant_total.abun,2) # Important plants account for 16% of the total plant abundance

# Create a matrix to create an altered proportion of important plants
multilayer_mat.herb.inc <- multilayer_mat

# Remove all plants that are important for herbivores
multilayer.herb_inc <- multilayer_mat[,imp.plant.herb_names] * 0 # Generates a loss of the important plants for herbivores
multilayer_mat.herb.inc[,imp.plant.herb_names] <- multilayer.herb_inc

# Check that a 10% increase in abundance has been achieved
multilayer_total.abun.inc <- sum(colSums(multilayer_mat.herb.inc))
multilayer_imp.abun.inc <- sum(colSums(multilayer_mat.herb.inc[,imp.plant.herb_names]))
round(multilayer_imp.abun.inc/multilayer_total.abun.inc,2) # Important plants account for 0% of the total plant abundance


## Simulate the robustness of the network to random species removals

# Split out the individual layers from the multilayer network
plant.poll_mat.herb.inc <- multilayer_mat.herb.inc[1:257,]
plant.herb_mat.herb.inc <- multilayer_mat.herb.inc[258:285,]
plant.lm.para_mat.herb.inc <- multilayer_mat.herb.inc[286:381,]

# Remove plants that are not shared between layers
plant.poll_mat.herb.inc <- plant.poll_mat.herb.inc[,commonplants]
plant.herb_mat.herb.inc <- plant.herb_mat.herb.inc[,commonplants]
plant.lm.para_mat.herb.inc <- plant.lm.para_mat.herb.inc[,commonplants]

# Choose number of simulations:
nperm = 1000

# Choose threshold value (0-1)
threshold = 0.5

# Set the starting values and storage objects
poll_herb.inc_Rvalues <- NULL
herb_herb.inc_Rvalues <- NULL
lmpara_herb.inc_Rvalues <- NULL

triggertally <- NULL
cascadelength <- NULL
c=1
trigmat<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
j=0
plantsurvivorsC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))
exttimesC<-matrix(nrow=nperm, ncol=ncol(plant.poll_mat))

# Start the permutations of the weighted model
for(k in 1:nperm){
  
  pl.p_mat<-plant.poll_mat.herb.inc      # make copy of original plant-pollinator matrix
  pl.h_mat<-plant.herb_mat.herb.inc      # make copy of original plant-herbivore matrix
  pl.lmp_mat<-plant.lm.para_mat.herb.inc # make copy of original plant-leaf miner parasitoid matrix
  
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
  
  poll_herb.inc_Rvalues[k]<-sum(poll_survivorsx)/(ncol(pl.p_mat)*nrow(pl.p_mat))         # Record the robustness value for each permutation (survivors)
  herb_herb.inc_Rvalues[k]<-sum(herb_survivorsx)/(ncol(pl.h_mat)*nrow(pl.h_mat))         # Record the robustness value for each permutation (survivors)
  lmpara_herb.inc_Rvalues[k]<-sum(lmpara_survivorsx)/(ncol(pl.lmp_mat)*nrow(pl.lmp_mat)) # Record the robustness value for each permutation (survivors)
  
}

# Median values for groups
median(poll_herb.inc_Rvalues)
median(herb_herb.inc_Rvalues)
median(lmpara_herb.inc_Rvalues)

# Percentage change in the robustness (R values)
1-(median(poll_Rvalues)/median(poll_herb.inc_Rvalues))
1-(median(herb_Rvalues)/median(herb_herb.inc_Rvalues))
1-(median(lmpara_Rvalues)/median(lmpara_herb.inc_Rvalues))

# Histograms of R values for groups
plot4a <- ggplot(aes(x = poll_herb.inc_Rvalues), data = data.frame(poll_herb.inc_Rvalues)) + geom_histogram(color = "black", fill = "yellowgreen") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(poll_herb.inc_Rvalues), color = "red") + ggtitle("")
plot4b <- ggplot(aes(x = herb_herb.inc_Rvalues), data = data.frame(herb_herb.inc_Rvalues)) + geom_histogram(color = "black", fill = "yellowgreen") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(herb_herb.inc_Rvalues), color = "red") + ggtitle("")
plot4c <- ggplot(aes(x = lmpara_herb.inc_Rvalues), data = data.frame(lmpara_herb.inc_Rvalues)) + geom_histogram(color = "black", fill = "yellowgreen") + theme_bw() + scale_x_continuous(limits = c(0,1)) + geom_vline(xintercept = median(lmpara_herb.inc_Rvalues), color = "red") +  ggtitle("")

# Arrange the plots
grid.arrange(plot4a, plot4b, plot4c, nrow = 1)


## Overall plots
grid.arrange(plot1a, plot1b, plot1c, plot2a, plot2b, plot2c, plot3a, plot3b, plot3c, plot4a, plot4b, plot4c, nrow = 4)




#### NETWORK VISUALISATION ####

# Convert the matrix to a graph object 
multilayer_net <- graph_from_incidence_matrix(multilayer_mat)

# Change some of the vertex attributes
V(multilayer_net)$code <- substr(names(V(multilayer_net)), 1, 4)
V(multilayer_net)$color <- V(multilayer_net)$code

# Set the colours to the appropriate values
V(multilayer_net)$color <- gsub("01PL", "darkgreen", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("03AP", "yellowgreen", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("06MI", "darkblue", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("04PR", "orange", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("05SE", "darkorange", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("02FV", "gold", V(multilayer_net)$color)
V(multilayer_net)$color <- gsub("12BF", "gold", V(multilayer_net)$color)

# Plot the network
plot(multilayer_net, vertex.label = NA, vertex.size = 6, vertex.color = V(multilayer_net)$color)
legend("topleft", legend = c("Plants", "Aphids", "Primary aphid parasitoids", "Secondary aphid parasitoids", "Flower visitors", "Leaf miner parasitoids"), 
       col = unique(V(multilayer_net)$color), bty = "n", pch = 20, pt.cex = 3, cex = 1, horiz = F, inset = c(0.005,0.005))
