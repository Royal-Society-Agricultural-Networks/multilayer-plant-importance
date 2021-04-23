#### Functions for the calculation of robustness for multilayer networks #### 

# Calculate standard robustness based on primary extinctions 
rob.fn <- function(myorder, myweb, myn.guilds, mycounts) { 
  
  ## MODEL SETUP ## 
  
  # Make arrays to store the results
  myres <- array(NA, c(mycounts[1,3], myn.guilds))
  myres.qcell <- array(NA, c(mycounts[1,3], myn.guilds))
  myres.qmax <- array(NA, c(mycounts[1,3], myn.guilds))
  
  # Initial setup for modelling
  tmp1 <- myweb                                                                # create a unique dataframe for modelling
  init.s.guild <- mycounts[,3]                                                 # record the number of species 
  init.nlinks.cols <- colSums(tmp1[mycounts[1,1]:mycounts[1,2],] > 0)          # record the number of unique links
  init.sumlinks.cols <- colSums(tmp1[mycounts[1,1]:mycounts[1,2],])            # record the total number of links
  init.sumlinks.cols.aphid <- colSums(tmp1[mycounts[3,1]:mycounts[3,2],])      # record the number of links to aphids
  init.sumlinks.cols.seedfeeder <- colSums(tmp1[mycounts[7,1]:mycounts[7,2],]) # record the number of links to seedfeeders
  init.sumlinks.cols.mammal <- colSums(tmp1[mycounts[9,1]:mycounts[9,2],])     # record the number of links to mammals/rodents
  
  # Create vectors for link data
  current.sumlinks.cols <- init.sumlinks.cols                       # create current version for the total number of links
  current.sumlinks.cols.aphid <- init.sumlinks.cols.aphid           # create current version for the number of links to aphids
  current.sumlinks.cols.seedfeeder <- init.sumlinks.cols.seedfeeder # create current version for the number of links to seedfeeders
  current.sumlinks.cols.mammal <- init.sumlinks.cols.mammal         # create current version for the number of links to mammals/rodents
  new.sumlinks.cols <- init.sumlinks.cols                           # set the new total links
  
  # Create arrays for the link data across guilds
  init.nlink.guild <- array(NA, myn.guilds)                                                                # create an array to store unique links for each guild
  for (g in 1:myn.guilds) {init.nlink.guild[g] <- sum(init.nlinks.cols[mycounts[g,1]:mycounts[g,2]])}      # store unique links for each guild
  init.sumlinks.guild <- array(NA, myn.guilds)                                                             # create an array to store total links for each guild
  for (g in 1:myn.guilds) {init.sumlinks.guild[g] <- sum(init.sumlinks.cols[mycounts[g,1]:mycounts[g,2]])} # store total links for each guild
  
  
  ## EXTINCTION PROCESS ## 
  
  # Start with the extinction process (based on specified order - myorder)
  for (i in 1:length(myorder)) { # for each plant in list 'myorder'
    tmp1[myorder[i],] <- 0       # make it extinct (by setting the weighted or binary interactions to zero)
    
    # PRIMARY EXTINCTIONS #
    
    # All guilds directly interacting with plants (not primary or secondary aphid parasitoids, or seed feeder parasitoids/ectoparasites)
    prev.sumlinks.cols <- new.sumlinks.cols                                       # set the previous links to the new links
    new.sumlinks.cols <- colSums(tmp1[mycounts[1,1]:mycounts[1,2],])              # calculate the new links 
    sumlinks.differences <- prev.sumlinks.cols - new.sumlinks.cols                # check the differences between previous and new links
    secondary.extinctions <- array(0, c(length(init.sumlinks.cols)))              # create an array to store secondary extinctions
    secondary.extinctions[new.sumlinks.cols == 0 & sumlinks.differences > 0] <- 1 # store information on secondary extinctions (for removal below)
    current.sumlinks.cols[secondary.extinctions > 0] <- 0                         # set the abundance of disconnected species to zero
    
    # Store the results  
    for (j in c(1:3, 6:12, myn.guilds)) {                                                                               # for all guilds directly interacting with plants
      myres[i,j] <- sum(colSums(tmp1[mycounts[1,1]:mycounts[1,2],mycounts[j,1]:mycounts[j,2]]) > 0) / (init.s.guild[j]) # store count data for standard results
      myres.qmax[i,j] <- sum(current.sumlinks.cols[mycounts[j,1]:mycounts[j,2]]) / init.sumlinks.guild[j]               # store count data for quantitative results
    }
    for (j in c(1:myn.guilds)) {
      myres.qcell[i,j]<-sum(colSums(tmp1[mycounts[1,1]:mycounts[1,2],mycounts[j,1]:mycounts[j,2]])) / init.sumlinks.guild[j] # store count data for qualitative results
    }    
    
    # SECONDARY EXTINCTIONS #
    
    # Aphid parasites (primary and secondary)
    if(sum(secondary.extinctions[mycounts[3,1]:mycounts[3,2]]) > 0) {                            # if there are recorded aphid extinctions
      for (k in which(secondary.extinctions[mycounts[3,1]:mycounts[3,2]] > 0)) {                 # for each aphid extinction
        tmp1[mycounts[3,1] - 1 + k,] <- 0                                                        # set the identified species counts to zero
        prev.sumlinks.cols.aphid <- init.sumlinks.cols.aphid                                     # set the initial total links to the previous total links
        new.sumlinks.cols.aphid <- colSums(tmp1[mycounts[3,1]:mycounts[3,2],])                   # calculate the new links
        sumlinks.differences.aphid <- prev.sumlinks.cols.aphid - new.sumlinks.cols.aphid         # check the differences between previous and new links
        tertiary.extinctions <- array(0, c(length(init.sumlinks.cols)))                          # create an array to store tertiary extinctions
        tertiary.extinctions[new.sumlinks.cols.aphid == 0 & sumlinks.differences.aphid > 0] <- 1 # generate values for tertiary extinctions
        current.sumlinks.cols.aphid[tertiary.extinctions > 0] <- 0                               # set the abundance of disconnected species to zero
      }        
    }
    for (k in c(4)) { # store the results for primary aphid parasites
      myres[i,k] <- sum(colSums(tmp1[mycounts[3,1]:mycounts[3,2],mycounts[k,1]:mycounts[k,2]])>0) / (init.s.guild[k]) # store count data for standard results
      myres.qmax[i,k] <- sum(current.sumlinks.cols.aphid[mycounts[k,1]:mycounts[k,2]]) / init.sumlinks.guild[k]       # store count data for quantitative results
    }
    for (k in c(5)) { # store the results for secondary aphid parasites
      myres[i,k] <- sum(colSums(tmp1[mycounts[3,1]:mycounts[3,2],mycounts[k,1]:mycounts[k,2]]) > 0) / (init.s.guild[k]) # store count data for standard results
      myres.qmax[i,k] <- sum(current.sumlinks.cols.aphid[mycounts[k,1]:mycounts[k,2]]) / init.sumlinks.guild[k]         # store count data for quantitative results
    }
    
    # Seed-feeder parasites
    if(sum(secondary.extinctions[mycounts[7,1]:mycounts[7,2]]) > 0) {                                      # if there are recorded seed-feeder extinctions                                     
      for (k in which(secondary.extinctions[mycounts[7,1]:mycounts[7,2]] > 0)) {                            # for each seed-feeder exinction
        tmp1[mycounts[7,1] - 1 + k,] <- 0                                                                     # set the identified species counts to zero 
        prev.sumlinks.cols.seedfeeder <- init.sumlinks.cols.seedfeeder                                      # set the initial total links to the previous total links
        new.sumlinks.cols.seedfeeder <- colSums(tmp1[mycounts[7,1]:mycounts[7,2],])                         # calculate the new links
        sumlinks.differences.seedfeeder <- prev.sumlinks.cols.seedfeeder - new.sumlinks.cols.seedfeeder     # check the differences between previous and new links
        tertiary.extinctions <- array(0, c(length(init.sumlinks.cols)))                                     # create an array to store tertiary extinctions
        tertiary.extinctions[new.sumlinks.cols.seedfeeder == 0 & sumlinks.differences.seedfeeder > 0] <- 1  # generate values for tertiary extinctions
        current.sumlinks.cols.seedfeeder[tertiary.extinctions > 0] <- 0                                     # set the abundance of disconnected species to zero
      }        
    }
    for (k in c(13)) { # store the results for seed-feeder parasites 
      myres[i,k] <- sum(colSums(tmp1[mycounts[7,1]:mycounts[7,2],mycounts[k,1]:mycounts[k,2]]) > 0) / (init.s.guild[k]) # store count data for standard results
      myres.qmax[i,k] <- sum(current.sumlinks.cols.seedfeeder[mycounts[k,1]:mycounts[k,2]]) / init.sumlinks.guild[k]    # store count data for quantitative results
    }
    
    # Mammal ectoparasites (fleas)
    if(sum(secondary.extinctions[mycounts[9,1]:mycounts[9,2]]) > 0) {                               # if there are recorded mammal extinctions
      for (k in which(secondary.extinctions[mycounts[9,1]:mycounts[9,2]] > 0)) {                    # for each mammal exinction
        tmp1[mycounts[9,1] - 1 + k,] <- 0                                                             # set the identified species counts to zero 
        prev.sumlinks.cols.mammal <- init.sumlinks.cols.mammal                                      # set the initial total links to the previous total links
        new.sumlinks.cols.mammal <- colSums(tmp1[mycounts[9,1]:mycounts[9,2],])                     # calculate the new links
        sumlinks.differences.mammal <- prev.sumlinks.cols.mammal - new.sumlinks.cols.mammal         # check the differences between previous and new links
        tertiary.extinctions <- array(0, c(length(init.sumlinks.cols)))                             # create an array to store tertiary extinctions
        tertiary.extinctions[new.sumlinks.cols.mammal == 0 & sumlinks.differences.mammal > 0] <- 1  # generate values for tertiary extinctions
        current.sumlinks.cols.mammal[tertiary.extinctions > 0] <- 0                                 # set the abundance of disconnected species to zero
      }        
    }
    for (k in c(14)) { # store the results for mammal ectoparasites
      myres[i,k] <- sum(colSums(tmp1[mycounts[9,1]:mycounts[9,2],mycounts[k,1]:mycounts[k,2]]) > 0) / (init.s.guild[k]) # store count data for standard results
      myres.qmax[i,k] <- sum(current.sumlinks.cols.mammal[mycounts[k,1]:mycounts[k,2]]) / init.sumlinks.guild[k]        # store count data for quantitative results
    }
  }
  list(myres, myres.qmax, myres.qcell, init.s.guild, mycounts[,3]) # list the results 
}   


# Creating xtab count data
get.counts<-function(myxtab) {
  summarynames<-substr(colnames(myxtab),1,4)
  richnessofeachguild<-table(summarynames)
  cumulativerichnessofeachguild<-cumsum(richnessofeachguild)
  startcounts<-c(1,cumulativerichnessofeachguild[1:length(cumulativerichnessofeachguild)-1]+1)
  endcounts<-c(cumulativerichnessofeachguild)
  lengthcounts<-endcounts-startcounts+1
  myguilds<-unique(substr(colnames(myxtab),1,4))
  counts<-array(c(startcounts,endcounts,lengthcounts),dim=c(length(myguilds),3))
  rownames(counts)<-myguilds
  colnames(counts)<-c("first","last","number")
  counts
}

# Calculate area under the curve (R) for species extinction scenarios
robustness.auc <- function(myarray, pl.or.fv, mycounts) {
  if (pl.or.fv == "pl") mylength <- mycounts[1,3] - 1     # for plant extinction processes
  if (pl.or.fv == "fv") mylength <- mycounts[2,3] - 1     # for flower visitor extinction processes
  if (length(dim(myarray))==2) {                          # if there are more than 2 dimensions
    res.auc <- array(NA, c(1, n.guilds))                  # create the area under the curve results array
    res.auc <- apply(myarray / mylength, 2, sum)          # store the data for each iteration
  } else {                                                # if there are not 2 dimensions 
    res.auc <- array(NA, c(n.iter, n.guilds))             # create the area under the curve results array
    for (guild in 1:n.guilds){                            # for each guild in the network
      m2 <- myarray[,guild,]                              # store the results array from each guild in the robustness model
      res.auc[,guild] <- apply(m2[,] / mylength, 2, sum)  # store the data for each iteration
    }
  }
  res.auc
}

# Calculate median values for numerous iterations 
robustness.median<-function(myarray,myguildnames) {
  res.med<-array(NA,c(dim(myarray)[1],length(myguildnames)))
  colnames(res.med)<-myguildnames
  for (guild in 1:length(myguildnames)){
    m2<-myarray[,guild,]
    res.med[,guild]<-apply(m2[,],1,median)
  }
  res.med
}

# Calculate the upper 90% values for numerous iterations
robustness.upper90<-function(myarray,myguildnames) {
  funct.upper90<-function(x) quantile(x,probs=0.95)
  res.med<-array(NA,c(dim(myarray)[1],length(myguildnames)))
  colnames(res.med)<-myguildnames
  for (guild in 1:length(myguildnames)){
    m2<-myarray[,guild,]
    res.med[,guild]<-apply(m2[,],1,funct.upper90)
  }
  res.med
}

# Calculate the lower 90% values for numerous iteractions
robustness.lower90<-function(myarray,myguildnames) {
  funct.lower90<-function(x) quantile(x,probs=0.05)
  res.med<-array(NA,c(dim(myarray)[1],length(myguildnames)))
  colnames(res.med)<-myguildnames
  for (guild in 1:length(myguildnames)){
    m2<-myarray[,guild,]
    res.med[,guild]<-apply(m2[,],1,funct.lower90)
  }
  res.med
}

# Calculate the upper 95% quantile
funct.upper95<-function(x) quantile(x,probs=0.975, na.rm = T)

# Calculate the lower 95% quantile
funct.lower95<-function(x) quantile(x,probs=0.025, na.rm = T)