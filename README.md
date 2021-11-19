# multilayer-plant-importance
Code and data for calculating the relative importance of plant species in plant-focused (bottom-up) multilayer networks

I hope this is interesting and useful for people! 



Thanks to @tribustale who spotted a mismatch in the code and the data we had supplied. For those who don't want to go to the issues section I will briefly describe the problem and the fix here: 

* The edgelist2matrix() function was selecting the last record for the abundance and eveness based analyses (reported in the supp. mat. of the manuscript) which was therefore not generating an interaction matrix at the farm scale (rather it was selecting an assortment of different interactions from across habitats in Norwood Farm)

* The dataset originally supplied in this repository (and that is still there: nore2aggregated.csv) had the interactions separated out for each habitat, rather than aggregated at the farm-scale 

* I have now uploaded the farm-scale data (see the Dryad page for more information on how the weighted average interaction strengths were derived: https://datadryad.org/stash/dataset/doi:10.5061/dryad.3s36r118)

* There are no changes to the results (i.e., optimal plant mixes) from the species richness-based analyses that we present in the main manuscript

* There are some small changes to the abundance- and eveness-based results that we report in Table S1 of the manuscript. In a handful of the mixes, there is 1 plant substituted from the 5 present in the original results provided in the supp. mat. of the manuscript. 

Hopefully no-one else has any more problems. But if you do have any issues, or want to chat about the code/data, then please get in contact with me either through: fredric.windsor@newcastle.ac.uk or fmwindsor@gmail.com 
