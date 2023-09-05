### This code was written by Tinka Bakker and Ameen Abu-Hanna as part of the SIMPLIFY project at Amsterdam UMC NL
### This code is used for a research paper by T. Bakker & J.E. Klopotowska et al. 
### The research paper has been submitted to a journal and has been shared on github as part of the Data Sharing Agreement stated in the paper


############### reading packages ################# 
# packages used 
install.packages('tidyverse')
install.packages('plyr')
install.packages('data.table')
library(tidyverse)
library(plyr)
library(data.table)

############### reading data ################# 

# set your working directory to the pathway where the data is stored, use setwd()

# two datasets should be imported:
# one with the drug administrations, called drug_data
# and one with all pDDI combinations, called ddi_data

# to read the example data run:
drug_data = read.csv2('drug_data.csv')
drug_data$time = as.POSIXct(drug_data$time, format="%d-%m-%Y %H:%M") # make date variable into a POSIXct format. Mind the format of your variable and change format in the code if necessary

# and run:
ddi_data = read.csv2('ddi_data.csv')

# to optimize calculations put it in this form
ddi_data = as.data.table(ddi_data)
setkey(ddi_data, drug_a)

# if you want to use the algorithm on other data than the provided example data, make sure that it has the same form. 
# in drug_data, each row should represent either the start of a drug administration, or the stop of a drug administration

# each row should contain 4 variables: 
# 1: admission_id; this should be a number unique for an admission (numeric)
# 2: drug_id; this number should indicate a specific drug and should correspond to your pDDI dataset (numeric)
# 3: time; the date and time of the drug administration, make sure it's of the type POSIXct (POSIXct)
# 4: type; indicating whether it is the stop time or start time of a drug (character)

# ddi_data should contain a pDDI database, stating which combinations of drug_id's form an interaction

# each row should contain 3 variables:
# 1: drug_a; this should indicate a drug which interacts with drug_b
# 2: drug_b; this should indicate a drug which interacts with drug_a
# 3: interaction; this should indicate the interaction number of the interaction between drug_a and drug_b

# BEWARE: make sure that interactions between drug_a and drug_b are represented twice in ddi_data, 
# namely by indicating a row where the first drug is drug_a and the second drug is drug_b, AND a row where the first drug is drug_b and the second drug is drug_a
# for example, if APPLES interact with ORANGES, make sure there are two rows in ddi_data:

# drug_a  ; drug_b  ; interaction
# APPLES  ; ORANGES ; 30
# ORANGES ; APPLES  ; 30

############### preparing the data for pDDI detection #################  

# We defined a pDDI as the administration of two drugs known to interact, administered within 24 hours of each other

# set parameter gap time
gap_time = 24   # set the gap time, which was 24 hours in our study

# make a vector with all admission ids
admissions = unique(drug_data$admission_id)

# make sure the data is ordered by admission first, and by time second
drug_data = drug_data[order(drug_data$admission_id, drug_data$time),]

# split the dataset into starts and stops
drug_data_starts = subset(drug_data, type == "start")
drug_data_stops = subset(drug_data, type == "stop")

# add the gaptime (24h) to the stop times
drug_data_stops$time = drug_data_stops$time + (gap_time*60*60) # we do *60*60 because it adds up in minutes

# bind the starts and stops together again
drug_data_2 = rbind(drug_data_starts, drug_data_stops)

# order the data by time again
drug_data_2 = drug_data_2[order(drug_data_2$time),]

# If the time interval between administrations of the same drug did not exceed 24 hours, 
# the separate administrations were merged into one drug administration record. 
# The resulting record was given the start time of the first administration and the stop time of the last administration.

# This function mergres separate administration into one drug administration record
process_1adm = function(df, admission){
  this_admission = subset(df, df$admission_id == admission) # select drugs for this admission
  active_drugs = vector() # here the active drugs will be stored. drugs can be added and removed over time
  index = vector() # here the indices of the rows we want to keep are stored
  
  for(i in 1:nrow(this_admission)){ # iterate through all drug stops and starts of this admission
    if(this_admission$type[i] == 'start'){ # if this row is a drug start
      if(this_admission$drug_id[i] %in% active_drugs){ # check if it is in active drugs already
      }else{
        index = append(index,i) # if not, add this row's index to index, this will be the start time of the resulting drug administration record
      }
      active_drugs = append(active_drugs, this_admission$drug_id[i]) # and add the drug to active drugs
    }
    if(this_admission$type[i] == 'stop'){ # if this row is a drug stop
      if(sum(active_drugs == this_admission$drug_id[i]) == 1){ # if this drug is on the active drugs list once
        index = append(index, i) # add this row's index to index, this will be the stop time of the resulting drug administration record
        active_drugs = active_drugs[-match(this_admission$drug_id[i], active_drugs)] # remove the drug from active drugs, as it stopped
      }
      if(sum(active_drugs == this_admission$drug_id[i]) > 1){ # if its on the active drug list more than once, 
        active_drugs = active_drugs[-match(this_admission$drug_id[i], active_drugs)] # remove from the active drug list once and loop will be continued
      }
    }
  }
  return(this_admission[index,]) # return the rows containing the resulting start and stop time(s)
}


df_grow = data.frame() # create a new data frame

# apply the function process_1adm to all admissions in the dataset
for(o in admissions){
  p = process_1adm(drug_data_2, o) # apply function to specific admission
  df_grow = rbind(p, df_grow) # add the results for admission o to the growing dataframe
}

# this is the dataset we will detect pDDIs is, it is now prepared for the detection algorithm
drug_data_3 = df_grow

# make sure its ordered by time again
drug_data_3 = drug_data_3[order(drug_data_3$time),]

# you can check out what it looks like now
View(drug_data_3)

############### pDDI detection ################# 

# ddi_check function, this function checks if the drug indicated by drug_id forms any interactions with the active drugs, according to pDDI database 'ddis'
ddi_check = function(drug_id, ddis, active_drugs){ 
  p = ddis[drug_a == drug_id] # select rows containing a drug combination with drug_a involved / which interactions are possible with drug_a?
  p = p[p$drug_b %in% active_drugs] # further restriction: select drug combinations involving drug_a of which drug_b is on the active list / which of the possible interactions with drug_a also have drug_b on the active list?
  p = p$interaction # select the interaction number of the detected pDDIs
  return(p) # return interaction numbers
}

# this function uses ddi_check to detect pDDIs in the drug data
ddi_detect = function(drugs,admission){ # drugs is a dataframe containing all drugs of all admissions
  adm = admission
  sub = subset(drugs, drugs$admission_id == adm) # select the drugs for admission 'adm'
  active <- NULL # start with empty active list
  interactions  <- NULL # start with empty interactions list
  
  if(length(sub$admission_id) > 0 ){ # if there is any medication for this admission
    for(i in 1:nrow(sub)){ # iterate through drug administration starts and stops
      if(sub$type[i] == 'start'){ # if of type start
        ints = ddi_check(sub$drug_id[i], ddi_data, active) # check if there are any interactions between the drug started and the active medications
        interactions = append(interactions, ints) # add those to the interaction list 
        active = append(active,sub$drug_id[i]) # put the started drug on the active list
      }else{
        active = active[-match(sub$drug_id[i], active)] # if the row is of type stop, remove this medication from the active list
      }
    }
  }
  return(unique(interactions)) # return the list of interactions
}

# make an empty list to store the results in
results = list(admissions, list())

# make an empty data.frame to store results in per admission
results_clean = data.frame(admission_id = admissions)

# number of admissions to be analysed
n_admissions = length(admissions)

# iterate through all admissions and apply the ddi_detect function on each admission
i = 1
for(a in admissions[1:n_admissions]){
  results[[2]][[i]] = ddi_detect(drug_data_3,a) # apply function to specific admission and store result in the results list
  i = i+1
}

admissions_ids = results[[1]] # get admission ids out of list structure
interactions = results[[2]] # get results out of list structure

for(j in 1:n_admissions){ # iterate through admissions
  n = length(interactions[[j]]) # count number of interactions for admission j
  if(n > 0){ # if there are any
    for(k in 1:n){ # iterate through interactions
      results_clean[j,k+1] = interactions[[j]][k] # put them in results_clean for a better viewable results dataframe
    }
  }
}

# store the final results in a variable called final_results
final_results = results_clean

# final results contains the results of the pDDI detection algorithm, the first column indicates the admission ID, the columns following indicate the pDDIs that were present in that admission.
# if V2 is empty for a specific admission, it means that that admission had no pDDIs. 

# write the results into a CSV file
write.csv(final_results, 'final_results.csv', row.names = F)





