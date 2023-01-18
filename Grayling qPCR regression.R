################################################################################
################################ Will Samuel ###################################
####################### 2022 Grayling eDNA Modeling ############################
################################################################################


library(tidyverse)


setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis/2022 Summer eDNA")


qpcr_dat <- read.csv(file = 'eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = 'eDNA Data (2).csv')

effort_dat <- read.csv(file = 'C:\Users\npwil\OneDrive\Desktop\School\Grad School\Thesis\Data and Analysis\2022 Summer Fish Sampling\Data.csv')

sonde_WS_dat <- read.csv(file = "First, Second, and third Stint of intensive eDNA Sonde  AND broad sites data_export (1).csv")

sonde_CPC_dat <- read.csv(file = "CARIBOU_POKER_CREEK_EDNA_SAMPLING_FOR_WS_2022.csv")



# Data Manipulation -------------------------------------------------------

##Clean up the Sonde Data

#Remove funk measurements


##Join all the datasets by sample_num or barcode_num






# Prepare parameters for modeling -----------------------------------------

##Standardize the qPCR data by L filtered 



# VIF ---------------------------------------------------------------------





# Models ------------------------------------------------------------------


#GLM with Abundnace

#GLM with Biomass

#GLM with Abundance Density

#GLM with Biomass Density


#Dredge the best model of each? 


#Compare alternative models and pick the best one to predict, 
#ORRRRR Predict with all of them and average the results... 




# Model Evaluation --------------------------------------------------------




# Predicting to other sites -----------------------------------------------




