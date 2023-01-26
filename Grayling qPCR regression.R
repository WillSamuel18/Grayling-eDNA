################################################################################
################################ Will Samuel ###################################
####################### 2022 Grayling eDNA Modeling ############################
################################################################################


library(tidyverse)


setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis")


qpcr_dat <- read.csv(file = '2022 Summer eDNA/eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = '2022 Summer eDNA/eDNA Data (2).csv')

effort_dat <- read.csv(file = '2022 Summer eDNA/W_Samuel_Fishing_Effort_Datasheet.csv') #Need to fill this out and create a .csv

grayling_dat <- read.csv(file = '2022 Summer Fish Sampling/Data/W_Samuel_Grayling_Datasheet2.csv') 

sonde_WS_dat <- read.csv(file = "2022 Summer eDNA/First, Second, and third Stint of intensive eDNA Sonde  AND broad sites data_export (1).csv")

sonde_CPC_dat <- read.csv(file = "2022 Summer eDNA/CARIBOU_POKER_CREEK_EDNA_SAMPLING_FOR_WS_2022.csv")



# Data Manipulation -------------------------------------------------------

##Clean up the Sonde Data

#Remove funky measurements


##Join all the datasets by sample_ID
edna_dat <- edna_dat %>% rename(sample_ID = barcode_num)
qpcr_dat <- qpcr_dat %>% rename(sample_ID = SampleId)

qpcr_dat <- merge(x = edna_dat, y = qpcr_dat, by = "sample_ID", all.x = TRUE)
str(qpcr.dat)
head(qpcr.dat)


#Standardize the qPCR measurements by volume
qpcr_dat$L_filtered_vial <- as.numeric(qpcr.dat$L_filtered_vial, na.rm = TRUE)
qpcr.dat$L_filtered_datasheet <- as.numeric(qpcr.dat$L_filtered_datasheet, na.rm = TRUE)


qpcr_dat <- qpcr.dat %>% 
  mutate("liters_filtered_avg" = 
           ifelse(is.na(L_filtered_datasheet),L_filtered_vial,
                  ifelse(is.na(L_filtered_vial), L_filtered_datasheet,((L_filtered_vial+L_filtered_datasheet)/2)))) %>% #This averages the volume measurements when there are both, defaults to the other when there are NAs
  mutate("copies_per_L" = AvgCopyNum/liters_filtered_avg)
head(qpcr.dat)


#Check for significant contamination
qpcr_contaminated <- qpcr_dat %>% filter(str_detect(notes, 'contamin')) 
qpcr_clean <- qpcr_dat %>% filter(!str_detect(notes, 'contamin')) 
summary(qpcr_clean$copies_per_L)
summary(qpcr_contaminated$copies_per_L)
#Looks good. None of the values are way higher than the rest of the data, revist this later when looking at averaging replicates

#Repeat this for "problem" filters (e.g., torn, etc. )
qpcr_problem <- qpcr_dat %>% filter(str_detect(problem., 'Y')) 
qpcr-clean2 <- qpcr_dat %>% filter(!str_detect(problem., 'Y')) 
summary(qpcr_clean2$copies_per_L)
summary(qpcr_problem$copies_per_L)
#This looks okay too. I can decide later whether to include these or not. 





?rep()
#Summarize fish data
grayling_dat$fish <- rep(1, times = 283)
grayling_dat$Fork_Length <- as.numeric(grayling_dat$Fork_Length, na.rm = TRUE)

grayling_sums <- grayling_dat %>% 
  filter(Reach_Type == "Control") %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date, Sampling_Method) %>% 
  #group_by(Reach_Type) %>% 
  summarize("Abundance" = sum(fish), 
            "Biomass_g" = sum(Weight, na.rm = TRUE),
            "Length_total_mm" = sum(Fork_Length, na.rm = TRUE)) 


#Calculate the combined sampling effort

#effort_dat <- effort_dat %>% 
#  filter(Reach_Type == "Control") %>% 
#  mutate("Effort" = (Angling_total+Efish_total)) #Need to come up with an equation to combine the sampling strategies


#Using Multigear Mean Standardization (MGMS) see Gibson-Reinemer et al. 2016 for more details

#Combine the effort and catch data
effort_dat <- effort_dat %>% 
  filter(Reach_Type == "Control")

grayling_sums <- merge(x = grayling_sums, y = effort_dat, by = "Date", all.x = TRUE)

#Step 1: calculate CPUE for each survey method
angling_effort_dat <- grayling_sums %>% 
  filter(Sampling_Method == Angling) %>% 
  mutate("Angling_CPUE_abun" = (Abundance/Angling_total)) %>% 
  mutate("Angling_CPUE_biom" = (Biomass/Angling_total)) 

#Step 2: calculate the mean total CPUE
mean.st.angling.abun <- mean(angling_effort_dat$Angling_CPUE_abun)
mean.st.angling.biom <- mean(angling_effort_dat$Angling_CPUE_biom)

#Step 3: use the mean total CPUE to create the standardized effort
angling_effort_dat <- angling_effort_dat %>% 
  mutate("Angling_CPUE_abun_STD" = (Angling_CPUE_abun/mean.st.angling.abun)) %>% 
  mutate("Angling_CPUE_biom_STD" = (Angling_CPUE_biom/mean.st.angling.biom)) 


#Repeat for E-fishing
Efish_effort_dat <- grayling_sums %>% 
  filter(Sampling_Method == E-Fishing) %>% 
  mutate("Efish_CPUE_abun" = (Abundance/Efish_total)) %>% #Maybe needs to be Time_on_Efish
  mutate("Efish_CPUE_biom" = (Biomass/Efish_total)) #Maybe needs to be Time_on_Efish

mean.st.Efish.abun <- mean(angling_effort_dat$Efish_CPUE_abun)
mean.st.Efish.biom <- mean(angling_effort_dat$Efish_CPUE_binom)

Efish_effort_dat <- Efish_effort_dat %>% 
  mutate("Efish_CPUE_abun_STD" = (Efish_CPUE_abun/mean.st.Efish.abun)) %>% 
  mutate("Efish_CPUE_biom_STD" = (Efish_CPUE_binom/mean.st.Efish.biom)) 


#Step 4: Combine the datasets back into one
st_effort_dat <- cbind(Angling_effort_dat, Efish_effort_dat) #Might be rbind not cbind


#Step 5: assign this standardized effort back to the dataset, calculate the MGMS_effort
st_effort_dat <- st_effort_dat %>% 
  mutate("MGMS_CPUE_abun" = (Efish_CPUE_abun_STD+Angling_CPUE_abun_STD)) %>% 
  mutate("MGMS_CPUE_biom" = (Efish_CPUE_biom_STD+Angling_CPUE_biom_STD)) 



  












# Prepare parameters for modeling -----------------------------------------











# Exploratory Analysis ----------------------------------------------------



qpcr_dat <- qpcr_dat %>% rename(Date = date)


qpcr_sums <- qpcr_dat %>% 
  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date) %>% 
  summarize("copies_per_L" = mean(copies_per_L)) 
 

grayling_sums <- merge(x = grayling_sums, y = qpcr_sums, by = "Date", all.x = TRUE)

plot(copies_per_L ~ Abundance, data = grayling_sums)
summary(lm(copies_per_L ~ Abundance, data = grayling_sums))

plot(copies_per_L ~ Biomass_g, data = grayling_sums)
summary(lm(copies_per_L ~ Biomass_g, data = grayling_sums))

#Hopefully this outlook improves once I standardize by effort


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




