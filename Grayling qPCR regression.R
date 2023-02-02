################################################################################
################################ Will Samuel ###################################
####################### 2022 Grayling eDNA Modeling ############################
################################################################################


library(tidyverse)  #for data manipulation
library(MuMIn)      #For model dredging



setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis")


qpcr_dat <- read.csv(file = '2022 Summer eDNA/eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = '2022 Summer eDNA/eDNA Data (2).csv')

effort_dat <- read.csv(file = '2022 Summer Fish Sampling/Data/W_Samuel_Fishing_Effort_Datasheet.csv') #Need to update the area of the CPC stream reaches

grayling_dat <- read.csv(file = '2022 Summer eDNA/W_Samuel_Grayling_Datasheet2+CPC.csv') 

#sonde_WS_dat <- read.csv(file = "2022 Summer eDNA/First, Second, and third Stint of intensive eDNA Sonde  AND broad sites data_export (1).csv")

#sonde_CPC_dat <- read.csv(file = "2022 Summer eDNA/CARIBOU_POKER_CREEK_EDNA_SAMPLING_FOR_WS_2022.csv")



# Data Manipulation -------------------------------------------------------

#Remove funky measurements#





##Join the datasets by sample_ID
edna_dat <- edna_dat %>% rename(sample_ID = barcode_num)
qpcr_dat <- qpcr_dat %>% rename(sample_ID = SampleId)


edna_dat <- merge(x = edna_dat, y = qpcr_dat, by = "sample_ID", all.x = TRUE)
str(edna_dat)
head(edna_dat)


#Filter out top of reaches for now
edna_dat <- edna_dat %>% 
  filter(!str_detect(sample_num, 'TR')) %>% 
  filter(!str_detect(sample_num, 'TR')) %>% 
  filter(!str_detect(sample_num, 'B')) %>% 
  filter(!str_detect(sample_num, 'C')) %>% 
    filter(!str_detect(sample_num, '.4')) #remove negative controls



#Standardize the qPCR measurements by volume
edna_dat$L_filtered_vial <- as.numeric(edna_dat$L_filtered_vial, na.rm = TRUE)
edna_dat$L_filtered_datasheet <- as.numeric(edna_dat$L_filtered_datasheet, na.rm = TRUE)


edna_dat <- edna_dat %>% 
  mutate("liters_filtered_avg" = 
           ifelse(is.na(L_filtered_datasheet),L_filtered_vial,
                  ifelse(is.na(L_filtered_vial), L_filtered_datasheet, ((L_filtered_vial+L_filtered_datasheet)/2)))) %>% #This averages the volume measurements when there are both, defaults to the other when there are NAs
  mutate("copies_per_L" = AvgCopyNum/liters_filtered_avg) %>% 
  mutate("temp_log" = ifelse(is.na(temp_log), temp_ds, temp_log)) %>% #Takes the datasheet sonde measurements when there NAs from the log. Does the same below
  mutate("ph_log" = ifelse(is.na(ph_log), ph_ds, ph_log)) %>% 
  mutate("sc_log" = ifelse(is.na(sc_log), sc_ds, sc_log)) %>% 
  mutate("hdo_ml.L_log" = ifelse(is.na(hdo_ml.L_log), hdo_ml.L_ds, hdo_ml.L_log)) %>% 
  mutate("hdo_perc_sat_log" = ifelse(is.na(hdo_perc_sat_log), hdo_perc_sat_ds, hdo_perc_sat_log)) %>% 
  mutate("turb_log" = ifelse(is.na(turb_log), turb_ds, turb_log)) 
head(edna_dat)







#Check for significant contamination
#qpcr_contaminated <- edna_dat %>% filter(str_detect(notes, 'contamin')) 
#qpcr_clean <- edna_dat %>% filter(!str_detect(notes, 'contamin')) 
#summary(qpcr_clean$copies_per_L)
#summary(qpcr_contaminated$copies_per_L)
#Looks good. None of the values are way higher than the rest of the data, revist this later when looking at averaging replicates

#Repeat this for "problem" filters (e.g., torn, etc. )
#qpcr_problem <- edna_dat %>% filter(str_detect(problem., 'Y')) 
#qpcr_clean2 <- edna_dat %>% filter(!str_detect(problem., 'Y')) 
#summary(qpcr_clean2$copies_per_L)
#summary(qpcr_problem$copies_per_L)
#This looks okay too. I can decide later whether to include these or not. 





#Summarize fish data
grayling_dat$fish <- rep(1, times = 429)
grayling_dat$Fork_Length <- as.numeric(grayling_dat$Fork_Length, na.rm = TRUE)

grayling_sums <- grayling_dat %>% 
  filter(Reach_Type == "Control") %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date, Sampling_Method, Site) %>% 
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

effort_dat <- merge(x = effort_dat , y = grayling_sums, by = "Date", all.x = TRUE)

#Step 1: calculate CPUE for each survey method
angling_effort_dat <- effort_dat %>% 
  filter(Sampling_Method == "Angling") %>% 
  mutate("Angling_CPUE_abun" = (Abundance/(Angling_total*Num_anglers))) %>% 
  mutate("Angling_CPUE_biom" = (Biomass_g/(Angling_total*Num_anglers))) 

#Step 2: calculate the mean total CPUE
mean.st.angling.abun <- mean(angling_effort_dat$Angling_CPUE_abun, na.rm = TRUE)
mean.st.angling.biom <- mean(angling_effort_dat$Angling_CPUE_biom, na.rm = TRUE)

#Step 3: use the mean total CPUE to create the standardized effort
angling_effort_dat <- angling_effort_dat %>% 
  mutate("Angling_CPUE_abun_STD" = (Angling_CPUE_abun/mean.st.angling.abun)) %>% 
  mutate("Angling_CPUE_biom_STD" = (Angling_CPUE_biom/mean.st.angling.biom)) 


#Repeat for E-fishing
Efish_effort_dat <- effort_dat %>% 
  filter(Sampling_Method == "E-Fishing") %>% 
  mutate("Efish_CPUE_abun" = (Abundance/Efish_total)) %>% #Maybe needs to be Time_on_Efish
  mutate("Efish_CPUE_biom" = (Biomass_g/Efish_total)) #Maybe needs to be Time_on_Efish

mean.st.Efish.abun <- mean(Efish_effort_dat$Efish_CPUE_abun, na.rm = TRUE)
mean.st.Efish.biom <- mean(Efish_effort_dat$Efish_CPUE_biom, na.rm = TRUE)

Efish_effort_dat <- Efish_effort_dat %>% 
  mutate("Efish_CPUE_abun_STD" = (Efish_CPUE_abun/mean.st.Efish.abun)) %>% 
  mutate("Efish_CPUE_biom_STD" = (Efish_CPUE_biom/mean.st.Efish.biom)) 


#Step 4: Combine the datasets back into one
str(angling_effort_dat)
str(Efish_effort_dat)


st_effort_dat <- merge(x = angling_effort_dat , y = Efish_effort_dat, by = c("Date", "Site_Num"), all.x = TRUE)



#library(plyr)       #For the rbind.fill function
#st_effort_dat <- rbind.fill(angling_effort_dat, Efish_effort_dat)
#detach(package:plyr,unload=TRUE)


#Step 5: assign this standardized effort back to the dataset, calculate the MGMS_effort
st_effort_dat <- st_effort_dat %>% 
  mutate("Angling_CPUE_abun_STD" = ifelse(is.na(Angling_CPUE_abun_STD), 0, Angling_CPUE_abun_STD)) %>% 
  mutate("Angling_CPUE_biom_STD" = ifelse(is.na(Angling_CPUE_biom_STD), 0, Angling_CPUE_biom_STD)) %>% 
  mutate("Efish_CPUE_abun_STD" = ifelse(is.na(Efish_CPUE_abun_STD), 0, Efish_CPUE_abun_STD)) %>% 
  mutate("Efish_CPUE_biom_STD" = ifelse(is.na(Efish_CPUE_biom_STD), 0, Efish_CPUE_biom_STD)) %>% 
  
  mutate("MGMS_CPUE_abun" = (Efish_CPUE_abun_STD+Angling_CPUE_abun_STD)) %>% 
  mutate("MGMS_CPUE_biom" = (Efish_CPUE_biom_STD+Angling_CPUE_biom_STD)) 


str(st_effort_dat)
  




# Prepare parameters for modeling -----------------------------------------





#NA VALUES???
#Remove or correct negative values? 
#(9 > pH < 5)???? pH on the Salcha trip

plot(hdo_ml.L_log ~ hdo_perc_sat_log, data = edna_dat)
summary(lm(hdo_ml.L_log ~ hdo_perc_sat_log, data = edna_dat))



# Exploratory Analysis ----------------------------------------------------




#qpcr_dat <- qpcr_dat %>% rename(Date = date)


#qpcr_sums <- qpcr_dat %>% 
#  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
#  group_by(Date) %>% 
#  summarize("copies_per_L" = mean(copies_per_L)) 
 

#grayling_sums <- merge(x = grayling_sums, y = qpcr_sums, by = "Date", all.x = TRUE)



#Using MGMS standardized effort
### MAKE SURE TO REMOVE ANGLE CREEK SAMPLS ON 6/22/22 AND 7/29/22
st_effort_dat <- st_effort_dat[-c(5),]

edna_dat <- edna_dat %>% rename(Date = date)

edna_sums <- edna_dat %>% 
  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date) %>% 
  summarize("copies_per_L" = mean(copies_per_L),
            "Vel_ms" = mean(c(flow_1, flow_2, flow_3), na.rm = TRUE), 
            "water_temp" = mean(temp_log, na.rm = TRUE), 
            "pH" = mean(ph_log, na.rm = TRUE), 
            "SC" = mean(sc_log, na.rm = TRUE), 
            "HDO" = mean(hdo_ml.L_log, na.rm = TRUE), 
            "HDO_perc" = mean(hdo_perc_sat_log, na.rm = TRUE), 
            "Turb" = mean(turb_log, na.rm = TRUE))



st_effort_dat <- merge(x = st_effort_dat, y = edna_sums, by = "Date", all.x = TRUE)

st_effort_dat <- st_effort_dat[-c(1),] #remove Belle creek on 6/13/22, its an outlier





###Should I log transform this or not?????????????

plot(log(copies_per_L) ~ MGMS_CPUE_abun, data = st_effort_dat)
summary(lm(copies_per_L ~ MGMS_CPUE_abun+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat))

plot(copies_per_L ~ MGMS_CPUE_biom, data = st_effort_dat)
#text(log(copies_per_L)~MGMS_CPUE_biom, labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 2)
summary(lm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat))

hist(grayling_dat$Fork_Length) #make this a multi panel plot for each site/sampling event (Date). 





#Using each replicate as a datapoint, instead of the means of the replicates
#st_effort_dat_all <- merge(x = edna_dat, y = st_effort_dat, by = "Date", all.x = TRUE)
#st_effort_dat_all <- merge(x = st_effort_dat_all, y = qpcr_dat, by = "sample_ID", all.x = TRUE)

#plot(log(copies_per_L.x) ~ MGMS_CPUE_biom, data = st_effort_dat_all)
#text(log(copies_per_L)~MGMS_CPUE_biom, labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 2)
#summary(lm(copies_per_L.x ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat_all))

#Since R^2 is much bigger than adjusted R^2, its indicating model overfit. So I think its best to use the means for everything. 


# VIF ---------------------------------------------------------------------





# Models ------------------------------------------------------------------


### MAKE SURE TO REMOVE ANGLE CREEK SAMPLS ON 6/24/22 AND 7/29/22


#st_effort_dat <- st_effort_dat[-c(5),]

#qpcr_dat <- qpcr_dat %>% rename(Date = date)

#qpcr_sums <- qpcr_dat %>% 
#  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
#  group_by(Date) %>% 
#  summarize("copies_per_L" = mean(copies_per_L)) 

#st_effort_dat <- merge(x = st_effort_dat, y = qpcr_sums, by = "Date", all.x = TRUE)

#st_effort_dat <- st_effort_dat[-c(1),] #remove Belle creek on 6/13/22, its an outlier


glm(log(copies_per_L) ~ MGMS_CPUE_biom+ )


#VIF/corr plot

#GLM with Abundnace

#GLM with Biomass

#GLM with Abundance Density

#GLM with Biomass Density


#Dredge the best model of each? 


#Compare alternative models and pick the best one to predict, 
#ORRRRR Predict with all of them and average the results... 




# Model Evaluation --------------------------------------------------------




# Predicting to other sites -----------------------------------------------






