################################################################################
################################ Will Samuel ###################################
################ 2022 Grayling eDNA N-Mixture Modeling #########################
################################################################################


library(tidyverse)  #for data manipulation
library(lubridate)
library(MuMIn)      #For model dredging
library(ggcorrplot) #to make a correlation plot
library(rcompanion) #for pseudo R-squared
library(stringr)    #To split a column in 2 (used for efish time)
library(cowplot)
library(unmarked)   #For the N-Mixture model
library(patchwork)
library(gridExtra)
library(AICcmodavg)	#For model averaging which can help you find the most important predictors



setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis")




qpcr_dat <- read.csv(file = '2022 Summer eDNA/eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = '2022 Summer eDNA/eDNA Data (2).csv')

#effort_dat <- read.csv(file = '2022 Summer Fish Sampling/Data/W_Samuel_Fishing_Effort_Datasheet.csv') 

#grayling_dat <- read.csv(file = '2022 Summer eDNA/W_Samuel_Grayling_Datasheet2+CPC.csv') 

# Data Manipulation -------------------------------------------------------




##Join the datasets by sample_ID
edna_dat <- edna_dat %>% rename(sample_ID = barcode_num)
qpcr_dat <- qpcr_dat %>% rename(sample_ID = SampleId)


edna_dat <- merge(x = edna_dat, y = qpcr_dat, by = "sample_ID", all.x = TRUE)
str(edna_dat)
head(edna_dat)

edna_dat <- edna_dat %>% rename(sample_num = ï..sample_num)



#Filter out top of reaches on intensive sites
edna_dat <- edna_dat %>% 
  #filter(!str_detect(sample_num, 'TR')) %>% 
  filter(transect == 'A') %>% 
  filter(!str_detect(sample_num, '\\.4')) #remove negative controls



#Standardize the qPCR measurements by volume
edna_dat$L_filtered_vial <- as.numeric(edna_dat$L_filtered_vial, na.rm = TRUE)
edna_dat$L_filtered_datasheet <- as.numeric(edna_dat$L_filtered_datasheet, na.rm = TRUE)


edna_dat <- edna_dat %>% 
  mutate("liters_filtered_avg" = 
           ifelse(is.na(L_filtered_datasheet),L_filtered_vial,
                  ifelse(is.na(L_filtered_vial), L_filtered_datasheet, ((L_filtered_vial+L_filtered_datasheet)/2)))) %>% #This averages the volume measurements when there are both, defaults to the other when there are NAs
  mutate("copies_per_L" = AvgCopyNum/liters_filtered_avg) %>% 
  mutate("temp_log" = ifelse(is.na(temp_log), temp_ds, temp_log)) %>% #Takes the datasheet sonde measurements when there NAs from the log. Does the same below
  mutate("ph_log" = ifelse(ph_log == "NA", NA, ph_log)) %>%
  mutate("ph_log" = ifelse(is.na(ph_log), ph_ds, ph_log)) %>%
  #mutate("ph_log" = ifelse(11 < ph_log < 3, ph_log, "NA")) %>%  #Remove faulty pH measurements
  mutate("sc_log" = ifelse(is.na(sc_log), sc_ds, sc_log)) %>% 
  mutate("hdo_ml.L_log" = ifelse(is.na(hdo_ml.L_log), hdo_ml.L_ds, hdo_ml.L_log)) %>% 
  mutate("hdo_perc_sat_log" = ifelse(is.na(hdo_perc_sat_log), hdo_perc_sat_ds, hdo_perc_sat_log)) %>% 
  mutate("turb_log" = ifelse(is.na(turb_log), turb_ds, turb_log)) 


edna_dat <- edna_dat %>% 
  mutate(ph_log = ifelse(ph_log < 3 | ph_log > 11, NA, ph_log)) #Remove faulty pH measurements


#View(edna_dat)







#Check for significant contamination
qpcr_contaminated <- edna_dat %>% filter(str_detect(notes, 'contamin')) 
qpcr_clean <- edna_dat %>% filter(!str_detect(notes, 'contamin')) 
summary(qpcr_clean$copies_per_L)
summary(qpcr_contaminated$copies_per_L)
#Looks good. None of the values are way higher than the rest of the data, revist this later when looking at averaging replicates

#Repeat this for "problem" filters (e.g., torn, etc. )
qpcr_problem <- edna_dat %>% filter(str_detect(problem., 'Y')) 
qpcr_clean2 <- edna_dat %>% filter(!str_detect(problem., 'Y')) 
summary(qpcr_clean2$copies_per_L)
summary(qpcr_problem$copies_per_L)
#This looks okay too. 


#Remove the clearly contaminated filters
edna_dat <- edna_dat %>% 
  filter(!str_detect(notes, 'contaminated|contamination|torn|hole|bubble')) #remove contaminated filters





#Format data by site (n = 66????)

subset_dat <- edna_dat %>% 
  select(sample_num, Site_Num, copies_per_L, liters_filtered_avg) %>% 
  mutate(Site_Num = ifelse(str_detect(sample_num , '1-1.'), "1--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '1-2.'), "1--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '1-3.'), "1--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '2-1.'), "2--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '2-2.'), "2--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '2-3.'), "2--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '3-1.'), "3--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '3-2.'), "3--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '3-3.'), "3--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '4-1.'), "4--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '4-2.'), "4--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '4-3.'), "4--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '5-1.'), "5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '5-2.'), "5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '5-3.'), "5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '6-1.'), "5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '6-2.'), "5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '6-3.'), "5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '7-1.'), "6--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '7-2.'), "6--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '7-3.'), "6--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R5-1.'), "R5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R5-2.'), "R5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R5-3.'), "R5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R20-1.'), "R20--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R20-2.'), "R20--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R20-3.'), "R20--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R33-1.'), "R33--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R33-2.'), "R33--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R33-3.'), "R33--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R65-1.'), "R65--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R65-2.'), "R65--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R65-3.'), "R65--3", Site_Num))
head(subset_dat)

site_num <- subset_dat %>%
  select(Site_Num) %>% 
  distinct(Site_Num)

replicate_1 <- subset_dat %>%
  filter(str_detect(sample_num, '\\.1')) %>% 
  select(Site_Num, copies_per_L)
  
replicate_2 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.2')) %>% 
  select(Site_Num, copies_per_L)

replicate_3 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.3')) %>% 
  select(Site_Num, copies_per_L)

replicate_4 <- subset_dat %>% ###This is blank since all the .4 filters were negative controls
  filter(str_detect(sample_num, '\\.4')) %>% 
  select(Site_Num, copies_per_L)

replicate_5 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.5')) %>% 
  select(Site_Num, copies_per_L)

replicate_6 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.6')) %>% 
  select(Site_Num, copies_per_L)

replicate_7 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.7')) %>% 
  select(Site_Num, copies_per_L)

replicate_8 <- subset_dat %>% 
  filter(str_detect(sample_num, '\\.8')) %>% 
  select(Site_Num, copies_per_L)


liters_filtered <- subset_dat %>%
  group_by(Site_Num) %>%
  summarize(mean_liters_filtered_avg = mean(liters_filtered_avg, na.rm = TRUE))



replicate_dat <- merge(x = site_num, y = replicate_1, by = "Site_Num", all.x = TRUE)
replicate_dat <- replicate_dat %>% rename(rep_1_copies_per_L = copies_per_L)

replicate_dat <- merge(x = replicate_dat, y = replicate_2, by = "Site_Num", all.x = TRUE)
replicate_dat <- replicate_dat %>% rename(rep_2_copies_per_L = copies_per_L)

replicate_dat <- merge(x = replicate_dat, y = replicate_3, by = "Site_Num", all.x = TRUE)
replicate_dat <- replicate_dat %>% rename(rep_3_copies_per_L = copies_per_L)

replicate_dat <- merge(x = replicate_dat, y = replicate_5, by = "Site_Num", all.x = TRUE)
replicate_dat <- replicate_dat %>% rename(rep_4_copies_per_L = copies_per_L)

replicate_dat <- merge(x = replicate_dat, y = replicate_6, by = "Site_Num", all.x = TRUE)
replicate_dat <- replicate_dat %>% rename(rep_5_copies_per_L = copies_per_L)


replicate_dat <- merge(x = replicate_dat, y = liters_filtered, by = "Site_Num", all.x = TRUE)

View(replicate_dat)



#Look at the replicates 
ggplot(replicate_dat, aes(x=Site_Num)) +
  geom_point(aes(y=rep_1_copies_per_L))+
  geom_point(aes(y=rep_2_copies_per_L))+
  geom_point(aes(y=rep_3_copies_per_L))+
  geom_point(aes(y=rep_4_copies_per_L))+
  geom_point(aes(y=rep_5_copies_per_L))


ggplot(replicate_dat, aes(x=Site_Num)) +
  geom_histogram(aes(y=mean(rep_1_copies_per_L, rep_2_copies_per_L, rep_3_copies_per_L, na.rm = T)))



intensive_sites <- replicate_dat %>% 
  filter(str_detect(Site_Num, '-'))

ggplot(intensive_sites, aes(x=Site_Num)) +
  geom_point(aes(y=rep_1_copies_per_L))+
  geom_point(aes(y=rep_2_copies_per_L))+
  geom_point(aes(y=rep_3_copies_per_L))+
  geom_point(aes(y=rep_4_copies_per_L))+
  geom_point(aes(y=rep_5_copies_per_L))+
  scale_y_continuous(n.breaks = 5, limits = c(0,11500))




broad_sites <- replicate_dat %>% 
  filter(!str_detect(Site_Num, '-'))

ggplot(broad_sites, aes(x=Site_Num)) +
  geom_point(aes(y=rep_1_copies_per_L))+
  geom_point(aes(y=rep_2_copies_per_L))+
  geom_point(aes(y=rep_3_copies_per_L))+
  geom_point(aes(y=rep_4_copies_per_L))+
  geom_point(aes(y=rep_5_copies_per_L))




#Join to the covariates:
covariate_dat <- edna_dat %>% 
  select(sample_num, Site_Num, date, drainage, site, lat.dd, long.dd, weather, flow_1, flow_2, flow_3, temp_log, ph_log, sc_log, hdo_ml.L_log, hdo_perc_sat_log, turb_log) %>% 
  mutate(Site_Num = ifelse(str_detect(sample_num , '1-1.'), "1--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '1-2.'), "1--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '1-3.'), "1--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '2-1.'), "2--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '2-2.'), "2--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '2-3.'), "2--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '3-1.'), "3--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '3-2.'), "3--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '3-3.'), "3--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '4-1.'), "4--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '4-2.'), "4--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '4-3.'), "4--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '5-1.'), "5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '5-2.'), "5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '5-3.'), "5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '6-1.'), "5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '6-2.'), "5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '6-3.'), "5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , '7-1.'), "6--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '7-2.'), "6--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , '7-3.'), "6--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R5-1.'), "R5--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R5-2.'), "R5--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R5-3.'), "R5--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R20-1.'), "R20--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R20-2.'), "R20--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R20-3.'), "R20--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R33-1.'), "R33--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R33-2.'), "R33--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R33-3.'), "R33--3", Site_Num),
         
         Site_Num = ifelse(str_detect(sample_num , 'R65-1.'), "R65--1", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R65-2.'), "R65--2", Site_Num),
         Site_Num = ifelse(str_detect(sample_num , 'R65-3.'), "R65--3", Site_Num)) %>% 
  group_by(Site_Num) %>% 
  summarize(date = first(date), 
            drainage = first(drainage), 
            site = first(site), 
            lat.dd = first(lat.dd), 
            long.dd = first(long.dd), 
            weather = first(weather), 
            flow_1 = mean(flow_1, na.rm = TRUE), 
            flow_2 = mean(flow_2, na.rm = TRUE), 
            flow_3 = mean(flow_3, na.rm = TRUE), 
            temp_log = mean(temp_log, na.rm = TRUE), 
            ph_log = ph_log, 
            sc_log = mean(sc_log, na.rm = TRUE), 
            hdo_ml.L_log = mean(hdo_ml.L_log, na.rm = TRUE), 
            hdo_perc_sat_log = mean(hdo_perc_sat_log, na.rm = TRUE), 
            turb_log = mean(turb_log, na.rm = TRUE)) %>% 
    mutate(flow = ((flow_1+flow_2+flow_3)/3)) %>% 
    select(-flow_1, -flow_2, -flow_3) %>% 
    mutate()







date <- mdy(covariate_dat$date)
doy <- as.vector(yday(date))
doy

covariate_dat <- cbind(covariate_dat, doy)

covariate_dat <- covariate_dat %>% rename(doy = "...15")



covariate_dat



#Remove duplicate rows, not sure why they are there...
covariate_dat <- covariate_dat %>%
  distinct(Site_Num, .keep_all = TRUE)

# Print the result
print(covariate_dat)




r_dat <- merge(x = replicate_dat, y = covariate_dat, by = "Site_Num", all.x = TRUE, all.y = F)

#r_dat <- r_dat[,-12]
#r_dat <- distinct(r_dat)


write.csv(r_dat, "2022 Summer eDNA/r_dat.csv")


##Need to add other covariates, like beaver density and UPDATED (2022) wildfire history. Summarize these by HUC in ArcGIS, use the eDNA lat/long to relate them to site numbers. 

#VB.beaver.fire.dat <- read.csv("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis/Beaver Project/Beaver_Project_2/Theoretical Passage Data Files/VB_with_beaver_fire_data.csv")

#str(VB.beaver.fire.dat)

#DONT USE FIRE DATA FROM THIS ^^ SINCE IT ONLY GOES TO 2017. The beaver data is still good though



spatial_dat <- read.csv("2022 Summer eDNA/eDNA_VB_FIRE_NETMAP_eDNApoints_Table.csv")
str(spatial_dat)


spatial_dat <- spatial_dat %>% 
  mutate(Surveyed = ifelse(is.na(Surveyed), 0, Surveyed),
         Burned_Numerical = ifelse(Percent_burned > 0, 1, 0),
         )
         





###Need to associate all the eDNA points with the data from the HUCs!!!

#First, filter out the early repeat sampling events
r_dat <- read.csv("2022 Summer eDNA/r_dat.csv")
str(r_dat)
#View(r_dat)

r_dat <- r_dat %>%
  filter(!str_detect(Site_Num, "--1|--2"),
         !str_detect(Site_Num, "R20|R33|R65"))
  



#Remove Site 162 from eDNA data, little Chena by Nordale 
r_dat <- r_dat[-c(27),]


#Remove Site 62, duplicate polygon from the same issue
which(spatial_dat$Site_Num == 85)
spatial_dat[c(42),]

which(spatial_dat$Site_Num == 62)
spatial_dat[c(39),]

which(spatial_dat$Site_Num == "")
spatial_dat[c(64),]


#Remove mystery site (85) and repeat value (62)
spatial_dat <- spatial_dat[-c(42),]
spatial_dat <- spatial_dat[-c(39),]
#spatial_dat <- spatial_dat[-c(55),]
#spatial_dat <- spatial_dat[-c(64:69),]#Remove empty rows from excel
#spatial_dat <- spatial_dat[-c(65),]#Remove empty rows from excel

spatial_dat <- spatial_dat %>% 
  filter(!Site_Num == "")





#I can't find eDNA data to match up with site 85 on the spatial dataset (McKay Creek HUC. It looks like that might have been a mistaken point



r_dat_sitenums <- data.frame(sort(r_dat$Site_Num))

#vect <- c(10000)
#r_dat_sitenums <- rbind(r_dat_sitenums, vect)
#r_dat_sitenums <- rbind(r_dat_sitenums, vect)


spatial_sitenums <- data.frame(sort(spatial_dat$Site_Num))

sitenums <- cbind(r_dat_sitenums, spatial_sitenums)
sitenums
View(sitenums)








combined_dat <- merge(x = r_dat, y = spatial_dat, by = "Site_Num", all.x = T, all.y = T)



View(combined_dat)


combined_dat <- combined_dat %>% 
  mutate("Alldam_density" = Alldam_Count/VB_AreaSqKm,
         "Ydam_density" = Ydam_Count/VB_AreaSqKm,
         "Time_Since_Last_Burn" = 2022-Year,
         "Time_Since_Last_Burn" = ifelse(is.na(Time_Since_Last_Burn), 100, Time_Since_Last_Burn))



combined_dat <- combined_dat %>% 
  mutate("mean_abundance" = rowMeans(.[c("rep_1_copies_per_L", "rep_2_copies_per_L", "rep_3_copies_per_L", "rep_4_copies_per_L", "rep_5_copies_per_L")], na.rm = TRUE))
         
         
         #(rep_1_copies_per_L + rep_2_copies_per_L + rep_3_copies_per_L+rep_4_copies_per_L+rep_5_copies_per_L)/5)

combined_dat$mean_abundance




write.csv(combined_dat, "2022 Summer eDNA/combined_dat.csv")






# Contamination analysis --------------------------------------------------



qpcr_dat <- read.csv(file = '2022 Summer eDNA/eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = '2022 Summer eDNA/eDNA Data (2).csv')


##Join the datasets by sample_ID
edna_dat <- edna_dat %>% rename(sample_ID = barcode_num)
qpcr_dat <- qpcr_dat %>% rename(sample_ID = SampleId)


edna_dat <- merge(x = edna_dat, y = qpcr_dat, by = "sample_ID", all.x = TRUE)
str(edna_dat)
head(edna_dat)

edna_dat <- edna_dat %>% rename(sample_num = ï..sample_num)



#Filter out top of reaches on intensive sites
edna_dat <- edna_dat %>% 
  #filter(!str_detect(sample_num, 'TR')) %>% 
  filter(transect == 'A') #%>% 
  #filter(!str_detect(sample_num, '\\.4')) #remove negative controls



#Standardize the qPCR measurements by volume
edna_dat$L_filtered_vial <- as.numeric(edna_dat$L_filtered_vial, na.rm = TRUE)
edna_dat$L_filtered_datasheet <- as.numeric(edna_dat$L_filtered_datasheet, na.rm = TRUE)


edna_dat <- edna_dat %>% 
  mutate("liters_filtered_avg" = 
           ifelse(is.na(L_filtered_datasheet),L_filtered_vial,
                  ifelse(is.na(L_filtered_vial), L_filtered_datasheet, ((L_filtered_vial+L_filtered_datasheet)/2)))) %>% #This averages the volume measurements when there are both, defaults to the other when there are NAs
  mutate("copies_per_L" = AvgCopyNum/liters_filtered_avg) %>% 
  mutate("temp_log" = ifelse(is.na(temp_log), temp_ds, temp_log)) %>% #Takes the datasheet sonde measurements when there NAs from the log. Does the same below
  mutate("ph_log" = ifelse(ph_log == "NA", NA, ph_log)) %>%
  mutate("ph_log" = ifelse(is.na(ph_log), ph_ds, ph_log)) %>%
  #mutate("ph_log" = ifelse(11 < ph_log < 3, ph_log, "NA")) %>%  #Remove faulty pH measurements
  mutate("sc_log" = ifelse(is.na(sc_log), sc_ds, sc_log)) %>% 
  mutate("hdo_ml.L_log" = ifelse(is.na(hdo_ml.L_log), hdo_ml.L_ds, hdo_ml.L_log)) %>% 
  mutate("hdo_perc_sat_log" = ifelse(is.na(hdo_perc_sat_log), hdo_perc_sat_ds, hdo_perc_sat_log)) %>% 
  mutate("turb_log" = ifelse(is.na(turb_log), turb_ds, turb_log)) 


edna_dat <- edna_dat %>% 
  mutate(ph_log = ifelse(ph_log < 3 | ph_log > 11, NA, ph_log)) #Remove faulty pH measurements


#View(edna_dat)





#Remove the clearly contaminated filters
edna_dat <- edna_dat %>% 
  filter(!str_detect(notes, 'contam|contaminated|contamination|torn|hole|bubble')) #remove contaminated filters


negative_controls <- edna_dat %>% 
  filter(str_detect(sample_num, '\\.4'))

eDNA_samples <- edna_dat %>% 
  filter(!str_detect(sample_num, '\\.4')) %>%  #remove negative controls
  select(copies_per_L,date, Site_Num) %>% 
  group_by(date, Site_Num) %>% 
  summarize(copies_per_L = mean(copies_per_L, na.rm=T))




p_samples <- ggplot(eDNA_samples, aes(x = copies_per_L))+
  geom_histogram(fill = "darkgreen", col = "black", alpha = 0.7, bins = 40)+
  labs(title = "A", x = "eDNA Concentration (Copies/L)", y = "Count")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold"))  # Make the title bold
  #coord_cartesian(ylim = c(0, 70))


summary(eDNA_samples$copies_per_L) #1379.7
sd(eDNA_samples$copies_per_L) #1705.66
max(eDNA_samples$copies_per_L) #1705.66



p_samples




c_samples <- ggplot(negative_controls, aes(x = copies_per_L))+
  geom_histogram(fill = "darkblue", col = "black", alpha = 0.7, bins = 40)+
  labs(title = "B", x = "eDNA Concentration (Copies/L)", y = element_blank())+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold"))  # Make the title bold
  #coord_cartesian(ylim = c(0, 70))



contamination_plot <- p_samples+c_samples
contamination_plot


summary(negative_controls$copies_per_L)
cbind(negative_controls$sample_num, negative_controls$copies_per_L)
View(negative_controls)
#9 samples had contamination, 7/9 occurred during fish sampling days



ggsave(plot= contamination_plot,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA_contamination.jpeg",
       dpi = 1000, 
       height = 3,
       width = 7,
       units = "in")






# Summary statistics ------------------------------------------------------



combined_dat <- read.csv("2022 Summer eDNA/combined_dat.csv")
#combined_dat <- r_dat[,-c(1,3)]


#Turn the quantitative values (copies per liter) into presence/absence

combined_dat <- combined_dat %>% 
  mutate(rep_1_copies_per_L = ifelse(rep_1_copies_per_L > 0, 1, 0),
         rep_2_copies_per_L = ifelse(rep_2_copies_per_L > 0, 1, 0),
         rep_3_copies_per_L = ifelse(rep_3_copies_per_L > 0, 1, 0),
         rep_4_copies_per_L = ifelse(rep_4_copies_per_L > 0, 1, 0),
         rep_5_copies_per_L = ifelse(rep_5_copies_per_L > 0, 1, 0)
  )




combined_dat <- combined_dat %>%
  mutate(
    temp_log = ifelse(is.na(temp_log), mean(temp_log, na.rm=T), temp_log),
    #ph_log = ifelse(ph_log == "NA", NA, ph_log), #I did this earlier on edna_dat
    ph_log = ifelse(is.na(ph_log), mean(ph_log, na.rm=T), ph_log),
    turb_log = ifelse(is.na(turb_log), mean(turb_log, na.rm=T), turb_log),
    flow = ifelse(is.na(flow), mean(flow, na.rm=T), flow)
    
  )





ggplot(combined_dat, aes(temp_log))+
  geom_histogram(aes(alpha = ))+
  theme_cowplot()


combined_dat$mean_MEANANNCMS



ggplot(combined_dat, aes(flow, mean_MEANANNCMS))+
  geom_point()+
  theme_cowplot()






# Occupancy model ---------------------------------------------------------




#?unmarked

#occu fits occurrence models with no linkage between abundance and detection (MacKenzie et al. 2002).

#occuRN fits abundance models to presence/absence data by exploiting the link between detection probability and abundance (Royle and Nichols 2003).

#occuFP fits occupancy models to data characterized by false negatives and false positive detections (e.g., Royle and Link [2006] and Miller et al. [2011]).

#pcount fits N-mixture models (aka binomial mixture models) to repeated count data (Royle 2004a, Kery et al 2005).

#gpcount fits the generalized N-mixture model described by Chandler et al. (2011) to repeated count data collected using the robust design.

#gmultmix fits a generalized form of the multinomial-mixture model of Royle (2004b) that allows for estimating availability and detection probability.

#?occu()

#?pcount


####______________________________________________
#Example from unmarked, copy this
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
plot(pferUMF, panels=4)
# add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))

# observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

(fm <- occu(~ obsvar1 ~ 1, pferUMF))

confint(fm, type='det', method = 'normal')
confint(fm, type='det', method = 'profile')

# estimate detection effect at obsvars=0.5
(lc <- linearComb(fm['det'],c(1,0.5)))

# transform this to probability (0 to 1) scale and get confidence limits
(btlc <- backTransform(lc))
confint(btlc, level = 0.9)

# Empirical Bayes estimates of proportion of sites occupied
re <- ranef(fm)
sum(bup(re, stat="mode"))







combined_dat <- read.csv("2022 Summer eDNA/combined_dat.csv")
#combined_dat <- r_dat[,-c(1,3)]


#Turn the quantitative values (copies per liter) into presence/absence

combined_dat <- combined_dat %>% 
  mutate(rep_1_copies_per_L = ifelse(rep_1_copies_per_L > 0, 1, 0),
         rep_2_copies_per_L = ifelse(rep_2_copies_per_L > 0, 1, 0),
         rep_3_copies_per_L = ifelse(rep_3_copies_per_L > 0, 1, 0),
         rep_4_copies_per_L = ifelse(rep_4_copies_per_L > 0, 1, 0),
         rep_5_copies_per_L = ifelse(rep_5_copies_per_L > 0, 1, 0)
  )




combined_dat <- combined_dat %>%
  mutate(
    temp_log = ifelse(is.na(temp_log), mean(temp_log, na.rm=T), temp_log),
    #ph_log = ifelse(ph_log == "NA", NA, ph_log), #I did this earlier on edna_dat
    ph_log = ifelse(is.na(ph_log), mean(ph_log, na.rm=T), ph_log),
    turb_log = ifelse(is.na(turb_log), mean(turb_log, na.rm=T), turb_log),
    flow = ifelse(is.na(flow), mean(flow, na.rm=T), flow)
    
  )



#An "unmarked" data frame for the presence-count model is created using the unmarkedFramePCount function. It takes three main arguments:

# y: The matrix of observed count data.
#siteCovs: The site-level covariates.
#obsCovs: The observation-level covariates.

#?unmarkedFramePCount

y <- data.frame(combined_dat[,4:8], digits = 0)


y <- y %>% 
  mutate(rep_1_copies_per_L = round(rep_1_copies_per_L/10),
         rep_2_copies_per_L = round(rep_2_copies_per_L/10),
         rep_3_copies_per_L = round(rep_3_copies_per_L/10),
  )




sitecovs <- data.frame(combined_dat[,c(10, 14:22, 67, 73:74, 77:107, 121:123)])  #May want to include the L filtered for each individual filters in the future

str(sitecovs)
#'data.frame':	62 obs. of  44 variables:
#  $ drainage        : chr  "Salcha" "Chena" "Salcha" "Salcha" ...
#$ weather         : chr  "Sunny" "Partially Cloudy" "Sunny" "Cloudy" ...
#$ temp_log        : num  7.21 7.31 5.86 8.08 8.73 ...
#$ ph_log          : num  7.45 7.73 6.74 7.45 3.2 ...
#$ sc_log          : num  NA 137.6 283.1 181.4 90.7 ...
#$ hdo_ml.L_log    : num  NA 12.3 12.2 10.5 11.3 ...
#$ hdo_perc_sat_log: num  NA 104.3 99.6 90.6 99.5 ...
#$ turb_log        : num  4.02 1.16 0.09 3.2 39.5 ...
#$ flow            : num  0.03 0.51 0.373 0.243 0.397 ...
#$ doy             : int  223 206 224 226 224 226 226 222 222 222 ...
#$ Year            : int  1997 1967 NA 1968 NA 1966 2016 1966 NA 1969 ...
#$ Burned_Numerical: int  1 1 0 1 0 1 1 1 0 1 ...
#$ Percent_burned  : int  34 70 0 1 0 81 34 88 0 68 ...
#$ mean_FitElev    : num  363 0 612 335 796 ...
#$ mean_ELEV_M     : num  361 341 606 334 791 ...
#$ mean_GRADIENT   : num  0.0437 0.0478 0.1173 0.0249 0.0871 ...
#$ max_MAX_GRAD_D  : num  0.184 0.278 0.241 0.15 0.258 ...
#$ mean_MNANPRC_M  : num  0.481 0.46 0.51 0.452 0.643 ...
#$ mean_MEANANNCMS : num  0.2757 0.129 0.0315 0.1937 0.0822 ...
#$ max_MEANANNCMS  : num  35.158 0.243 0.108 1.856 0.286 ...
#$ mean_WIDTH_M    : num  0.891 1.363 0.603 1.18 0.949 ...
#$ mean_DEPTH_M    : num  0.0529 0.2227 0.0457 0.0623 0.0587 ...
#$ mean_STRM_ORDER : num  2.45 2.33 1.53 2.57 1.89 ...
#$ max_STRM_ORDER  : num  5.55 5 3 5 3 ...
#$ mean_AZIMTH_DEG : num  173 165 136 188 286 ...
#$ mean_SINUOSITY  : num  0.99 0 0.982 1.013 1.005 ...
#$ mean_VAL_WIDTH  : num  46.2 94.6 14.4 70.4 16.3 ...
#$ mean_FP_WIDTH   : num  33.5 67.8 12.6 54.1 14 ...
#$ mean_VWI_Floor  : num  54 62.7 25.5 74.6 23.4 ...
#$ mean_ValCnstrnt : num  42.6 46.1 22.3 62.3 19.8 ...
#$ max_DROPMAX     : num  8.49 23.08 95.91 13.68 35.09 ...
#$ mean_IP_Chinook : num  0.1678 0.1956 0.0585 0.2519 0.089 ...
#$ mean_IP_Chum    : num  0.2947 0.3067 0.0346 0.4023 0.0442 ...
#$ max_Fish        : int  1 1 1 1 1 1 1 1 1 1 ...
#$ max_FISH_Anadr  : int  0 1 0 0 0 0 0 0 0 0 ...
#$ mean_FlowVel    : num  4.43e-03 1.66 8.54e-01 5.01e-07 1.02 ...
#$ max_FlowVel     : num  0.7775 2.7538 1.3544 0.0204 1.7924 ...
#$ mean_BFQ        : num  4.57e-02 5.55e-01 3.20e-02 5.98e-06 8.66e-02 ...
#$ max_BFQ         : num  8.015 2.383 0.109 0.243 0.361 ...
#$ mean_StrmPow    : num  1.18 2.02e+02 2.45e+01 5.86e-05 5.15e+01 ...
#$ max_StrmPow     : num  235.63 1139.67 60.53 2.38 251.09 ...
#$ max_BeavHab     : int  1 1 1 1 1 1 1 1 1 1 ...
#$ max_GEP_Cum     : num  0.353 0.346 0.509 0.353 0.568 ...
#$ mean_GEP        : num  0.181 0.166 0.381 0.146 0.317 ...
#$ Alldam_density  : num  0.237 0.107 0 0 0 ...
#$ Ydam_density    : num  0.237 0.107 0 0 0 ...
#$ Time_Since_Last_Burn: int  25 55 100 54 100 56 6 56 100 53 ...




sum(combined_dat$rep_1_copies_per_L, na.rm = T) #54 87.1%
sum(combined_dat$rep_2_copies_per_L, na.rm = T) #54 
sum(combined_dat$rep_3_copies_per_L, na.rm = T) #54 
sum(combined_dat$rep_4_copies_per_L, na.rm = T)
sum(combined_dat$rep_5_copies_per_L, na.rm = T)


#Right now we we'll just make these nothing, but we might want to include filter specific covariates (e.g., L filtered)
#m <- 62
#n <- 1
#obscovs <- list(x2 = round(matrix(runif(m * n), m, n)))

covar = rep(1, 62)
obscovs <- data.frame(covar)  


edna_matrix <- unmarkedFrameOccu(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)



#starting_values <- c(sigma = 6.721, turb_log = 3.9983, flow = 0.4935, covar = 1)

#Null model
m <- occu(formula = ~ 1 ~ 1, data = edna_matrix, se = T)

summary(m)





global <- occu(formula = ~ temp_log + turb_log + flow + doy ~ 1, data = edna_matrix, se = T)
summary(global) #AIC 322.0

m1  <- occu(formula = ~temp_log ~ 1 , data = edna_matrix, se = TRUE)
m1
summary( m1 )
#AIC 318.914860679275 

m2  <- occu(formula = ~turb_log ~ 1 , data = edna_matrix, se = TRUE)
summary( m2 )
#AIC 317.949065145275 

m3  <- occu(formula = ~flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m3 )
#AIC 318.631011619057 

m4  <- occu(formula = ~doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m4 )
#AIC 327.764892397531 

m5  <- occu(formula = ~temp_log + turb_log ~ 1 , data = edna_matrix, se = TRUE)
summary( m5 )
#AIC 319.346743480782 

m6  <- occu(formula = ~temp_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m6 )
#AIC 320.602494020965 

m7  <- occu(formula = ~temp_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m7 )
#AIC 320.775809424329 

m8  <- occu(formula = ~turb_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m8 )
#AIC 319.900496525724 

m9  <- occu(formula = ~turb_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m9 )
#AIC 319.244151678218 

m10  <- occu(formula = ~flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m10 )
#AIC 328.650639273185 

m11  <- occu(formula = ~temp_log + turb_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m11 )
#AIC 322.16139765398 

m12  <- occu(formula = ~temp_log + turb_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m12 )
#AIC 320.094118914255 

m13  <- occu(formula = ~temp_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m13 )
#AIC 322.317051370482 

m14  <- occu(formula = ~turb_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m14 )
#AIC 322.421279638573 

m15  <- occu(formula = ~temp_log + turb_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m15 )
#AIC 322.008857495693 



#You can also do that in a loop function

# Define the predictor variables
predictors <- c("temp_log", "turb_log", "flow", "doy")

# Create an empty list to store the models, AIC values, and predictors
models <- list()
aic_values <- numeric()
predictor_sets <- character()

# Create the global model
global_formula <- as.formula("~ temp_log + turb_log + flow + doy ~ 1")
global_model <- occu(formula = global_formula, data = edna_matrix, se = TRUE)
global_aic <- AIC(global_model)
cat("global <- occu(formula =", deparse(global_formula), ", data = edna_matrix, se = TRUE)\n")
cat("summary(global)\n")
cat(paste("#AIC", global_aic, "\n\n"))

# Loop through all possible combinations of predictor variables
for (i in 1:length(predictors)) {
  comb <- combn(predictors, i)
  for (j in 1:ncol(comb)) {
    formula <- as.formula(paste("~", paste(comb[, j], collapse = " + "), "~ 1"))
    model_name <- paste("m", length(models) + 1, sep = "")
    model <- occu(formula = formula, data = edna_matrix, se = TRUE)
    models[[model_name]] <- model
    aic_value <- AIC(model)
    aic_values <- c(aic_values, aic_value)
    predictor_set <- paste(comb[, j], collapse = " + ")
    predictor_sets <- c(predictor_sets, predictor_set)
    cat(paste(model_name, " <- occu(formula =", deparse(formula), ", data = edna_matrix, se = TRUE)\n"))
    cat(paste("summary(", model_name, ")\n"))
    cat(paste("#AIC", aic_value, "\n\n"))
  }
}

# Create a table summarizing model names, AIC values, and predictors
model_summary <- data.frame(Model = c("global", names(models)), AIC = c(global_aic, aic_values), Predictors = c("All", predictor_sets))

# Order the table by AIC values from smallest to largest
model_summary <- model_summary[order(model_summary$AIC), ]

# Print the sorted table
print(model_summary)


#Top selected model was m2 (turbidity)
#m2  <- occu(formula = ~turb_log ~ 1 , data = edna_matrix, se = TRUE)
#summary( m2 )
#AIC 317.949065145275 

#Model      AIC                       Predictors
#3      m2 317.9491                         turb_log
#4      m3 318.6310                             flow
#2      m1 318.9149                         temp_log
#10     m9 319.2442                   turb_log + doy
#6      m5 319.3467              temp_log + turb_log
#9      m8 319.9005                  turb_log + flow
#13    m12 320.0941        temp_log + turb_log + doy
#7      m6 320.6025                  temp_log + flow
#8      m7 320.7758                   temp_log + doy
#1  global 322.0089                              All
#16    m15 322.0089 temp_log + turb_log + flow + doy
#12    m11 322.1614       temp_log + turb_log + flow
#14    m13 322.3171            temp_log + flow + doy
#15    m14 322.4213            turb_log + flow + doy
#5      m4 327.7649                              doy
#11    m10 328.6506                       flow + doy




logLik(m2)
logLik(m3)
logLik(m1)
logLik(m9)
logLik(m5)
logLik(m8)
logLik(m12)
logLik(m6)
logLik(m7)
logLik(global)
logLik(m11)
logLik(m13)
logLik(m14)
logLik(m4)
logLik(m10)





#Turbidity and temperature are the most important things
















ggplot(combined_dat, aes(x=doy, y=temp_log))+
  geom_smooth()






#This will estimate the proportion of sites occupied 

INSERTTOPMODELHERE
ts2 <- occu(~1 ~cat,ts)
ts2


m2  <- occu(formula = ~turb_log ~ 1 , data = edna_matrix, se = TRUE)



backTransform(m2, 'state') #This will tell you the proportion of sites occupied
#97.8% of sites are occupied, +/- 2.53% SE







m2 <- occu(formula = ~ turb_log  ~ 1, data = edna_matrix, se = T)
summary(m2)
plogis(coef(m2, type="det")) # Should be close to p
confint(m2, type='det', method = 'normal')
confint(m2, type='det', method = 'profile')

backTransform(m2, 'state')
backTransform(linearComb(m2, coefficients = c(1,0), type = 'det'))
#Overall detectability was 0.742 +/- 0.0301 SE, with turb set to the mean











newData <- data.frame(temp_log = rep(seq(from = 2, to = 15, length.out = 10)))

z <- predict(m1,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

temp_plot <- ggplot(z, aes(temp_log, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Water Temperature (C)", y = "Detection probability")

temp_plot








newData <- data.frame(turb_log = rep(seq(from = 0, to = 61, length.out = 10)))

z <- predict(m2,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

turb_plot <- ggplot(z, aes(turb_log, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Turbidity (NTU)", y = "Detection probability")

turb_plot






newData <- data.frame(flow = rep(seq(from = 0, to = 1.136, length.out = 10)))

z <- predict(m3,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

flow_plot <- ggplot(z, aes(flow, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Flow Velocity (m/s)", y = "Detection probability")

flow_plot






newData <- data.frame(doy = rep(seq(from = 192, to = 226, length.out = 10)))

z <- predict(m4,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

doy_plot <- ggplot(z, aes(doy, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Day of Year", y = "Detection probability")

doy_plot





detection_occupancy_plot <- turb_plot | temp_plot | flow_plot | doy_plot
detection_occupancy_plot



ggsave(plot= detection_occupancy_plot,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/occupancy_detectability_rep1-3.jpeg",
       dpi = 1000, 
       height = 5,
       width = 20,
       units = "in")























m2 <- occu(formula = ~ temp_log + turb_log + ph_log + flow + doy + drainage ~ 1, data = edna_matrix, se = T)

#???????m2 <- occu(formula = ~ 1 ~ temp_log + turb_log + ph_log + flow + doy + drainage, data = edna_matrix, se = T)


summary(m2)

plogis(coef(m, type="det")) # Should be close to p

confint(m2, type='det', method = 'normal')
confint(m2, type='det', method = 'profile')



m3 <- occu(formula = ~ Burned_Numerical + Percent_burned + Time_Since_Last_Burn ~ 1, data = edna_matrix, se = T)


summary(m3)

plogis(coef(m3, type="det")) 

predicted_occu <- predict(m3, type = "state")
plot(predicted_occu)


confint(m3, type='det', method = 'normal')
confint(m3, type='det', method = 'profile')



m3 <- occu(formula = ~  ~ 1, data = edna_matrix, se = T)



m4 <- occu(formula = ~ temp_log + turb_log + flow + doy ~ 1, data = edna_matrix, se = T)


m4 <- occu(formula = ~ temp_log + turb_log + flow + doy ~ Burned_Numerical, data = edna_matrix, se = T) #+ Percent_burned + Time_Since_Last_Burn

summary(m4)

plogis(coef(m4, type="det")) 






# Create a data frame with covariate values for prediction
new_data <- data.frame(
  Burned_Numerical = combined_dat$Burned_Numerical,
  Percent_burned = combined_dat$Percent_burned,
  Time_Since_Last_Burn = combined_dat$Time_Since_Last_Burn)

# Predict occupancy probabilities using the new data
predicted_occupancy <- predict(m3, type = "occu", newdata = new_data)

predicted_occupancy <- fitted(m3, newdata = new_data)$Occu


predicted_occu <- predict(m3, type = "occu", newdata = edna_matrix)
plot(predicted_occu)





















# N Mixture Model (with abundance ranges) ---------------------------------



# Read in data and prep it for the model ----------------------------------





combined_dat <- read.csv("2022 Summer eDNA/combined_dat.csv")
#combined_dat <- r_dat[,-c(1,3)]

str(combined_dat, list.len = ncol(combined_dat))


#Round quantitative values 
#combined_dat <- combined_dat %>% 
#  mutate(rep_1_copies_per_L = round(rep_1_copies_per_L),
#         rep_2_copies_per_L = round(rep_2_copies_per_L),
#         rep_3_copies_per_L = round(rep_3_copies_per_L),
#         rep_4_copies_per_L = round(rep_4_copies_per_L),
#         rep_5_copies_per_L = round(rep_5_copies_per_L),
#         )




combined_dat <- combined_dat %>%
  mutate(
    temp_log = ifelse(is.na(temp_log), mean(temp_log, na.rm=T), temp_log),
    #ph_log = ifelse(ph_log == "NA", NA, ph_log), #I did this earlier on edna_dat
    ph_log = ifelse(is.na(ph_log), mean(ph_log, na.rm=T), ph_log),
    turb_log = ifelse(is.na(turb_log), mean(turb_log, na.rm=T), turb_log),
    flow = ifelse(is.na(flow), mean(flow, na.rm=T), flow)
    
  )



combined_dat <- combined_dat %>%
  mutate("log_Ydam" = (log(Ydam_density+1)),
         "exp_Ydam" = Ydam_density^2)

combined_dat$log_Ydam
combined_dat$Ydam_density




summary(combined_dat$mean_abundance) #1491.4
sd(combined_dat$mean_abundance) #1902.95
max(combined_dat$mean_abundance) #8771.18

combined_dat$mean_abundance



y <- data.frame(combined_dat[,4:6]) #, digits = 0 #Only including the first 3 replicates

y <- y %>% 
  mutate(rep_1_copies_per_L = round(rep_1_copies_per_L/10),
         rep_2_copies_per_L = round(rep_2_copies_per_L/10),
         rep_3_copies_per_L = round(rep_3_copies_per_L/10),
         )


#y <- y %>% 
#  mutate(rep_1_copies_per_L = ifelse(rep_1_copies_per_L > 0, 1, 0),
#         rep_2_copies_per_L = ifelse(rep_2_copies_per_L > 0, 1, 0),
#         rep_3_copies_per_L = ifelse(rep_3_copies_per_L > 0, 1, 0),
#         rep_4_copies_per_L = ifelse(rep_4_copies_per_L > 0, 1, 0),
#         rep_5_copies_per_L = ifelse(rep_5_copies_per_L > 0, 1, 0)
#         )

sitecovs <- data.frame(combined_dat[,c(9, 11, 16:23, 28:29, 68, 74:75, 78:108, 122:127)])  #May want to include the L filtered for each individual filters in the future



covar = rep(1, 62)
obscovs <- data.frame(covar)  


edna_matrix <- unmarkedFramePCount(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)



m.occu <- pcount(formula = ~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 1, data = edna_matrix, se = T)
summary(m.occu) #AIC 






# Select the predictors for DETECTION using AIC ---------------------------



# Define the predictor variables
predictors <- c("temp_log", "turb_log", "flow", "doy", "mean_liters_filtered_avg")

# Create an empty list to store the models, AIC values, predictors, and log-likelihood values
models <- list()
aic_values <- numeric()
predictor_sets <- character()
log_likelihood_values <- numeric()  # Added to store log-likelihood values

# Create the global model
global_formula <- as.formula("~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 1")
global_model <- pcount(formula = global_formula, data = edna_matrix, se = TRUE, K = 1165)
global_aic <- AIC(global_model)
cat("global <- pcount(formula =", deparse(global_formula), ", data = edna_matrix, se = TRUE)\n")
cat("summary(global)\n")
cat(paste("#AIC", global_aic, "\n\n"))

# Loop through all possible combinations of predictor variables
for (i in 1:length(predictors)) {
  comb <- combn(predictors, i)
  for (j in 1:ncol(comb)) {
    formula <- as.formula(paste("~", paste(comb[, j], collapse = " + "), "~ 1"))
    model_name <- paste("m", length(models) + 1, sep = "")
    model <- pcount(formula = formula, data = edna_matrix, se = TRUE)
    models[[model_name]] <- model
    aic_value <- AIC(model)
    aic_values <- c(aic_values, aic_value)
    predictor_set <- paste(comb[, j], collapse = " + ")
    predictor_sets <- c(predictor_sets, predictor_set)
    cat(paste(model_name, " <- pcount(formula =", deparse(formula), ", data = edna_matrix, se = TRUE)\n"))
    cat(paste("summary(", model_name, ")\n"))
    cat(paste("#AIC", aic_value, "\n"))
    
    # Calculate and store the log-likelihood
    log_likelihood <- logLik(model)
    log_likelihood_values <- c(log_likelihood_values, log_likelihood)
    cat(paste("Log Likelihood for", model_name, ":", log_likelihood, "\n"))
    
    cat("\n")
  }
}

# Create a table summarizing model names, AIC values, predictors, and log-likelihood values
model_summary <- data.frame(Model = c("global", names(models)), AIC = c(global_aic, aic_values), Predictors = c("All", predictor_sets), LogLikelihood = c(logLik(global_model), log_likelihood_values))

# Order the table by AIC values from smallest to largest
model_summary <- model_summary[order(model_summary$AIC), ]

# Print the sorted table
print(model_summary)


#Model      AIC                                                  Predictors LogLikelihood
#1  global 17984.49                                                         All     -8985.246
#32    m31 17984.49 temp_log + turb_log + flow + doy + mean_liters_filtered_avg     -8985.246
#28    m27 18035.03       temp_log + turb_log + flow + mean_liters_filtered_avg     -9011.516
#19    m18 18047.46              temp_log + turb_log + mean_liters_filtered_avg     -9018.731
#25    m24 18066.35                   turb_log + doy + mean_liters_filtered_avg     -9028.176
#24    m23 18067.14                  turb_log + flow + mean_liters_filtered_avg     -9028.572
#13    m12 18067.75                         turb_log + mean_liters_filtered_avg     -9029.876
#31    m30 18067.75            turb_log + flow + doy + mean_liters_filtered_avg     -9027.876
#30    m29 18285.66            temp_log + flow + doy + mean_liters_filtered_avg     -9136.830
#22    m21 18326.57                   temp_log + doy + mean_liters_filtered_avg     -9158.283
#21    m20 18399.36                  temp_log + flow + mean_liters_filtered_avg     -9194.682
#20    m19 18428.39                                       temp_log + flow + doy     -9209.195
#10     m9 18430.23                         temp_log + mean_liters_filtered_avg     -9211.117
#9      m8 18432.54                                              temp_log + doy     -9212.269
#16    m15 18462.86                              doy + mean_liters_filtered_avg     -9227.431
#26    m25 18464.92                       flow + doy + mean_liters_filtered_avg     -9227.459
#23    m22 18475.06                                       turb_log + flow + doy     -9232.528
#15    m14 18480.68                             flow + mean_liters_filtered_avg     -9236.340
#6      m5 18486.13                                    mean_liters_filtered_avg     -9240.065
#27    m26 18497.73                            temp_log + turb_log + flow + doy     -9242.866
#18    m17 18499.12                                   temp_log + turb_log + doy     -9244.562
#12    m11 18516.15                                              turb_log + doy     -9254.075
#17    m16 18517.50                                  temp_log + turb_log + flow     -9253.750
#14    m13 18532.80                                                  flow + doy     -9262.398
#5      m4 18533.23                                                         doy     -9263.613
#7      m6 18533.45                                         temp_log + turb_log     -9262.727
#2      m1 18584.50                                                    temp_log     -9289.248
#8      m7 18586.14                                             temp_log + flow     -9289.068
#3      m2 18589.30                                                    turb_log     -9291.649
#11    m10 18590.88                                             turb_log + flow     -9291.438
#4      m3 18604.16                                                        flow     -9299.079
#29    m28 20512.81        temp_log + turb_log + doy + mean_liters_filtered_avg    -10250.405


write.csv(model_summary, file = "2022 Summer eDNA/detection_model_summary.csv")



###Plot the relationships between detection and environmental covariates

#^^ moved down to plotting section



# Explore predictors using AIC --------------------------------------------






str(combined_dat, list.len = ncol(combined_dat))
#'data.frame':	62 obs. of  123 variables:
#  $ X                    : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Site_Num             : chr  "1" "1--3" "10" "109" ...
#$ X.x                  : int  1 4 5 6 7 8 9 10 11 12 ...
#$ rep_1_copies_per_L   : num  1 1 1 1 1 0 1 1 NA 1 ...
#$ rep_2_copies_per_L   : num  1 1 1 1 1 0 1 1 1 1 ...
#$ rep_3_copies_per_L   : num  1 1 1 1 1 0 1 1 1 1 ...
#$ rep_4_copies_per_L   : num  NA NA NA NA NA 1 1 NA 1 NA ...
#$ rep_5_copies_per_L   : num  NA NA NA NA NA NA 1 NA NA NA ...
#$ date                 : chr  "8/11/2022" "7/25/2022" "8/12/2022" "8/14/2022" ...
#$ drainage             : chr  "Salcha" "Chena" "Salcha" "Salcha" ...
#$ site                 : chr  "Cache Creek" "Colorado Creek" "Maiden Creek" "Redmond Creek" ...
#$ lat.dd               : num  64.6 64.9 64.7 64.7 64.7 ...
#$ long.dd              : num  146 147 145 147 145 ...
#$ weather              : chr  "Sunny" "Partially Cloudy" "Sunny" "Cloudy" ...
#$ temp_log             : num  7.21 7.31 5.86 8.08 8.73 ...
#$ ph_log               : num  7.45 7.73 6.74 7.45 3.2 ...
#$ sc_log               : num  NA 137.6 283.1 181.4 90.7 ...
#$ hdo_ml.L_log         : num  NA 12.3 12.2 10.5 11.3 ...
#$ hdo_perc_sat_log     : num  NA 104.3 99.6 90.6 99.5 ...
#$ turb_log             : num  4.02 1.16 0.09 3.2 39.5 ...
#$ flow                 : num  0.03 0.51 0.373 0.243 0.397 ...
#$ doy                  : int  223 206 224 226 224 226 226 222 222 222 ...
#$ OBJECTID             : int  24 61 19 10 16 27 11 26 25 12 ...
#$ Join_Count           : int  3 2 0 1 0 1 3 1 0 1 ...
#$ Stream_Name_HUC      : chr  "190803050608-Salcha River" "Fourmile Creek-Chena River" "Maiden Creek-Salcha River" "190803050806-Salcha River" ...
#$ Stream_Name_2        : chr  "190803050608-Salcha River" "Colorado Creek" "" "" ...
#$ VB_AreaSqKm          : num  4.223 9.352 0.689 22.606 1.588 ...
#$ Surveyed             : int  1 1 1 1 1 1 1 1 1 1 ...
#$ Ydam_Count           : int  1 1 0 0 0 0 2 22 0 12 ...
#$ Alldam_Count         : int  1 1 0 0 0 0 2 24 0 12 ...
#$ mean_alldam_Area_m2  : num  4059 1491 NA NA NA ...
#$ mean_Ydam_Area_m2    : num  4059 1491 NA NA NA ...
#$ AreaSqKm_1           : num  157 102 NA NA NA ...
#$ HUC12_1              : num  1.91e+11 1.91e+11 NA NA NA ...
#$ X.y                  : logi  NA NA NA NA NA NA ...
#$ HUType_1             : chr  "S" "M" "" "" ...
#$ HUMod_1              : chr  "NM" "IT" "" "" ...
#$ ToHUC_1              : num  1.91e+11 1.91e+11 NA NA NA ...
#$ Surveyed_1           : int  1 1 NA NA NA NA 1 1 NA 1 ...
#$ Shape_Area_12        : num  1.57e+08 1.02e+08 NA NA NA ...
#$ BURNED_AreaSqKm      : num  1.433 6.525 0 0.287 0 ...
#$ FIRE_NAME            : chr  "Butte Creek" "ANACONDA CREEK" "" "REDMOND CREEK" ...
#$ RECORDNUMB           : int  NA NA NA NA NA NA 468 NA NA NA ...
#$ ACRES                : num  3659 1492 NA 1242 NA ...
#$ AFSNUMBER            : chr  "B462" "Y58" "" "Z82" ...
#$ DOFNUMBER            : int  711462 NA NA NA NA NA 611468 NA NA NA ...
#$ USFSNUMBER           : chr  "P00462" " " "" " " ...
#$ ADDNUMBER            : chr  "4222" " " "" " " ...
#$ PERIMETERD           : chr  "2001-09-15" "1967-11-01" "" "1968-11-01" ...
#$ LATESTPERI           : chr  "Yes" "Yes" "" "Yes" ...
#$ SOURCE               : chr  "Image Interpretation" "Not Recorded" "" "Digitized" ...
#$ SOURCEMETH           : chr  " " " " "" "Image" ...
#$ SOURCECLAS           : chr  " " " " "" " " ...
#$ AGENCYACRE           : num  0 0 NA 0 NA 0 0 0 NA 0 ...
#$ COMMENTS             : chr  "Fire perimeter modified using Landsat Image Service." "Digitized from 1978 Landsat image P75 R14." "" "Digitized from 1976 Landsat image P73 R15." ...
#$ FIREID               : chr  "19892" "26150" "" "26135" ...
#$ FIREYEAR             : int  1997 1967 NA 1968 NA 1966 2016 1966 NA 1969 ...
#$ UPDATETIME           : chr  "2019-03-14" "2008-03-12" "" "2008-03-12" ...
#$ UPDATEUSER           : chr  "jljenkins" "jhrobak" "" "jhrobak" ...
#$ USEDONFINA           : chr  " " "No" "" "No" ...
#$ FPOUTDATE            : chr  "1997-08-29" "1958-07-26" "" "" ...
#$ FPMERGEDDA           : chr  "" "" "" "" ...
#$ IRWINID              : chr  " " " " "" " " ...
# PRESCRIBED           : chr  "N" "N" "" "N" ...
#$ FIRESEASON           : chr  "1990 - 1999" "1960 - 1969" "" "1960 - 1969" ...
#$ Shape_Leng           : num  21570 10060 NA 10627 NA ...
#$ Year                 : int  1997 1967 NA 1968 NA 1966 2016 1966 NA 1969 ...
#$ Years                : int  1990 1960 NA 1960 NA 1960 2010 1960 NA 1960 ...
#$ Shape_Le_1           : num  21570 10060 NA 10627 NA ...
#$ AltMode              : chr  "0" "" "" "" ...
#$ Shape_Length_1       : chr  "21570.22454" "10059.62716" "" "10626.93625" ...
#$ Shape_Area_1         : chr  "14807261.84" "6037881.849" "" "5025774.838" ...
#$ Burned_Numerical     : int  1 1 0 1 0 1 1 1 0 1 ...
#$ Percent_burned       : int  34 70 0 1 0 81 34 88 0 68 ...
#$ sum_LENGTH_M         : num  28646 84255 11411 150229 26864 ...
#$ sum_AREA_SQKM        : num  13129 10437 549 44806 2366 ...
#$ mean_FitElev         : num  363 0 612 335 796 ...
#$ mean_ELEV_M          : num  361 341 606 334 791 ...
#$ mean_GRADIENT        : num  0.0437 0.0478 0.1173 0.0249 0.0871 ...
#$ max_MAX_GRAD_D       : num  0.184 0.278 0.241 0.15 0.258 ...
#$ mean_MNANPRC_M       : num  0.481 0.46 0.51 0.452 0.643 ...
#$ mean_MEANANNCMS      : num  0.2757 0.129 0.0315 0.1937 0.0822 ...
#$ max_MEANANNCMS       : num  35.158 0.243 0.108 1.856 0.286 ...
#$ mean_WIDTH_M         : num  0.891 1.363 0.603 1.18 0.949 ...
#$ mean_DEPTH_M         : num  0.0529 0.2227 0.0457 0.0623 0.0587 ...
#$ mean_STRM_ORDER      : num  2.45 2.33 1.53 2.57 1.89 ...
#$ max_STRM_ORDER       : num  5.55 5 3 5 3 ...
#$ mean_AZIMTH_DEG      : num  173 165 136 188 286 ...
#$ mean_SINUOSITY       : num  0.99 0 0.982 1.013 1.005 ...
#$ mean_VAL_WIDTH       : num  46.2 94.6 14.4 70.4 16.3 ...
#$ mean_FP_WIDTH        : num  33.5 67.8 12.6 54.1 14 ...
#$ mean_VWI_Floor       : num  54 62.7 25.5 74.6 23.4 ...
#$ mean_ValCnstrnt      : num  42.6 46.1 22.3 62.3 19.8 ...
#$ max_DROPMAX          : num  8.49 23.08 95.91 13.68 35.09 ...
#$ mean_IP_Chinook      : num  0.1678 0.1956 0.0585 0.2519 0.089 ...
#$ mean_IP_Chum         : num  0.2947 0.3067 0.0346 0.4023 0.0442 ...
#$ max_Fish             : int  1 1 1 1 1 1 1 1 1 1 ...
#$ max_FISH_Anadr       : int  0 1 0 0 0 0 0 0 0 0 ...
#$ mean_FlowVel         : num  4.43e-03 1.66 8.54e-01 5.01e-07 1.02 ...
#$ max_FlowVel          : num  0.7775 2.7538 1.3544 0.0204 1.7924 ...
#$ mean_BFQ             : num  4.57e-02 5.55e-01 3.20e-02 5.98e-06 8.66e-02 ...
#$ max_BFQ              : num  8.015 2.383 0.109 0.243 0.361 ...
#$ mean_StrmPow         : num  1.18 2.02e+02 2.45e+01 5.86e-05 5.15e+01 ...
#$ max_StrmPow          : num  235.63 1139.67 60.53 2.38 251.09 ...
#$ max_BeavHab          : int  1 1 1 1 1 1 1 1 1 1 ...
#$ max_GEP_Cum          : num  0.353 0.346 0.509 0.353 0.568 ...
#$ mean_GEP             : num  0.181 0.166 0.381 0.146 0.317 ...
#$ sum_Length_KILOMETERS: num  27.4 79 10.8 142 25.4 ...
#$ Polyline_Count       : int  294 838 115 1498 263 253 1993 1135 210 2391 ...
#$ Join_Count_12        : chr  "4" "3" "2" "1" ...
#$ eDNA_Site_Name       : chr  "WS EDNA 1" "WS EDNA SITE 1A" "WS EDNA 10" "WS EDNA 109" ...
#$ Sampletype           : chr  "Test" "Train" "Test" "Test" ...
#$ Elevation            : num  288 0 431 209 435 ...
#$ DateTime             : chr  "2022-08-11" "" "2022-08-13" "2022-08-14" ...
#$ POINT_X              : num  -146 -147 -145 -147 -145 ...
#$ POINT_Y              : num  64.6 64.9 64.7 64.5 64.7 ...
#$ POINT_Z              : num  2.88e+02 -4.70e-05 4.31e+02 2.09e+02 4.35e+02 ...
#$ HUC12                : num  1.91e+11 1.91e+11 1.91e+11 1.91e+11 1.91e+11 ...
#$ Shape_Length         : num  67517 194688 32612 354693 76461 ...
#$ Shape_Area           : num  4222942 9352427 688862 22605899 1587973 ...
#$ Alldam_density       : num  0.237 0.107 0 0 0 ...
#$ Ydam_density         : num  0.237 0.107 0 0 0 ...
#$ Time_Since_Last_Burn : int  25 55 100 54 100 56 6 56 100 53 ...





#Plot all  the predictor relationships to start
combined_dat$rep_1_copies_per_L

combined_dat <- combined_dat %>% 
  mutate("mean_abundance2" = (rep_1_copies_per_L + rep_2_copies_per_L + rep_3_copies_per_L)/3)


combined_dat$mean_abundance2




predictors <- c("temp_log", "flow", "doy", "turb_log", "mean_liters_filtered_avg", "Surveyed", "Ydam_density", 
                "Alldam_density", "Burned_Numerical", "Percent_burned", "Time_Since_Last_Burn",
                "mean_ELEV_M", "mean_GRADIENT", "mean_MEANANNCMS", "max_MAX_GRAD_D", "mean_WIDTH_M", 
                "mean_DEPTH_M", "max_STRM_ORDER", "mean_SINUOSITY", "VB_AreaSqKm", 
                "mean_StrmPow")


predictor_labels <- c("Temperature", "Flow", "DOY", "Turbidity", "mean_liters_filtered_avg", "Surveyed", 
                      "Ydam_density", "Alldam_density", "Burned_Numerical", "Percent_burned",
                      "Time_Since_Last_Burn", 
                      "mean_ELEV_M", "mean_GRADIENT", "mean_MEANANNCMS", "max_MAX_GRAD_D", "mean_WIDTH_M",
                      "mean_DEPTH_M", "max_STRM_ORDER", "mean_SINUOSITY", "VB_AreaSqKm",
                      "mean_StrmPow")

# Create an empty list to store the plots
plot_list <- list()

# Loop through each predictor and create a scatterplot
#for (i in 1:length(predictors)) {
#  predictor <- predictors[i]
# label <- predictor_labels[i]

#  plot <- ggplot(combined_dat, aes_string(x = predictor, y = "mean_abundance")) +
#    geom_point() +
#    geom_smooth()+
#    labs(x = label, y = "Mean Abundance") +
#    theme_bw()
#
#  plot_list[[predictor]] <- plot
#}

# Arrange the plots in a grid
#grid.arrange(grobs = plot_list, ncol = 4)  #


# Loop through each predictor and create a scatterplot
for (i in 1:length(predictors)) {
  predictor <- predictors[i]
  label <- predictor_labels[i]
  
  plot <- ggplot(combined_dat, aes_string(x = predictor, y = "mean_abundance2")) +
    geom_point() +
    geom_smooth(method = "lm")+
    labs(x = label, y = "Mean Abundance") +
    theme_bw()
  
  plot_list[[predictor]] <- plot
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, ncol = 4) 



combined_dat$Ydam_density

ggplot(combined_dat, aes(x = Ydam_density, y = mean_abundance)) +
  geom_point() +
  geom_smooth()+
  labs(x = "Ydam_density", y = "Mean Abundance") +
  theme_bw()


ggplot(combined_dat, aes(x = log(Ydam_density), y = mean_abundance)) +
  geom_point() +
  geom_smooth()+
  labs(x = "Ydam_density", y = "Mean Abundance") +
  theme_bw()



x_values <- seq(min(combined_dat$Ydam_density), max(combined_dat$Ydam_density), length.out = 100)
y_values <- 8000 * exp(-1.5 * x_values)

ggplot(combined_dat, aes(x = Ydam_density, y = mean_abundance)) +
  geom_point() +
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "darkblue", lwd = 1.5, alpha = 0.6) +
  labs(x = "Ydam_density", y = "Mean Abundance") +
  theme_bw()








#NEED TO DO VIF!!---------------------------------------------------------------

VIF.dat <- combined_dat %>% 
  select(#Wildfire
          Burned_Numerical, Percent_burned, Time_Since_Last_Burn,
           #Beaver
           Ydam_density, # Ydam_density*Surveyed,
           #Geomorphology
           max_STRM_ORDER, mean_SINUOSITY, VB_AreaSqKm, max_MAX_GRAD_D, #mean_DEPTH_M,  mean_WIDTH_M,
           #Hydrology
           mean_MEANANNCMS, mean_StrmPow)




vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

postVIF_func <- vif_func(VIF.dat)
postVIF_func
#"Burned_Numerical"     
#"Percent_burned"       
#"Time_Since_Last_Burn" 
#"Ydam_density"        
#"mean_WIDTH_M"         
#"max_STRM_ORDER"       
#"mean_SINUOSITY"       
#"VB_AreaSqKm"         
#"max_MAX_GRAD_D"       
#"mean_StrmPow"


VIF<-function(X) {
  #Computes Variance Inflation Factors for a Predictor Matrix
  #INPUTS:
  #X is a matrix (or data frame) of the predictors (no column of ones).
  cat("REMINDER: Your input matrix should not include the response\n")
  a<-1/sqrt(dim(X)[1]-1)*scale(X)
  b<-cbind(diag(solve(t(a)%*%a)))
  dimnames(b)<-list(dimnames(X)[[2]],"VIF")
  return(b)
}

VIF(VIF.dat)

#VIF
#Burned_Numerical     5.261004
#Percent_burned       1.936812
#Time_Since_Last_Burn 4.875727
#Ydam_density         1.142997
#max_STRM_ORDER       1.922846
#mean_SINUOSITY       6.132428
#VB_AreaSqKm          1.983526
#max_MAX_GRAD_D       3.148857
#mean_MEANANNCMS      1.268992
#mean_StrmPow         8.586885

#Did not include stream width and depth due to collinearity with MEANANNCMS and SINUOSITY












#Include temp, flow, and doy in all models

#m <- pcount(~ 1 ~ 1, data = edna_matrix, K = 1165, se = T)
#summary(m)


#Consider interaction terms between beaver density and the basins surveyed
#beaver <- pcount(~ temp_log + flow + doy ~ Ydam_density + Ydam_density:Surveyed
#                    , data = edna_matrix, K = 1165, se = T)
#summary(beaver)





m.occu <- pcount(~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg  
                   #+ turb_log*mean_liters_filtered_avg 
                   ~ 1
                    , data = edna_matrix, K = 1165, se = T)
summary(m.occu)


# Extract coefficients
coefficients <- coef(m.occu)

# Extract standard errors
standard_errors <- sqrt(diag(vcov(m.occu)))

# Z-Score for a 90% confidence interval
z_score <- qnorm(0.95)  # 0.95 corresponds to 95% confidence level

# Calculate lower and upper bounds of the confidence intervals
lower_bounds <- coefficients - z_score * standard_errors
upper_bounds <- coefficients + z_score * standard_errors

# Create a data frame to store the results
results <- data.frame(Parameter = names(coefficients),
                      Estimate = coefficients,
                      Lower_CI = lower_bounds,
                      Upper_CI = upper_bounds)

# Print the results
print(results)

















m.global <- pcount(
                  #Occupancy factors
                  ~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg  ~ 
                  #Fire
                  Burned_Numerical + Percent_burned + Time_Since_Last_Burn +
                  #Beaver
                  Ydam_density + # Ydam_density*Surveyed+
                  #Geomorphology
                  max_STRM_ORDER + mean_SINUOSITY + VB_AreaSqKm + max_MAX_GRAD_D + # mean_WIDTH_M + mean_DEPTH_M +  
                  #Hydrology
                  mean_MEANANNCMS + mean_StrmPow, 
            data = edna_matrix, K = 1165, se = T)

summary(m.global)


#Abundance (log-scale):
#Estimate       SE      z  P(>|z|)
#(Intercept)           3.64959 0.095616  38.17 0.00e+00
#Burned_Numerical      0.80291 0.063982  12.55 4.03e-36
#Percent_burned       -0.00207 0.000370  -5.59 2.21e-08
#Time_Since_Last_Burn  0.00739 0.000672  11.01 3.59e-28
#Ydam_density          0.11825 0.010905  10.84 2.14e-27
#max_STRM_ORDER        0.24957 0.015167  16.46 7.66e-61
#mean_SINUOSITY       -0.85814 0.064170 -13.37 8.70e-41
#VB_AreaSqKm           0.01357 0.000648  20.95 1.80e-97
#max_MAX_GRAD_D        2.05754 0.277065   7.43 1.12e-13
#mean_MEANANNCMS       0.51628 0.038294  13.48 2.00e-41
#mean_StrmPow         -0.00514 0.000335 -15.35 3.76e-53

#Detection (logit-scale):
#  Estimate      SE      z   P(>|z|)
#(Intercept)               14.9241 0.55140  27.07 2.49e-161
#temp_log                  -0.2017 0.01188 -16.98  1.22e-64
#turb_log                   0.0315 0.00858   3.67  2.43e-04
#flow                       0.0821 0.08102   1.01  3.11e-01
#doy                       -0.0607 0.00222 -27.35 1.07e-164
#mean_liters_filtered_avg  -0.1180 0.03195  -3.69  2.23e-04

#AIC: 15588.74 
#Number of sites: 62
#optim convergence code: 0
#optim iterations: 319 
#Bootstrap iterations: 0 






# Extract coefficients
coefficients <- coef(m.global)

# Extract standard errors
standard_errors <- sqrt(diag(vcov(m.global)))

# Z-Score for a 90% confidence interval
z_score <- qnorm(0.95)  # 0.95 corresponds to 95% confidence level

# Calculate lower and upper bounds of the confidence intervals
lower_bounds <- coefficients - z_score * standard_errors
upper_bounds <- coefficients + z_score * standard_errors

# Create a data frame to store the results
results <- data.frame(Parameter = names(coefficients),
                      Estimate = coefficients,
                      Lower_CI = lower_bounds,
                      Upper_CI = upper_bounds)

# Print the results
print(results)

#Parameter     Estimate     Lower_CI     Upper_CI
#lam(Int)                                   lam(Int)  3.953319089  3.800181325  4.106456854
#lam(Burned_Numerical)         lam(Burned_Numerical)  0.601177482  0.497206334  0.705148631
#lam(Percent_burned)             lam(Percent_burned) -0.001622954 -0.002229070 -0.001016838
#lam(Time_Since_Last_Burn) lam(Time_Since_Last_Burn)  0.005824518  0.004730243  0.006918794
#lam(Ydam_density)                 lam(Ydam_density)  0.108111752  0.090098901  0.126124602
#lam(max_STRM_ORDER)             lam(max_STRM_ORDER)  0.233744458  0.208822089  0.258666827
#lam(mean_SINUOSITY)             lam(mean_SINUOSITY) -0.851131877 -0.955393565 -0.746870188
#lam(VB_AreaSqKm)                   lam(VB_AreaSqKm)  0.014047427  0.012981338  0.015113516
#lam(max_MAX_GRAD_D)             lam(max_MAX_GRAD_D)  1.458712532  1.011997837  1.905427228
#lam(mean_MEANANNCMS)           lam(mean_MEANANNCMS)  0.587531731  0.524583319  0.650480142
#lam(mean_StrmPow)                 lam(mean_StrmPow) -0.004742126 -0.005283588 -0.004200664
#p(Int)                                       p(Int) 11.834322961 10.982497378 12.686148544
#p(temp_log)                             p(temp_log) -0.129811738 -0.151289741 -0.108333736
#p(flow)                                     p(flow) -0.284755375 -0.409205844 -0.160304905
#sp(doy)                                       p(doy) -0.049361736 -0.053051247 -0.045672225













m.fire <- pcount(~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 
                   #Fire
                   Burned_Numerical + Percent_burned + Time_Since_Last_Burn
            , data = edna_matrix, K = 1165, se = T)
summary(m.fire)


m.geo <- pcount(~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 
                #Geomorphology
                max_STRM_ORDER + mean_SINUOSITY + VB_AreaSqKm + max_MAX_GRAD_D, #  + mean_WIDTH_M + mean_DEPTH_M +  
                data = edna_matrix, K = 1165, se = T)
summary(m.geo)



m.hydro <- pcount(~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 
                #Hydrology
                mean_MEANANNCMS + mean_StrmPow, 
                data = edna_matrix, K = 1165, se = T)
summary(m.hydro)



m.beaver <- pcount(~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 
                  #Beaver
                  Ydam_density, 
                  data = edna_matrix, K = 1165, se = T)
summary(m.beaver)

# Extract coefficients
coefficients <- coef(m.beaver)
# Extract standard errors
standard_errors <- sqrt(diag(vcov(m.beaver)))
# Z-Score for a 90% confidence interval
z_score <- qnorm(0.95)  # 0.95 corresponds to 95% confidence level
# Calculate lower and upper bounds of the confidence intervals
lower_bounds <- coefficients - z_score * standard_errors
upper_bounds <- coefficients + z_score * standard_errors
# Create a data frame to store the results
results <- data.frame(Parameter = names(coefficients),
                      Estimate = coefficients,
                      Lower_CI = lower_bounds,
                      Upper_CI = upper_bounds)

# Print the results
print(results)






#predict(m.beaver)



#log_Ydam

#m.beaver2 <- pcount(~ temp_log + flow + doy ~ 
                     #Beaver
#                      log_Ydam, 
#                   data = edna_matrix, K = 1165, se = T)
#summary(m.beaver2)



#m.beaver2 <- pcount(~ temp_log + flow + doy ~ 
#                      #Beaver
#                      exp_Ydam, 
#                    data = edna_matrix, K = 1165, se = T)
#summary(m.beaver2)







# Model averaging (decided not to do this because of the AIC weight ------------

  

Cand.mod <- list(m.occu, m.global, m.fire, m.geo, m.hydro, m.beaver) 

Modnames <- c("m.occu", "m.global", "m.fire", "m.geo", "m.hydro", "m.beaver") 

aictab(Cand.mod, Modnames)



#First look at detectability covariates

AICcmodavg::modavg(parm = "temp_log", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "detect")
#Model-averaged estimate: -0.13
#Unconditional SE: 0.01
#90% Unconditional confidence interval: -0.15, -0.11


AICcmodavg::modavg(parm = "flow", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "detect")
#Model-averaged estimate: -0.28
#Unconditional SE: 0.08
#90% Unconditional confidence interval: -0.41, -0.16


AICcmodavg::modavg(parm = "doy", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "detect")
#Model-averaged estimate: -0.05
#Unconditional SE: 0
#90% Unconditional confidence interval: -0.05, -0.05



#Then look at abundance covariates

AICcmodavg::modavg(parm = "max_MAX_GRAD_D", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 1.46
#Unconditional SE: 0.27
#90% Unconditional confidence interval: 1.01, 1.91



AICcmodavg::modavg(parm = "VB_AreaSqKm", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0.01
#Unconditional SE: 0
#90% Unconditional confidence interval: 0.01, 0.02



AICcmodavg::modavg(parm = "mean_SINUOSITY", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: -0.85
#Unconditional SE: 0.06
#90% Unconditional confidence interval: -0.96, -0.75



AICcmodavg::modavg(parm = "max_STRM_ORDER", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0.23
#Unconditional SE: 0.02
#90% Unconditional confidence interval: 0.21, 0.26



AICcmodavg::modavg(parm = "mean_StrmPow", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0
#Unconditional SE: 0
#90% Unconditional confidence interval: -0.01, 0



AICcmodavg::modavg(parm = "mean_MEANANNCMS", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0.59
#Unconditional SE: 0.04
#90% Unconditional confidence interval: 0.52



AICcmodavg::modavg(parm = "Time_Since_Last_Burn", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0.01
#Unconditional SE: 0
#90% Unconditional confidence interval: 0, 0.01



AICcmodavg::modavg(parm = "Percent_burned", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0
#Unconditional SE: 0
#90% Unconditional confidence interval: 0, 0



AICcmodavg::modavg(parm = "Ydam_density", cand.set = Cand.mod, modnames = Modnames, conf.level = 0.90, parm.type = "lambda")
#Model-averaged estimate: 0.11
#Unconditional SE: 0.01
#90% Unconditional confidence interval: 0.09, 0.13


















#Summary information about the best predicting model

summary(m.occu)
plogis(coef(m.occu, type="det")) # Should be close to p
confint(m.occu, type='det', method = 'normal')
confint(m.occu, type='det', method = 'profile')
#               0.025       0.975
#p(Int)       5.86697400  8.48411663
#p(temp_log) -0.17319633 -0.10220529
#p(flow)     -0.25050556  0.13288192
#p(doy)      -0.03238909 -0.02122046

backTransform(m.occu, 'state')
backTransform(linearComb(m.occu, coefficients = c(1,0,0,0,0), type = 'det'))
#Overall detectability was 0.742 +/- 0.0301 SE, with turb set to the mean



summary(m.global)
plogis(coef(m.global, type="det")) # Should be close to p
confint(m.global, type='det', method = 'normal')
confint(m.global, type='det', method = 'profile')

backTransform(m.global, 'state')
backTransform(linearComb(m.global, coefficients = c(1,0), type = 'det'))
#Overall detectability was 0.742 +/- 0.0301 SE, with turb set to the mean


modSel(), crossVal()

predict(), ranef(), posteriorSamples()


#####detection probabilities by observer 1=Nate, 2=Brian plus 95% CI's
timeinvar <- backTransform(linearComb(lwd_nmix_pois_t, c(1, 1), type="det"))
timevar <- backTransform(linearComb(lwd_nmix_pois_t, c(1, 2), type="det"))

timeinvar_CI <- confint(backTransform(linearComb(lwd_nmix_pois_t, c(1, 1), type="det")))
timevar_CI <- confint(backTransform(linearComb(lwd_nmix_pois_t, c(1, 2), type="det")))











# Plotting ----------------------------------------------------------------






combined_dat <- read.csv("2022 Summer eDNA/combined_dat.csv")
#combined_dat <- r_dat[,-c(1,3)]

str(combined_dat, list.len = ncol(combined_dat))


combined_dat <- combined_dat %>%
  mutate(
    temp_log = ifelse(is.na(temp_log), mean(temp_log, na.rm=T), temp_log),
    #ph_log = ifelse(ph_log == "NA", NA, ph_log), #I did this earlier on edna_dat
    ph_log = ifelse(is.na(ph_log), mean(ph_log, na.rm=T), ph_log),
    turb_log = ifelse(is.na(turb_log), mean(turb_log, na.rm=T), turb_log),
    flow = ifelse(is.na(flow), mean(flow, na.rm=T), flow)
    
  )



combined_dat <- combined_dat %>%
  mutate("log_Ydam" = (log(Ydam_density+1)),
         "exp_Ydam" = Ydam_density^2)

combined_dat$log_Ydam
combined_dat$Ydam_density


y <- data.frame(combined_dat[,4:6]) #, digits = 0 #Only including the first 3 replicates

y <- y %>% 
  mutate(rep_1_copies_per_L = round(rep_1_copies_per_L/10),
         rep_2_copies_per_L = round(rep_2_copies_per_L/10),
         rep_3_copies_per_L = round(rep_3_copies_per_L/10),
  )


sitecovs <- data.frame(combined_dat[,c(9, 11, 16:23, 28:29, 68, 74:75, 78:108, 122:127)])  #May want to include the L filtered for each individual filters in the future


covar = rep(1, 62)
obscovs <- data.frame(covar)  


edna_matrix <- unmarkedFramePCount(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)













m.occu  <- pcount(formula = ~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg
                  ~ 1 , data = edna_matrix, se = TRUE)
summary( m.occu )
#AIC 17984.49


###Plot the relationships between detection and environmental covariates


min(combined_dat$turb_log) #0
max(combined_dat$turb_log) #61.03
min(combined_dat$temp_log) #2.88
max(combined_dat$temp_log) #15.1
min(combined_dat$flow) #0
max(combined_dat$flow) #1.136
min(combined_dat$doy) #192
max(combined_dat$doy) #226
min(combined_dat$mean_liters_filtered_avg) #0.5
max(combined_dat$mean_liters_filtered_avg) #5.166

#newData <- data.frame(turb_log = rep(seq(from = 0, to = 61, length.out = 10)),
#                      temp_log = rep(seq(from = 0, to = 2, length.out = 10)),
#                      flow = rep(seq(from = 0, to = 1.15, length.out = 10)),
#                      doy = rep(seq(from = 192, to = 226, length.out = 10)))





newData <- data.frame(temp_log = rep(seq(from = 2, to = 15, length.out = 10)),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10))




z <- predict(m.occu,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

temp_plot <- ggplot(z, aes(temp_log, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  #scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, by = 0.1)) +
  #scale_x_continuous(limits = c(2, 14), breaks = seq(2, 14, by = 2)) +
    labs(x = "Water Temperature (C)", y = "Detection Probability")

temp_plot







newData <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(seq(from = 0, to = 61, length.out = 10)),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10))


z <- predict(m.occu,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

turb_plot <- ggplot(z, aes(turb_log, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Turbidity (NTU)", y = element_blank())

turb_plot






newData <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(seq(from = 0, to = 1.136, length.out = 10)),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10))


z <- predict(m.occu,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

flow_plot <- ggplot(z, aes(flow, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Flow Velocity (m/s)", y = element_blank())

flow_plot






newData <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(seq(from = 192, to = 226, length.out = 10)),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10))



z <- predict(m.occu,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

doy_plot <- ggplot(z, aes(doy, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Day of Year", y = element_blank())

doy_plot







newData <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(seq(from = 0.5, to = 5.1, length.out = 10)))


z <- predict(m.occu,type='det',newdata=newData,appendData=T)
z

theme_set(theme_bw(base_size = 24)) 

L_plot <- ggplot(z, aes(mean_liters_filtered_avg, Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size=2) +
  labs(x = "Liters Filtered", y = element_blank())

L_plot







detection_plot <-  temp_plot | turb_plot | flow_plot | doy_plot | L_plot
detection_plot



ggsave(plot= detection_plot,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/N-mixture_detectability2.jpeg",
       dpi = 1000, 
       height = 5,
       width = 23,
       units = "in")






combined_dat$mean_liters_filtered_avg


#ggplot(combined_dat, aes(turb_log, mean_liters_filtered_avg))+
#  geom_point()



ggplot(combined_dat, aes(doy, temp_log))+
  geom_point()


summary(lm(temp_log~doy, data= combined_dat))


ggplot(combined_dat, aes(doy, temp_log))+
  geom_point()












#Panel plot VV none of these worked, so I do it manually below
new_data <- list()


for (col_name in names(sitecovs)) {
  # Check if the column is numeric
  if (is.numeric(sitecovs[[col_name]])) {
    # Calculate the minimum and maximum values
    min_val <- min(sitecovs[[col_name]], na.rm = TRUE)
    max_val <- max(sitecovs[[col_name]], na.rm = TRUE)
    
    # Generate num_values values between the minimum and maximum
    new_values <- seq(from = min_val, to = max_val, length.out = 10)
  }
    new_data[[col_name]] <- new_values
}

new_data <- data.frame(new_data)
new_data





# Define the list of predictors you want to plot
predictors <- c("temp_log", "turb_log", "flow", "doy", "mean_liters_filtered_avg",
                "max_MAX_GRAD_D", "VB_AreaSqKm", "mean_SINUOSITY", "max_STRM_ORDER",
                "mean_StrmPow", "mean_MEANANNCMS", 
                "Time_Since_Last_Burn",  "Percent_burned", "Burned_Numerical", 
                "Ydam_density"   
)

predictor_labels <- c("Temperature", "Turbidity", "Flow", "DOY", "Liters Filtered", 
                      "Gradient (%)", "Valley Bottom Area", "Sinuosity", "Stream Order",
                      "Stream Power", "Mean Annual Discharge",
                      "Time_Since_Last_Burn (years)",  "Percent Burned", "Burned?",
                      "Beaver Dam Density"
)

# Create an empty list to store the plots
plot_list <- list()

# Predict the response variable using the model for all predictors
for (predictor in predictors) {
  
  #Create the new data
  if (is.numeric(sitecovs[[col_name]])) {
    # Calculate the minimum and maximum values
    min_val <- min(sitecovs[[col_name]], na.rm = TRUE)
    max_val <- max(sitecovs[[col_name]], na.rm = TRUE)
    
    # Generate num_values values between the minimum and maximum
    new_values <- seq(from = min_val, to = max_val, length.out = 10)
  }
  
  # Predict the response variable using the model
  z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)
  
  # Create the plot for the current predictor and store it in the list
  plot <- ggplot(z, aes(x = .data[[predictor]], y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
    geom_line(size = 1) +
    labs(x = predictor_labels[predictor == predictors], y = "eDNA concentration") +
    theme_bw()
  
  plot_list[[predictor]] <- plot
}

# Combine the plots into a single panel
combined_plot <- wrap_plots(plotlist = plot_list)

# Print the combined plot
combined_plot












#Create the new dataset
new_data <- list()


for (col_name in names(sitecovs)) {
  # Check if the column is numeric
  if (is.numeric(sitecovs[[col_name]])) {
    # Calculate the minimum and maximum values
    min_val <- min(sitecovs[[col_name]], na.rm = TRUE)
    max_val <- max(sitecovs[[col_name]], na.rm = TRUE)
    
    # Generate num_values values between the minimum and maximum
    new_values <- seq(from = min_val, to = max_val, length.out = 10)
  }
  new_data[[col_name]] <- new_values
}

new_data <- data.frame(new_data)
new_data











# Define the list of predictors you want to plot
predictors <- c("temp_log", "turb_log", "flow", "doy", "mean_liters_filtered_avg",
                "max_MAX_GRAD_D", "VB_AreaSqKm", "mean_SINUOSITY", "max_STRM_ORDER",
                "mean_StrmPow", "mean_MEANANNCMS", 
                "Time_Since_Last_Burn",  "Percent_burned", "Burned_Numerical", 
                "Ydam_density"   
)

predictor_labels <- c("Temperature", "Turbidity", "Flow", "DOY", "Liters Filtered", 
                      "Gradient (%)", "Valley Bottom Area", "Sinuosity", "Stream Order",
                      "Stream Power", "Mean Annual Discharge",
                      "Time_Since_Last_Burn (years)",  "Percent Burned", "Burned?",
                      "Beaver Dam Density"
)





# Create an empty list to store the plots
plot_list <- list()

# Predict the response variable using the model for all predictors
for (predictor in predictors) {
  
  #Fit model
  model <- pcount(formula = ~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg ~ 
                    [[predictor]], data = edna_matrix, se = T, K = 1165)
  
  #Create the new data
  if (is.numeric(sitecovs[[predictor]])) 
    # Calculate the minimum and maximum values
    min_val <- min(sitecovs[[predictor]], na.rm = TRUE)
    max_val <- max(sitecovs[[predictor]], na.rm = TRUE)
    
    # Generate num_values values between the minimum and maximum
    new_values <- seq(from = min_val, to = max_val, length.out = 62)
  
  
  new_data <- cbind(combined_dat[,c(9,16,21,22,23)], new_values)
  
  # Predict the response variable using the model
  z <- predict(model, type = 'det', newdata = new_data, appendData = TRUE)
  
  # Create the plot for the current predictor and store it in the list
  plot <- ggplot(z, aes(x = .data[[predictor]], y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
    geom_line(size = 1) +
    labs(x = predictor_labels[predictor == predictors], y = "eDNA concentration") +
    theme_bw()
  
  plot_list[[predictor]] <- plot
}

# Combine the plots into a single panel
combined_plot <- wrap_plots(plotlist = plot_list)

# Print the combined plot
combined_plot


















####Doing it manually because I can't figure out the loop


new_data <- list()


for (col_name in names(sitecovs)) {
  # Check if the column is numeric
  if (is.numeric(sitecovs[[col_name]])) {
    # Calculate the minimum and maximum values
    min_val <- min(sitecovs[[col_name]], na.rm = TRUE)
    max_val <- max(sitecovs[[col_name]], na.rm = TRUE)
    
    # Generate num_values values between the minimum and maximum
    new_values <- seq(from = min_val, to = max_val, length.out = 10)
  }
  new_data[[col_name]] <- new_values
}

new_data <- data.frame(new_data)

new_data <- new_data %>% 
  mutate(Burned_Numerical = c(0,0,0,0,0,1,1,1,1,1))

new_data






model <- pcount(formula = ~ 1 ~ 
                  max_MAX_GRAD_D
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p1 <- ggplot(z, aes(x = max_MAX_GRAD_D, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Gradient", y = "eDNA Concentration (Copies/L)") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  

p1



model <- pcount(formula = ~ 1 ~ 
                  VB_AreaSqKm
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p2 <- ggplot(z, aes(x = VB_AreaSqKm, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Valley Bottom Area", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely



model <- pcount(formula = ~ 1 ~ 
                  mean_SINUOSITY
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p3 <- ggplot(z, aes(x = mean_SINUOSITY, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Sinuosity", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely




model <- pcount(formula = ~ 1 ~ 
                  max_STRM_ORDER
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p4 <- ggplot(z, aes(x = max_STRM_ORDER, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Stream Order", y = "eDNA concentration") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely




model <- pcount(formula = ~ 1 ~ 
                  mean_StrmPow
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p5 <- ggplot(z, aes(x = mean_StrmPow, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Stream Power", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely





model <- pcount(formula = ~ 1 ~ 
                  mean_MEANANNCMS
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p6 <- ggplot(z, aes(x = mean_MEANANNCMS, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Mean Annual Discharge", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely




model <- pcount(formula = ~ 1 ~ 
                  Time_Since_Last_Burn
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p7 <- ggplot(z, aes(x = Time_Since_Last_Burn, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Time Since Last Burn", y = "eDNA concentration") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely






model <- pcount(formula = ~ 1 ~ 
                  Percent_burned
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p8<- ggplot(z, aes(x = Percent_burned, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Percent Burned", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely




model <- pcount(formula = ~ 1 ~ 
                  Burned_Numerical
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p9 <- ggplot(z, aes(x = Burned_Numerical, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Burned?", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely






model <- pcount(formula = ~ 1 ~ 
                  Ydam_density
                , data = edna_matrix, se = T, K = 1165)

z <- predict(model, type = 'state', newdata = new_data, appendData = TRUE)

p10 <- ggplot(z, aes(x = Ydam_density, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(x = "Beaver Dam Density", y = "eDNA concentration") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely








panel <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                   nrow = 4  
                   #labels = "AUTO"  # Automatically label the plots
)
panel


ggsave(plot= panel,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA_predictor_panel.jpeg",
       dpi = 1000, 
       height = 8,
       width = 6.5,
       units = "in")














combined_dat$mean_abundance



x_values <- seq(min(combined_dat$Ydam_density), max(combined_dat$Ydam_density), length.out = 100)
y_values <- (8000 * exp(-2.2 * x_values))+800

p11 <- ggplot(combined_dat, aes(x = Ydam_density, y = mean_abundance)) +
  geom_point(size = 2.5, alpha = 0.7) +
  #geom_smooth(formula = y ~ log(x), se = FALSE) + # method = "lm", # Logarithmic line of best fit
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "darkblue", lwd = 1.5, alpha = 0.75) +
  labs(x = "Beaver Dam Density", y = element_blank()) +
  theme_bw()#+
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely

p11





# Arrange the plots into a grid
grid <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3  # Set the number of columns to 3
)

# Arrange p10 and p11 in a separate row
row4 <- plot_grid(p10, p11, ncol = 2)

# Combine the two rows into a final grid
final_grid <- plot_grid(
  grid,
  row4,
  ncol = 1,  # Set the number of columns to 1 for a single column layout
  rel_heights = c(3, 1.5)  # Adjust the relative heights of the rows
)

# Print the final grid
print(final_grid)






ggsave(plot= final_grid,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA_predictor_panel2.jpeg",
       dpi = 1000, 
       height = 9,
       width = 6.5,
       units = "in")

















#Final attempt!!



####Doing it manually because I can't figure out the loop






#Format the basic new data

new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                      
                      max_MAX_GRAD_D = rep(seq(from = min(combined_dat$max_MAX_GRAD_D), to = max(combined_dat$max_MAX_GRAD_D), length.out = 10)),
                      VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                      mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                      max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                      mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                      mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                      Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                      Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                      Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                      Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)

newData




#Global model:



m.global <- pcount(
  #Occupancy factors
  ~ temp_log + turb_log + flow + doy + mean_liters_filtered_avg  ~ 
    #Fire
    Burned_Numerical + Percent_burned + Time_Since_Last_Burn +
    #Beaver
    Ydam_density + # Ydam_density*Surveyed+
    #Geomorphology
    max_STRM_ORDER + mean_SINUOSITY + VB_AreaSqKm + max_MAX_GRAD_D + # mean_WIDTH_M + mean_DEPTH_M +  
    #Hydrology
    mean_MEANANNCMS + mean_StrmPow, 
  data = edna_matrix, K = 1165, se = T)

summary(m.global)











#Predictions and Plots




new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                      turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                      flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                      doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                      mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                      
                      max_MAX_GRAD_D = rep(seq(from = min(combined_dat$max_MAX_GRAD_D), to = max(combined_dat$max_MAX_GRAD_D), length.out = 10)),
                      VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                      mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                      max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                      mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                      mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                      Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                      Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                      Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                      Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)

new_data 

z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p1 <- ggplot(z, aes(x = max_MAX_GRAD_D, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "A", x = "Gradient (%)", y = "eDNA Concentration (Copies/L)") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely.
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p1






new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(seq(from = min(combined_dat$VB_AreaSqKm), to = max(combined_dat$VB_AreaSqKm), length.out = 10)),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)



z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p2 <- ggplot(z, aes(x = VB_AreaSqKm, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "B", x = "Valley Bottom Area (km^2)", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold

p2







new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),

                        mean_SINUOSITY = rep(seq(from = min(combined_dat$mean_SINUOSITY), to = max(combined_dat$mean_SINUOSITY), length.out = 10)),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)



z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p3 <- ggplot(z, aes(x = mean_SINUOSITY, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "C", x = "Sinuosity", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p3







new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        
                        max_STRM_ORDER = rep(seq(from = min(combined_dat$max_STRM_ORDER), to = max(combined_dat$max_STRM_ORDER), length.out = 10)),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)


z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p4 <- ggplot(z, aes(x = max_STRM_ORDER, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "D", x = "Stream Order", y = "eDNA Concentration (Copies/L)") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p4







new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        
                        mean_StrmPow = rep(seq(from = min(combined_dat$mean_StrmPow), to = max(combined_dat$mean_StrmPow), length.out = 10)),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)



z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p5 <- ggplot(z, aes(x = mean_StrmPow, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "E", x = "Stream Power (N/s)", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p5










new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        
                        mean_MEANANNCMS = rep(seq(from = min(combined_dat$mean_MEANANNCMS), to = max(combined_dat$mean_MEANANNCMS), length.out = 10)),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)


z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p6 <- ggplot(z, aes(x = mean_MEANANNCMS, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "F", x = "Mean Annual Discharge (cms)", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p6








new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        
                        Time_Since_Last_Burn = rep(seq(from = min(combined_dat$Time_Since_Last_Burn), to = max(combined_dat$Time_Since_Last_Burn), length.out = 10)),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)



z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p7 <- ggplot(z, aes(x = Time_Since_Last_Burn, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "G", x = "Time Since Last Burn (Years)", y = "eDNA Concentration (Copies/L)") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold

  
p7









new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        
                        Percent_burned = rep(seq(from = min(combined_dat$Percent_burned), to = max(combined_dat$Percent_burned), length.out = 10)),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)



z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p8<- ggplot(z, aes(x = Percent_burned, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "H", x = "Percent Burned", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold


p8







new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
                        Ydam_density = rep(mean(combined_dat$Ydam_density, na.rm = T), length.out = 10)
)




z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p9 <- ggplot(z, aes(x = Burned_Numerical, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "I", x = "Burned?", y = element_blank()) +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold

p9



#p9 <- ggplot(z, aes(x = as.factor(Burned_Numerical), y = Predicted)) +
#  geom_point(size = 5) +
#  labs(x = "Burned?", y = element_blank()) +
#  scale_y_continuous(labels = function(x) x * 10) +
#  theme_bw()

#p9







new_data  <- data.frame(temp_log = rep(mean(combined_dat$temp_log, na.rm = T), length.out = 10),
                        turb_log = rep(mean(combined_dat$turb_log, na.rm = T), length.out = 10),
                        flow = rep(mean(combined_dat$flow, na.rm = T), length.out = 10),
                        doy = rep(mean(combined_dat$doy, na.rm = T), length.out = 10),
                        mean_liters_filtered_avg = rep(mean(combined_dat$mean_liters_filtered_avg, na.rm = T), length.out = 10),
                        
                        max_MAX_GRAD_D = rep(mean(combined_dat$max_MAX_GRAD_D, na.rm = T), length.out = 10), 
                        VB_AreaSqKm = rep(mean(combined_dat$VB_AreaSqKm, na.rm = T), length.out = 10),
                        mean_SINUOSITY = rep(mean(combined_dat$mean_SINUOSITY, na.rm = T), length.out = 10),
                        max_STRM_ORDER = rep(mean(combined_dat$max_STRM_ORDER, na.rm = T), length.out = 10),
                        mean_StrmPow = rep(mean(combined_dat$mean_StrmPow, na.rm = T), length.out = 10),
                        mean_MEANANNCMS = rep(mean(combined_dat$mean_MEANANNCMS, na.rm = T), length.out = 10),
                        Time_Since_Last_Burn = rep(mean(combined_dat$Time_Since_Last_Burn, na.rm = T), length.out = 10),
                        Percent_burned = rep(mean(combined_dat$Percent_burned, na.rm = T), length.out = 10),
                        Burned_Numerical = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        Ydam_density = rep(seq(from = min(combined_dat$Ydam_density), to = max(combined_dat$Ydam_density), length.out = 10))
                        
)


z <- predict(m.global, type = 'state', newdata = new_data, appendData = TRUE)

p10 <- ggplot(z, aes(x = Ydam_density, y = Predicted)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_line(size = 1) +
  labs(title = "J", x = "Beaver Dam Density (Dams/km^2", y = "eDNA Concentration (Copies/L)") +
  scale_y_continuous(labels = function(x) x * 10)+ #back transform the eDNA concentration so it matches the original numbers OR
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold

p10







panel <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                   nrow = 4  
                   #labels = "AUTO"  # Automatically label the plots
)
panel


ggsave(plot= panel,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA_predictor_panel_mean.jpeg",
       dpi = 1000, 
       height = 9,
       width = 6.5,
       units = "in")














combined_dat$mean_abundance



x_values <- seq(min(combined_dat$Ydam_density), max(combined_dat$Ydam_density), length.out = 100)
y_values <- (8000 * exp(-2.2 * x_values))+800

p11 <- ggplot(combined_dat, aes(x = Ydam_density, y = mean_abundance)) +
  geom_point(size = 2.5, alpha = 0.7) +
  #geom_smooth(formula = y ~ log(x), se = FALSE) + # method = "lm", # Logarithmic line of best fit
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "darkblue", lwd = 1.5, alpha = 0.75) +
  labs(title = "K", x = "Beaver Dam Density (dams/km^2)", y = element_blank()) +
  theme_bw()+
#theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())  #Remove the eDNA values completely
  theme(plot.title = element_text(face = "bold"))  # Make the title bold

p11





# Arrange the plots into a grid
grid <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3  # Set the number of columns to 3
)

# Arrange p10 and p11 in a separate row
row4 <- plot_grid(p10, p11, ncol = 2)

# Combine the two rows into a final grid
final_grid <- plot_grid(
  grid,
  row4,
  ncol = 1,  # Set the number of columns to 1 for a single column layout
  rel_heights = c(3, 1.5)  # Adjust the relative heights of the rows
)

# Print the final grid
print(final_grid)






ggsave(plot= final_grid,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA_predictor_panel2_mean.jpeg",
       dpi = 1000, 
       height = 14,
       width = 9,
       units = "in")















