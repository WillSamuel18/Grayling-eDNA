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
  select(sample_num, Site_Num, copies_per_L) %>% 
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

#View(replicate_dat)



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
print(distinct_covariate_dat)




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



write.csv(combined_dat, "2022 Summer eDNA/combined_dat.csv")






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

#y <- y %>% 
#  mutate(rep_1_copies_per_L = ifelse(rep_1_copies_per_L > 0, 1, 0),
#         rep_2_copies_per_L = ifelse(rep_2_copies_per_L > 0, 1, 0),
#         rep_3_copies_per_L = ifelse(rep_3_copies_per_L > 0, 1, 0),
#         rep_4_copies_per_L = ifelse(rep_4_copies_per_L > 0, 1, 0),
#         rep_5_copies_per_L = ifelse(rep_5_copies_per_L > 0, 1, 0)
#         )

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
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/occupancy_detectability.jpeg",
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





# Real data
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
(fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))
fm.mallard
(fm.mallard.nb <- pcount(~ date + I(date^2) ~ length + elev, mixture = "NB", mallardUMF, K=30))
fm.mallard.nb










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

sitecovs <- data.frame(combined_dat[,c(10, 14:22, 27, 67, 73:74, 77:107, 121:123)])  #May want to include the L filtered for each individual filters in the future



covar = rep(1, 62)
obscovs <- data.frame(covar)  


edna_matrix <- unmarkedFramePCount(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)



test <- pcount(formula = ~ 1 ~ 1, data = edna_matrix, se = T, K = 1165)
summary(test) #AIC 
#still need to back transform


global <- pcount(formula = ~ temp_log + turb_log + flow + doy ~ 1, data = edna_matrix, se = T, K = 1165)
summary(global) #AIC 

turb <- pcount(formula = ~ turb_log ~ 1, data = edna_matrix, se = T, K = 1165)
summary(turb) #AIC 


edna_matrix <- unmarkedFramePCount(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)








global <- occuRN(formula = ~ temp_log + turb_log + flow + doy ~ 1, data = edna_matrix, se = T)
summary(global) #AIC 322.0








m1  <- pcount(formula = ~temp_log ~ 1 , data = edna_matrix, se = TRUE)
m1
summary( m1 )
#AIC 318.914860679275 

m2  <- pcount(formula = ~turb_log ~ 1 , data = edna_matrix, se = TRUE)
summary( m2 )
#AIC 317.949065145275 

m3  <- pcount(formula = ~flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m3 )
#AIC 318.631011619057 

m4  <- pcount(formula = ~doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m4 )
#AIC 327.764892397531 

m5  <- pcount(formula = ~temp_log + turb_log ~ 1 , data = edna_matrix, se = TRUE)
summary( m5 )
#AIC 319.346743480782 

m6  <- pcount(formula = ~temp_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m6 )
#AIC 320.602494020965 

m7  <- pcount(formula = ~temp_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m7 )
#AIC 320.775809424329 

m8  <- pcount(formula = ~turb_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m8 )
#AIC 319.900496525724 

m9  <- pcount(formula = ~turb_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m9 )
#AIC 319.244151678218 

m10  <- pcount(formula = ~flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m10 )
#AIC 328.650639273185 

m11  <- pcount(formula = ~temp_log + turb_log + flow ~ 1 , data = edna_matrix, se = TRUE)
summary( m11 )
#AIC 322.16139765398 

m12  <- pcount(formula = ~temp_log + turb_log + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m12 )
#AIC 320.094118914255 

m13  <- pcount(formula = ~temp_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m13 )
#AIC 322.317051370482 

m14  <- pcount(formula = ~turb_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
summary( m14 )
#AIC 322.421279638573 

m15  <- pcount(formula = ~temp_log + turb_log + flow + doy ~ 1 , data = edna_matrix, se = TRUE)
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

#    Model      AIC                       Predictors
#14    m13 18428.39            temp_log + flow + doy
#8      m7 18432.54                   temp_log + doy
#15    m14 18475.06            turb_log + flow + doy
#1  global 18497.73                              All
#16    m15 18497.73 temp_log + turb_log + flow + doy
#13    m12 18499.12        temp_log + turb_log + doy
#10     m9 18516.15                   turb_log + doy
#12    m11 18517.50       temp_log + turb_log + flow
#11    m10 18532.80                       flow + doy
#5      m4 18533.23                              doy
#6      m5 18533.45              temp_log + turb_log
#2      m1 18584.50                         temp_log
#7      m6 18586.14                  temp_log + flow
#3      m2 18589.30                         turb_log
#9      m8 18590.88                  turb_log + flow
#4      m3 18604.16                             flow





logLik(m13)
logLik(m7)
logLik(m14)
logLik(global)
logLik(m12)
logLik(m9)
logLik(m11)
logLik(m10)
logLik(m4)
logLik(m5)
logLik(m1)
logLik(m6)
logLik(m2)
logLik(m8)
logLik(m3)











###Plot the relationships between detection and environmental covariates


min(combined_dat$turb_log) #0
max(combined_dat$turb_log) #61.03
min(combined_dat$temp_log) #2.88
max(combined_dat$temp_log) #15.1
min(combined_dat$flow) #0
max(combined_dat$flow) #1.136
min(combined_dat$doy) #192
max(combined_dat$doy) #226

#newData <- data.frame(turb_log = rep(seq(from = 0, to = 61, length.out = 10)),
#                      temp_log = rep(seq(from = 0, to = 2, length.out = 10)),
#                      flow = rep(seq(from = 0, to = 1.15, length.out = 10)),
#                      doy = rep(seq(from = 192, to = 226, length.out = 10)))



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





detection_plot <- turb_plot | temp_plot | flow_plot | doy_plot
detection_plot



ggsave(plot= detection_plot,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/N-mixture_detectability.jpeg",
       dpi = 1000, 
       height = 5,
       width = 20,
       units = "in")















m <- pcount(~ 1 ~ 1, data = edna_matrix, K = 1165, se = T)
summary(m)


m <- pcount(~ temp_log + turb_log + flow + doy ~ 1
            , data = edna_matrix, K = 1165, se = T)
summary(m)



fire <- pcount(~ turb_log ~ Burned_Numerical + Percent_burned + Time_Since_Last_Burn
            , data = edna_matrix, K = 1165, se = T)
summary(fire)


netmap <- pcount(~ turb_log ~ max_MEANANNCMS + VB_AreaSqKm + mean_GRADIENT, data = edna_matrix, K = 1165, se = T)
summary(netmap)


beaver <- pcount(~ turb_log ~ Alldam_density, data = edna_matrix, K = 1165, se = T)
summary(beaver)


beaver2 <- pcount(~ turb_log ~ Ydam_density, data = edna_matrix, K = 1165, se = T)
summary(beaver2)


all <- pcount(~ turb_log ~ Burned_Numerical +  max_MEANANNCMS + mean_GRADIENT + Alldam_density, data = edna_matrix, K = 1165, se = T)
summary(all)


#Max of 6 abundance replicates, and already using 1-4 detection covariates


AIC(m) #421.8
AIC(fire) #422.5
AIC(netmap) #422.5
AIC(beaver) #418.5
AIC(beaver2) #418
AIC(all) #424.4







reanf
























































r_dat <- read.csv("2022 Summer eDNA/r_dat.csv")
r_dat <- r_dat[,-1]


str(r_dat)


#Need to replace missing values

r_dat <- r_dat %>%
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

y <- data.frame(round(r_dat[,2:6], digits = 0))



y <- y %>% 
  mutate(rep_1_copies_per_L = round(rep_1_copies_per_L),
         rep_2_copies_per_L = round(rep_2_copies_per_L),
         rep_3_copies_per_L = round(rep_3_copies_per_L),
         rep_4_copies_per_L = round(rep_4_copies_per_L),
         rep_5_copies_per_L = round(rep_5_copies_per_L))



y <- y %>% 
  mutate(rep_1_copies_per_L = rep_1_copies_per_L/10,
         rep_2_copies_per_L = rep_2_copies_per_L/10,
         rep_3_copies_per_L = rep_3_copies_per_L/10,
         rep_4_copies_per_L = rep_4_copies_per_L/10,
         rep_5_copies_per_L = rep_5_copies_per_L/10)


sitecovs <- data.frame(r_dat[,c(7,8,12:20)])  #May want to include the L for each individual filters in the future

#cat <- sitecovs[,c(1:3,11)]
#sitecovs <- round(sitecovs[,4:10])

#sitecovs <- cbind(cat, sitecovs)

#Right now we we'll just make these nothing, but we might want to include filter specific covariates (e.g., L filtered)
#m <- 84
#n <- 1
#obscovs <- list(x2 = round(matrix(runif(m * n), m, n)))

covar = rep(1, 84)
obscovs <- data.frame(covar)  


edna_matrix <- unmarkedFramePCount(y=y, siteCovs=sitecovs)   #, obsCovs = obscovs
edna_matrix
summary(edna_matrix)



#?pcount()

m <- pcount(~ 1 ~ 1
            , data = edna_matrix, K = 1165, se = T)



m <- pcount(~ ph_log + temp_log + turb_log + flow + doy + drainage ~ 1
       , data = edna_matrix, K = 1165, se = T)

# + 
#date + drainage + weather + 


summary(m)

plogis(coef(m, type="det")) # Should be close to p
