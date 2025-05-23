################################################################################
################################ Will Samuel ###################################
####################### 2022 Grayling eDNA Modeling ############################
################################################################################


library(tidyverse)  #for data manipulation
library(MuMIn)      #For model dredging
library(ggcorrplot) #to make a correlation plot
library(rcompanion) #for pseudo R-squared
library(stringr)    #To split a column in 2 (used for efish time)
library(cowplot)


setwd("C:/Users/npwil/OneDrive/Desktop/School/Grad School/Thesis/Data and Analysis")


qpcr_dat <- read.csv(file = '2022 Summer eDNA/eDNA results (JV)/JVB1924-qpcr-tabulated-data.csv')

edna_dat <- read.csv(file = '2022 Summer eDNA/eDNA Data (2).csv')

effort_dat <- read.csv(file = '2022 Summer Fish Sampling/Data/W_Samuel_Fishing_Effort_Datasheet.csv') 

grayling_dat <- read.csv(file = '2022 Summer eDNA/W_Samuel_Grayling_Datasheet2+CPC.csv') 

#sonde_WS_dat <- read.csv(file = "2022 Summer eDNA/First, Second, and third Stint of intensive eDNA Sonde  AND broad sites data_export (1).csv")

#sonde_CPC_dat <- read.csv(file = "2022 Summer eDNA/CARIBOU_POKER_CREEK_EDNA_SAMPLING_FOR_WS_2022.csv")



# Data Manipulation -------------------------------------------------------




##Join the datasets by sample_ID
edna_dat <- edna_dat %>% rename(sample_ID = barcode_num)
qpcr_dat <- qpcr_dat %>% rename(sample_ID = SampleId)


edna_dat <- merge(x = edna_dat, y = qpcr_dat, by = "sample_ID", all.x = TRUE)
str(edna_dat)
head(edna_dat)


#Filter out top of reaches for now
edna_dat <- edna_dat %>% 
  #filter(!str_detect(sample_num, 'TR')) %>% 
  filter(transect == 'A') %>% 
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
  mutate("ph_log" = ifelse(ph_log > 1, ph_log, "NA")) %>%  #Remove faulty pH measurements
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
  filter(Reach_Type == "Control", Sampling_Method == "Angling") %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date) %>% 
  #group_by(Reach_Type) %>% 
  summarize("Abundance_angl" = sum(fish), 
            "Biomass_g_angl" = sum(Weight, na.rm = TRUE),
            "Length_total_mm_angl" = sum(Fork_Length, na.rm = TRUE))
            
grayling_sums2 <- grayling_dat %>% 
  filter(Reach_Type == "Control", Sampling_Method == "E-Fishing") %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date, Site_Num) %>% 
  #group_by(Reach_Type) %>%    
  summarize("Abundance_efish" = sum(fish), 
            "Biomass_g_efish" = sum(Weight, na.rm = TRUE),
            "Length_total_mm_efish" = sum(Fork_Length, na.rm = TRUE))
            

grayling_sums <- merge(x = grayling_sums , y = grayling_sums2, 
                       by = c("Date"), all.x = TRUE, all.y = TRUE)



#Calculate the combined sampling effort

#effort_dat <- effort_dat %>% 
#  filter(Reach_Type == "Control") %>% 
#  mutate("Effort" = (Angling_total+Efish_total)) #Need to come up with an equation to combine the sampling strategies


#Using Multigear Mean Standardization (MGMS) see Gibson-Reinemer et al. 2016 for more details

#Combine the effort and catch data
effort_dat <- effort_dat %>% 
  filter(Reach_Type == "Control")

effort_dat <- merge(x = effort_dat , y = grayling_sums, by = c("Date", "Site_Num"), all.x = TRUE)

#Calculate the efishing time in seconds instead of minutes:seconds
effort_dat[c('min', 'sec')] <- as.numeric(str_split_fixed(effort_dat$Time_on_Efish, ':', 2))

effort_dat <- effort_dat %>% 
  mutate("Time_on_Efish_sec" = ((min*60)+sec))


#Step 1: calculate CPUE for each survey method
st_effort_dat <- effort_dat %>% 
  mutate("Angling_CPUE_abun" = (Abundance_angl/((Angling_total*Num_anglers)*Area_m2))) %>% 
  mutate("Angling_CPUE_biom" = (Biomass_g_angl/((Angling_total*Num_anglers)*Area_m2))) %>% 
  mutate("Efish_CPUE_abun" = (Abundance_efish/(Time_on_Efish_sec*Area_m2))) %>% #Maybe needs to be Efish_total
  mutate("Efish_CPUE_biom" = (Biomass_g_efish/(Time_on_Efish_sec*Area_m2))) #Maybe needs to be Efish_total

#Step 2: calculate the mean total CPUE
mean.st.angling.abun <- mean(st_effort_dat$Angling_CPUE_abun, na.rm = TRUE)
mean.st.angling.biom <- mean(st_effort_dat$Angling_CPUE_biom, na.rm = TRUE)
mean.st.Efish.abun <- mean(st_effort_dat$Efish_CPUE_abun, na.rm = TRUE)
mean.st.Efish.biom <- mean(st_effort_dat$Efish_CPUE_biom, na.rm = TRUE)

#Step 3: use the mean total CPUE to create the standardized effort
st_effort_dat <- st_effort_dat %>% 
  mutate("Angling_CPUE_abun_STD" = (Angling_CPUE_abun/mean.st.angling.abun)) %>% 
  mutate("Angling_CPUE_biom_STD" = (Angling_CPUE_biom/mean.st.angling.biom)) %>% 
  mutate("Efish_CPUE_abun_STD" = (Efish_CPUE_abun/mean.st.Efish.abun)) %>% 
  mutate("Efish_CPUE_biom_STD" = (Efish_CPUE_biom/mean.st.Efish.biom)) 



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




#Summarize to include all replicates
#st_effort_dat_all <- merge(x = edna_dat, y = st_effort_dat, 
#                          by = c("Site_resample"), all.x = T, all.y = F)

#st_effort_dat_all <- st_effort_dat_all %>% 
#  mutate("Vel_ms" = mean(c(flow_1, flow_2, flow_3), na.rm = TRUE),
#            "water_temp" = temp_log, 
#            "pH" = ph_log, na.rm = TRUE, 
#            "SC" = sc_log, 
#            "HDO" = hdo_ml.L_log, 
#            "HDO_perc" = hdo_perc_sat_log, 
#            "Turb" = turb_log) %>% 
#  mutate("Vel_ms" = ifelse(is.na(Vel_ms), mean(Vel_ms, na.rm = T), Vel_ms), 
#       "water_temp" = ifelse(is.na(water_temp), mean(water_temp, na.rm = T), water_temp),
#       #"pH" = ifelse(is.na(pH), mean(pH, na.rm = T), pH),
#       "SC" = ifelse(is.na(SC), mean(SC, na.rm = T), SC),
#       "HDO" = ifelse(is.na(HDO), mean(HDO, na.rm = T), HDO),
#       "HDO_perc" = ifelse(is.na(HDO_perc), mean(HDO_perc, na.rm = T), HDO_perc),
#       "Turb" = ifelse(is.na(Turb), mean(Turb, na.rm = T), Turb)) %>% 
#  mutate("log_copies_per_L" = copies_per_L,
#         "log_MGMS_CPUE_abun" = log(MGMS_CPUE_abun),
#         "log_MGMS_CPUE_biom" = log(MGMS_CPUE_biom))
         
#st_effort_dat_all <- st_effort_dat_all %>% 
#  mutate("copies_per_L_log" = ifelse(copies_per_L>0,log(copies_per_L),0),
#         "MGMS_CPUE_abun_log" = ifelse(MGMS_CPUE_abun>0,log(MGMS_CPUE_abun),0),
#         "MGMS_CPUE_biom_log" = ifelse(MGMS_CPUE_biom>0,log(MGMS_CPUE_biom),0))



str(st_effort_dat_all) #93 obs


#st_effort_dat <- st_effort_dat[-c(1,23),] #remove outliers: Belle creek on 6/13/22, and CPC reach 33 on  8/10/22
#st_effort_dat <- st_effort_dat[-c(4,17),] #remove Angel creek on 6/22/22 and 7/29/22, when there was no water

#View(st_effort_dat) #24 obs to use for modeling




#Average replicates
edna_dat <- edna_dat %>% rename(Date = date)

st_effort_dat <- merge(x = st_effort_dat, y = edna_dat, 
                       by = c("Date", "Site_Num"), all.x = F, all.y = F)

st_effort_dat <- st_effort_dat %>% 
  filter(transect == 'A') %>% 
  group_by(Date, Site_Num) %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  mutate(ph_log = ifelse(ph_log > 1, ph_log, "NA")) %>%  #Remove faulty pH measurements
  summarize("copies_per_L" = mean(copies_per_L, na.rm = T),
            "mean_L_filtered" = mean(liters_filtered_avg, na.rm = T),
            "Vel_ms" = mean(c(flow_1, flow_2, flow_3), na.rm = TRUE), 
            "water_temp" = mean(temp_log, na.rm = TRUE), 
            "pH" = mean(ph_log, na.rm = TRUE), 
            "SC" = mean(sc_log, na.rm = TRUE), 
            "HDO" = mean(hdo_ml.L_log, na.rm = TRUE), 
            "HDO_perc" = mean(hdo_perc_sat_log, na.rm = TRUE), 
            "Turb" = mean(turb_log, na.rm = TRUE),
            "MGMS_CPUE_abun" = mean(MGMS_CPUE_abun, na.rm = T),
            "MGMS_CPUE_biom" = mean(MGMS_CPUE_biom, na.rm = T)) %>% 
  mutate("Vel_ms" = ifelse(is.na(Vel_ms), mean(Vel_ms, na.rm = T), Vel_ms), 
            "mean_L_filtered" = ifelse(is.na(mean_L_filtered), mean(mean_L_filtered, na.rm = T), mean_L_filtered),
            "water_temp" = ifelse(is.na(water_temp), mean(water_temp, na.rm = T), water_temp),
            "pH" = ifelse(is.na(pH), mean(pH, na.rm = T), pH),
            "SC" = ifelse(is.na(SC), mean(SC, na.rm = T), SC),
            "HDO" = ifelse(is.na(HDO), mean(HDO, na.rm = T), HDO),
            "HDO_perc" = ifelse(is.na(HDO_perc), mean(HDO_perc, na.rm = T), HDO_perc),
            "Turb" = ifelse(is.na(Turb), mean(Turb, na.rm = T), Turb)) %>% 
  mutate("copies_per_L_log" = log(copies_per_L),
         "MGMS_CPUE_abun_log" = log(MGMS_CPUE_abun),
         "MGMS_CPUE_biom_log" = log(MGMS_CPUE_biom))


st_effort_dat <- st_effort_dat %>% 
  mutate("copies_per_L_log" = ifelse(copies_per_L>0,log(copies_per_L),0),
         "MGMS_CPUE_abun_log" = ifelse(MGMS_CPUE_abun>0,log(MGMS_CPUE_abun),0),
         "MGMS_CPUE_biom_log" = ifelse(MGMS_CPUE_biom>0,log(MGMS_CPUE_biom),0))




str(st_effort_dat) #28 obs

st_effort_dat <- st_effort_dat[-c(5,18),] #remove Angel creek on 6/22/22 and 7/29/22, when there was no water
#st_effort_dat <- st_effort_dat[-c(1,19,23),] #remove outliers: Belle creek on 6/13/22, and CPC reach 33 on 8/10/22


#View(st_effort_dat) #24 obs to use for modeling





# Exploratory Analysis ----------------------------------------------------




###Should I log transform this or not?????????????

plot(log(MGMS_CPUE_abun) ~ log(copies_per_L), data = st_effort_dat)
text(log(MGMS_CPUE_abun) ~ log(copies_per_L), labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 3)
#text(log(copies_per_L)~MGMS_CPUE_abun, labels=Date,data=st_effort_dat, cex=0.9, font=2, pos = 3)


#st_effort_dat <- st_effort_dat %>% subset(!Date == "6/9/2022")

plot(MGMS_CPUE_abun_log ~ copies_per_L_log, data = st_effort_dat)
text(MGMS_CPUE_abun_log ~ copies_per_L_log, labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 3)
text(MGMS_CPUE_abun_log ~ copies_per_L_log, labels=Date,data=st_effort_dat, cex=0.9, font=2, pos = 4)

summary(lm(copies_per_L_log ~ MGMS_CPUE_abun_log+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))
plot(copies_per_L_log ~ MGMS_CPUE_abun_log, data = st_effort_dat)


#summary(lm(copies_per_L_log ~ MGMS_CPUE_biom_log+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))





plot(log(copies_per_L) ~ MGMS_CPUE_biom, data = st_effort_dat)
text(log(copies_per_L) ~ MGMS_CPUE_biom, labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 3)
#text(log(copies_per_L)~MGMS_CPUE_biom, labels=Date,data=st_effort_dat, cex=0.9, font=2, pos = 3)
summary(lm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat_all))

hist(grayling_dat$Fork_Length) #make this a multi panel plot for each site/sampling event (Date). 


plot(log(MGMS_CPUE_abun) ~ log(MGMS_CPUE_biom), data = st_effort_dat)
cor(st_effort_dat$MGMS_CPUE_abun,st_effort_dat$MGMS_CPUE_biom)




st_effort_dat$cop

ggplot(st_effort_dat, aes(copies_per_L))+
  geom_histogram(fill = "lightblue", col = "black")+
  theme_cowplot()






# VIF ---------------------------------------------------------------------

plot(hdo_ml.L_log ~ hdo_perc_sat_log, data = edna_dat)
summary(lm(hdo_ml.L_log ~ hdo_perc_sat_log, data = edna_dat))



forVIF <- st_effort_dat %>% 
  select("Vel_ms", "water_temp", "pH", "SC", "HDO", "HDO_perc", "Turb")  %>% 
  mutate(Vel_ms = ifelse(is.na(Vel_ms), mean(Vel_ms, na.rm = T), Vel_ms)) %>% 
  mutate(water_temp = ifelse(is.na(water_temp), mean(water_temp, na.rm = T), water_temp)) %>%
  mutate(pH = ifelse(is.na(pH), mean(pH, na.rm = T), pH)) %>% 
  mutate(SC = ifelse(is.na(SC), mean(SC, na.rm = T), SC)) %>% 
  mutate(HDO = ifelse(is.na(HDO), mean(HDO, na.rm = T), HDO)) %>% 
  mutate(HDO_perc = ifelse(is.na(HDO_perc), mean(HDO_perc, na.rm = T), HDO_perc)) %>% 
  mutate(Turb = ifelse(is.na(Turb), mean(Turb, na.rm = T), Turb)) 
  


corr <- round(cor(forVIF), 1)
corr
ggcorrplot(corr, type = "lower", lab = TRUE, method = "circle") 
#HDO and HDO_perc
#HDO and Turb
#water temp and HDO
#Water temp and Turb 


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


fromVIF <- vif_func(forVIF)
#So, remove HDO or HDO_percent. Duh... I think I want to leave the others, and check for over prediction later on


# Models ------------------------------------------------------------------


### MAKE SURE TO REMOVE ANGLE CREEK SAMPLS ON 6/24/22 AND 7/29/22


#st_effort_dat <- st_effort_dat[-c(5),]

#qpcr_dat <- qpcr_dat %>% rename(Date = date)

#qpcr_sums <- qpcr_dat %>% 
#  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
#  group_by(Date) %>% 
#  summarize("copies_per_L" = mean(copies_per_L)) 

#st_effort_dat <- merge(x = st_effort_dat, y = qpcr_sums, by = "Date", all.x = TRUE)





date <- mdy(st_effort_dat$Date)
doy <- data.frame(yday(date))
doy


st_effort_dat2 <- c(st_effort_dat, doy)
str(st_effort_dat2)





m.global.abun <- glm(copies_per_L ~ MGMS_CPUE_abun+Vel_ms+water_temp+Turb+yday.date.+mean_L_filtered, data = st_effort_dat2)
options(scipen = 1000)
summary(m.global.abun)

nagelkerke(m.global.abun)
accuracy(list(m.global.abun))
#Pseudo.R.squared                        NA -> mean
#McFadden                                 0.126
#Cox and Snell (ML)                       0.911
#Nagelkerke (Cragg and Uhler)             0.911
#Efron.r.squared                          0.269


(0.911+0.911+0.269+0.126)/4 #R2 = 0.554





#switch eDNA and fish abundance
m.global.abun <- glm(MGMS_CPUE_abun ~ copies_per_L+Vel_ms+water_temp+Turb+yday.date.+mean_L_filtered, data = st_effort_dat2)
options(scipen = 1000)
summary(m.global.abun)

nagelkerke(m.global.abun)
accuracy(list(m.global.abun))
#Pseudo.R.squared                        NA -> mean
#McFadden                                 0.167
#Cox and Snell (ML)                       0.633
#Nagelkerke (Cragg and Uhler)             0.635
#Efron.r.squared                          0.346


(0.633+0.635+0.346+0.167)/4 #R2 = 0.445




options(na.action = "na.fail")
global_dredge_abun <- dredge(m.global.abun, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge_abun

#Dredge prefers ~MGMS_CPUE_abun+Vel_ms




m.global.biom <- lm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
summary(m.global.biom)

nagelkerke(m.global.biom)
accuracy(list(m.global.biom))
#Pseudo.R.squared                   With NA     NA -> mean
#McFadden                            0.139        0.139
#Cox and Snell (ML)                  0.933        0.933
#Nagelkerke (Cragg and Uhler)        0.933        0.933
#Efron.r.squared                     0.326        0.326


options(na.action = "na.fail")
global_dredge_biom <- dredge(m.global.biom, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge
#Dredge prefers ~HDO+pH






#Same models, but using all the datapoints instead of the averaged data
m.global.abun.all <- glm(copies_per_L ~ MGMS_CPUE_abun+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat_all)
summary(m.global.abun.all)

nagelkerke(m.global.abun.all)
accuracy(list(m.global.abun.all))
#Pseudo.R.squared                   With NA     NA -> mean
#McFadden                            0.102        0.01
#Cox and Snell (ML)                  0.852        0.16
#Nagelkerke (Cragg and Uhler)        0.852        0.16
#Efron.r.squared                     0.162        0.16

options(na.action = "na.fail")
global_dredge_abun_all <- dredge(m.global.abun.all, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge_abun

#Dredge prefers ~MGMS_CPUE_abun+Vel_ms




m.global.biom.all <- glm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat_all)
summary(m.global.biom.all)

nagelkerke(m.global.biom.all)
accuracy(list(m.global.biom.all))
#Pseudo.R.squared                   With NA     NA -> mean
#McFadden                            0.102        0.01
#Cox and Snell (ML)                  0.853        0.167
#Nagelkerke (Cragg and Uhler)        0.853        0.167
#Efron.r.squared                     0.169        0.167


options(na.action = "na.fail")
global_dredge_biom_all <- dredge(m.global.biom.all, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge
#Dredge prefers ~HDO+pH




#Compare alternative models and pick the best one to predict, 
#ORRRRR Predict with all of them and average the results... 




# Model Evaluation --------------------------------------------------------

#Testing the global abundance model
global_pred_abun <- vector()

for(i in 1:24){
  validate.dat <- st_effort_dat[i,] #Single out one value
  training.dat <- st_effort_dat[-i,] #Use the training data using the data -i
  #m.global <- glm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
  global_pred_abun[i] <- predict(m.global.abun, newdata = validate.dat, type = "response")
  
}


obs <- as.numeric(st_effort_dat$copies_per_L)
obs
comp_abun <- vector()
comp_abun$global_pred_abun <- as.numeric(global_pred_abun)


comp_abun <- data.frame(obs, global_pred_abun)


comp_abun <- comp_abun %>%
  mutate(global_pred_abun = ifelse(is.na(global_pred_abun), NA, global_pred_abun)) %>% 
  mutate("diff" = obs-global_pred_abun) %>% 
  mutate("diff_abs" = ifelse(diff < 0, -diff, diff)) 
  
sum(comp_abun$diff_abs, na.rm=T)#Difference between the observed and predicted values 15648.66



#Testing the global abundance model
global_pred_biom <- vector()

for(i in 1:93){
  validate.dat <- st_effort_dat[i,] #Single out one value
  training.dat <- st_effort_dat[-i,] #Use the training data using the data -i
  #m.global <- glm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
  global_pred_biom[i] <- predict(m.global.biom, newdata = validate.dat, type = "response")
  
}


comp_biom$global_pred_biom <- as.numeric(global_pred_biom)

comp_biom <- vector()
comp_biom <- data.frame(obs, global_pred_biom)


comp_biom <- comp_biom %>%
  mutate("diff" = obs-global_pred_abun) %>% 
  mutate("diff_abs" = ifelse(diff < 0, -diff, diff))

sum(comp_biom$diff_abs, na.rm=T)#Difference between the observed and predicted values 76473.62


# Predicting to other sites -----------------------------------------------







# Caribou Poker Creek -----------------------------------------------------


plot(copies_per_L ~ Site_Num, data = st_effort_dat)

CPCfigure.dat <- st_effort_dat %>% 
  filter(Site_Num > 6) %>% 
  select(Site_Num, copies_per_L, MGMS_CPUE_abun) %>% 
  group_by(Site_Num) %>% 
  mutate(copies_per_L = mean(copies_per_L, na.rm = T),
         MGMS_CPUE_abun = mean(MGMS_CPUE_abun, na.rm =T))


ggplot(CPCfigure.dat, aes(Site_Num, copies_per_L))+
  geom_smooth(method = "lm")+
  geom_point(size = 3)+
  labs(x = "Site", y = "eDNA Concentration (Copies/L)")+
  theme_cowplot()






#Mean Annual eDNA per reach
st_effort_dat %>%
  filter(Site_Num > 6) %>%
  filter(!Date == c("6/28/2022")) %>%
  group_by(Site_Num) %>% 
  summarize(mean(copies_per_L))

#       Site_Num `mean(copies_per_L)`
#1        7                3448.
#2        8                2205.
#3        9                1087.
#4       10                 604.


(604+1087+2205+3448)/4
#Mean eDNA conc. is 1836

dat <- c(604, 1087, 2205, 3448)
sd(dat) #SD = 1266.693

#Mean flow in CPC
st_effort_dat %>% 
  filter(Site_Num >6) %>% 
  summarize(Vel_ms = mean(Vel_ms))


#Length of CPC
3528+4723 # = 8251
#8.2km 

#So, over 8251 m at 0.34 m/s, the DNA could make it through the whole stream in...
8251/0.34 # = 24267.65 seconds
24267/60 # = 404.45 minutes
404.45/60 # ~6.74 HOURS!!!




#Fish Data

dat2 <- c(0.812, 1.36, 1.31, 1.36, 5.43, 0.442, 0.987, 0.89)
mean(dat2) #mean = 1.574
sd(dat2) #SD = 1.59




cpc_dat <- read.csv(file = '2022 Summer eDNA/2022_CPCRW_Reach_Fish_Data-WS.csv')

str(cpc_dat)


#Summarize fish data
cpc_dat$fish <- rep(1, times = 349)
#grayling_dat$Fork_Length <- as.numeric(grayling_dat$Fork_Length, na.rm = TRUE)

cpc_sums <- cpc_dat %>% 
  filter(method == "Angling") %>% 
  group_by(Date, Reach_ID) %>% 
  summarize("Abundance_angl" = sum(fish), 
            "Biomass_g_angl" = sum(weight..g., na.rm = TRUE),
            "Length_total_mm_angl" = sum(length..mm., na.rm = TRUE),
            "Angling_Effort" = mean(Angling.minutes*Number.of.Anglers))


cpc_dat[c('min', 'sec')] <- as.numeric(str_split_fixed(cpc_dat$Electrofishing.minutes, ':', 2))

cpc_dat <- cpc_dat %>% 
  mutate("Time_on_Efish_sec" = ((min*60)+sec))

cpc_sums2 <- cpc_dat %>% 
  filter(method == "Electrofishing") %>% 
  group_by(Date, Reach_ID) %>% 
  summarize("Abundance_efish" = sum(fish), 
            "Biomass_g_efish" = sum(weight..g., na.rm = TRUE),
            "Length_total_mm_efish" = sum(length..mm., na.rm = TRUE),
            "Efish_Effort" = mean(Time_on_Efish_sec))




cpc_sums <- merge(x = cpc_sums , y = cpc_sums2, 
                       by = c("Date","Reach_ID"), all.x = TRUE, all.y = TRUE)



#Calculate the combined sampling effort
#Using Multigear Mean Standardization (MGMS) see Gibson-Reinemer et al. 2016 for more details

#Step 1: calculate CPUE for each survey method
cpc_sums <- cpc_sums %>% 
  mutate("Angling_CPUE_abun" = (Abundance_angl/Angling_Effort), 
        "Angling_CPUE_biom" = (Biomass_g_angl/Angling_Effort), 
        "Efish_CPUE_abun" = (Abundance_efish/Efish_Effort),
        "Efish_CPUE_biom" = (Biomass_g_efish/Efish_Effort)) 

#Step 2: calculate the mean total CPUE
mean.st.angling.abun <- mean(cpc_sums$Angling_CPUE_abun, na.rm = TRUE)
mean.st.angling.biom <- mean(cpc_sums$Angling_CPUE_biom, na.rm = TRUE)
mean.st.Efish.abun <- mean(cpc_sums$Efish_CPUE_abun, na.rm = TRUE)
mean.st.Efish.biom <- mean(cpc_sums$Efish_CPUE_biom, na.rm = TRUE)

#Step 3: use the mean total CPUE to create the standardized effort
cpc_sums <- cpc_sums %>% 
  mutate("Angling_CPUE_abun_STD" = (Angling_CPUE_abun/mean.st.angling.abun)) %>% 
  mutate("Angling_CPUE_biom_STD" = (Angling_CPUE_biom/mean.st.angling.biom)) %>% 
  mutate("Efish_CPUE_abun_STD" = (Efish_CPUE_abun/mean.st.Efish.abun)) %>% 
  mutate("Efish_CPUE_biom_STD" = (Efish_CPUE_biom/mean.st.Efish.biom)) 



#Step 5: assign this standardized effort back to the dataset, calculate the MGMS_effort
cpc_sums <- cpc_sums %>% 
  mutate("Angling_CPUE_abun_STD" = ifelse(is.na(Angling_CPUE_abun_STD), 0, Angling_CPUE_abun_STD)) %>% 
  mutate("Angling_CPUE_biom_STD" = ifelse(is.na(Angling_CPUE_biom_STD), 0, Angling_CPUE_biom_STD)) %>% 
  mutate("Efish_CPUE_abun_STD" = ifelse(is.na(Efish_CPUE_abun_STD), 0, Efish_CPUE_abun_STD)) %>% 
  mutate("Efish_CPUE_biom_STD" = ifelse(is.na(Efish_CPUE_biom_STD), 0, Efish_CPUE_biom_STD)) %>% 
  
  mutate("MGMS_CPUE_abun" = (Efish_CPUE_abun_STD+Angling_CPUE_abun_STD)) %>% 
  mutate("MGMS_CPUE_biom" = (Efish_CPUE_biom_STD+Angling_CPUE_biom_STD)) 



cpc_sums



cpc_sums %>%
  filter(!Date == c("6/28/2022")) %>%
  group_by(Reach_ID) %>% 
  summarize(MGMS_CPUE_abun = sum(MGMS_CPUE_abun),
            MGMS_CPUE_biom = sum(MGMS_CPUE_biom))




ggplot(cpc_sums, aes(Reach_ID, MGMS_CPUE_abun)) + 
  geom_point() +
  theme_cowplot()


summary(cpc_sums$MGMS_CPUE_abun)


summary_info_fish <- cpc_sums %>% 
  select(Reach_ID, MGMS_CPUE_abun) %>% 
  group_by(Reach_ID) %>% 
  summarize(MGMS_CPUE_abun = mean(MGMS_CPUE_abun))


ggplot(summary_info_fish, aes(Reach_ID, MGMS_CPUE_abun)) + 
  geom_point() +
  theme_cowplot()


summary(summary_info_fish$MGMS_CPUE_abun)

#Fish CPUE, averaged across the year for each reach


# Plotting ----------------------------------------------------------------




summary(lm(MGMS_CPUE_abun ~ copies_per_L+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))
#Why wont this work???
st_effort_dat


plot_dat <- st_effort_dat %>% 
  select(MGMS_CPUE_abun_log, copies_per_L_log, Date) %>% 
  group_by(Date) %>% 
  summarize(MGMS_CPUE_abun_log = MGMS_CPUE_abun_log, 
            copies_per_L_log = mean(copies_per_L_log, na.rm =T))






p <- ggplot(plot_dat, aes(x=MGMS_CPUE_abun_log, y=copies_per_L_log))+
  geom_smooth(method = lm, alpha = 0.4, linewidth = 1.5)+
  #stat_smooth(method = lm, formula = y ~ ifelse(x>0,log(x),))+
  geom_point()+
  labs(x = "log(CPUE Fish Abundance)", y = "log(eDNA Concentration) (Copies/L)")+
  #scale_x_continuous(trans = "exp")+
  theme_cowplot()


p


ggsave(plot= p,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA model FINAL.jpeg",
       dpi = 1000, 
       height = 4,
       width = 4,
       units = "in")









ggplot(data=comp_abun, aes(x=obs, y = global_pred_biom,))+
  geom_point()+
  theme_cowplot()






#Abundance Model
summary(lm(MGMS_CPUE_abun ~ copies_per_L+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))


p1 <- ggplot(st_effort_dat, aes(x=copies_per_L, y=MGMS_CPUE_abun))+
  geom_smooth(method = lm, alpha = 0.2, linewidth = 1.5)+
  geom_point()+
  labs(x = "eDNA Concentration (Copies/L)", y = "CPUE Fish Abundance")+
  #geom_text(aes(label = label), size = 3, hjust = 0, vjust = 0)+
  #geom_text(aes(label = label2), size = 3, hjust = 0, vjust = 0)+
  theme_cowplot()
p1

ggsave(plot= p1,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA model 1.jpeg",
       dpi = 1000, 
       height = 4,
       width = 4,
       units = "in")


summary(lm(MGMS_CPUE_abun_log ~ copies_per_L_log+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))



p2 <- ggplot(st_effort_dat, aes(x=copies_per_L_log, y= MGMS_CPUE_abun_log))+
  geom_smooth(method = lm, alpha = 0.2, linewidth = 1.5)+
  #stat_smooth(method = lm, formula = y ~ ifelse(x>0,log(x),))+
  geom_point()+
  labs(x = "log(eDNA Concentration) (Copies/L)", y = "log(CPUE Fish Abundnace)")+
  #geom_text(aes(label = label), size = 3, hjust = 0, vjust = 0)+
  #geom_text(aes(label = label2), size = 3, hjust = 0, vjust = 0)+
  theme_cowplot()

p2

ggsave(plot= p2,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA model 2.jpeg",
       dpi = 1000, 
       height = 4,
       width = 4,
       units = "in")


plot(copies_per_L_log~MGMS_CPUE_abun_log,data = st_effort_dat)
text(copies_per_L_log ~ MGMS_CPUE_abun_log, labels=Site_Num,data=st_effort_dat, cex=0.9, font=2, pos = 3)

       
       
       
summary(lm(copies_per_L_log ~ MGMS_CPUE_abun_log+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))


p3 <- ggplot(st_effort_dat, aes(x=log(MGMS_CPUE_abun), y=log(copies_per_L)))+
  geom_smooth(method = lm, alpha = 0.6, linewidth = 1.5)+
  geom_point()+
  labs(x = "CPUE Fish Abundance", y = "eDNA Concentration (Copies/L)")+
  #geom_text(aes(label = label), size = 3, hjust = 0, vjust = 0)+
  #geom_text(aes(label = label2), size = 3, hjust = 0, vjust = 0)+
  theme_cowplot()
p3

ggsave(plot= p3,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA model 3.jpeg",
       dpi = 1000, 
       height = 4,
       width = 4,
       units = "in")


summary(lm(copies_per_L_log ~ MGMS_CPUE_abun_log+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat))


lable <- c("Multiple R^2 = 0.53")
lable2 <- c("Adjusted R^2 = 0.33")



p4 <- ggplot(st_effort_dat, aes(x=MGMS_CPUE_abun_log, y=copies_per_L_log))+
  geom_smooth(method = lm, alpha = 0.6, linewidth = 1.5)+
  #stat_smooth(method = lm, formula = y ~ ifelse(x>0,log(x),))+
  geom_point()+
  labs(x = "log(CPUE Fish Abundance)", y = "log(eDNA Concentration) (Copies/L)")+
  #scale_x_continuous(trans = "exp")+
  theme_cowplot()


p4

ggsave(plot= p4,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/eDNA model 4.jpeg",
       dpi = 1000, 
       height = 4,
       width = 4,
       units = "in")

       
       
       

summary(edna_dat)


intensive <- edna_dat %>%
  filter(transect %in% c("A", "B", "C")) %>% 
  mutate("copies_per_L_intensive" = copies_per_L) %>% 
  mutate("copies_per_L_broad" = NA)

  

broad <- edna_dat %>%
  filter(!transect %in% c("A", "B", "C")) %>% 
  mutate("copies_per_L_broad" = copies_per_L) %>% 
  mutate("copies_per_L_intensive" = NA)



summary(intensive$copies_per_L)
#Min.   1st Qu. Median  Mean   3rd Qu.  Max. 
#0.00    0.00   92.33  587.66  704.59 5303.47 
summary(broad$copies_per_L)
#Min.   1st Qu.Median  Mean   3rd Qu.   Max. 
#0.0     0.0   332.7  1089.7  1346.7 11644.3        



curves <- rbind(intensive, broad)

curves$copies_per_L_broad <- as.numeric(curves$copies_per_L_broad)
curves$copies_per_L_intensive <- as.numeric(curves$copies_per_L_intensive)


#the distribution of the intensive vs broad eDNA data
###RUNNING INTO ISSUES WITH THIS PLOT...?

ggplot(curves, aes(x=copies_per_L))+
  geom_histogram(aes(x=copies_per_L_intensive))+
  geom_histogram(aes(x=copies_per_L_broad, color = "red", fill = "darkred"), alpha = 0.5)+
  #geom_density(aes(x=copies_per_L_intensive), lwd = 2)+
  #geom_density(aes(x=copies_per_L_broad, color = "red"), lwd = 2)+
  scale_x_continuous(limits = c(0,12000))+
  scale_y_continuous(limits = c(0, 50)) +
  theme_cowplot()
  


#Without 0's in the data
ggplot(curves, aes(x=copies_per_L))+
  geom_histogram(aes(x=copies_per_L_intensive))+
  geom_histogram(aes(x=copies_per_L_broad, color = "red", fill = "darkred"), alpha = 0.5)+
  #geom_density(aes(x=copies_per_L_intensive), lwd = 2)+
  #geom_density(aes(x=copies_per_L_broad, color = "red"), lwd = 2)+
  theme_cowplot()+
  scale_x_continuous(limits = c(0,12000))+
  scale_y_continuous(limits = c(0,50))

curves <- curves %>%
  #filter(copies_per_L_intensive > 0, 
  #      copies_per_L_broad > 0)




#Example Regression plot
CPUE <- c(1,2,3,5,6,7,8,9,10,11,11,12,14,16,18)
eDNA <- c(60,100,115,180,600,380,600,1100,900,1000, 1350, 1500,1750, 1900, 2100)
  
#rnorm(n = 50, st)

example <- data.frame(cbind(CPUE, eDNA))
example


example_plot <- ggplot(data = example, aes(x=CPUE, y =eDNA))+
  geom_smooth(method = "lm")+
  geom_point(size = 2)+
  theme_cowplot()+
  labs(x = "Catch Per Unit Effort", y= "eDNA Concentration")
example_plot



ggsave(plot= example_plot,
       filename = "2022 Summer eDNA/Grayling-eDNA R/Figures/example_plot.jpeg",
       dpi = 1000, 
       height = 3.5,
       width = 4,
       units = "in")



