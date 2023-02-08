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
  mutate("Angling_CPUE_abun" = (Abundance_angl/(Angling_total*Num_anglers))) %>% 
  mutate("Angling_CPUE_biom" = (Biomass_g_angl/(Angling_total*Num_anglers))) %>% 
  mutate("Efish_CPUE_abun" = (Abundance_efish/Time_on_Efish_sec)) %>% #Maybe needs to be Efish_total
  mutate("Efish_CPUE_biom" = (Biomass_g_efish/Time_on_Efish_sec)) #Maybe needs to be Efish_total

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
st_effort_dat_all <- merge(x = edna_dat, y = st_effort_dat, 
                           by = c("Date", "Site_Num"), all.x = T, all.y = F)
st_effort_dat_all <- st_effort_dat_all %>% 
  mutate("Vel_ms" = mean(c(flow_1, flow_2, flow_3), na.rm = TRUE),
            "water_temp" = temp_log, 
            "pH" = ph_log, na.rm = TRUE, 
            "SC" = sc_log, 
            "HDO" = hdo_ml.L_log, 
            "HDO_perc" = hdo_perc_sat_log, 
            "Turb" = turb_log) %>% 
  mutate("Vel_ms" = ifelse(is.na(Vel_ms), mean(Vel_ms, na.rm = T), Vel_ms), 
       "water_temp" = ifelse(is.na(water_temp), mean(water_temp, na.rm = T), water_temp),
       "pH" = ifelse(is.na(pH), mean(pH, na.rm = T), pH),
       "SC" = ifelse(is.na(SC), mean(SC, na.rm = T), SC),
       "HDO" = ifelse(is.na(HDO), mean(HDO, na.rm = T), HDO),
       "HDO_perc" = ifelse(is.na(HDO_perc), mean(HDO_perc, na.rm = T), HDO_perc),
       "Turb" = ifelse(is.na(Turb), mean(Turb, na.rm = T), Turb)) 
         
str(st_effort_dat_all) #93 obs


#st_effort_dat <- st_effort_dat[-c(1,23),] #remove outliers: Belle creek on 6/13/22, and CPC reach 33 on  8/10/22
#st_effort_dat <- st_effort_dat[-c(4,17),] #remove Angel creek on 6/22/22 and 7/29/22, when there was no water

#View(st_effort_dat) #24 obs to use for modeling





#Average replicates
edna_dat <- edna_dat %>% rename(Date = date)

edna_sums <- edna_dat %>% 
  filter(transect == 'A') %>% #You can group by control reach for better data clarity, but this reduces the number of points in your regression. 
  group_by(Date, Site_Num) %>% 
  mutate(ph_log == ifelse(ph_log > 1, ph_log, "NA")) %>%  #Remove faulty pH measurements
  summarize("copies_per_L" = mean(copies_per_L),
            "Vel_ms" = mean(c(flow_1, flow_2, flow_3), na.rm = TRUE), 
            "water_temp" = mean(temp_log, na.rm = TRUE), 
            "pH" = mean(ph_log, na.rm = TRUE), 
            "SC" = mean(sc_log, na.rm = TRUE), 
            "HDO" = mean(hdo_ml.L_log, na.rm = TRUE), 
            "HDO_perc" = mean(hdo_perc_sat_log, na.rm = TRUE), 
            "Turb" = mean(turb_log, na.rm = TRUE)) %>% 
  mutate("Vel_ms" = ifelse(is.na(Vel_ms), mean(Vel_ms, na.rm = T), Vel_ms), 
  "water_temp" = ifelse(is.na(water_temp), mean(water_temp, na.rm = T), water_temp),
  "pH" = ifelse(is.na(pH), mean(pH, na.rm = T), pH),
  "SC" = ifelse(is.na(SC), mean(SC, na.rm = T), SC),
  "HDO" = ifelse(is.na(HDO), mean(HDO, na.rm = T), HDO),
  "HDO_perc" = ifelse(is.na(HDO_perc), mean(HDO_perc, na.rm = T), HDO_perc),
  "Turb" = ifelse(is.na(Turb), mean(Turb, na.rm = T), Turb)) 




st_effort_dat <- merge(x = st_effort_dat, y = edna_sums, 
                       by = c("Date", "Site_Num"), all.x = F, all.y = F)
str(st_effort_dat) #28 obs


st_effort_dat <- st_effort_dat[-c(1,23),] #remove outliers: Belle creek on 6/13/22, and CPC reach 33 on  8/10/22
st_effort_dat <- st_effort_dat[-c(4,17),] #remove Angel creek on 6/22/22 and 7/29/22, when there was no water

#View(st_effort_dat) #24 obs to use for modeling





# Exploratory Analysis ----------------------------------------------------




###Should I log transform this or not?????????????

plot(log(copies_per_L) ~ MGMS_CPUE_abun, data = st_effort_dat_all)
summary(lm(copies_per_L ~ MGMS_CPUE_abun+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat_all))


plot(log(copies_per_L) ~ MGMS_CPUE_biom, data = st_effort_dat_all)
text(log(copies_per_L) ~ MGMS_CPUE_biom, labels=Site_Num,data=st_effort_dat_all, cex=0.9, font=2, pos = 3)
#text(log(copies_per_L)~MGMS_CPUE_biom, labels=Date,data=st_effort_dat, cex=0.9, font=2, pos = 3)
summary(lm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO_perc+Turb, data = st_effort_dat_all))

hist(grayling_dat$Fork_Length) #make this a multi panel plot for each site/sampling event (Date). 






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



m.global.abun <- glm(copies_per_L ~ MGMS_CPUE_abun+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
summary(m.global.abun)

nagelkerke(m.global.abun)
accuracy(list(m.global.abun))
#Pseudo.R.squared                   With NA     NA -> mean
#McFadden                            0.237        0.146
#Cox and Snell (ML)                  0.995        0.944
#Nagelkerke (Cragg and Uhler)        0.965        0.944
#Efron.r.squared                     0.221        0.222

options(na.action = "na.fail")
global_dredge_abun <- dredge(m.global.abun, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge_abun

#Dredge prefers ~MGMS_CPUE_abun+Vel_ms




m.global.biom <- glm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
summary(m.global.biom)

nagelkerke(m.global.biom)
accuracy(list(m.global.biom))
#Pseudo.R.squared                   With NA     NA -> mean
#McFadden                            0.237        0.146
#Cox and Snell (ML)                  0.995        0.944
#Nagelkerke (Cragg and Uhler)        0.965        0.944
#Efron.r.squared                     0.212        0.208


options(na.action = "na.fail")
global_dredge_biom <- dredge(m.global.biom, trace = 2, evaluate = TRUE)
options(na.action = "na.omit")
global_dredge
#Dredge prefers ~HDO+pH



#GLM with Abundnace

#GLM with Biomass

#GLM with Abundance Density

#GLM with Biomass Density


#Dredge the best model of each? 


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
  
sum(comp_abun$diff_abs)#Difference between the observed and predicted values 76473.62
comp_abun <- comp_abun %>% mutate()



#Testing the global abundance model
global_pred_biom <- vector()

for(i in 1:93){
  validate.dat <- st_effort_dat_all[i,] #Single out one value
  training.dat <- st_effort_dat_all[-i,] #Use the training data using the data -i
  #m.global <- glm(copies_per_L ~ MGMS_CPUE_biom+Vel_ms+water_temp+pH+SC+HDO+Turb, data = st_effort_dat)
  global_pred_biom[i] <- predict(m.global.biom, newdata = validate.dat, type = "response")
  
}


comp_biom$global_pred_biom <- as.numeric(global_pred_biom)

comp_biom <- vector()
comp_biom <- data.frame(obs, global_pred_biom)


comp_biom <- comp_biom %>%
  mutate("diff" = obs-global_pred_abun) %>% 
  mutate("diff_abs" = ifelse(diff < 0, -diff, diff))

sum(comp_biom$diff_abs)#Difference between the observed and predicted values 76473.62


# Predicting to other sites -----------------------------------------------






# Plotting ----------------------------------------------------------------



ggplot(data=comp_abun, aes(x=obs, y = global_pred_biom,))+
  geom_point()+
  theme_cowplot()


