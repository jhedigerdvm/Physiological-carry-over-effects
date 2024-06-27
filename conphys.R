
#install.packages(c('tidyverse', 'lme4', 'lmerTest', 'AICcmodavg', 'car', 'mgcv', 'MuMIn'))


library(tidyverse) #help with workflow
library(lme4) #mixed effect models
library(lmerTest) #mixed effect models with p-values
library(AICcmodavg) #rank models
library(car) #check model errors
library(MuMIn) #dredging to try ever combination of model parameters
library(rmarkdown)
library(here)
library(ggplot2)

All <- read.csv('../raw/FinalFormat_EnvirCort_Females.csv', header = T) #Read in full dataset with Fecal Cort Levels, Environmental Data, and animal biological data

for(i in 1:nrow(All)){
  
  if(All$Age[i] == .5) {All$ageclass[i] <- "fawn"} else if(All$Age[i] %in% c(1.5, 2.5, 3.5)){All$ageclass[i] <- "imm"} else if(All$Age[i] %in% c(4.5, 5.5, 6.5))
  {All$ageclass[i] <- "mat"} 
}


All <- All %>% mutate_at(19:29, scale) %>% #Scale and center all predictor variables
  filter(!is.na(Lactation)) %>% #remove any missing entries for lactation
  mutate_at(c("LeftMaster", "Year", "ageclass", "Capt_Loc", "Lactation"), factor) #Factor some of the biological data and random effects

All$Log_Fecal <- log(All$Fecal) #Log transform the fecal cort levels to approximate normality in the model errors


Does <- All %>% filter(Age %in% c(2.5, 3.5, 4.5, 5.5, 6.5)) #filter dataset to just reproductively mature does 
Does1 <- Does[!duplicated(Does$LeftMaster, na.rm = TRUE) | duplicated(Does$LeftMaster, fromLast = TRUE, na.rm = TRUE), ] #this code removes any duplicates but keeps the first occasion
duplicated(Does1$LeftMaster)
Does <- Does1

#As the rut approaches, increase testosterone may influence buck cort levels. 
for(i in 1:nrow(Does)) { #For loop to cycle through each individual row
  
  #i <- 1
  Does$StudyJulian[i] <- as.Date(Does$Date[i]) - as.Date(paste0(Does$Year[i], '-10-01')) #Estimate the number of days since October 1 of that year
}

head(All)
#rain - milimeters
#biomass - kg/ha
#sand - %
#brush cover - %
#lactation - na(males), yes, no
#Gross BC (Boone & Crockett) - cm
#AdjWt - kg


#this creates a vector of the column names containing important rainperiods
rain.var <- c("rain_1month_Previous", "rain_2month_Previous", "rain_July_August")

#create empty list to store competitive models
f.model <- list()

#part of a progress bar
k <- 1

#This for loop cycles through each rain period variable and substitutes it the model formula
for (i in rain.var) {
  
  #paste the progress bar message
  cat('Working on', i, ':', k, 'of 3', '\n')
  
  #i <- rain.var[1] #debugger variable
  
  #the global model formula for the full data set and only sexually mature females
  f <- paste0('Log_Fecal ~ ageclass + Lactation*(', i, ' + Percent_Sand + WoodyCover + BCS + StudyJulian) + (1|Capt_Loc) + (1|Year)')
  
  #fit the global model with sexually mature females
  female <- lmer(data = Does, formula = as.formula(f), na.action = na.fail)
  
  #This creates every model for every parameter combination in the global model while specifying Ageclass as being in every model
  f.combo <- dredge(female)
  
  #This gets the top performing models after dredging with a delta AIC less than 5
  out3 <- get.models(f.combo, subset = delta < 5)
  
  for (x in 1:length(out3)){
    names(out3)[x] <- paste((as.character(formula(out3[[x]]))[3]), collapse = ' ')
  }
  
  #this saves outs to final list of top performing models
  f.model[[i]] <- out3
  
  #update the progress message
  k <- k + 1
}

f.model <- unlist(f.model)

test3 <- separate(data = data.frame(name = names(f.model)), col = name, into = c('var', 'model'), sep = '\\.')
test3$Dup <- !(duplicated(test3$model))

f.model <- f.model[test3$Dup]
names(f.model) <- test3$model[test3$Dup]
```
Evaluation of these models for sexually mature females shows that July-August rainfall is important influence on fall cort levels but the effect is dependent on lactation


#for the sexually mature females the model with an interaction between lactation status and July-Aug rainfall was the best fitting
aic<-aictab(f.model)
#model summary for this top performing model
summary(f.model[['Lactation + rain_July_August + (1 | Capt_Loc) + (1 | Year) + Lactation:rain_July_August']])
summary(f.model[['WoodyCover + (1 | Capt_Loc) + (1 | Year)']])

#coef from model summary
ef <- fixef(f.model[['Lactation + rain_July_August + (1 | Capt_Loc) + (1 | Year) + Lactation:rain_July_August']])

#Next calculate the slope of rain effect on cort of lactating female
slope <- ef["rain_July_August"] + ef["LactationYes - Lactating:rain_July_August"]

#This coef for rain is for a change in 1 SD of rain so we need to know this sd
ch.rain <- (1* attr(All$rain_July_August, 'scaled:scale') + attr(All$rain_July_August, 'scaled:center')) - (0* attr(All$rain_July_August, 'scaled:scale') + attr(All$rain_July_August, 'scaled:center'))

#so log(cort) is reduced -0.3743819 for every 52.24 mm increase in rain, now we want this change per cm of rain
crt.rain <- (slope)/(ch.rain/10)

#now to get an interpretable change in value in normal cort units 
final <- (exp(crt.rain) - 1)*100

#so this how much in % of how cort changes per cm increase in rain.
final

#lets get our coef from model summary
ef.brush <- fixef(f.model[['WoodyCover + (1 | Capt_Loc) + (1 | Year)']])

#Next calculate the slope of brush effect on cort of females
slope.brush <- ef.brush["WoodyCover"]

#This coef for brush is for a change in 1 SD of brush so we need to know this sd
ch.brush <- (1* attr(All$WoodyCover, 'scaled:scale') + attr(All$WoodyCover, 'scaled:center')) - (0* attr(All$WoodyCover, 'scaled:scale') + attr(All$WoodyCover, 'scaled:center'))

#so log(cort) increases 0.1821032 for every  12.06459% increase in brush, now we want this change per 1% increase in brush
crt.brush <- (slope.brush)/(ch.brush)

#now to get an interpretable change in value in normal cort units 
final.brush <- (exp(crt.brush) - 1)*100

#so this how much in % of how cort changes per 1% increase in brush cover.
final.brush


