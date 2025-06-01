#Code used to clean the PSSS and prepare it for analysis 

#Created by Danielle Ethier
#Last edited April 2025

in.PSSS <- read.csv("Data/PSSS.csv") # reads in back-up copy of database 

#Thayer's Gull missing SpeciesCode. Need to assign here species_id == 5190
in.PSSS$SpeciesCode[in.PSSS$CommonName == "Iceland (Thayer's) Gull"] <- "ICGU"

in.PSSS$SpeciesCode<-as.factor(in.PSSS$SpeciesCode)

# filter data by months October to April
in.PSSS <- subset(in.PSSS, MonthCollected %in% c(10:12, 1:4))

# if duplicate records in data file; this keeps only one. There currently are no duplicates. 
in.PSSS <- distinct(in.PSSS)

#create day of year column 
in.PSSS$doy <- as.numeric(format(as.Date(paste(in.PSSS$YearCollected, in.PSSS$MonthCollected, in.PSSS$DayCollected, sep="-")), "%j"))
#create a new winter survey month column starting in October = 1 to April = 8
in.PSSS$wmonth <-in.PSSS$wmonth <- ifelse(in.PSSS$MonthCollected >= 10, 
                                      in.PSSS$MonthCollected - 9, 
                                      in.PSSS$MonthCollected + 3)

# keep only one survey per month by selecting the first survey if there are duplicates. There are currently none. 
in.PSSS <- in.PSSS %>% group_by(SurveyAreaIdentifier, YearCollected, MonthCollected) %>% slice_min(doy) %>% ungroup()

# create a new column called wyear, which groups surveys by winter year
# for example, January-April 2005 would be wyear 2004
in.PSSS$wyear <- ifelse(in.PSSS$MonthCollected %in% c(1:4), in.PSSS$YearCollected-1, 
                         in.PSSS$YearCollected)

# filter data by years 
in.PSSS <- subset(in.PSSS, wyear >= Y1 & wyear <= Y2)

#Remove missing DecimalLatitude or DecimalLongitude
in.PSSS <- in.PSSS %>% filter(!is.na(DecimalLatitude) & !is.na(DecimalLongitude))
#Filter events to 45.06N to 50.64N latitude and 125.07W to 115.15W longitude
in.PSSS <- in.PSSS %>% filter(DecimalLatitude >= 45.06 & DecimalLatitude <= 50.64 & DecimalLongitude >= -125.07 & DecimalLongitude <= -115.15)

# Because there are errors in the Duration in Hours column, I will adjust the TimeObservationsEnded and TimeObservationStarted column. 
# Specifically, any value that is < 6 AM I will add 12 to bring it into 24 hour time. 

# Remove missing TimeObservationsStarted and TimeObservationsEnded
in.PSSS <- in.PSSS %>% filter(!is.na(TimeObservationsStarted) & !is.na(TimeObservationsEnded))
# Remove times with less than 3 nchar
in.PSSS <- in.PSSS %>% filter(nchar(TimeObservationsStarted) > 2 & nchar(TimeObservationsEnded) > 2)
# REmove times with greater tahn 4 nchar
in.PSSS <- in.PSSS %>% filter(nchar(TimeObservationsStarted) < 5 & nchar(TimeObservationsEnded) < 5)

#convert time into decimal hours
in.PSSS <- in.PSSS %>% mutate(DecimalTimeObservationsStarted = (TimeObservationsStarted %/% 100)+((TimeObservationsStarted %% 100)/60), DecimalTimeObservationsEnded = (TimeObservationsEnded %/% 100)+((TimeObservationsEnded %% 100)/60))

# Calculate the duration
  in.PSSS$DurationInHours2 <- calculate_duration(in.PSSS$DecimalTimeObservationsStarted, in.PSSS$DecimalTimeObservationsEnded)                      
            
# if DurationinHours is negative, we want to swap the survey start time and survey end time as they seem to be reversed.
  in.PSSS<- in.PSSS %>% mutate(DecimalTimeObservationsStarted2= ifelse(in.PSSS$DurationInHours2 < 0, DecimalTimeObservationsEnded, DecimalTimeObservationsStarted), DecimalTimeObservationsEnded2= ifelse(in.PSSS$DurationInHours2 < 0, DecimalTimeObservationsStarted, DecimalTimeObservationsEnded))

# recalculate duration in hours using the start and end times
  in.PSSS$DurationInHours <- calculate_duration(in.PSSS$DecimalTimeObservationsStarted2, in.PSSS$DecimalTimeObservationsEnded2)                      
  
#remove remaining negative DurationInHours
  in.PSSS <- in.PSSS %>% filter(DurationInHours >= 0)  
#Filter Duration in hours greater than 0.3 and less than 3
  in.PSSS<-in.PSSS[in.PSSS$DurationInHours > 0.03 & in.PSSS$DurationInHours < 5,]  
  
  # before we do more cleaning, make your events layer to ensure we have captured all the survey events. 
  # create an events matrix for future zero filling
  
  #events for monthly analysis
  event.PSSS <- in.PSSS %>% dplyr::select(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected, DecimalLatitude, DecimalLongitude, DurationInHours) %>% distinct()
  
  # if there are multiple events in a single day (now caused by Duration in Hours), take the minimum
  event.PSSS <- event.PSSS %>% group_by(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected) %>% slice_min(DurationInHours) %>% ungroup()
  
  #Remove ObservationCounts that are NA
  in.PSSS <- in.PSSS %>% filter(!is.na(ObservationCount))
  #Remove missing SurveyAreaIdentifiers
  in.PSSS <- in.PSSS %>% filter(!is.na(SurveyAreaIdentifier))
  
  #for PSSS we need to combine observation counts from the same survey that have different bearings and distances, otherwise there are duplicates
  in.PSSS <- in.PSSS %>%
    group_by(ProjectCode, SurveyAreaIdentifier, SpeciesCode, CommonName, wyear, wmonth, YearCollected, MonthCollected, DayCollected) %>%
    summarise(ObservationCount = sum(ObservationCount), .groups = 'drop')
  
  # retain columns that are needed for the analysis
  in.PSSS <- in.PSSS %>% dplyr::select(ProjectCode, SurveyAreaIdentifier, SpeciesCode, CommonName, ObservationCount,  wyear, wmonth, YearCollected, MonthCollected, DayCollected)

  # clean up some of the species names and codes
  in.PSSS<-in.PSSS %>% filter(!is.na(CommonName)) %>% filter(!CommonName %in% c("alcid sp.", "grebe sp.", "gull (small)", "diving duck sp.", "gull sp.", "scoter sp.", "cormorant sp.", "goldeneye sp.", "dabbling duck sp.", "merganser sp.", "loon sp."))
  
  #remove species that are detected less than 10 times over all years
  sample<-in.PSSS %>% group_by(CommonName) %>% summarise(n_tot = sum(ObservationCount, na.rm=TRUE)) %>% filter(n_tot>10)
  sample<-sample %>% filter(CommonName != "") %>% select(-n_tot)
  list<-sample$CommonName
  in.PSSS<-in.PSSS %>% filter(CommonName %in% list)
  
  #Group some species that are hard to identity
  in.PSSS <- in.PSSS %>%
    mutate(
      CommonName = case_match(
        CommonName,
        c("gull (large)", "Glaucous Gull", "Glaucous-winged Gull", "Western Gull", "Herring Gull", "Iceland (Thayer's) Gull", "Iceland (Thayer's Gull)", "WEGU x GWGU hybrid", "California Gull") ~ "Large Gull",
        .default = CommonName), 
    CommonName = case_match(
      CommonName,
      c("scaup sp.", "Lesser Scaup", "Greater Scaup", "Greater-Lesser Scaup") ~ "Greater-Lesser Scaup",
        .default = CommonName
      ), 
    CommonName = case_match(
      CommonName,
      c("Eared Grebe", "Horned Grebe") ~ "Eared-Horned Grebe",
        .default = CommonName
      ), 
    CommonName = case_match(
      CommonName,
      c("Canada Goose", "Cackling Goose") ~ "Canada-Cackling Goose",
        .default = CommonName
      ), 
    CommonName = case_match(
      CommonName,
      c("Clark's Grebe", "Western Grebe") ~ "Western-Clark's Grebe",
        .default = CommonName
      ))
   
  # write index.data to file
  write.csv(in.PSSS, "Data/PSSS.clean.csv")
  write.csv(event.PSSS, "Data/PSSS.events.csv")
 