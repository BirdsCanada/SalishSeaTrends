clean_BCCWS <- function(Y1, Y2) {
  # Read raw data
  in.BCCWS <- read.csv("Data/BCCWS_BMDE.csv")
  
  # Adjust SpeciesCode for MEW GULL species_id 5140 and 5142
  in.BCCWS$SpeciesCode[in.BCCWS$species_id == 5140] <- "SBIG"
  in.BCCWS$SpeciesCode[in.BCCWS$species_id == 5142] <- "SBIG"
  in.BCCWS$SpeciesCode <- as.factor(in.BCCWS$SpeciesCode)
  
  # Use Nearshore observation counts only
  in.BCCWS$ObservationCount <- as.numeric(in.BCCWS$ObservationCount3)
  
  # Merge species info (assuming 'sp' dataframe is loaded in environment)
  in.BCCWS <- merge(in.BCCWS, sp, by.x = "CommonName", by.y = "english_name", all.x = TRUE)
  
  # Filter data by months October to April
  in.BCCWS <- subset(in.BCCWS, MonthCollected %in% c(10:12, 1:4))
  
  # Filter out known bad date record
  in.BCCWS <- in.BCCWS %>% filter(!(SurveyAreaIdentifier == "LMSQ-3" & SpeciesCode == "BAEA" & YearCollected == 1999 & MonthCollected == 12 & DayCollected == 12))
  
  # Parse form ID and remove bad forms
  in.BCCWS$form.id <- gsub("BCCWS-", "", in.BCCWS$SamplingEventIdentifier)
  in.BCCWS <- subset(in.BCCWS, !(form.id %in% c(3794, 5469, 5063, 6945)))
  
  #Remove route being done for joint venture research (primarily in marshes)
  remove <- c("LM-GR-ctrl-WS01", "LM-GR-ctrl-WS02", "LM-SA-WS01", "LM-SA-WS02", 
              "LM-SA-ctrl-WS01", "LM-BU-ctrl-WS01", "LM-BU-WS01", "LM-MS-WS01", 
              "LM-MS-WS02", "LM-MS-ctrl-WS01", "LM-MS-ctrl-WS02", "LM-WV-ctrl-WS01", 
              "LM-WV-WS01", "LM-WC-WS01", "LM-WC-ctrl-WS01", "LM-GE-WS01", 
              "LM-SW-WS01", "LM-SW-GE-ctrl-WS01", "LM-NC-WS01", "LM-NC-WS02", 
              "LM-NC-ctrl-WS01", "LM-BE-ctrl-WS01", "LM-BE-WS01", "LM-BE-WS02")
  
  in.BCCWS <- in.BCCWS %>% 
    filter(!SurveyAreaIdentifier %in% remove)
  
  
  # Remove duplicate records
  in.BCCWS <- distinct(in.BCCWS)
  
  # Create day of year column
  in.BCCWS$doy <- as.numeric(format(as.Date(paste(in.BCCWS$YearCollected, in.BCCWS$MonthCollected, in.BCCWS$DayCollected, sep = "-")), "%j"))
  
  # Create survey month column (October=1 to April=8)
  in.BCCWS$wmonth <- ifelse(in.BCCWS$MonthCollected >= 10, in.BCCWS$MonthCollected - 9, in.BCCWS$MonthCollected + 3)
  
  # Keep only one survey per month (first survey by doy)
  in.BCCWS <- in.BCCWS %>% group_by(SurveyAreaIdentifier, YearCollected, MonthCollected) %>% slice_min(doy) %>% ungroup()
  
  # Create winter year column (Jan-Apr belong to previous year)
  in.BCCWS$wyear <- ifelse(in.BCCWS$MonthCollected %in% 1:4, in.BCCWS$YearCollected - 1, in.BCCWS$YearCollected)
  
  # Filter by winter years input
  in.BCCWS <- subset(in.BCCWS, wyear >= Y1 & wyear <= Y2)
  
  # Fix TimeObservationsStarted and Ended columns (<6 add 12)
  in.BCCWS$TimeObservationsEnded <- as.numeric(in.BCCWS$TimeObservationsEnded)
  in.BCCWS$TimeObservationsStarted <- as.numeric(in.BCCWS$TimeObservationsStarted)
  in.BCCWS$TimeObservationsEnded <- ifelse(in.BCCWS$TimeObservationsEnded < 6, in.BCCWS$TimeObservationsEnded + 12, in.BCCWS$TimeObservationsEnded)
  in.BCCWS$TimeObservationsStarted <- ifelse(in.BCCWS$TimeObservationsStarted < 6, in.BCCWS$TimeObservationsStarted + 12, in.BCCWS$TimeObservationsStarted)
  
  # Remove invalid times and durations
  in.BCCWS <- in.BCCWS %>% filter(TimeObservationsEnded < 24.59, DurationInHours > 0)
  
  # Remove rows with missing DurationInHours, SurveyAreaIdentifier, Latitude or Longitude
  in.BCCWS <- in.BCCWS %>% filter(!is.na(DurationInHours), !is.na(SurveyAreaIdentifier), !is.na(DecimalLatitude), !is.na(DecimalLongitude))
  
  # Filter spatial bounds
  in.BCCWS <- in.BCCWS %>% filter(DecimalLatitude >= 45.06 & DecimalLatitude <= 50.64 & DecimalLongitude >= -125.07 & DecimalLongitude <= -115.15)
  in.BCCWS <- in.BCCWS %>% filter(!(DecimalLatitude < 48.5 & DecimalLongitude < -124))
  
  # Filter DurationInHours between 0.3 and 10
  in.BCCWS <- in.BCCWS %>% filter(DurationInHours > 0.3 & DurationInHours < 10)
  
  # Create event matrix for zero filling
  event.BCCWS <- in.BCCWS %>%
    select(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected, DecimalLatitude, DecimalLongitude, DurationInHours) %>%
    distinct()
  
  # If multiple events in a day, take minimum DurationInHours
  event.BCCWS <- event.BCCWS %>% group_by(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected) %>% slice_min(DurationInHours) %>% ungroup()
  
  # Ensure single lat/lon per SurveyAreaIdentifier (take first)
  event.BCCWS <- event.BCCWS %>% group_by(ProjectCode, SurveyAreaIdentifier) %>% slice_min(DecimalLatitude) %>% ungroup()
  
  # Retain needed columns for analysis
  in.BCCWS <- in.BCCWS %>% select(ProjectCode, SurveyAreaIdentifier, SpeciesCode, CommonName, ObservationCount, wyear, YearCollected, wmonth, MonthCollected, DayCollected)
  
  # Remove species detected less than 10 times over all years
  sample <- in.BCCWS %>% group_by(CommonName) %>% summarise(n_tot = sum(ObservationCount, na.rm = TRUE)) %>% filter(n_tot > 10)
  sample <- sample %>% filter(CommonName != "") %>% select(-n_tot)
  list <- sample$CommonName
  in.BCCWS <- in.BCCWS %>% filter(CommonName %in% list)
  
  # Group some species that are hard to identify
  in.BCCWS <- in.BCCWS %>%
    mutate(
      CommonName = case_when(
        CommonName %in% c("gull (large)", "Iceland Gull", "Iceland Gull (Thayer's)", "Western x Glaucous-winged Gull (hybrid)", "Glaucous Gull", "Glaucous-winged Gull", "Western Gull", "Herring Gull", "Iceland (Thayer's) Gull", "Iceland (Thayer's Gull)", "WEGU x GWGU hybrid", "California Gull") ~ "Large Gull",
        CommonName %in% c("scaup sp.", "Lesser Scaup", "Greater Scaup", "Greater/Lesser Scaup") ~ "Greater-Lesser Scaup",
        CommonName %in% c("Eared Grebe", "Horned Grebe") ~ "Eared-Horned Grebe",
        CommonName %in% c("Canada Goose", "Cackling Goose") ~ "Canada-Cackling Goose",
        CommonName %in% c("Clark's Grebe", "Western Grebe") ~ "Western-Clark's Grebe",
        TRUE ~ CommonName
      )
    )
  
  
  
  # Remove rows with NA CommonName or ObservationCount
  in.BCCWS <- in.BCCWS %>% filter(!is.na(CommonName), !is.na(ObservationCount))
  
  # Write cleaned data to files
  write.csv(in.BCCWS, "Data/BCCWS.clean.csv", row.names = FALSE)
  write.csv(event.BCCWS, "Data/BCCWS.events.csv", row.names = FALSE)
  
  # Return cleaned data frames invisibly if needed
  list(in.BCCWS = in.BCCWS, event.BCCWS = event.BCCWS)
}
