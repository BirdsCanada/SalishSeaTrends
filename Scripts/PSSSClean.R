clean_PSSS <- function(Y1, Y2) {
  # Read raw data
  in.PSSS <- read.csv("Data/PSSS.csv")
  
  # Fix Thayer's Gull species code
  in.PSSS$SpeciesCode[in.PSSS$CommonName == "Iceland (Thayer's) Gull"] <- "ICGU"
  in.PSSS$SpeciesCode <- as.factor(in.PSSS$SpeciesCode)
  
  # Filter months October-April
  in.PSSS <- subset(in.PSSS, MonthCollected %in% c(10:12, 1:4))
  
  # Remove duplicates
  in.PSSS <- distinct(in.PSSS)
  
  # Create time tracking columns
  in.PSSS$doy <- as.numeric(format(as.Date(paste(in.PSSS$YearCollected, in.PSSS$MonthCollected, in.PSSS$DayCollected, sep = "-")), "%j"))
  in.PSSS$wmonth <- ifelse(in.PSSS$MonthCollected >= 10, in.PSSS$MonthCollected - 9, in.PSSS$MonthCollected + 3)
  
  # Keep first survey per month
  in.PSSS <- in.PSSS %>% group_by(SurveyAreaIdentifier, YearCollected, MonthCollected) %>% slice_min(doy) %>% ungroup()
  
  # Create winter year column
  in.PSSS$wyear <- ifelse(in.PSSS$MonthCollected %in% 1:4, in.PSSS$YearCollected - 1, in.PSSS$YearCollected)
  
  # Filter by input years
  in.PSSS <- subset(in.PSSS, wyear >= Y1 & wyear <= Y2)
  
  # Clean spatial data
  in.PSSS <- in.PSSS %>% 
    filter(!is.na(DecimalLatitude), !is.na(DecimalLongitude)) %>%
    filter(DecimalLatitude >= 45.06 & DecimalLatitude <= 50.64 & 
             DecimalLongitude >= -125.07 & DecimalLongitude <= -115.15)
  
  # Process time data
  in.PSSS <- in.PSSS %>%
    filter(nchar(TimeObservationsStarted) %in% 3:4, nchar(TimeObservationsEnded) %in% 3:4) %>%
    mutate(
      DecimalTimeObservationsStarted = (TimeObservationsStarted %/% 100) + ((TimeObservationsStarted %% 100)/60),
      DecimalTimeObservationsEnded = (TimeObservationsEnded %/% 100) + ((TimeObservationsEnded %% 100)/60)
    )
  
  # Calculate duration (requires calculate_duration() function)
  in.PSSS$DurationInHours2 <- calculate_duration(in.PSSS$DecimalTimeObservationsStarted, in.PSSS$DecimalTimeObservationsEnded)
  
  # Fix reversed times
  in.PSSS <- in.PSSS %>%
    mutate(
      DecimalTimeObservationsStarted2 = ifelse(DurationInHours2 < 0, DecimalTimeObservationsEnded, DecimalTimeObservationsStarted),
      DecimalTimeObservationsEnded2 = ifelse(DurationInHours2 < 0, DecimalTimeObservationsStarted, DecimalTimeObservationsEnded)
    )
  
  # Recalculate duration
  in.PSSS$DurationInHours <- calculate_duration(in.PSSS$DecimalTimeObservationsStarted2, in.PSSS$DecimalTimeObservationsEnded2)
  
  # Final duration filtering
  in.PSSS <- in.PSSS %>% 
    filter(DurationInHours >= 0, DurationInHours > 0.03, DurationInHours < 5)
  
  # Create event matrix
  event.PSSS <- in.PSSS %>%
    select(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected, 
           DecimalLatitude, DecimalLongitude, DurationInHours) %>%
    distinct() %>%
    group_by(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected) %>%
    slice_min(DurationInHours) %>%
    ungroup()
  
  # Clean observation data
  in.PSSS <- in.PSSS %>%
    filter(!is.na(ObservationCount), !is.na(SurveyAreaIdentifier)) %>%
    group_by(ProjectCode, SurveyAreaIdentifier, SpeciesCode, CommonName, wyear, wmonth, YearCollected, MonthCollected, DayCollected) %>%
    summarise(ObservationCount = sum(ObservationCount), .groups = "drop") %>%
    select(ProjectCode, SurveyAreaIdentifier, SpeciesCode, CommonName, ObservationCount, wyear, wmonth, YearCollected, MonthCollected, DayCollected)
  
  # Filter and group species
  in.PSSS <- in.PSSS %>%
    filter(!CommonName %in% c("alcid sp.", "grebe sp.", "gull (small)", "diving duck sp.", "gull sp.", "scoter sp.", 
                              "cormorant sp.", "goldeneye sp.", "dabbling duck sp.", "merganser sp.", "loon sp.")) %>%
    group_by(CommonName) %>%
    filter(sum(ObservationCount) > 10) %>%
    ungroup() %>%
    mutate(
      CommonName = case_match(
        CommonName,
        c("gull (large)", "Glaucous Gull", "Glaucous-winged Gull", "Western Gull", "Herring Gull",
          "Iceland (Thayer's) Gull", "Iceland (Thayer's Gull)", "WEGU x GWGU hybrid", "California Gull") ~ "Large Gull",
        c("scaup sp.", "Lesser Scaup", "Greater Scaup", "Greater-Lesser Scaup") ~ "Greater-Lesser Scaup",
        c("Eared Grebe", "Horned Grebe") ~ "Eared-Horned Grebe",
        c("Canada Goose", "Cackling Goose") ~ "Canada-Cackling Goose",
        c("Clark's Grebe", "Western Grebe") ~ "Western-Clark's Grebe",
        .default = CommonName
      )
    )
  
  # Write outputs
  write.csv(in.PSSS, "Data/PSSS.clean.csv", row.names = FALSE)
  write.csv(event.PSSS, "Data/PSSS.events.csv", row.names = FALSE)
  
  list(clean_data = in.PSSS, event_data = event.PSSS)
}
