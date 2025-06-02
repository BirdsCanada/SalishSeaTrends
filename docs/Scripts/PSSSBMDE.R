psss_to_bmde <- function(PSSS) {

  # Read PSSS data
  PSSS <- read.csv("Data/PSSS.csv", stringsAsFactors = FALSE)
  
  # Parse latitude and longitude from position string
  PSSS$lat <- sub(" W.*", "", PSSS$position)
  PSSS$long <- sub(".*W", "", PSSS$position)
  PSSS$lat <- gsub('N', '', PSSS$lat)
  PSSS$long <- gsub('W', '', PSSS$long)
  PSSS$DecimalLatitude <- as.numeric(conv_unit(PSSS$lat, from = 'deg_dec_min', to = 'dec_deg'))
  PSSS$DecimalLongitude <- as.numeric(conv_unit(PSSS$long, from = 'deg_dec_min', to = 'dec_deg')) * -1
  
  # Handle sites with different lat/long format
  test <- PSSS %>%
    filter(is.na(DecimalLatitude)) %>%
    separate(lat, into = c("DecimalLatitude", "DecimalLongitude"), sep = ",")
  PSSS <- PSSS %>%
    filter(!is.na(DecimalLatitude)) %>%
    select(-lat)
  PSSS <- rbind(PSSS, test)
  
  # Split survey_date into Month, Day, Year
  PSSS <- PSSS %>%
    separate(survey_date, into = c("Date", "del"), sep = " ") %>%
    select(-del) %>%
    separate(Date, into = c("MonthCollected", "DayCollected", "YearCollected"), sep = "/")
  
  # Wrangle raptor data into long format
  raptor1 <- PSSS %>%
    filter(!is.na(raptor1)) %>%
    mutate(common_name = raptor1, bird_count = raptor1_count, notes = raptor1_affect) %>%
    select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) %>%
    group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>%
    mutate(bird_count = sum(bird_count)) %>%
    distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all = TRUE)
  raptor2 <- PSSS %>%
    filter(!is.na(raptor2)) %>%
    mutate(common_name = raptor2, bird_count = raptor2_count, notes = raptor2_affect) %>%
    select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) %>%
    group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>%
    mutate(bird_count = sum(bird_count)) %>%
    distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all = TRUE)
  raptor3 <- PSSS %>%
    filter(!is.na(raptor3)) %>%
    mutate(common_name = raptor3, bird_count = raptor3_count, notes = raptor3_affect) %>%
    select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) %>%
    group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>%
    mutate(bird_count = sum(bird_count)) %>%
    distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all = TRUE)
  
  raptor <- bind_rows(raptor1, raptor2, raptor3) %>%
    filter(common_name == "Bald Eagle")
  
  # Remove raptor columns and bind raptor data back
  PSSS <- PSSS %>%
    select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count,
           -raptor1_affect, -raptor2_affect, -raptor3_affect)
  PSSS <- bind_rows(PSSS, raptor)
  
  # If common name is missing, assign from bird_name
  PSSS <- PSSS %>%
    mutate(common_name = ifelse(common_name == "", bird_name, common_name))
  
  # Remove bearing and distance
  PSSS <- PSSS %>% select(-bearing, -dist)
  
  # Replace Thayer's Gull with Ivory Gull
  PSSS <- PSSS %>%
    mutate(common_name = ifelse(common_name == "Thayer's Gull", "Ivory Gull", common_name))
  
  # Merge with species table
  PSSS <- merge(PSSS, sp, by.x = "common_name", by.y = "english_name", all.x = TRUE)
  
  # Rename columns for BMDE
  PSSS <- PSSS %>%
    rename(
      CommonName = common_name,
      SurveyAreaIdentifier = site_code,
      Locality = site_name,
      MinimumElevationInMeters = elevation,
      MaximumElevationInMeters = elevation,
      TimeCollected = start_time,
      TimeObservationsEnded = end_time,
      ObservationCount = bird_count,
      ObservationCount2 = large_flock_best,
      ObsCountAtLeast = large_flock_min,
      ObsCountAtMost = large_flock_max,
      FieldNotes = notes,
      Collector = name,
      ScientificName = scientific_name,
      SpeciesCode = species_code,
      AllSpeciesReported = is_complete
    )
  
  # Add/modify metadata columns
  PSSS$RouteIdentifier <- PSSS$SurveyAreaIdentifier
  PSSS$BasisOfRecord <- "Observation"
  PSSS$CollectionCode <- "PSSS"
  PSSS$Continent <- "North America"
  PSSS$Country <- "United States"
  PSSS$StateProvince <- "Washington"
  PSSS$ProtocolType <- "PointCount"
  PSSS$ProtocolSpeciesTargeted <- "Waterbirds"
  PSSS$ProtocolURL <- "https://seattleaudubon.org/wp-content/uploads/2021/01/PSSS_Protocol_2014-15.pdf"
  PSSS$SurveyAreaShape <- "300 m"
  PSSS$ObservationDescriptor <- "Total Count"
  PSSS$ObservationDescriptor2 <- "Large flock best estiamte"
  PSSS$TimeObservationsStarted <- PSSS$TimeCollected
  PSSS$ProjectCode <- "PSSS"
  PSSS$InstitutionCode <- "PSBO"
  PSSS$CatalogNumber <- PSSS$survey_id
  PSSS$GlobalUniqueIdentifier <- paste0("URN:catalog:", PSSS$InstitutionCode, ":", PSSS$CollectionCode, ":", PSSS$CatalogNumber)
  PSSS$ObservationDate <- paste0(PSSS$YearCollected, "-", PSSS$MonthCollected, "-", PSSS$DayCollected)
  PSSS$DecimalLatitude <- as.numeric(PSSS$DecimalLatitude)
  PSSS$DecimalLongitude <- as.numeric(PSSS$DecimalLongitude)
  
  # Add missing columns for BMDE
  BMDE_col <- unique(BMDE$local_name)
  missing <- setdiff(BMDE_col, names(PSSS))
  PSSS[missing] <- ""
  
  # Filter to BMDE columns and write output
  PSSS <- PSSS[BMDE_col]
  write.csv(PSSS, "Data/PSSS_BMDE.csv", row.names = FALSE, quote = FALSE)
  
  invisible(PSSS)
}
