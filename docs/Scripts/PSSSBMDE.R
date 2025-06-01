
sp.code<-sp.code %>% filter(authority=="BSCDATA") %>% dplyr::select(-authority, -species_id2, -rank) %>% distinct()
sp.tax<-sp.tax %>% dplyr::select(species_id, scientific_name, english_name) %>% distinct()

sp<-left_join(sp.code, sp.tax, by="species_id")
sp<-sp %>% distinct(english_name, .keep_all = TRUE)

#remove birds that are outside survey area (>300m)
# PSSS <- PSSS %>%
#   mutate(bearing = case_when(bearing == "" ~ NA, .default = bearing)) %>%
#   mutate(within_survey_area = case_when(is.na(bearing) == FALSE & is.na(dist) == FALSE & is.na(bird_count) == FALSE ~ TRUE,
#                                         horizon_obscured == 1 & is.na(bird_count) == FALSE ~ TRUE,
#                                         is.na(bearing) == FALSE & is.na(dist) == TRUE & is.na(bird_count) == FALSE ~ TRUE,
#                                         is.na(bearing) == TRUE & is.na(dist) == FALSE & is.na(bird_count) == FALSE ~ TRUE,
#                                         .default = FALSE)) %>% filter(is_complete == 1) %>% filter(within_survey_area==TRUE) %>% dplyr::select(-within_survey_area)

PSSS$lat<-sub(" W.*", "", PSSS$position)  
PSSS$long<-sub(".*W", "", PSSS$position)

PSSS$lat = gsub('N', '', PSSS$lat)
PSSS$long = gsub('W', '', PSSS$long)

PSSS$DecimalLatitude = measurements::conv_unit(PSSS$lat, from = 'deg_dec_min', to = 'dec_deg')
PSSS$DecimalLatitude<-as.numeric((PSSS$DecimalLatitude))
PSSS$DecimalLongitude = measurements::conv_unit(PSSS$long, from = 'deg_dec_min', to = 'dec_deg')
PSSS$DecimalLongitude<-as.numeric(PSSS$DecimalLongitude)
PSSS$DecimalLongitude=PSSS$DecimalLongitude*(-1)

##several sites have lat and long in a different format, which means they manually need fixed. 
#I have flagged this to Toby to fix in the underlying data. 

test<-PSSS %>% filter(is.na(DecimalLatitude)) %>% separate(lat, into=c("DecimalLatitude", "DecimalLongitude"), sep=",")
PSSS<-PSSS %>% filter(!is.na(DecimalLatitude)) %>% select(-lat)
PSSS<-rbind(PSSS, test)

# #manually assign the decimal lat and long from the location table fields DecimalLatitude and DecimalLongitude and DecminLatitude and DecminLongitude
# PSSS$DecimalLatitude[PSSS$site_code==60]<-47.12093
# PSSS$DecimalLongitude[PSSS$site_code==60]<--122.77609
# PSSS$DecimalLatitude[PSSS$site_code==131]<-47.596750
# PSSS$DecimalLongitude[PSSS$site_code==131]<--122.541700
# PSSS$DecimalLatitude[PSSS$site_code==194]<-48.61269
# PSSS$DecimalLongitude[PSSS$site_code==194]<--123.09778
# PSSS$DecimalLatitude[PSSS$site_code==108]<-47.16545
# PSSS$DecimalLongitude[PSSS$site_code==108]<--122.80804
# PSSS$DecimalLatitude[PSSS$site_code==171]<-48.61269
# PSSS$DecimalLongitude[PSSS$site_code==171]<--123.09778

# #break apart survey_date and reform into day, month, year
# PSSS<-PSSS %>% separate(survey_date, into=c("Date", "del"), sep=" ") %>% dplyr::select(-del) %>% separate(Date, into=c("MonthCollected", "DayCollected", "YearCollected"), sep="/") 

#break apart survey_date and reform into day, month, year
PSSS<-PSSS %>% separate(survey_date, into=c("Date", "del"), sep=" ") %>% dplyr::select(-del) %>% 
  separate(Date, into=c("MonthCollected", "DayCollected", "YearCollected"), sep="/") 


#wrangle raptor data into the long format since each species identification should be in a unique row. 
raptor1<-PSSS %>% filter(!is.na(raptor1)) %>% mutate(common_name = raptor1, bird_count = raptor1_count, notes= raptor1_affect)%>%  dplyr::select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) 
raptor1<-raptor1 %>% group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>% mutate(bird_count=sum(bird_count)) %>% distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all=TRUE)
raptor2<-PSSS %>% filter(!is.na(raptor2)) %>% mutate(common_name = raptor2, bird_count = raptor2_count, notes= raptor2_affect)%>%  dplyr::select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) 
raptor2<-raptor2 %>% group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>% mutate(bird_count=sum(bird_count)) %>% distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all=TRUE)
raptor3<-PSSS %>% filter(!is.na(raptor3)) %>% mutate(common_name = raptor3, bird_count = raptor3_count, notes= raptor3_affect) %>%  dplyr::select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) 
raptor3<-raptor3 %>% group_by(site_name, common_name, YearCollected, MonthCollected, DayCollected) %>% mutate(bird_count=sum(bird_count)) %>% distinct(common_name, site_name, YearCollected, MonthCollected, DayCollected, .keep_all=TRUE)

#bind raptor data back with PSSS data
raptor<-rbind(raptor1, raptor2)
raptor<-rbind(raptor, raptor3)
raptor<-raptor %>% filter(common_name =="Bald Eagle")

#now bind back with PSSS data
PSSS<-PSSS %>%  dplyr::select(-raptor1, -raptor2, -raptor3, -raptor1_count, -raptor2_count, -raptor3_count, -raptor1_affect, -raptor2_affect, -raptor3_affect) 
PSSS<-rbind(PSSS, raptor)

#If common name is misisng, assign the "bird_name" hybrid field
PSSS<-PSSS %>% mutate(common_name = ifelse(common_name =="", bird_name, common_name))

#remove rows with missing common name
#PSSS<-PSSS %>% filter(common_name !="")

#remove bearing and distance because we want each species/ site/ date to be a single row in the data set similar to BBCWS
PSSS<-PSSS %>% dplyr::select(-bearing, -dist)

#replace Thayer's Gull with Ivory Gull
PSSS<-PSSS %>% mutate(common_name = ifelse(common_name == "Thayer's Gull", "Ivory Gull", common_name))

#Merge with species ID
PSSS<-merge(PSSS, sp, by.x=c("common_name"), by.y= ("english_name"), all.x=TRUE)

#rename data columns to match BMDE

PSSS<-PSSS %>% dplyr::rename(CommonName = common_name, SurveyAreaIdentifier= site_code, Locality = site_name, MinimumElevationInMeters=elevation, MaximumElevationInMeters=elevation, TimeCollected = start_time, TimeObservationsEnded=end_time, ObservationCount = bird_count, ObservationCount2=large_flock_best, ObsCountAtLeast = large_flock_min, ObsCountAtMost = large_flock_max, FieldNotes=notes, Collector = name, ScientificName=scientific_name, SpeciesCode=species_code, AllSpeciesReported=is_complete)

PSSS$RouteIdentifier<-PSSS$SurveyAreaIdentifier
PSSS$BasisOfRecord <- "Observation"
PSSS$CollectionCode <- "PSSS"
PSSS$Continent <-"North America"
PSSS$Country<-"United States"
PSSS$StateProvince<-"Washington"
PSSS$ProtocolType <- "PointCount"
PSSS$ProtocolSpeciesTargeted <- "Waterbirds"
PSSS$ProtocolURL= "https://seattleaudubon.org/wp-content/uploads/2021/01/PSSS_Protocol_2014-15.pdf"
PSSS$SurveyAreaShape = "300 m"
#PSSS$EffortUnit1 = "Party-hours"
PSSS$ObservationDescriptor = "Total Count"
PSSS$ObservationDescriptor2 = "Large flock best estiamte" 
PSSS$TimeObservationsStarted=PSSS$TimeCollected
PSSS$ProjectCode = "PSSS"
PSSS$InstitutionCode<-"PSBO"
PSSS$CatalogNumber<-PSSS$survey_id
PSSS<-PSSS %>% mutate(GlobalUniqueIdentifier = paste0("URN:catalog:", InstitutionCode, ":", CollectionCode, ":", CatalogNumber, sep=""))
PSSS<-PSSS %>% mutate(ObservationDate = paste0(YearCollected, "-", MonthCollected, "-", DayCollected, sep=""))
PSSS$DecimalLatitude<-as.numeric(PSSS$DecimalLatitude)
PSSS$DecimalLongitude<-as.numeric(PSSS$DecimalLongitude)

#Now that we have specified all the data columns we can, we will create the BMDE standardized data table. 

#Identify the missing columns of data
BMDE_col<-unique(BMDE$local_name)
missing<-setdiff(BMDE_col, names(PSSS))
PSSS[missing]<-""


#filter to BMDE columns
PSSS<-PSSS[BMDE_col]

#PSSS<-PSSS %>% dplyr::select(-position, -zero_ref_point, -weather, -precipitation, #-sea_state, -tide_movement, -visibility_distance, -poor_visibility_reason, -poor_visibility_reason_other, -equipment, -walker_count, -dog_count, -power_boat_count, -unpowered_boat_count, -other_activities_name, -other_activities_count, -SurveyComments, -horizon_obscured, -bird_name, -survey_id, -SiteComments, -is_verifier, -eye_height, -arm_length, -binocular_magnification, -scope_magnification, -lat, -long, -species_id, -no_measurements, -FieldNotes)

write.csv(PSSS, "Data/PSSS_BMDE.csv", row.names = FALSE, quote=FALSE)


