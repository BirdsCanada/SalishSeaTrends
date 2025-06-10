#Analysis scripts

if(species.list == "all"){
  sp.list<-unique(sp.data$CommonName)
}else{
  sp.list<-species.list
}

area<- toupper(area)

if(area=="BCCWS"){
  sp.dat<-sp.data %>% filter(ProjectCode=="BCCWS")
  event<-events %>% filter(ProjectCode=="BCCWS")
}  

if(area=="PSSS"){
  sp.dat<-sp.data %>% filter(ProjectCode=="PSSS")
  event<-events %>% filter(ProjectCode=="PSSS")
  
  #This is fixed in the data cleaning script and can be removed during next round. 
  event<-event %>% filter(DurationInHours<=3)
} 

#Specify the spatial extent of the analysis
if(area=="SALISHSEA"){
   sp.dat<-sp.data 
  event<-events
}


##remove COVID data 
sp.dat<-sp.dat %>% filter(wyear != 2020)
event<-event %>% filter(wyear != 2020)

guild <- tolower(guild)

#If guild is set to "Yes" this will override the species list above.   
if(guild=="yes"){
  type <- tolower(type)
  if(type == "migration"){
  sp.list<-unique(sp.data$Migration)
  colnames(sp.dat)[colnames(sp.dat) == "Migration"] <- "Guild"
    }
  if(type=="diet"){
  sp.list<-unique(sp.data$Diet)
  colnames(sp.dat)[colnames(sp.dat) == "Diet"] <- "Guild"
  }
  if(type=="family"){
  sp.list<-unique(sp.data$family_name)  
  colnames(sp.dat)[colnames(sp.dat) == "family_name"] <- "Guild"
  }
}  

#Create a loop for the species list
  for(i in 1:length(sp.list)){
  # for(i in 30:length(sp.list)){ #restarting loop if it breaks
   dat<-NULL 
    #i<-1 #for testing
   
   print(paste("Currently analyzing species ", i, "/", sp.list[i], sep = "")) 
   
    if(guild =="yes"){
    dat <- sp.dat %>% filter(Guild==sp.list[i])
    dat<-dat %>% distinct(ProjectCode, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected, .keep_all = TRUE)
    sp.code<-sp.list[i]
    dat$SpeciesCode<-sp.list[i]
    species_name<- type
    species_sci_name<- " "
    sp.id<- " "
 
    }else{
    
    #Subset the data for the species
    dat <- sp.dat %>% filter(CommonName==sp.list[i])
    dat<-dat %>% distinct(ProjectCode, SpeciesCode, CommonName, species_id, SurveyAreaIdentifier, wyear, YearCollected, wmonth, MonthCollected, DayCollected, .keep_all = TRUE)
   
    sp.code <- unique(dat$SpeciesCode)
    if (length(sp.code) == 0) {
      sp.code <- unique(dat$CommonName) #for group species
    }
    
    if(length(sp.code)>1){
      sp.code<-unique(dat$CommonName)
    }
    
    sp.id<- unique(dat$species_id)
    if (length(sp.id) == 0) {
      sp.id <- unique(dat$CommonName) #for group species
    }
    
    species_name<-unique(dat$CommonName)
    
    species_sci_name <- unique(dat$scientific_name) #for group species
    if (length(species_sci_name) == 0) {
      species_sci_name <- unique(dat$CommonName) #for group species
    }
  }
    
##zero-fill the dat using the events dataframe##
    dat<-left_join(event, dat, by= c("ProjectCode", "SurveyAreaIdentifier", "wyear", "wmonth", "YearCollected", "MonthCollected", "DayCollected"))
#Observation Counts will be backfilled with a 0 whenever it is NA
    dat$ObservationCount[is.na(dat$ObservationCount)]<-0
    dat$CommonName[is.na(dat$CommonName)]<-species_name
    dat$SpeciesCode[is.na(dat$SpeciesCode)]<-sp.code
  
#remove extreme outliers 
    outlier<-(quantile(dat$ObservationCount, probs = c(0.99)))*3
    dat<-dat %>% filter(ObservationCount<outlier)

if(nrow(dat)>0){        
    
  # Count the number of unique years each site has data (regardless of detection)
  site_years_sampled <- aggregate(wyear ~ SurveyAreaIdentifier, data = dat, FUN = function(x) length(unique(x)))
  
  # Filter for sites with at least 10 years of data
  sites_to_keep <- subset(site_years_sampled, wyear >= 10)$SurveyAreaIdentifier
  
  # Filter the main dataset
  dat <- subset(dat, SurveyAreaIdentifier %in% sites_to_keep)  

#Remove SurveyAreaIdentifier from the data on where the sum of ObersevationCount is 0 across all years
#If a species was never detected on a route, we will not include that route in the species specific analysis
#This is considered out of range or in unsuitable habitat
    dat<-dat %>% group_by(SurveyAreaIdentifier) %>% filter(sum(ObservationCount)>0) %>% ungroup()
    routes<-n_distinct(dat$SurveyAreaIdentifier)
    
 #Minimum Data Requirements##
    
    #Now we will check that the minimum data requirements are met. 
    #We will want to ensure that the species was detected a minimum of X times in a year
    #That the species was detected in at least 1/2 of the survey years
    #And that greater than a certain % of sites have non-zero counts
    
    SpeciesMean<- dat %>% group_by(wyear) %>% summarize(YearMean = sum(ObservationCount)) %>% ungroup() %>% summarize(MeanMean = mean(YearMean))
    SpeciesMean$NumYears <- n_distinct(dat$wyear)
    
    #Now cheek the SpeciesMean object to see if the species meets the minimum data requirements 
    #all the variable must be large than the values min.abundance, min.years, zero.count, if TRUE continue with the analysis
    
    if(SpeciesMean$MeanMean>=min.abundance & SpeciesMean$NumYears>=min.years & routes>nsites){
  min.data <- TRUE 
      }else{
  min.data <- FALSE
    }

    print(paste("Did", sp.list[i], "meet minimum data requirements:", min.data))
    
    #only continue if the species meets the minimum data requirements      
if(min.data==TRUE){

###Model without spatial effect on abundance 
 
  wyears <- unique(dat$wyear)
  mean_wyear <- max(wyears)

#Create index variables
    dat <- dat %>% mutate( 
      std_yr = wyear - Y2,
      protocol = factor(ProjectCode), 
      kappa = as.integer(factor(dat$SurveyAreaIdentifier)),
      year_idx = as.integer(wyear - mean_wyear), #intercept is the expected count during the most recent year of data collection. 
      wmonth_idx = as.factor(wmonth), 
      sp_idx = as.integer(factor(CommonName)))%>%
      st_as_sf(coords = c("DecimalLongitude", "DecimalLatitude"), crs = 4326, remove = FALSE) %>% 
      st_transform(epsg6703km) %>%
        mutate(
          easting = st_coordinates(.)[, 1],
          northing = st_coordinates(.)[, 2]) %>%
        arrange(SurveyAreaIdentifier, wyear)
    
    # Aggregate data to annual max count
    dat <- dat %>%
      group_by(SurveyAreaIdentifier, wyear) %>%
      slice_max(ObservationCount, n = 1, with_ties = FALSE) %>%
      ungroup()

  
     #Plot the data for visual inspection

     #Convert the data to a spatial object
     dat_sf <- st_as_sf(dat, coords = c("DecimalLongitude", "DecimalLatitude"), crs = 4326)
     #remove zero ObservationCounts from the output plot
     dat_sf<-dat_sf %>% group_by(SurveyAreaIdentifier, geometry) %>% summarise(ObservationCount = mean(ObservationCount))
     map<-st_transform(map, crs = 4326)

     #write plot to Plots folder in Output
       y<- ggplot()+
       geom_sf(data = map, fill = "white", color = "black") +
       # Plot the points, size by the sum of ObservationCount within a year
       geom_sf(data=dat_sf, aes(size=ObservationCount)) +
       # Add a theme with a minimal design and change the font styles, to your preference
       theme_minimal() +
       #theme(legend.position = "bottom") +
       # To make the points in the legend larger without affecting map points
       guides(color = guide_legend(override.aes = list(size = 3))) +
       #make the text on the x-axis vertical
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       # Define the title and axis names
       labs(title = sp.list[i], x = "Longitude", y = "Latitude")+
       # Change legend lable
       scale_size_continuous(name = "Mean Observation Count")

       ggsave(
         filename = file.path(plot.dir, paste0(name, sp.list[i], "_Map_SPDE.jpeg")),
         plot = y,
         width = 8,
         height = 6,
         dpi = 150
       )

#create the datframe of covariates
N<-nrow(dat)

    if(guild=="yes"){
      
      dat$sp_idx[is.na(dat$sp_idx)] <- 999 #replace NA which are zero counts with generic sp_idx
      
      Covariates<- data.frame(
        Intercept=rep(1, N),
        kappa = dat$kappa, 
        sp_idx = as.integer(dat$sp_idx),
        year_idx=dat$year_idx
        )
      
    }else{

    Covariates<- data.frame(
      Intercept=rep(1, N),
       kappa = dat$kappa, 
       year_idx=dat$year_idx 
       )
    }

    #Create the Mesh
    #Make a set of distinct study sites for mapping kappa
    site_map <- dat %>%
      dplyr::select(SurveyAreaIdentifier, easting, northing) %>% distinct()
    
    #make unique locations
    Loc_unique<-dat %>% dplyr::select(SurveyAreaIdentifier, easting, northing) %>%
      distinct() %>%
      dplyr::select(easting, northing) %>% distinct() %>%
      st_drop_geometry() %>%
      as.matrix()

    sample_size<-nrow(Loc_unique)
    
    #make the mesh this way so that the point fall on the vertices of the lattice
    Loc_all<-dat %>% dplyr::select(easting, northing) %>% st_drop_geometry() %>% as.matrix()
    
     #change the crs of map to epsg6703km
    map <- st_transform(map, crs = epsg6703km)
    
    # Convert sf polygon to SpatialPolygons (sp format)
    map_sp <- as(map, "Spatial")
    
    # Convert to inla.mesh.segment
    boundary_segment <- inla.sp2segment(map_sp)
    
    mesh2<- inla.mesh.2d(Loc_unique,
      boundary = boundary_segment,  # Your polygon boundary
      max.edge = c(25, 50),       # Inner/outer resolution (in CRS units)
      cutoff = 2,                  # Minimum distance between vertices
      crs = fm_crs(dat)             # Use the same CRS as your data
    )
    
    spde <- inla.spde2.pcmatern(
      mesh = mesh2,
      prior.range = prior.range,
      prior.sigma = prior.sigma
    )
    
    #Spatial Fields
    # make index sets for the spatial model
    alpha_idx <- inla.spde.make.index(name = "alpha", n.spde = spde$n.spde) #n.repl is used to account for repeated measure at the same location over time.
    
    #Projection matrix A using all locations
    A_alph <- inla.spde.make.A(mesh=mesh2, loc=Loc_all)  # pg 218 this formula should include an alpha, default is 2 if the model does not include times.
   
    #Create Stack Object for INLA
    Stack <- inla.stack (
      tag="FitGAM",
      data=list(count=as.vector(dat$ObservationCount)), #response from the dataframe
      effects = list(Covariates=Covariates, alpha = alpha_idx),  #covariate list #removed tau=tau_idx
      A = list(
        1, #value of 1 to non-spatial terms
        A_alph
       )
    )
    
    if(guild=="yes"){
      
      formula.sp<- count ~ -1 + Intercept + 
      f(year_idx, model = "iid", hyper = hyper.iid) +  
      f(kappa, model="iid", hyper=hyper.iid) + #site effect replaced with spatail 
      f(sp_idx, model="iid", hyper=hyper.iid) + 
      f(alpha, model =spde)

    
      }else{
    
     formula.sp<- count ~ -1+ Intercept + 
      f(year_idx, model = "iid", hyper = hyper.iid) + 
      f(kappa, model="iid", hyper=hyper.iid) + #site effect replaced with spatial
      f(alpha, model =spde) 
    }   


    #fit the spatial model using INLA
    M1 <- inla(
      formula.sp, 
      family = fam, 
      data = inla.stack.data(Stack), 
      offset = log(dat$DurationInHours),
      control.predictor = list(
        A = inla.stack.A(Stack),  # Projector matrix
        compute = TRUE,            # Store latent predictor
        link = 1                   # Use the response-family link (log for NB)
      ),
      control.compute = list(
        dic = TRUE,                # For model comparison
        waic = TRUE,               # For model comparison
        config = TRUE,             # Enable posterior sampling
        cpo = TRUE                 # Optional: cross-validated PIT
      ),
      control.fixed = list(
        mean = 0,                  # Prior mean for fixed effects
        prec = 0.001               # Prior precision for fixed effects
      ), 
      control.family = list(
        hyper = list(theta = list(prior = "loggamma", param = c(3, 0.1)))
      )
    )
    
    # Calculate disperson statistic
    
    disp <- compute_dispersion_SPDE(M1 = M1, Stack = Stack, fam = fam)
    print(paste(sp.list[i], " Dispersions Statistic = ", disp, sep = ""))
    
    # Append to dispersion file
 
    dispersion_entry <- data.frame(
      area_code = area,
      SpeciesCode = sp.list[i],
      dispersion = disp
    )
    
    write.table(dispersion_entry,
                file = paste0(out.dir, name, "_DispersionStat_SPDE.csv"),
                append = TRUE,
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
  
    
    # Get estimation data indices (non-NA counts)
    stack_data <- inla.stack.data(Stack)
    est_indices <- which(!is.na(stack_data$count))
    
    # Observed values
    observed <- stack_data$count[est_indices]
    
    # Fitted values (posterior means)
    fitted <- M1$summary.fitted.values$mean[est_indices]
    
    df <- data.frame(Observed = observed, Fitted = fitted)
    
    p <- ggplot(df, aes(x = Observed, y = Fitted)) +
      geom_point(color = "dodgerblue") +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
      labs(
        x = "Observed Count",
        y = "Fitted Value",
        title = paste("Observed vs Fitted Values (SPDE Model)", sp.list[i])
      ) +
      theme_minimal(base_size = 14)
    
    print(p)
    
    ggsave(
      filename = file.path(plot.dir, paste0(name, sp.list[i], "_FitPlot_SPDE.jpeg")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 150
    )
    
    
#Calculate Posterior estimate of abundance
nsamples<- 100
post.sample1 <-NULL #clear previous
post.sample1<-inla.posterior.sample(nsamples, M1)

tmp0<-NULL
calculate_indices <- function(sample) {
  effects <- exp(sample$latent[1:nrow(dat)])  # Direct latent field access
  aggregate(effects ~ SurveyAreaIdentifier + wyear, data = dat, FUN = mean)
}

tmp0 <- bind_rows(
  lapply(post.sample1, calculate_indices), 
  .id = "sample"
)

tmp1 <- tmp0 %>%
  group_by(SurveyAreaIdentifier, wyear) %>%
  summarise(
    index = mean(effects),
    stdev = sd(effects),
    stderr = stdev / sqrt(n()),  # Uses group-wise sample count
    lower_ci = quantile(effects, 0.025),
    upper_ci = quantile(effects, 0.975),
    .groups = 'drop'
  ) %>% 
  mutate(area_code=SurveyAreaIdentifier)


#Assign data to output table 

  indices.csv<-tmp1 %>% dplyr::select(area_code, wyear, index, lower_ci, upper_ci, stdev, stderr) %>%
    ungroup() %>% mutate(
    species_code = sp.code,
    years = paste(min(dat$wyear), "-", max(dat$wyear), sep = ""),
    year = wyear,
    period ="all years",
    season = "winter",
    model_type = "ALPHA SPATIAL SPDE",
    species_id=sp.id,
    species_name=species_name,
    species_sci_name=species_sci_name,
    error="",
    #Assing missing data fields 
    upload_id="",
    trend_id="",
    smooth_upper_ci="",
    smooth_lower_ci="",
    upload_dt="",
    family=fam,
    results_code = "BCCWS/PSSS",
    version = Sys.Date(),
    season="Winter",
    trend_index="") #look at CMMN code to generate in next round of analysis
  
  # Run LOESS 
  indices.csv <- indices.csv %>%
    group_by(area_code) %>%
    arrange(year, .by_group = TRUE) %>%
    tidyr::drop_na(index) %>%  # Remove NA/NaN/Inf within groups
    mutate(
      LOESS_index = if (n() >= 10) {  # Require ≥3 data points for LOESS
        predict(
          loess(index ~ year, span = 0.55, na.action = na.exclude),
          newdata = data.frame(year = year)
        )
      } else {
        NA_real_
      }
    ) %>%
    ungroup()
  
  # Order output before printing to table
  indices.csv<-indices.csv %>% dplyr::select(results_code, version, area_code, season, period, species_code, species_id, year, index, stderr, stdev, upper_ci, lower_ci, LOESS_index, trend_index)
  
  # Write data to table
  write.table(indices.csv, 
              file = paste(out.dir,	name, "_AnnualIndices_SPDE.csv", sep = ""),
              row.names = FALSE, 
              append = TRUE, 
              quote = FALSE, 
              sep = ",", 
              col.names = FALSE)
  
  
  ##############################################################################
  ##END POINT TRENDS ##
  
  # Calculate endpoint trends per region and posterior sample
  trend_samples <- tmp0 %>%
    group_by(SurveyAreaIdentifier, sample) %>%
    summarise(
      n1 = effects[wyear == min(wyear)],
      nT = effects[wyear == max(wyear)],
      y1 = min(wyear),
      yT = max(wyear),
      years = yT - y1
    ) %>%
    mutate(
      # Geometric mean annual percent change
      annual_trend = (nT / n1)^(1 / years) - 1,
      percent_annual_change = 100 * annual_trend,
      # Total percent change over the period
      percent_change_total = 100 * (nT / n1 - 1)
    )
  
  # Summarise across posterior samples for credible intervals and uncertainty
  trend_summary <- trend_samples %>%
    group_by(SurveyAreaIdentifier) %>%
    summarise(
      trnd = mean(percent_annual_change),
      stdev = sd(percent_annual_change),
      stderr = stdev / sqrt(n()),
      lower_ci = quantile(percent_annual_change, 0.025),
      upper_ci = quantile(percent_annual_change, 0.975),
      percent_change = mean(percent_change_total)) %>% 
    mutate(
      index_type = "Endpoint Trend",
      Width_of_Credible_Interval = upper_ci - lower_ci,
      precision_cat = case_when(
        Width_of_Credible_Interval < 3.5 ~ "High",
        between(Width_of_Credible_Interval, 3.5, 6.7) ~ "Medium",
        TRUE ~ "Low"
      )) 
  
  #write output to table
  trend.out<-NULL
  trend.out <- trend_summary %>%
    mutate(model_type="ALPHA SPATIAL SPDE", 
           model_family = fam,
           year_start=Y1, 
           year_end=Y2,
           years = paste(Y1, "-", Y2, sep = ""),
           period ="all years",
           season = "winter",
           results_code = "BCCWS/PSSS",
           area_code = SurveyAreaIdentifier,
           version=2025, 
           species_code = sp.code,
           species_id=sp.id, 
           species_name=species_name,
           species_sci_name=species_sci_name,
           stderr = "",
           index_type= "Endpoint trend", 
           model_fit = "", 	
           percent_change_low ="", 
           percent_change_high = "",
           prob_decrease_0 = "",
           prob_decrease_25 = "",
           prob_decrease_30 = "",
           prob_decrease_50 = "",
           prob_increase_0 = "",
           prob_increase_33 = "",	
           prob_increase_100 = "",
           confidence = "",
           precision_num = "",
           suitability="",
           coverage_num = "",
           coverage_cat = "",
           goal = "",
           goal_lower = "",
           sample_size = sample_size,
           sample_size_units="Number of Sites",
           sample_total = "",
           subtitle = "",
           pval = "",
           pval_str = "",
           post_prob = "",
           trnd_order = "",
           dq = "",
           prob_LD = "",
           prob_MD = "",
           prob_LC = "",
           prob_MI = "",
           prob_LI = "",
           quantile_050 = "",
           quantile_165 = "",
           quantile_835 = "",
           quantile_950 = "",
           trend_id = "",
           upload_dt = "")
  
  write.trend<-trend.out %>% dplyr::select(results_code,	version,	area_code,	season,	period, species_code,	species_id,	years,year_start,	year_end,	trnd,	lower_ci, upper_ci, index_type, stderr,	model_type,	model_fit,	percent_change,	percent_change_low,	percent_change_high,	prob_decrease_0,	prob_decrease_25,	prob_decrease_30,	prob_decrease_50,	prob_increase_0,	prob_increase_33,	prob_increase_100, suitability, precision_num,	precision_cat,	coverage_num,	coverage_cat,	sample_size, sample_size_units, prob_LD, prob_MD, prob_LC, prob_MI, prob_LI)
  
  write.table(write.trend, 
              file = paste(out.dir, name, "_TrendsEndpoint_SPDE.csv", sep = ""), 
              row.names = FALSE, 
              append = TRUE, 
              quote = FALSE, 
              sep = ",", 
              col.names = FALSE)  

  ##############################################################################
  ##SLOPE TRENDS ##
  
  calc_slope <- function(df) {
    mod <- lm(log(effects) ~ wyear, data = df)
    slope <- coef(mod)[2]
    years <- max(df$wyear) - min(df$wyear)
    percent_annual_change <- 100 * (exp(slope) - 1)
    total_percent_change <- 100 * (exp(slope * years) - 1)
    tibble(
      slope = slope,
      percent_annual_change = percent_annual_change,
      total_percent_change = total_percent_change
    )
  }
  
  # Apply to each region and posterior sample
  slope_trends <- tmp0 %>%
    group_by(SurveyAreaIdentifier, sample) %>%
    group_modify(~ calc_slope(.x)) %>%
    ungroup()
  
  # Summarise across posterior samples for each region
  slope_summary <- slope_trends %>%
    group_by(SurveyAreaIdentifier) %>%
    summarise(
      trnd = mean(percent_annual_change),
      stdev = sd(percent_annual_change),
      stderr = stdev / sqrt(n()),
      lower_ci = quantile(percent_annual_change, 0.025),
      upper_ci = quantile(percent_annual_change, 0.975),
      percent_change = mean(total_percent_change),
    ) %>% 
    mutate(
      index_type = "Slope Trend",
      Width_of_Credible_Interval = upper_ci - lower_ci,
      precision_cat = case_when(
        Width_of_Credible_Interval < 3.5 ~ "High",
        between(Width_of_Credible_Interval, 3.5, 6.7) ~ "Medium",
        TRUE ~ "Low"
      ))
  
  #write output to table
  trend.out<-NULL
  trend.out <- slope_summary %>%
    mutate(model_type="ALPHA SPATIAL SPDE", 
           model_family = fam,
           years = paste(Y1, "-", Y2, sep = ""),
           year_start=Y1, 
           year_end = Y2,
           period ="all years",
           season = "winter",
           results_code = "BCCWS/PSSS",
           area_code = SurveyAreaIdentifier,
           version=2025, 
           species_code = sp.code,
           species_id=sp.id, 
           species_name=species_name,
           species_sci_name=species_sci_name,
           stderr = "",
           model_fit = "", 	
           percent_change_low ="", 
           percent_change_high = "",
           prob_decrease_0 = "",
           prob_decrease_25 = "",
           prob_decrease_30 = "",
           prob_decrease_50 = "",
           prob_increase_0 = "",
           prob_increase_33 = "",	
           prob_increase_100 = "",
           confidence = "",
           precision_num = "",
           suitability="",
           coverage_num = "",
           coverage_cat = "",
           goal = "",
           goal_lower = "",
           sample_size = sample_size,
           sample_size_units="Number of Sites",
           sample_total = "",
           subtitle = "",
           pval = "",
           pval_str = "",
           post_prob = "",
           trnd_order = "",
           dq = "",
           prob_LD = "",
           prob_MD = "",
           prob_LC = "",
           prob_MI = "",
           prob_LI = "",
           quantile_050 = "",
           quantile_165 = "",
           quantile_835 = "",
           quantile_950 = "",
           trend_id = "",
           upload_dt = "")
  
  write.trend<-trend.out %>% dplyr::select(results_code,	version,	area_code,	season,	period, species_code,	species_id,	years,year_start,	year_end,	trnd,	lower_ci, upper_ci, index_type, stderr,	model_type,	model_fit,	percent_change,	percent_change_low,	percent_change_high,	prob_decrease_0,	prob_decrease_25,	prob_decrease_30,	prob_decrease_50,	prob_increase_0,	prob_increase_33,	prob_increase_100, suitability, precision_num,	precision_cat,	coverage_num,	coverage_cat,	sample_size, sample_size_units, prob_LD, prob_MD, prob_LC, prob_MI, prob_LI)
  
  write.table(write.trend, 
              file = paste(out.dir, name, "_TrendsSlope_SPDE.csv", sep = ""), 
              row.names = FALSE, 
              append = TRUE, 
              quote = FALSE, 
              sep = ",", 
              col.names = FALSE)  
  
  
  
  } #end nrow data 
  } #end min.data  
  } #end SpeciesLoop


