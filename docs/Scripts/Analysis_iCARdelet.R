#Analysis iCAR

if(species.list == "all"){
  sp.list<-unique(sp.data$CommonName)
}else{
  sp.list<-species.list
}

  sp.dat<-sp.data 
  event<-events
  
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

#grid data are only those cells containing data. This layers is created in ArcGIS. 
nb1 <- spdep::poly2nb(poly, row.names=poly$data, queen=TRUE); nb1
is.symmetric.nb(nb1, verbose = FALSE, force = TRUE)
nb2INLA("nb1.graph", nb1)
nb1.adj <- paste(getwd(),"/nb1.graph", sep="")
g1 <- inla.read.graph("nb1.graph")


#Create a loop for the species list
for(i in 1:length(sp.list)){
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
    
    #Only retains routes on which a species was detected >1 years   
    # Count the number of years each site had detection
    site_years_detected <- aggregate(ObservationCount ~ SurveyAreaIdentifier + wyear, data = dat, FUN = sum)
    site_years_detected$Detected <- site_years_detected$ObservationCount > 0
    
    # Sum the number of years with detection for each site
    site_detect_years <- aggregate(Detected ~ SurveyAreaIdentifier, data = site_years_detected, FUN = sum)
    
    # Filter for sites with detection in more than 1 years
    sites_to_keep <- subset(site_detect_years, Detected > 1)$SurveyAreaIdentifier
    
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
      
      #Prepare the parameters
      wyears <- unique(dat$wyear)
      mean_wyear <- max(wyears)
      
      #Create index variables
      dat <- dat %>% mutate( 
        std_yr = wyear - Y2,
        alpha_i = as.integer(factor(dat$Name)),
        protocol = as.factor(dat$ProjectCode), 
        kappa = as.integer(factor(dat$SurveyAreaIdentifier)),
        year_idx = as.integer(wyear - mean_wyear), #intercept is the expected count during the most recent year of data collection. 
        wmonth_idx = as.factor(wmonth), 
        sp_idx = as.integer(factor(CommonName)))%>%
        arrange(SurveyAreaIdentifier, wyear)
      
      # Aggregate data to annual max count
      dat <- dat %>%
        group_by(SurveyAreaIdentifier, wyear) %>%
        slice_max(ObservationCount, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      # #Last, remove alpha_i with incomplete sampling over all years. 
      # #Indicates that the polygon is not well sampled for a given species and end point trend may not work
      # period_num = Y2-Y1 
      # 
      #  complete_sites <- dat %>%
      #   group_by(alpha_i) %>%
      #   summarise(has_all_years = n_distinct(wyear)) %>%
      #   filter(has_all_years == period_num) %>%
      #   select(alpha_i)
      # 
      # Filter original data frame to keep only complete sites
      # dat <- dat %>%
      #   filter(alpha_i %in% complete_sites$alpha_i)
      
      #Model Formula
      if(guild=="yes"){
        
        dat$sp_idx[is.na(dat$sp_idx)] <- 999 #replace NA which are zero counts with generic sp_idx
        
        formula<- ObservationCount ~ -1 + 
          f(year_idx, model = "iid", hyper = hyper.iid) + 
          f(sp_idx, model="iid", hyper=hyper.iid) +
          # cell ICAR random intercepts
          f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE, hyper = hyper.iid) +
          # random route intercepts
          f(kappa, model = "iid", hyper = hyper.iid)
       
      }else{
        
        formula<- ObservationCount ~ -1 + 
          f(year_idx, model = "iid", hyper = hyper.iid) + 
          # cell ICAR random intercepts
          f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE, hyper = hyper.iid) +
          # random route intercepts
          f(kappa, model="iid", hyper = hyper.iid)
        
      }
      
      M2<-try(inla(formula, family = fam, data = dat, offset = log(dat$DurationInHours),
                   control.predictor = list(compute = TRUE),
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
      ))
  
      if (inherits(M2, "try-error")) {
        cat("The iCAR model failed to run. Please check your input or model specifications.\n")
      } else {
        # Continue with normal processing
        print("Model iCAR model ran successfully.")
        
     #Dispersion Statistic
       Dispersion1 <- calculate_dispersion_iCAR(M2, dat$ObservationCount)
       print(paste(sp.list[i], " Dispersions Statistic = ", Dispersion1, sep = ""))
      
       # Append to dispersion file
       
       dispersion_entry <- data.frame(
         area_code = name,
         SpeciesCode = sp.list[i],
         dispersion = Dispersion1
       )
       
       write.table(dispersion_entry,
                   file = paste0(out.dir, name, "_DispersionStat_iCAR.csv"),
                   append = TRUE,
                   sep = ",",
                   row.names = FALSE,
                   col.names = FALSE)
      
      mu1 <- M2$summary.fitted.values$mean  
      dat$mu1<-mu1
      
      df <- data.frame(Observed = dat$ObservationCount, Fitted = dat$mu1)
      
      d <- ggplot(df, aes(x = Observed, y = Fitted)) +
        geom_point(color = "dodgerblue") +
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
        labs(
          x = "Observed Count",
          y = "Fitted Value",
          title = paste("Observed vs Fitted Values (iCAR Model)", sp.list[i])
        ) +
        theme_minimal(base_size = 14)
      
      ggsave(
        filename = file.path(plot.dir, paste0(name, sp.list[i], "_FitPlot_iCAR.jpeg")),
        plot = d,
        width = 8,
        height = 6,
        dpi = 150
      )
      
      print(d)
  
      ##Remove polygons with no survey sites
      cells_with_counts <- unique(dat$alpha_i[which(!is.na(dat$ObservationCount))]) 
     
      # Filter dat to only those alpha_i, then calculate sample size per polygon
      sample_size <- dat %>%
        filter(alpha_i %in% cells_with_counts) %>%
        group_by(alpha_i) %>%
        summarise(sample_size = n_distinct(SurveyAreaIdentifier), .groups = "drop")
      
      #Calculate Posterior estimate of abundance
      nsamples<- 1000
      post.sample1 <-NULL #clear previous
      post.sample1<-inla.posterior.sample(nsamples, M2)
      
      tmp0<-NULL
              calculate_indices <- function(sample) {
              effects <- exp(sample$latent[1:nrow(dat)])  # Direct latent field access
              aggregate(effects ~ alpha_i + wyear, data = dat, FUN = mean)
            }
           
             tmp0 <- bind_rows(
              lapply(post.sample1, calculate_indices), 
              .id = "sample"
            )
        
             tmp1 <- tmp0 %>%
               group_by(alpha_i, wyear) %>%
               summarise(
                 index = mean(effects),
                 stdev = sd(effects),
                 stderr = stdev / sqrt(n()),  # Uses group-wise sample count
                 lower_ci = quantile(effects, 0.025),
                 upper_ci = quantile(effects, 0.975),
                 .groups = 'drop'
               )
        

      #link back this the site name using the grid_key
      grid3<-grid %>% st_drop_geometry() %>% 
        select(alpha_i, Name, Area) %>% 
        mutate(area_code = Name) %>% 
        distinct()
      
      tmp1<-left_join(tmp1, grid3, by="alpha_i")
      
      #Assign data to output table 
      indices.csv<-tmp1 %>% dplyr::select(wyear, index, lower_ci, upper_ci, stdev, stderr, area_code) %>% mutate(
        species_code = sp.code,
        years = paste(min(dat$wyear), "-", max(dat$wyear), sep = ""),
        year = wyear,
        period ="all years",
        season = "winter",
        model_type = "iCAR ALPHA SPATIAL",
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
        season="Winter"
        #trend_index=""
        )
      
      #LOESS
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
      
      #will print after computing end point trends. 
      
      ##############################################################################
      ##END POINT TRENDS ##
      
      # Calculate endpoint trends per region and posterior sample
      trend_samples <- tmp0 %>%
        group_by(alpha_i, sample) %>%
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
        group_by(alpha_i) %>%
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
          )) %>% 
        left_join(grid3, by = "alpha_i")
      
      
      #write output to table   
      trend.out<-NULL
      trend.out <- trend_summary %>%
        mutate(model_type="iCAR ALPHA SPATIAL", 
               model_family = fam,
               years = paste(Y1, "-", Y2, sep = ""),
               year_start=Y1, 
               year_end=Y2,
               period ="all years",
               season = "winter",
               results_code = "BCCWS/PSSS",
               version=Sys.Date(), 
               area_code=trend_summary$area_code,
               species_code = sp.code,
               species_id=sp.id, 
               species_name=species_name,
               species_sci_name=species_sci_name,
               index_type= "Endpoint Trend", 
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
      
      trend.out<-left_join(trend.out, sample_size, by="alpha_i")
      
      write.trend<-trend.out %>% dplyr::select(results_code,	version,	area_code,	season,	period, species_code,	species_id,	years,year_start,	year_end,	trnd,	lower_ci, upper_ci, index_type, stderr,	model_type,	model_fit,	percent_change,	percent_change_low,	percent_change_high,	prob_decrease_0,	prob_decrease_25,	prob_decrease_30,	prob_decrease_50,	prob_increase_0,	prob_increase_33,	prob_increase_100, suitability, precision_num,	precision_cat,	coverage_num,	coverage_cat,	sample_size, sample_size_units, prob_LD, prob_MD, prob_LC, prob_MI, prob_LI)
      
      write.table(write.trend, 
                  file = paste(out.dir, name, "_TrendsEndpoint_iCAR.csv", sep = ""), 
                  row.names = FALSE, 
                  append = TRUE, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = FALSE)  
      
      ##############################################################################
      ##SLOPE TRENDS ##
      
      calc_slope <- function(df) {
        mod <- lm(log(effects) ~ wyear, data = df)
        intercept <- coef(mod)[1]
        slope <- coef(mod)[2]
        years <- max(df$wyear) - min(df$wyear)
        percent_annual_change <- 100 * (exp(slope) - 1)
        total_percent_change <- 100 * (exp(slope * years) - 1)
        tibble(
          intercept = intercept,
          slope = slope,
          percent_annual_change = percent_annual_change,
          total_percent_change = total_percent_change
        )
      }
      
      #Apply to each region and poterior sample
      slope_trends<-tmp0 %>% 
        group_by(alpha_i, sample) %>% 
        group_modify(~calc_slope(.x)) %>% 
        ungroup()
      
      
     # Summarise across posterior samples for each region
       slope_summary <- slope_trends %>%
         group_by(alpha_i) %>%
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
                )) %>% left_join(grid3, by = "alpha_i")
         
       #write output to table
      trend.out<-NULL
      trend.out <- slope_summary %>%
        mutate(model_type="ALPHA SPATIAL iCAR", 
               model_family = fam,
               years = paste(Y1, "-", Y2, sep = ""),
               year_start=Y1, 
               year_end = Y2,
               period ="all years",
               season = "winter",
               results_code = "BCCWS/PSSS",
               version= Sys.Date(), 
               species_code = sp.code,
               species_id=sp.id, 
               species_name=species_name,
               species_sci_name=species_sci_name,
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
      
      trend.out<-left_join(trend.out, sample_size, by="alpha_i")
      
      write.trend<-trend.out %>% dplyr::select(results_code,	version,	area_code,	season,	period, species_code,	species_id,	years,year_start,	year_end,	trnd,	lower_ci, upper_ci, index_type, stderr,	model_type,	model_fit,	percent_change,	percent_change_low,	percent_change_high,	prob_decrease_0,	prob_decrease_25,	prob_decrease_30,	prob_decrease_50,	prob_increase_0,	prob_increase_33,	prob_increase_100, suitability, precision_num,	precision_cat,	coverage_num,	coverage_cat,	sample_size, sample_size_units, prob_LD, prob_MD, prob_LC, prob_MI, prob_LI)
      
      write.table(write.trend, 
                  file = paste(out.dir, name, "_TrendsSlope_iCAR.csv", sep = ""), 
                  row.names = FALSE, 
                  append = TRUE, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = FALSE)  
     
     
       ############################################################################
      ##Trend Slope estimates to index file
     
       trendlines <- slope_trends %>%
        left_join(tmp0 %>% select(alpha_i, wyear) %>% distinct(), by = "alpha_i", relationship = "many-to-many") %>%
        group_by(alpha_i, sample) %>%
        mutate(fitted = exp(intercept + slope * wyear)) %>%
        ungroup()
      
       trend_summary <- trendlines %>%
         group_by(alpha_i, wyear) %>%
         summarise(
           trend_index = mean(fitted)
         ) %>%
         ungroup() %>%
         left_join(grid3, by="alpha_i") %>% 
         select(trend_index, area_code, wyear)
        
       
        full_area_index <- trendlines %>%
         group_by(wyear) %>%
         summarise(
           trend_index = mean(fitted)
         ) %>%
         mutate(area_code = "Full Study Area")
       
       trend_summary <- bind_rows(trend_summary, full_area_index)
       names(trend_summary)[names(trend_summary) == "wyear"] <- "year"
       
       indices.csv<-indices.csv %>% left_join(trend_summary, by=c("year", "area_code"))
       
       # Order output before printing to table
       indices.csv<-indices.csv %>% ungroup %>% dplyr::select(results_code, version, area_code, season, period, species_code, species_id, year, index, stderr, stdev, upper_ci, lower_ci, LOESS_index, trend_index)
       
       # Write data to table
       write.table(indices.csv, 
                   file = paste(out.dir,	name, "_AnnualIndices_iCAR.csv", sep = ""),
                   row.names = FALSE, 
                   append = TRUE, 
                   quote = FALSE, 
                   sep = ",", 
                   col.names = FALSE)
       
      
      ############################################################################
      ### FULL STUDY AREA INDEX ###
   
      #Join region areas to posterior samples
      grid_unique <- grid %>% select(alpha_i, Area) %>% distinct()
      
      tmp2 <- tmp0 %>% 
        left_join(grid_unique, by = "alpha_i")
      
      #Calculate area-weighted indices per simulation
      area_weighted_indices <- tmp2 %>%
        group_by(wyear, sample) %>%
        summarise(
          # Area-weighted mean calculation
          weighted_index = sum(effects * Area) / sum(Area),
          .groups = "drop"
        )
      
      #Summarize across simulation
      tmp2_area <- area_weighted_indices %>%
        group_by(wyear) %>%
        summarise(
          # Mean instead of median for proper expectation
          index = mean(weighted_index),
          # Quantile-based credible intervals
          lower_ci = quantile(weighted_index, 0.025),
          upper_ci = quantile(weighted_index, 0.975),
          # Optional variance metrics
          stdev = sd(weighted_index),
          stderr = stdev / sqrt(n()),
          .groups = "drop"
        )
      
      #Assign data to output table 
      indices.csv<-NULL 
      indices.csv<-tmp2_area %>% dplyr::select(wyear, index, lower_ci, upper_ci, stdev, stderr) %>% mutate(
        species_code = sp.code,
        index_type= "Endpoint Trend",
        area_code = "Full Study Area",
        years = paste(min(dat$wyear), "-", max(dat$wyear), sep = ""),
        year = wyear,
        period ="all years",
        season = "winter",
        model_type = "iCAR ALPHA SPATIAL",
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
        trend_index="")
      
      #LOESS
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
      indices.csv<-indices.csv %>% ungroup %>% dplyr::select(results_code, version, area_code, season, period, species_code, species_id, year, index, stderr, stdev, upper_ci, lower_ci, LOESS_index, trend_index)
      
      # Write data to table
      write.table(indices.csv, 
                  file = paste(out.dir,	name, "_AnnualIndices_iCAR.csv", sep = ""),
                  row.names = FALSE, 
                  append = TRUE, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = FALSE)
      
      ##############################################################################
      ##END POINT TRENDS FULL AREA##
      
      # Calculate endpoint trends per region and posterior sample
      trend_samples2 <- area_weighted_indices %>%
        group_by(sample) %>%
        summarise(
          n1 = weighted_index[wyear == min(wyear)],
          nT = weighted_index[wyear == max(wyear)],
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
      trend_summary2 <- trend_samples2 %>%
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
    
      sample_size_all<-sum(sample_size$sample_size)
      
      #write output to table
      trend.out<-NULL
      trend.out <- trend_summary2 %>%
        mutate(model_type="iCAR ALPHA SPATIAL", 
               model_family = fam,
               years = paste(Y1, "-", Y2, sep = ""),
               year_start=Y1, 
               year_end = Y2,
               period ="all years",
               season = "winter",
               results_code = "BCCWS/PSSS",
               area_code = "Full Study Area",
               version=Sys.Date(), 
               species_code = sp.code,
               species_id=sp.id, 
               species_name=species_name,
               species_sci_name=species_sci_name,
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
               sample_size = sample_size_all,
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
                  file = paste(out.dir, name, "_TrendsEndpoint_iCAR.csv", sep = ""), 
                  row.names = FALSE, 
                  append = TRUE, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = FALSE)  
      
      ##############################################################################
      ##SLOPE TRENDS FULL AREA##
      
      calc_slope2 <- function(df) {
        mod <- lm(log(weighted_index) ~ wyear, data = df)
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
      slope_trends2 <- area_weighted_indices %>%
        group_by(sample) %>%
        group_modify(~ calc_slope2(.x)) %>%
        ungroup()
      
      # Summarise across posterior samples for each region
      slope_summary2 <- slope_trends2 %>%
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
            TRUE ~ "Low"))
      
      #write output to table
      trend.out<-NULL
      trend.out <- slope_summary2 %>%
        mutate(model_type="ALPHA SPATIAL iCAR", 
               model_family = fam,
               years = paste(Y1, "-", Y2, sep = ""),
               area_code = "Full Study Area",
               year_start=Y1, 
               year_end = Y2,
               period ="all years",
               season = "winter",
               results_code = "BCCWS/PSSS",
               version= Sys.Date(), 
               species_code = sp.code,
               species_id=sp.id, 
               species_name=species_name,
               species_sci_name=species_sci_name,
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
               sample_size = sample_size_all,
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
                  file = paste(out.dir, name, "_TrendsSlope_iCAR.csv", sep = ""), 
                  row.names = FALSE, 
                  append = TRUE, 
                  quote = FALSE, 
                  sep = ",", 
                  col.names = FALSE)  
  
        } # end try error
      } #end min data
    } #end if nrows = 0
  } #end sp.list
