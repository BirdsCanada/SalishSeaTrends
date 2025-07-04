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

# # Identify all unique wyear values
# all_years <- unique(event$wyear)
# 
# # Filter basins with samples in every wyear
# event <- event %>%
#   group_by(Name) %>%
#   filter(all(all_years %in% wyear)) %>%
#   ungroup()

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
      years <- sort(unique(dat$wyear)) 
      last_year <- max(years)
      year_idx_vector <- years - last_year
      mean_year<-mean(years)
      sd_year<-sd(years)
      
      
      #Create index variables
      dat <- dat %>% mutate( 
        std_yr = wyear - Y2,
        alpha_i = as.integer(factor(dat$Name)),
        tau_i = as.integer(factor(dat$Name)),
        protocol = as.factor(dat$ProjectCode), 
        kappa = as.integer(factor(dat$SurveyAreaIdentifier)),
        year_idx = as.integer(wyear - last_year), #intercept is the expected count during the most recent year of data collection. 
        wmonth_idx = as.factor(wmonth), 
        sp_idx = as.integer(factor(CommonName)), 
        # Create area-year interaction term
        area_year = as.integer(factor(paste(Name, wyear, sep = "_"))))%>%
        arrange(SurveyAreaIdentifier, wyear)
      
      # Aggregate data to annual max count
      dat <- dat %>%
        group_by(SurveyAreaIdentifier, wyear) %>%
        slice_max(ObservationCount, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      #hyper.iid.simple <- list(prec = list(initial = 0, fixed = TRUE))
      
      #Model Formula
      if(guild=="yes"){
        
        dat$sp_idx[is.na(dat$sp_idx)] <- 999 #replace NA which are zero counts with generic sp_idx
        
        formula<- ObservationCount ~ -1 + 
          f(year_idx, model = "iid", hyper = hyper.iid) + 
          f(sp_idx, model="iid", hyper=hyper.iid) +
          # cell ICAR random intercepts
          f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE, hyper = hyper.iid) +
          # random route intercepts
          f(kappa, model = "iid", hyper = hyper.iid)+ 
          # cell ICAR random year slopes
          f(tau_i, std_yr, model="besag", graph=g1, constr=FALSE, scale.model=TRUE,
             hyper = hyper.iid) +
          f(area_year, model="iid", hyper=hyper.iid)
         
          
       
      }else{
        
        formula<- ObservationCount ~ -1 + 
          f(year_idx, model = "iid", hyper = hyper.iid) + 
          # cell ICAR random intercepts
          f(alpha_i, model="besag", graph=g1, constr=FALSE, scale.model=TRUE, hyper = hyper.iid) +
          # random route intercepts
          f(kappa, model="iid", hyper = hyper.iid) +
          # cell ICAR random year slopes
          f(tau_i, std_yr, model="besag", graph=g1, constr=FALSE, scale.model=TRUE, hyper = hyper.iid)+
          f(area_year, model="iid", hyper=hyper.iid)
        
      }
      
      if(fam == "nbinomial"){
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
      }
      
      if(fam =="poisson"){
        M2 <- try(INLA::inla(
          formula, family = fam,  data = dat, offset = log(dat$DurationInHours),
          control.predictor = list(compute = TRUE),
          control.compute = list(
            dic = TRUE, 
            waic = TRUE, 
            config = TRUE, 
            cpo = TRUE
          ),
          control.fixed = list(
            mean = 0,
            prec = 0.001
          )
          # Removed control.family block
        ))
        
      }
  
      if (inherits(M2, "try-error")) {
        cat("The iCAR model failed to run. Please check your input or model specifications.\n")
      } else {
        # Continue with normal processing
        print("Model iCAR model ran successfully.")
        
     #Dispersion Statistic
       Dispersion1 <- calculate_dispersion_iCAR(M2, dat$ObservationCount, fam)
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
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1.2) +
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
       # Extract effects for this posterior sample
       alpha_vals    <- sample$latent[grep("^alpha_i", rownames(sample$latent))]
       tau_vals      <- sample$latent[grep("^tau_i", rownames(sample$latent))]
       year_vals     <- sample$latent[grep("^year_idx", rownames(sample$latent))]
       area_year_vals <- sample$latent[grep("^area_year", rownames(sample$latent))]
       
       # Make a cell-year grid
       expand_grid(
         alpha_i = seq_along(alpha_vals),
         year_idx = year_idx_vector
       ) %>%
         mutate(
           # Lookup: year_vals should correspond to year_idx
           year_effect = year_vals[match(year_idx, sort(unique(year_idx_vector)))],
           # Create area_year index (must match how you coded it in your data/model)
           area_year = as.integer(factor(paste(alpha_i, year_idx, sep = "_"))),
           # Lookup area_year effect
           area_year_effect = area_year_vals[area_year],
           # Linear predictor with area-year interaction
           log_index = alpha_vals[alpha_i] + year_effect + tau_vals[alpha_i] * year_idx + area_year_effect,
           index1 = exp(log_index),
           wyear = last_year + year_idx
         )
     }
     
    
     # # 1000 posterior samples in post.sample1
     tmp0 <- bind_rows(
       lapply(post.sample1, calculate_indices),
       .id = "sample"
     )
     
     tmp1<-NULL
     tmp1 <- tmp0 %>%
       drop_na(index1) %>% 
       group_by(alpha_i, wyear) %>%
       summarise(
         index = mean(index1, na.rm = TRUE),
         stdev   = sd(index1, na.rm = TRUE),
         stderr   = sd(index1, na.rm = TRUE) / sqrt(n()),
         lower_ci   = quantile(index, 0.025, na.rm = TRUE),
         upper_ci   = quantile(index, 0.975, na.rm = TRUE),
         .groups = "drop"
       )
     

      #link back this the site name using the grid_key
      grid3<-grid %>% st_drop_geometry() %>% 
        select(alpha_i, Name, Area) %>% 
        mutate(area_code = Name) %>% 
        distinct()
      
      tmp1<-left_join(tmp1, grid3, by="alpha_i")
      
      ###FULL BASING INDEX VALUES###
      all_samples_with_area<-tmp0 %>% left_join(grid3 %>%  select(alpha_i, Area), by="alpha_i")
    
      area_weighted_indices <- all_samples_with_area%>%
        group_by(sample, wyear) %>%
        summarise(
          weighted_index = sum(index1 * Area, na.rm = TRUE) / sum(Area), 
          .groups = "drop"
        )
      
      tmp2_area<-NULL
      tmp2_area <- area_weighted_indices %>%
        group_by(wyear) %>%
        summarise(
          index = mean(weighted_index, na.rm = TRUE),
          stdev = sd(weighted_index, na.rm = TRUE),
          stderr = stdev / sqrt(n()),
          lower_ci = quantile(weighted_index, 0.025, na.rm = TRUE),
          upper_ci = quantile(weighted_index, 0.975, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          alpha_i =999,
          Area = 999,
          Name = "AreaWeighted", 
          area_code="BCR_5" # Identifier for merged output
         )
      
      tmp1 <-rbind(tmp1, tmp2_area)
      
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
        results_code = "BCCWS_PSSS",
        version = as.integer(format(Sys.Date(), "%Y")),
        season="Winter"
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
      
      
       ##This is written to table at the end once the slope trends are calculated
      
      ########################################################################
      ##SLOPE TRENDS ##
      
      years_span <- length(unique(years)) - 1 
      # For each posterior sample, extract tau_i values
      extract_tau <- function(sample) {
        tibble(
          alpha_i = seq_along(grep("^tau_i", rownames(sample$latent))),
          tau = sample$latent[grep("^tau_i", rownames(sample$latent))]
        )
      }

      tau_trends <- bind_rows(
        lapply(post.sample1, extract_tau),
        .id = "sample"
      ) %>%
        mutate(
          sample = as.integer(sample),
          percent_annual_change = 100 * (exp(tau) - 1),
          # If you want total percent change for a given period (years_span), add:
           total_percent_change = 100 * (exp(tau * years_span) - 1)
        )
      
      slope_summary <- tau_trends %>%
        group_by(alpha_i) %>%
        summarise(
          trnd     = mean(percent_annual_change),
          stdev    = sd(percent_annual_change),
          stderr   = stdev / sqrt(n()),
          lower_ci = quantile(percent_annual_change, 0.025),
          upper_ci = quantile(percent_annual_change, 0.975),
          .groups = "drop"
        ) %>%
        mutate(
          index_type = "Slope Trend",
          precision_num = upper_ci - lower_ci,
          precision_cat = case_when(
            precision_num < 3.5           ~ "High",
            between(precision_num, 3.5, 6.7) ~ "Medium",
            TRUE                                       ~ "Low"
          )
        ) %>%
        left_join(grid3, by = "alpha_i")  # merge with spatial/area info
      
      # If you want to calculate percent change over a span of years, do:
       # number of steps between min/max year
      slope_summary <- slope_summary %>%
        mutate(percent_change = 100 * (exp(trnd / 100 * years_span) - 1))
   
      #Full area
      calc_area_weighted_slope <- function(df) {
        # Check for at least 5 unique years
        if (length(unique(df$wyear)) < 5) {
          return(tibble(slope = NA_real_, 
                        percent_annual_change = NA_real_))
        }
        
        # Proceed with trend calculation
        fit <- lm(log(weighted_index) ~ wyear, data = df)
        slope <- coef(fit)[2]
        percent_annual_change <- 100 * (exp(slope) - 1)
        
        tibble(slope = slope, 
               percent_annual_change = percent_annual_change)
      }
      
      
      area_slope_trends <- area_weighted_indices %>%
        group_by(sample) %>%
        group_modify(~{
          fit <- lm(log(weighted_index) ~ wyear, data = .x)
          slope <- coef(fit)[2]
          percent_annual_change = 100 * (exp(slope) - 1)
          percent_change = 100 * (exp(slope * years_span) - 1)
          tibble(
            slope = slope,
            percent_annual_change = percent_annual_change,
            percent_change = percent_change
          )
        }) %>%
        ungroup()
      
      years_span <- length(unique(area_weighted_indices$wyear)) - 1
      
      area_slope_summary <- area_slope_trends %>% 
        summarise(
          trnd     = mean(percent_annual_change),
          stdev    = sd(percent_annual_change),
          stderr   = stdev / sqrt(n()),
          lower_ci = quantile(percent_annual_change, 0.025),
          upper_ci = quantile(percent_annual_change, 0.975),
          percent_change = mean(100 * (exp(slope * years_span) - 1))
        ) %>%
        mutate(
          aplha_i = 9999, 
          index_type = "Slope Trend",
          precision_num = upper_ci - lower_ci,
          precision_cat = case_when(
            precision_num < 3.5           ~ "High",
            between(precision_num, 3.5, 6.7) ~ "Medium",
            TRUE                                       ~ "Low"
          ),
          Name = "BCR_5", 
          Area = 9999, 
          area_code="BCR_5",
          percent_change = percent_change)
      
      all_cols <- union(names(slope_summary), names(area_slope_summary))
      
      # Add missing columns to slope_summary
      for(col in setdiff(all_cols, names(slope_summary))) {
        slope_summary[[col]] <- NA
      }
      
      # Add missing columns to area_slope_summary
      for(col in setdiff(all_cols, names(area_slope_summary))) {
        area_slope_summary[[col]] <- NA
      }
      
      
      slope_summary <- slope_summary[, all_cols]
      area_slope_summary <- area_slope_summary[, all_cols]
      
      
      final_slope_summary <- bind_rows(
        slope_summary,      # cell-specific slopes from your code above
        area_slope_summary  # area-weighted summary just created
      )
      
      #write output to table
      trend.out<-NULL
      trend.out <- final_slope_summary %>%
        mutate(model_type="ALPHA SPATIAL iCAR", 
               model_family = fam,
               years = paste(Y1, "-", Y2, sep = ""),
               year_start=Y1, 
               year_end = Y2,
               period ="all years",
               season = "winter",
               results_code = "BCCWS_PSSS",
               version = as.integer(format(Sys.Date(), "%Y")), 
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
               #precision_num = "",
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
     
     
    #############################################################################  
    ##Calculate trend_index for plotting on NatureCounts
      
      #For each alpha_i 
      std_yr_seq <- seq(min(dat$std_yr), max(dat$std_yr), by = 1)
      
      extract_tau_trend_points <- function(sample, std_yr_seq) {
        tau_vals <- sample$latent[grep("^tau_i", rownames(sample$latent))]
        alpha_vals <- sample$latent[grep("^alpha_i", rownames(sample$latent))]
        
        # For each area and each standardized year, calculate the trend point
        expand_grid(
          alpha_i = seq_along(tau_vals),
          std_yr = std_yr_seq
        ) %>%
          mutate(
            tau = tau_vals[alpha_i],
            alpha = alpha_vals[alpha_i],
            log_index = alpha + tau * std_yr,
            index = exp(log_index)
          )
      }
      
      trend_points <- bind_rows(
        lapply(post.sample1, extract_tau_trend_points, std_yr_seq = std_yr_seq),
        .id = "sample"
      ) %>%
        mutate(sample = as.integer(sample))
      
      trend_index <- trend_points %>%
        group_by(alpha_i, std_yr) %>%
        summarize(
          trend_index = mean(index)
        ) %>% 
        mutate(year = std_yr + Y2) %>% ungroup()
      
      trends_index<-left_join(trend_index, grid3, by="alpha_i")
      
      trends_index<-trends_index %>% ungroup() %>% 
        select(trend_index, area_code, year)
      
      
      #Calculate slope points per year for full study area using weighted indices
      area_slope_points <- area_weighted_indices %>%
        group_by(sample) %>%
        group_modify(~{
          fit <- lm(log(weighted_index) ~ wyear, data = .x)
          slope <- coef(fit)[2]
          intercept <- coef(fit)[1]
          
        # Generate predicted values for EACH YEAR
          tibble(
            wyear = unique(.x$wyear),  # Unique years in this sample
            log_pred = intercept + slope * wyear,
            index_pred = exp(log_pred),
            percent_annual_change = 100 * (exp(slope) - 1),
            percent_change = 100 * (exp(slope * years_span) - 1)
          )
        }) %>%
        ungroup()
      
      # Summarize across samples
      slope_points_summary <- area_slope_points %>%
        group_by(wyear) %>%
        summarise(
          trend_index = mean(index_pred)
        ) %>% 
        mutate(area_code = "BCR_5", year=wyear) %>% select(-wyear)
      
      #Combine output
       trends_index<-rbind(trends_index, slope_points_summary)
      
       indices.csv<-left_join(indices.csv, trends_index, by=c("area_code", "year"))
      
      # Order output before printing to table
      indices.csv<-indices.csv %>% ungroup %>% dplyr::select(results_code, version, area_code, season, period, species_code, species_id, year, index, stderr, stdev, upper_ci, lower_ci, LOESS_index, trend_index)
      
      # Write indices to table
      write.table(indices.csv,
                  file = paste(out.dir,	name, "_AnnualIndices_iCAR.csv", sep = ""),
                  row.names = FALSE,
                  append = TRUE,
                  quote = FALSE,
                  sep = ",",
                  col.names = FALSE)
      
    
#############################################################################

      
        } # end try error
      } #end min data
    } #end if nrows = 0
  } #end sp.list
