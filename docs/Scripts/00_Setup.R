#To update the project versions and dependancies, run the following command 
#renv::snapshot()

options(timeout = 1200)

#Create data folder and output folders in working directory
if(!dir.exists("Data")) dir.create("Data")
if(!dir.exists("Data/Spatial")) dir.create("Data/Spatial")
if(!dir.exists("Output")) dir.create("Output")
if(!dir.exists("Output/Plots")) dir.create("Output/Plots")

source("Scripts/PSSSBMDE.R")
source("Scripts/BCCWSClean.R")
source("Scripts/PSSSClean.R")
source("Scripts/OutputTables.R")

#Assign directories
out.dir <- "Output/"
data.dir <- "Data/"
spatial.dir <- "Data/Spatial/"
plot.dir <- "Output/Plots/"

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 

#Load required libraries
librarian::shelf("BirdsCanada/naturecounts", tidyverse, sf, mapview, sdmpredictors,
                 svMisc, terra, geojsonsf, leaflet, HelpersMG, gdalUtilities, ggplot2,
                 exactextractr, readxl, reshape, ggmap, gridExtra, ggspatial, prettymapr, 
                 rnaturalearth, mapview, INLA, mgcv, sn, fmesher, inlabru, splines, 
                 maps, splancs, spdep, igraph, ggspatial, terra, tidyterra, stringr, reshape2, 
                 measurements, ggOceanMaps, leaflet, readr, shiny, rsconnect, quarto)

BMDE<-meta_bmde_fields("core")

sp.c<-meta_species_codes()
sp.c<-sp.c %>% filter(authority=="BSCDATA") %>% distinct() %>% select(species_id, species_code)

sp.tax<-meta_species_taxonomy()

sp<-left_join(sp.c, sp.tax, by="species_id")
sp<-sp %>% distinct(english_name, .keep_all = TRUE)
sp<-sp %>% dplyr::select(species_code, scientific_name, english_name) %>% distinct()

run_analysis <- function(model = c("SPDE", "iCAR")) {
  model <- match.arg(model)
# Source the appropriate analysis script
  if (model == "SPDE") {
    source("Scripts/Analysis_SPDE.R", local = knitr::knit_global())
  } else if (model == "iCAR") {
    source("Scripts/Analysis_iCAR.R", local = knitr::knit_global())
  }
  message("Analysis for '", model, "'model has been run. Check Output folder for results.")
}

# Source the appropriate graph script
graph_results <- function(name, trend = c("slope", "endpoint")) {
  trend <- tolower(trend)
  trend <- match.arg(trend)
  assign("name", name, envir = knitr::knit_global())
  assign("trend", trend, envir = knitr::knit_global())
  suppressWarnings(source("Scripts/Graph_iCAR.R", local = knitr::knit_global()))
  message("Graph results for iCAR model (", name, ") have been run. Check Output/Plot folder for results.")
}
    
#Calculate Dispersion Statistic for SPDE 
compute_dispersion_SPDE <- function(M1, Stack, fam) {
  # Extract indices of estimation points from the stack
  stack_data <- inla.stack.data(Stack)
  est_indices <- which(!is.na(stack_data$count))  # Identify non-NA responses
  # Align observed and fitted values
  observed <- stack_data$count[est_indices]
  fitted_mean <- M1$summary.fitted.values$mean[est_indices]
  # Calculate effective parameters
  p_eff <- M1$dic$p.eff
  n <- length(observed)
  # Compute Pearson residuals
  if(fam == "nbinomial") {
    theta <- M1$summary.hyperpar$mean[1]
    pearson <- (observed - fitted_mean) / sqrt(fitted_mean * (1 + fitted_mean/theta))
  } else if(fam == "poisson") {
    pearson <- (observed - fitted_mean) / sqrt(fitted_mean)
  } else {
    stop("Family must be 'poisson' or 'nbinomial'")
  }
  # Dispersion statistic
  dispersion_stat <- sum(pearson^2) / (n - p_eff)
  
  return(dispersion_stat)
} 

#Calculate Dispersion Statistic for iCAR 
calculate_dispersion_iCAR <- function(inla_model, observed, fam = "nbinomial") {
  # Extract fitted values
  mu <- inla_model$summary.fitted.values$mean
    # Initialize theta (only used for negative binomial)
  theta <- NULL
    # Family-specific calculations
  if(fam == "nbinomial") {
    # Check for negative binomial family and extract theta (size)
    theta_name <- grep("size for the nbinomial", 
                       rownames(inla_model$summary.hyperpar), 
                       value = TRUE)
    if (length(theta_name) == 0) {
      stop("Negative binomial 'size' parameter not found in model. Check family specification.")
    }
    theta <- inla_model$summary.hyperpar[theta_name, "mean"]
    
    # Pearson residuals for negative binomial
    pearson_resid <- (observed - mu) / sqrt(mu + (mu^2)/theta)
  } else if(fam == "poisson") {
    # Pearson residuals for Poisson
    pearson_resid <- (observed - mu) / sqrt(mu)
  } else {
    stop("Unsupported family. Use 'nbinomial' or 'poisson'")
  }
  # Effective number of parameters
  if (!is.null(inla_model$dic$p.eff)) {
    p_eff <- inla_model$dic$p.eff
  } else {
    warning("Effective number of parameters (p.eff) not found; using number of fixed effects instead.")
    p_eff <- nrow(inla_model$summary.fixed)
  }
  # Dispersion statistic calculation
  N <- length(observed)
  dispersion <- sum(pearson_resid^2) / (N - p_eff)
  
  return(dispersion)
}

# Function to calculate duration in hours
calculate_duration <- function(start, end) {
  # Convert start and end times to total minutes
  start_minutes <- floor(start) * 60 + round((start - floor(start)) * 60)
  end_minutes <- floor(end) * 60 + round((end - floor(end)) * 60)
  
  # Calculate the duration in minutes
  duration_minutes <- end_minutes - start_minutes
  
  # Convert duration back to hours
  duration_hours <- duration_minutes / 60
  
  return(duration_hours)
}

#Guild default is "No" unless otherwise changed in the setup script
guild<-"no"


epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)
