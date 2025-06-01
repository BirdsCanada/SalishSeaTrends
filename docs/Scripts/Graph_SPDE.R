# Read trend data and clean

if (trend == "endpoint") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsEndPoint_SPDE.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd)
} else if (trend == "slope") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsSlope_SPDE.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd)
}

# Ensure xy and map are sf objects with correct CRS
if(!inherits(xy, "sf")) stop("xy must be an sf object")
if(!inherits(map, "sf")) stop("map must be an sf object")

# Create individual maps for each species and save them
species_list3 <- unique(all.trends$species_code)

for(current_sp in species_list3) {
  # Filter data for current species
  species_data <- all.trends %>% 
    filter(species_code == current_sp)
  
  # Join spatial data with trend data
  Grid_sp <- left_join(xy, species_data, by = c("SurveyAreaIdentifier" = "area_code"))
  
  # Calculate percentiles to exclude extreme outliers
  lower_lim <- quantile(Grid_sp$trnd, 0.025, na.rm = TRUE)
  upper_lim <- quantile(Grid_sp$trnd, 0.975, na.rm = TRUE)
  
  # Create diverging color scale centered at 0
  max_abs <- max(abs(c(lower_lim, upper_lim)))
  
  p <- ggplot() +
    geom_sf(data = Grid_sp, aes(fill = trnd), size = 2, shape = 21, color = "black") +
    scale_fill_gradient2(
      name = "Annual Trends (%)",
      low = "red", mid = "white", high = "blue",
      midpoint = 0,
      limits = c(-max_abs, max_abs),
      na.value = "grey90"
    ) +
    geom_sf(data = map, fill = NA, color = "grey30", size = 0.2) +
    labs(title = paste("Annual Trends for", current_sp)) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save map
  ggsave(filename = file.path(plot.dir, paste0(trend, "_SPDE_", name, gsub(" ", "_", current_sp), ".jpeg")), 
         plot = p, width = 8, height = 6, dpi = 300)
}
