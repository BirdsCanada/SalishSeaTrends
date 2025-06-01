# Read trend data and clean

trend<-tolower(trend)

if (trend == "endpoint") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsEndPoint_iCAR.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd)
} else if (trend == "slope") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsSlope_iCAR.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd)
}


# Ensure Grid and map are sf objects
if(!inherits(Grid, "sf")) stop("Grid must be an sf object")
if(!inherits(map, "sf")) stop("map must be an sf object")

species_list3 <- unique(all.trends$species_code)

for(current_sp in species_list3) {
  # Filter data for current species
  species_data <- all.trends %>% filter(species_code == current_sp)
  
  # Join and clip
  Grid_sp <- Grid %>% left_join(species_data, by = c("Name" = "area_code"))
  Grid_sp<- st_transform(Grid_sp, st_crs(mapCRS))
  Grid_sp_clipped <- sf::st_intersection(Grid_sp, map)
  
  
  # Calculate percentiles for robust color scaling
  lower_lim <- quantile(Grid_sp$trnd, 0.025, na.rm = TRUE)
  upper_lim <- quantile(Grid_sp$trnd, 0.975, na.rm = TRUE)
  max_abs <- max(abs(c(lower_lim, upper_lim)), na.rm = TRUE)
  
  # Create diverging color scale centered at 0
  p <- ggplot() +
    geom_sf(data = Grid_sp_clipped, aes(fill = trnd), color = "white", size = 0.1) +
    scale_fill_gradient2(
      name = "Annual Trends (%)",
      low = "red", mid = "white", high = "blue",
      midpoint = 0,
     # limits = c(-max_abs, max_abs),
      na.value = "grey90"
    ) +
    geom_sf(data = map, fill = NA, color = "grey30", size = 0.2) +
    labs(title = paste("Annual Trends for", current_sp)) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save map
  ggsave(
    filename = file.path(plot.dir, paste0("iCAR_", name, gsub(" ", "_", current_sp), ".jpeg")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}



