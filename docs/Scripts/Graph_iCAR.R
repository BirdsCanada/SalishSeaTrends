# Read trend data and clean

trend<-"slope"

if (trend == "endpoint") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsEndPoint_iCAR.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd, lower_ci, upper_ci)
} else if (trend == "slope") {
  all.trends <- read.csv(file.path(out.dir, paste0(name, "_TrendsSlope_iCAR.csv"))) %>% 
    drop_na(results_code) %>% 
    select(area_code, species_code, trnd, lower_ci, upper_ci)
}

#Remove Full Study Area Results
all.trends<-all.trends %>% filter(area_code != "Full Study Area")

# Column for significance and direction
all.trends <- all.trends %>%
  mutate(
    signif = ifelse(lower_ci > 0 & upper_ci > 0, "pos",
                    ifelse(lower_ci < 0 & upper_ci < 0, "neg", NA))
  )


# Ensure Grid and map are sf objects
if(!inherits(poly, "sf")) stop("poly must be an sf object")
if(!inherits(map, "sf")) stop("map must be an sf object")

species_list3 <- unique(all.trends$species_code)

for(current_sp in species_list3) {
  # Filter data for current species
  species_data <- all.trends %>% filter(species_code == current_sp)

# Join and clip
Grid_sp <- poly %>% left_join(species_data, by = c("Name" = "area_code"))
Grid_sp <- st_transform(Grid_sp, st_crs(mapCRS))
Grid_sp_clipped <- sf::st_intersection(Grid_sp, map)

# Calculate percentiles for robust color scaling
lower_lim <- quantile(Grid_sp$trnd, 0.025, na.rm = TRUE)
upper_lim <- quantile(Grid_sp$trnd, 0.975, na.rm = TRUE)
max_abs <- max(abs(c(lower_lim, upper_lim)), na.rm = TRUE)

# Prepare centroids for hash marks
Grid_sp_clipped$centroid <- st_centroid(Grid_sp_clipped$geometry)
centroids <- st_coordinates(Grid_sp_clipped$centroid)
Grid_sp_clipped$x <- centroids[,1]
Grid_sp_clipped$y <- centroids[,2]

# Filter for significant areas only
sig_points <- Grid_sp_clipped %>% filter(!is.na(signif))

# Create diverging color scale centered at 0 and add hash marks
p <- ggplot() +
  geom_sf(data = Grid_sp_clipped, aes(fill = trnd), color = "white", size = 0.1) +
  scale_fill_gradient2(
    name = "Annual Trends (%)",
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    na.value = "grey90"
  ) +
  geom_sf(data = map, fill = NA, color = "grey30", size = 0.2) +
  # Add hash marks for significant trends
  geom_point(
    data = sig_points,
    aes(x = x, y = y, shape = signif),
    size = 3,
    color = "black"
  ) +
  scale_shape_manual(
    values = c("neg" = 4, "pos" = 3),  # 4 = x, 3 = +
    labels = c("neg" = "Significant Negative", "pos" = "Significant Positive"),
    name = "Significance"
  ) +
  labs(title = paste("Annual Trends for", current_sp)) +
    theme_minimal() +  # This line removes all axes, grids, and backgrounds
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


