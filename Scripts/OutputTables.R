## Set up output tables

output_tables <- function(name, model = c("SPDE", "iCAR")) {
  model <- match.arg(model)
  
  # File suffix based on model
  suffix <- ifelse(model == "SPDE", "_SPDE.csv", "_iCAR.csv")
  
  # 1. Indices Table
  indices.csv <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 15))
  names(indices.csv) <- c("results_code", "version", "area_code", "season", "period",
                          "species_code", "species_id", "year", "index", "stderr", "stdev",
                          "upper_ci", "lower_ci", "LOESS_index", "trend_index")
  write.table(indices.csv,
              file = file.path(out.dir, paste0(name, "_AnnualIndices", suffix)),
              row.names = FALSE, append = FALSE, quote = FALSE, sep = ",")
  
  # 2. Trends Table
  trends.csv <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 39))
  names(trends.csv) <- c("results_code", "version", "area_code", "season", "period",
                         "species_code", "species_id", "years", "year_start", "year_end",
                         "trnd", "lower_ci", "upper_ci", "index_type", "stderr", "model_type",
                         "model_fit", "percent_change", "percent_change_low", "percent_change_high",
                         "prob_decrease_0", "prob_decrease_25", "prob_decrease_30", "prob_decrease_50",
                         "prob_increase_0", "prob_increase_33", "prob_increase_100", "suitability",
                         "precision_num", "precision_cat", "coverage_num", "coverage_cat",
                         "sample_size", "sample_size_units", "prob_LD", "prob_MD", "prob_LC",
                         "prob_MI", "prob_LI")
  # write.table(trends.csv,
  #             file = file.path(out.dir, paste0(name, "_TrendsEndpoint", suffix)),
  #             row.names = FALSE, append = FALSE, quote = FALSE, sep = ",")
  # 

  write.table(trends.csv,
              file = file.path(out.dir, paste0(name, "_TrendsSlope", suffix)),
              row.names = FALSE, append = FALSE, quote = FALSE, sep = ",")
  
  # 3. Dispersion Table
  dispersion.csv <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
  names(dispersion.csv) <- c("area_code", "SpeciesCode", "dispersion")
  write.table(dispersion.csv,
              file = file.path(out.dir, paste0(name, "_DispersionStat", suffix)),
              row.names = FALSE, append = FALSE, quote = FALSE, sep = ",")
  
  message("Output tables created for model: ", name, "_", model)
}

