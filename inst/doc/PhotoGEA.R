## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  fig.align = "center"
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages('PhotoGEA')

## ----example, eval = FALSE----------------------------------------------------
# # Load required packages
# library(PhotoGEA)
# library(lattice)
# 
# # Define a vector of paths to the files we wish to load; in this case, we are
# # loading example files included with the PhotoGEA package
# file_paths <- c(
#   PhotoGEA_example_file_path('c3_aci_1.xlsx'),
#   PhotoGEA_example_file_path('c3_aci_2.xlsx')
# )
# 
# # Load the data from each file
# licor_exdf_list <- lapply(file_paths, function(fpath) {
#   read_gasex_file(fpath, 'time')
# })
# 
# # Get the names of all columns that are present in all of the Licor files
# columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)
# 
# # Extract just these columns
# licor_exdf_list <- lapply(licor_exdf_list, function(x) {
#   x[ , columns_to_keep, TRUE]
# })
# 
# # Combine the data from all the files
# licor_data <- do.call(rbind, licor_exdf_list)
# 
# # Define a new column that uniquely identifies each curve
# licor_data[, 'curve_id'] <-
#   paste(licor_data[, 'species'], '-', licor_data[, 'plot'] )
# 
# # Organize the data
# licor_data <- organize_response_curve_data(
#     licor_data,
#     'curve_id',
#     c(9, 10, 16),
#     'CO2_r_sp'
# )
# 
# # Calculate the total pressure
# licor_data <- calculate_total_pressure(licor_data)
# 
# # Calculate temperature-dependent values of C3 photosynthetic parameters
# licor_data <- calculate_temperature_response(licor_data, c3_temperature_param_bernacchi)
# 
# # The default optimizer uses randomness, so we will set a seed to ensure the
# # results from this fit are always identical
# set.seed(1234)
# 
# # Fit all curves in the data set
# aci_results <- consolidate(by(
#   licor_data,
#   licor_data[, 'curve_id'],
#   fit_c3_aci,
#   Ca_atmospheric = 420
# ))

## ----echo = FALSE-------------------------------------------------------------
timing <- system.time({
# Load required packages
library(PhotoGEA)
library(lattice)

# Define a vector of paths to the files we wish to load; in this case, we are
# loading example files included with the PhotoGEA package
file_paths <- c(
  PhotoGEA_example_file_path('c3_aci_1.xlsx'),
  PhotoGEA_example_file_path('c3_aci_2.xlsx')
)

# Load the data from each file
licor_exdf_list <- lapply(file_paths, function(fpath) {
  read_gasex_file(fpath, 'time')
})

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
  x[ , columns_to_keep, TRUE]
})

# Combine the data from all the files
licor_data <- do.call(rbind, licor_exdf_list)

# Define a new column that uniquely identifies each curve
licor_data[, 'curve_id'] <-
  paste(licor_data[, 'species'], '-', licor_data[, 'plot'] )

# Organize the data
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_id',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Calculate the total pressure
licor_data <- calculate_total_pressure(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_temperature_response(licor_data, c3_temperature_param_bernacchi)

# The default optimizer uses randomness, so we will set a seed to ensure the
# results from this fit are always identical
set.seed(1234)

# Fit all curves in the data set
aci_results <- consolidate(by(
  licor_data,
  licor_data[, 'curve_id'],
  fit_c3_aci,
  Ca_atmospheric = 420
))
})

## ----echo = FALSE-------------------------------------------------------------
timing

## -----------------------------------------------------------------------------
plot_c3_aci_fit(aci_results, 'curve_id', 'Ci', ylim = c(-10, 80))

## -----------------------------------------------------------------------------
xyplot(
  A_residuals ~ Ci | curve_id,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', aci_results$fits$units$Ci, ')'),
  ylab = paste0('Assimilation rate residual (measured - fitted)\n(', aci_results$fits$units$A, ')'),
)

## -----------------------------------------------------------------------------
aci_results$parameters[, c('curve_id', 'Tp_at_25_lower', 'Tp_at_25', 'Tp_at_25_upper')]

## -----------------------------------------------------------------------------
barchart_with_errorbars(
  aci_results$parameters[, 'Vcmax_at_25'],
  aci_results$parameters[, 'species'],
  xlab = 'Species',
  ylab = paste0('Vcmax at 25 degrees C (', aci_results$parameters$units$Vcmax_at_25, ')'),
  ylim = c(0, 200)
)

