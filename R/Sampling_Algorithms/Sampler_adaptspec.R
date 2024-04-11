Sampler_adaptspec = function(timeseries, B, iter = 1000, tmin = 30){
  # Extract n from timeseries
  n = length(timeseries)
  # Checks for inputs
  # Check Smax = n/tmin >= 1;
  Smax = n/tmin
  try(if(Smax == 0) stop("No segments can be created with given timeseries and tmin value"))
  # if 1 then give a warning.
  try(if(Smax == 1) warning("Only 1 segment can be created with given timeseries and tmin value"))
  
  
}