# Side-by-side comparison of CAPTAIN's predicted priority (no cost reward vs.
# cost reward on), both already exported to georeferenced rasters by
# captain_priority_to_raster.R.

library(here)
library(terra)

captain_dir <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain/captain2-main/ugf_data"

no_cost <- rast(file.path(captain_dir, "predictions_20250507", "captain_priority.tif"))
with_cost <- rast(file.path(captain_dir, "predictions_20250507", "captain_priority_cost.tif"))
protected_areas_captain <- rast(file.path(captain_dir, "environmental_layers", "protected_areas.tif"))

# same sentinel-bin + gradient scheme as captain_priority_to_raster.R
display_layer <- function(priority) ifel(protected_areas_captain == 1, -0.2, priority)

gradient_cols <- hcl.colors(6, "YlGnBu", rev = TRUE)
captain_cols <- c("grey40", gradient_cols)
captain_breaks <- c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1)

png(here("outputs/captain_priority_compare.png"), width = 4000, height = 1300, res = 150)
par(mfrow = c(1, 2))

plot(display_layer(no_cost),
     col = captain_cols,
     breaks = captain_breaks,
     main = "No cost reward",
     plg = list(legend = c("Existing protected", "0", "0.2", "0.4", "0.6", "0.8", "1")))

plot(display_layer(with_cost),
     col = captain_cols,
     breaks = captain_breaks,
     main = "Cost reward on (0.3)",
     plg = list(legend = c("Existing protected", "0", "0.2", "0.4", "0.6", "0.8", "1")))

par(mfrow = c(1, 1))
dev.off()
