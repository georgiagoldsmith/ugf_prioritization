# Convert CAPTAIN's predicted priority grid (exported as CSV from
# predict_ugf_model.py) into a georeferenced raster, using cocoa_ugf as the
# template. Every CAPTAIN input layer (species tifs, costs.tif, disturbance.tif)
# was written from r_3km <- project(r, cocoa_ugf, ...), so the CSV shares
# cocoa_ugf's exact extent/CRS/resolution/dimensions with no flip needed --
# rxr.open_rasterio(...).to_numpy() and terra both read rasters top-row-first.

library(here)
library(terra)
library(sf)

# rebuild cocoa_ugf's final 3km grid independently (rather than assuming it's
# already in the session) -- matches prioritzr.R lines 8-44 exactly, since a
# stale/partial workspace here would silently produce the wrong template
ugf_boundary <- st_read(here("data/UGF_gp.shp", "UGF_gp.shp"))
ugf_vect_3857 <- vect(st_transform(ugf_boundary, 3857))

cocoa_ugf <- rast(here("data/cocoa suitability", "cocoa_ugf.tif"))
cocoa_ugf <- extend(cocoa_ugf, ext(ugf_vect_3857), snap = "out")
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)
ugf_boundary_rast <- rasterize(ugf_vect_3857, cocoa_ugf, field = 1, touches = TRUE)
cocoa_ugf <- ifel(is.na(cocoa_ugf) & !is.na(ugf_boundary_rast), 0, cocoa_ugf)
cocoa_template_3km <- rast(ext(cocoa_ugf), resolution = 3000, crs = crs(cocoa_ugf))
cocoa_ugf <- project(cocoa_ugf, cocoa_template_3km, method = "bilinear")
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)

captain_dir <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain/captain2-main/ugf_data"
priority_csv <- file.path(captain_dir, "predictions_20250507", "ugf_pred_totalarea_pa_fix20250507_priority_grid.csv")

priority_mat <- as.matrix(read.csv(priority_csv, header = FALSE))

stopifnot(nrow(priority_mat) == nrow(cocoa_ugf), ncol(priority_mat) == ncol(cocoa_ugf))

captain_priority <- rast(priority_mat, extent = ext(cocoa_ugf), crs = crs(cocoa_ugf))
names(captain_priority) <- "captain_priority"

# graph_to_grid() fills non-planning-unit cells (outside the study area) with
# 0 by default, not NA -- indistinguishable on-disk from "eligible but never
# selected". Re-mask with cocoa_ugf's boundary so the true study-area edge
# shows up correctly, matching how the prioritizr rasters are already masked.
captain_priority <- mask(captain_priority, cocoa_ugf)

png(here("outputs/captain_priority.png"), width = 2000, height = 1200, res = 150)
plot(captain_priority, main = "CAPTAIN predicted protection priority")
dev.off()

writeRaster(captain_priority,
           file.path(captain_dir, "predictions_20250507", "captain_priority.tif"),
           overwrite = TRUE)

# classified version, but keeping the actual 0-1 selection-frequency gradient
# for newly-selected cells rather than flattening them to one flat color --
# existing protected areas still get their own single distinct color so they
# read separately from CAPTAIN's own gradient. protected_areas.tif is the same
# raster exported earlier for CAPTAIN's add_to_existing_protected_areas lock-in.
protected_areas_captain <- rast(file.path(captain_dir, "environmental_layers", "protected_areas.tif"))

# push protected cells below the 0-1 priority range into their own sentinel
# bin (-0.2) so they get a distinct flat color instead of joining the gradient
captain_display <- ifel(protected_areas_captain == 1, -0.2, captain_priority)

# priority is only ever 0, 0.2, 0.4, 0.6, 0.8, 1 (mean of 5 binary runs), so
# 6 gradient bins + 1 flat sentinel bin for protected areas
gradient_cols <- hcl.colors(6, "YlGnBu", rev = TRUE)
captain_cols <- c("grey40", gradient_cols)
captain_breaks <- c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1)

png(here("outputs/captain_priority_classified.png"), width = 2000, height = 1200, res = 150)
plot(captain_display,
     col = captain_cols,
     breaks = captain_breaks,
     main = "CAPTAIN Conservation Prioritization Solution",
     plg = list(legend = c("Existing protected", "0", "0.2", "0.4", "0.6", "0.8", "1")))
dev.off()
