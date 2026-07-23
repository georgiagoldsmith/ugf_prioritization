
# export layers for CAPTAIN

captain_dir <- "/Users/georgiagoldsmith/Documents/Bren/CI-internship/captain/captain2-main/ugf_data"
sdm_dir <- file.path(captain_dir, "present_habitat_suitability")
env_dir <- file.path(captain_dir, "environmental_layers")
dir.create(sdm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)

# one .tif per species, 0-1 suitability
for (i in 1:nlyr(species_stack_norm)) {
  sp <- names(species_stack_norm)[i]
  fname <- paste0(gsub(" ", "_", sp), ".tif")
  writeRaster(species_stack_norm[[i]], file.path(sdm_dir, fname), overwrite = TRUE)
}

# reproject the full LCCS land cover classification (not just the binary
# urban/cropland layers) onto the 3km cocoa grid -- the disturbance tiering
# below needs the finer class detail that urban_ugf/cropland_ugf collapse away
lccs_ugf_raw <- rast(here("data/LC/lccs_ugf.tif"))
lccs_3km <- project(lccs_ugf_raw, cocoa_ugf, method = "near")
lccs_3km <- mask(lccs_3km, ugf_vect_3857, touches = TRUE)

# classes with no conservation value: urban, rainfed/irrigated cropland,
# bare areas, water -- excluded from selection entirely (cost far above any
# achievable budget) rather than just scored as "disturbed"
exclude_classes <- c(190, 10, 11, 20, 200, 201, 210)

# disturbance tiers for everything else, rescaled 0-1
tier_high <- c(30, 110, 120, 122, 130)     # 0.67: mosaic cropland, herbaceous mosaic, shrubland/grassland regrowth
tier_mid  <- c(40, 100)                    # 0.33: mosaic natural veg w/ cropland encroachment
tier_low  <- c(50, 60, 62, 160, 170, 180)  # 0:    intact/flooded forest (incl. mangrove)

disturbance <- ifel(lccs_3km %in% exclude_classes, 1,
                    ifel(lccs_3km %in% tier_high, 0.67,
                         ifel(lccs_3km %in% tier_mid, 0.33,
                              ifel(lccs_3km %in% tier_low, 0, NA))))

# cost: cocoa suitability everywhere, except excluded classes get a cost far
# above any achievable budget so they can never be selected
cost_final <- ifel(lccs_3km %in% exclude_classes, 1e6, cocoa_ugf_norm)

# cost and disturbance layers (NA cells are converted to 0 by CAPTAIN's
# own loader via nan_to_zero=True, so no need to gap-fill before export)
writeRaster(cost_final, file.path(env_dir, "costs.tif"), overwrite = TRUE)
writeRaster(disturbance, file.path(env_dir, "disturbance.tif"), overwrite = TRUE)
