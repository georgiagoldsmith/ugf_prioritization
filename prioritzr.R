

library(here)
library(terra)
library(sf)
library(prioritizr)

# Load data
protected_areas <- st_read(here("data/protected areas", "protected_areas_ugf.shp"))
cocoa_ugf <- rast(here("data/cocoa suitability", "cocoa_ugf.tif"))
bird_priority_raster <- rast(here("outputs/bird_priority_raster.tif"))
urban <- rast(here("data/LC/urban_ugf.tif"))
cropland_ugf <- rast(here("data/LC/cropland_ugf.tif"))
ugf_boundary <- st_read(here("data/UGF_gp.shp", "UGF_gp.shp"))
ugf_vect_3857 <- vect(st_transform(ugf_boundary, 3857))

# cocoa_ugf's grid is the template every other layer gets reprojected onto,
# but it falls slightly short of the true boundary extent on the west and
# east edges (source data grid doesn't quite reach the boundary bbox).
# Extend it now so nothing downstream inherits that shortfall; the new NA
# cells get picked up by the gap-fill step below.
cocoa_ugf <- extend(cocoa_ugf, ext(ugf_vect_3857), snap = "out")
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)

# AHAD cocoa suitability has gaps within the UGF boundary (areas outside its
# prediction envelope). Fill those with 0 suitability rather than leaving them
# NA, so they stay in the planning-unit set instead of being silently dropped
# from the whole study area. This has to happen BEFORE downscaling below --
# bilinear resampling needs valid neighbors on all sides, so leftover NA gaps
# here would bleed into neighboring cells during interpolation and erode the
# valid-data footprint at the finer resolution.
ugf_boundary_rast <- rasterize(ugf_vect_3857, cocoa_ugf, field = 1, touches = TRUE)
cocoa_ugf <- ifel(is.na(cocoa_ugf) & !is.na(ugf_boundary_rast), 0, cocoa_ugf)

# downscale from ~9.6km to 3km resolution (closer to the other inputs'
# native resolution) while keeping the planning-unit count tractable for
# the solver (~48k cells instead of ~4.5M at the other layers' full 310m)
cocoa_template_3km <- rast(ext(cocoa_ugf), resolution = 3000, crs = crs(cocoa_ugf))
cocoa_ugf <- project(cocoa_ugf, cocoa_template_3km, method = "bilinear")

# re-mask at the new resolution for a precise boundary trace (the 9.6km
# mask's edges are too coarse to be accurate once upsampled to 3km)
ugf_boundary_rast <- rasterize(ugf_vect_3857, cocoa_ugf, field = 1, touches = TRUE)
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)

# check cocoa raster
print(cocoa_ugf)
res(cocoa_ugf)
crs(cocoa_ugf)
png(here("outputs/cocoa_suitability.png"), width = 2000, height = 1200, res = 150)
plot(cocoa_ugf)
dev.off()

# check priority raster
print(bird_priority_raster)

# check protected areas
print(protected_areas)

# reproject bird priority raster to match cocoa
bird_priority_3km <- project(bird_priority_raster, cocoa_ugf, method = "bilinear")

# rasterize protected areas to match cocoa
protected_areas_v <- vect(protected_areas)
# touches=TRUE is deliberately omitted here: it's appropriate for tracing one
# large smooth boundary polygon, but protected areas are many small, scattered
# reserves with a high perimeter-to-area ratio, so touches=TRUE would inflate
# their footprint far beyond their true extent (verified: ~39% vs a true ~22%
# at 9.6km resolution). Default center-based rasterization is more accurate.
protected_areas_rast <- rasterize(protected_areas_v, cocoa_ugf, field = 1, background = 0)

# reproject and resample urban raster to match cocoa template
urban_3km <- project(urban, cocoa_ugf, method = "near")
urban_3km <- mask(urban_3km, ugf_vect_3857, touches = TRUE)
cropland_3km <- project(cropland_ugf, cocoa_ugf, method = "near")
cropland_3km <- mask(cropland_3km, ugf_vect_3857, touches = TRUE)

# check
plot(urban_3km)
png(here("outputs/urban_9km.png"), width = 2000, height = 1200, res = 150)
plot(urban_3km)
dev.off()

plot(cropland_3km)
png(here("outputs/cropland_9km.png"), width = 2000, height = 1200, res = 150)
plot(cropland_3km)
dev.off()

# check all match
ext(cocoa_ugf)
ext(bird_priority_3km)
ext(protected_areas_rast)

res(cocoa_ugf)
res(bird_priority_3km)
res(protected_areas_rast)

bird_priority_3km <- mask(bird_priority_3km, ugf_vect_3857, touches = TRUE)
protected_areas_rast <- mask(protected_areas_rast, ugf_vect_3857, touches = TRUE)

# only lock out urban/cropland cells that are NOT already protected
urban_not_protected <- ifel(protected_areas_rast == 1, NA, urban_3km)
crop_not_protected <- ifel(protected_areas_rast == 1, NA, cropland_3km)

# add_locked_out_constraints() takes a single layer, so combine urban and
# cropland into one locked-out mask instead of passing them separately
locked_out <- ifel(urban_not_protected == 1 | crop_not_protected == 1, 1, 0)

png(here("outputs/data_check.png"), width = 4000, height = 1000, res = 150)
par(mfrow = c(1, 4))
plot(cocoa_ugf, main = "Cocoa Suitability (Cost)")
plot(bird_priority_3km, main = "Bird Priority (Feature)")
plot(protected_areas_rast, main = "Protected Areas (Locked In)")
plot(locked_out, main = "Locked Out (Urban + Cropland)")
par(mfrow = c(1, 1))
dev.off()

# build species raster stack
species_stack <- NULL

for (i in 1:nrow(birds_star_t)) {
  sp <- birds_star_t$scientific_name[i]
  sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
  if (!file.exists(sp_file)) {
    sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
  }
  if (file.exists(sp_file)) {
    r <- rast(sp_file)
    r_3km <- project(r, cocoa_ugf, method = "near")
    r_3km <- mask(r_3km, ugf_vect_3857, touches = TRUE)
    names(r_3km) <- sp
    if (is.null(species_stack)) species_stack <- r_3km
    else species_stack <- c(species_stack, r_3km)
  }
}

# normalize each species raster to 0-1
species_stack_norm <- species_stack / global(species_stack, "max", na.rm = TRUE)$max

# rescale cost to 0-1
cocoa_ugf_norm <- cocoa_ugf / global(cocoa_ugf, "max", na.rm = TRUE)$max

############################
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
# bare areas, water. These get flagged as fully disturbed (1.0) below;
# actual exclusion from selection happens on the Python side via a cost
# override AFTER budget is computed from eligible-cell costs only (mirrors
# Horn et al.'s own `graph_cost[...] = budget + 1` pattern). Do NOT bake an
# exclusion sentinel into costs.tif itself -- an earlier version did this
# with a flat 1e6, which inflated CAPTAIN's budget calculation since
# budget = sum(cost) * protection_target is computed from this file's sum.
exclude_classes <- c(190, 10, 11, 20, 200, 201, 210)

# disturbance tiers for everything else, rescaled 0-1
tier_high <- c(30, 110, 120, 122, 130)     # 0.67: mosaic cropland, herbaceous mosaic, shrubland/grassland regrowth
tier_mid  <- c(40, 100)                    # 0.33: mosaic natural veg w/ cropland encroachment
tier_low  <- c(50, 60, 62, 160, 170, 180)  # 0:    intact/flooded forest (incl. mangrove)

disturbance <- ifel(lccs_3km %in% exclude_classes, 1,
              ifel(lccs_3km %in% tier_high, 0.67,
              ifel(lccs_3km %in% tier_mid, 0.33,
              ifel(lccs_3km %in% tier_low, 0, NA))))

# costs.tif is real cocoa suitability everywhere; excluded cells are
# identified downstream by disturbance == 1, not by an inflated cost here
writeRaster(cocoa_ugf_norm, file.path(env_dir, "costs.tif"), overwrite = TRUE)
writeRaster(disturbance, file.path(env_dir, "disturbance.tif"), overwrite = TRUE)

############################
# Using min_set_objective()
#############################

# get star_t scores in same order as species stack
weights <- birds_star_t$star_t[match(names(species_stack), birds_star_t$scientific_name)]

#problem
# per-species targets scaled by star_t score (0-0.3 range)
weights_norm <- weights / max(weights, na.rm = TRUE) # rescale star-t scores to 0-1
species_targets <- weights_norm * 0.3

# require at least 30% of total area selected, as an explicit linear
# constraint rather than a pseudo-feature with its own relative target
total_cells <- global(!is.na(cocoa_ugf), "sum", na.rm = TRUE)$sum
area_threshold <- round(total_cells * 0.3)
area_layer <- ifel(!is.na(cocoa_ugf), 1, NA)

# cost = cocoa suitability, so the solver minimizes total cocoa exposure
# instead of treating every cell as equally "expensive"
p <- problem(cocoa_ugf_norm, features = species_stack_norm) |>
  add_min_set_objective() |>
  add_relative_targets(species_targets) |>
  add_locked_in_constraints(protected_areas_rast) |>
  add_locked_out_constraints(locked_out) |>
  add_binary_decisions() |>
  add_gurobi_solver()

s <- solve(p)
cat("% selected:", global(s, "sum", na.rm = TRUE)$sum / total_cells * 100, "%\n")
plot(s)

# classify: 0 = not selected, 1 = newly selected, 2 = existing protected (selected via lock-in)
solution_classified <- ifel(s == 1 & protected_areas_rast == 1, 2,
                            ifel(s == 1, 1, 0))

png(here("outputs/prioritization_solution.png"), width = 2000, height = 1200, res = 150)

plot(solution_classified,
     col = c("grey80", "darkgreen", "lightgreen"),
     main = "Conservation Prioritization Solution",
     legend = FALSE)

legend("bottomleft",
       legend = c("Not selected", "Newly selected", "Existing protected areas"),
       fill = c("grey80", "darkgreen", "lightgreen"),
       cex = 1.2)

dev.off()

eval_target_coverage_summary(p, s)

################
# p2 (linear constraints)
################
# get star_t scores in same order as species stack
weights <- birds_star_t$star_t[match(names(species_stack), birds_star_t$scientific_name)]

#problem
# per-species targets scaled by star_t score (0-0.3 range)
weights_norm <- weights / max(weights, na.rm = TRUE) # rescale star-t scores to 0-1
species_targets <- weights_norm * 0.3

# require at least 30% of total area selected, as an explicit linear
# constraint rather than a pseudo-feature with its own relative target
total_cells <- global(!is.na(cocoa_ugf), "sum", na.rm = TRUE)$sum
area_threshold <- round(total_cells * 0.3)
area_layer <- ifel(!is.na(cocoa_ugf), 1, NA)

# cost = cocoa suitability, so the solver minimizes total cocoa exposure
# instead of treating every cell as equally "expensive"
p2 <- problem(cocoa_ugf_norm, features = species_stack_norm) |>
  add_min_set_objective() |>
  add_relative_targets(species_targets) |>
  add_linear_constraints(threshold = area_threshold, sense = "=", data = area_layer) |> 
  add_locked_in_constraints(protected_areas_rast) |>
  add_locked_out_constraints(locked_out) |>
  add_binary_decisions() |>
  add_gurobi_solver()

s2 <- solve(p2)
cat("% selected:", global(s, "sum", na.rm = TRUE)$sum / total_cells * 100, "%\n")
plot(s2)

# classify: 0 = not selected, 1 = newly selected, 2 = existing protected (selected via lock-in)
solution_classified <- ifel(s2 == 1 & protected_areas_rast == 1, 2,
                            ifel(s2 == 1, 1, 0))

png(here("outputs/prioritization_solution2.png"), width = 2000, height = 1200, res = 150)

plot(solution_classified,
     col = c("grey80", "darkgreen", "lightgreen"),
     main = "Conservation Prioritization Solution",
     legend = FALSE)

legend("bottomleft",
       legend = c("Not selected", "Newly selected", "Existing protected areas"),
       fill = c("grey80", "darkgreen", "lightgreen"),
       cex = 1.2)

dev.off()

eval_target_coverage_summary(p2, s2)


##################################
# Using min_shortfall_objective()
#################################

# uniform cost so budget is a literal area cap (min_shortfall's budget is
# expressed in "cost units" -- needs cost=1 per cell for it to mean area)
cost_uniform <- ifel(!is.na(cocoa_ugf), 1, NA)

# get star_t scores in same order as species stack
weights <- birds_star_t$star_t[match(names(species_stack), birds_star_t$scientific_name)]

#problem
# flat 30% target for every species -- priority is now handled entirely by
# add_feature_weights(weights) below, not by varying the target itself
species_targets_1 <- 0.3

total_cells <- global(!is.na(cocoa_ugf), "sum", na.rm = TRUE)$sum
budget <- round(total_cells * 0.3)

# downweight each species' suitability by cocoa suitability spatially (cells
# attractive to cocoa contribute less to representation), on top of the
# per-species star_t feature weights below (add_feature_weights() only takes
# one scalar per species, so it can't itself vary by cocoa -- that has to
# happen in the feature values instead)
species_stack_downweighted <- species_stack_norm * (1 - cocoa_ugf_norm)

p1 <- problem(cost_uniform, features = species_stack_downweighted) |>
  add_min_shortfall_objective(budget) |>
  add_relative_targets(species_targets_1) |>
  add_feature_weights(weights) |>
  add_locked_out_constraints(locked_out) |>
  add_locked_in_constraints(protected_areas_rast) |>
  add_binary_decisions() |>
  add_gurobi_solver()

s1 <- solve(p1)
cat("% selected:", global(s1, "sum", na.rm = TRUE)$sum / total_cells * 100, "%\n")
plot(s1)

# classify: 0 = not selected, 1 = newly selected, 2 = existing protected (selected via lock-in)
solution_classified <- ifel(s1 == 1 & protected_areas_rast == 1, 2,
                      ifel(s1 == 1, 1, 0))

png(here("outputs/prioritization_solution_1.png"), width = 2000, height = 1200, res = 150)

plot(solution_classified,
     col = c("grey80", "darkgreen", "lightgreen"),
     main = "Conservation Prioritization Solution",
     legend = FALSE)

legend("bottomleft",
       legend = c("Not selected", "Newly selected", "Existing protected areas"),
       fill = c("grey80", "darkgreen", "lightgreen"),
       cex = 1.2)

dev.off()

eval_target_coverage_summary(p1, s1)

# check % selected
selected_cells <- global(s, "sum", na.rm = TRUE)$sum
protected_cells <- global(protected_areas_rast, "sum", na.rm = TRUE)$sum
cat("% selected:", selected_cells / total_cells * 100, "%\n")
cat("% from protected areas:", protected_cells / total_cells * 100, "%\n")
cat("% from new selections:", (selected_cells - protected_cells) / total_cells * 100, "%\n")
plot(s)

#########################
# Mimics Horn et al. scenario 2b (prioritizR_optimisation.Rmd) using UGF data:
#########################

# min_set objective, protected areas locked in, economic costs as the cost
# layer, exact 30% area constraint, species feature weights + relative
# targets. Boundary penalties (spatial connectivity) deliberately omitted.
#
# Assumes prioritzr.R has already been run through at least the point where
# cocoa_ugf_norm, species_stack_norm, protected_areas_rast, weights,
# area_layer, area_threshold, and total_cells exist in the session.

library(prioritizr)

# reuse the same 30% area constraint already built for the main analysis
p_2b <- problem(x = cocoa_ugf_norm, features = species_stack_norm) |>
  add_min_set_objective() |>
  add_binary_decisions() |>
  add_locked_in_constraints(protected_areas_rast) |>
  add_locked_out_constraints(locked_out) |>
  add_linear_constraints(threshold = area_threshold, sense = "=", data = area_layer) |>
  add_feature_weights(weights) |>
  add_relative_targets(species_targets) |>
  add_gurobi_solver(gap = 0.1)

s_2b <- solve(p_2b)

cat("% selected (scenario 2b):", global(s_2b, "sum", na.rm = TRUE)$sum / total_cells * 100, "%\n")

# classify and save separately from the main prioritization outputs
solution_classified_2b <- ifel(s_2b == 1 & protected_areas_rast == 1, 2,
                               ifel(s_2b == 1, 1, 0))

png(here("outputs/prioritization_solution_scenario2b.png"), width = 2000, height = 1200, res = 150)
plot(solution_classified_2b,
     col = c("grey80", "darkgreen", "lightgreen"),
     main = "Conservation Prioritization Solution (Scenario 2b, no boundary penalties)",
     legend = FALSE)
legend("bottomleft",
       legend = c("Not selected", "Newly selected", "Existing protected areas"),
       fill = c("grey80", "darkgreen", "lightgreen"), cex = 1.2)
dev.off()

writeRaster(s_2b, here("outputs/prioritizr_scenario2b.tif"), overwrite = TRUE)


############################

# tests
urban_cells <- global(urban_not_protected, "sum", na.rm = TRUE)$sum

cat("Total cells:", total_cells, "\n")
cat("Protected (locked in):", protected_cells, "\n")
cat("Urban (locked out):", urban_cells, "\n")
cat("Need to reach 30%:", total_cells * 0.3, "\n")
cat("Already protected:", protected_cells, "\n")
cat("Additional cells needed:", total_cells * 0.3 - protected_cells, "\n")
cat("Available to select:", total_cells - protected_cells - urban_cells, "\n")









