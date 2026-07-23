library(terra)
library(here)
library(sf)

clipped_dir <- here("data/birds_aoh/clipped")
birds_star_t <- read.csv(here("outputs/birds_star_t.csv"))
ugf_boundary <- st_read(here("data/UGF_gp.shp/UGF_gp.shp"))
ugf_vect <- vect(st_transform(ugf_boundary, 4326))
aoh_table <- read.csv(here("outputs/aoh_table.csv"))

list.files(here("data/birds_aoh/clipped")) |> head(5)
r <- rast(file.path(clipped_dir, list.files(clipped_dir, pattern = "\\.tif$")[1]))
print(r)
res(r)

# use first raster as template
# create template covering full UGF extent
template <- rast(ext(ugf_vect), resolution = 0.0009920722, crs = "EPSG:4326")
values(template) <- 0
priority_raster <- template

for (i in 1:nrow(birds_star_t)) {
  sp <- birds_star_t$scientific_name[i]
  star_score <- birds_star_t$star_t[i]
  
  sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
  if (!file.exists(sp_file)) {
    sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
  }
  
  if (file.exists(sp_file)) {
    r <- rast(sp_file)
    r_proj <- project(r, template, method = "near")
    r_proj <- ifel(is.na(r_proj), 0, 1)
    priority_raster <- priority_raster + (r_proj * star_score)
    gc()
  }
}

priority_raster <- mask(priority_raster, ugf_vect)
plot(priority_raster)

#normalize values
priority_norm <- (priority_raster - global(priority_raster, "min", na.rm=TRUE)$min) / 
  (global(priority_raster, "max", na.rm=TRUE)$max - global(priority_raster, "min", na.rm=TRUE)$min)

plot(priority_norm)
writeRaster(priority_norm, here("outputs/bird_priority_normalized.tif"), overwrite=TRUE)
png(here("outputs/species_aoh_plots.png"), width = 3000, height = 1500, res = 150)

par(mfrow = c(3, 6))

for (i in 1:nrow(birds_star_t)) {
  sp <- birds_star_t$scientific_name[i]
  
  sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), "_R.tif"))
  if (!file.exists(sp_file)) {
    sp_file <- file.path(clipped_dir, paste0(gsub(" ", "_", sp), ".tif"))
  }
  
  if (file.exists(sp_file)) {
    r <- rast(sp_file)
    coltab(r) <- NULL  # remove embedded color table
    plot(r, main = sp, cex.main = 0.6, legend = FALSE, col="purple")
  }
}

par(mfrow = c(1, 1))
dev.off()
png(here("outputs/bird_priority_raster_18.png"), width = 2000, height = 1200, res = 150)
plot(priority_raster)
dev.off()


r <- rast(file.path(clipped_dir, list.files(clipped_dir, pattern = "\\.tif$")[1]))
coltab(r) <- NULL
NAflag(r)
global(r, "max", na.rm = TRUE)
values(r) |> table()



print(template)
r_criniger <- rast(file.path(clipped_dir, "Criniger_olivaceus.tif"))
print(r_criniger)
ext(template)
ext(r_criniger)

template <- rast(ext(ugf_vect), resolution = 0.0009920722, crs = "EPSG:4326")
values(template) <- 0
test_raster <- template

r <- rast(file.path(clipped_dir, "Criniger_olivaceus.tif"))
r_proj <- project(r, template, method = "near")
r_proj[is.na(r_proj)] <- 0
test_raster <- test_raster + (r_proj * 1.966692)  # its star_t score

plot(test_raster)





