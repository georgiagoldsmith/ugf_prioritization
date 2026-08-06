#################
#Updated from cocoa master grid script
################
library(terra)
library(sf)
library(tidyterra)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript CAPTAIN_plot_protection.R <priority_csv_path> <use_cost> <add_pa> <run_label> <diagnostics_txt_path>\n",
       "Got ", length(args), " argument(s) — this script must be run via Rscript, not sourced interactively.")
}
priority_csv     <- args[1]
use_cost         <- as.logical(args[2])
add_pa           <- as.logical(args[3])
run_label <- args[4]
diagnostics_path <- args[5]

if (!file.exists(priority_csv)) {
  stop("priority_csv path does not exist: ", priority_csv)
}

diagnostics_text <- if (file.exists(diagnostics_path)) {
  paste(readLines(diagnostics_path), collapse = "\n")
} else {
  "Diagnostics file not found."
}

#-----from cocoa_master_grid.R------###

# rebuild cocoa_ugf's final 3km grid independently (rather than assuming it's
# already in the session) -- matches prioritzr.R lines 8-44 exactly, since a
# stale/partial workspace here would silently produce the wrong template
ugf_boundary <- st_read("~/captain_test/captain2-main/captaincode/data_ugf/UGF_gp.shp") 
ugf_vect_3857 <- vect(st_transform(ugf_boundary, 3857))

cocoa_ugf <- rast("~/captain_test/captain2-main/captaincode/data_ugf/cocoa_ugf.tif")
cocoa_ugf <- extend(cocoa_ugf, ext(ugf_vect_3857), snap = "out")
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)
ugf_boundary_rast <- rasterize(ugf_vect_3857, cocoa_ugf, field = 1, touches = TRUE)
cocoa_ugf <- ifel(is.na(cocoa_ugf) & !is.na(ugf_boundary_rast), 0, cocoa_ugf)
cocoa_template_3km <- rast(ext(cocoa_ugf), resolution = 3000, crs = crs(cocoa_ugf))
cocoa_ugf <- project(cocoa_ugf, cocoa_template_3km, method = "bilinear")
cocoa_ugf <- mask(cocoa_ugf, ugf_vect_3857, touches = TRUE)

priority_mat <- as.matrix(read.csv(priority_csv, header = FALSE))
stopifnot(nrow(priority_mat) == nrow(cocoa_ugf), ncol(priority_mat) == ncol(cocoa_ugf))
captain_priority <- rast(priority_mat, extent = ext(cocoa_ugf), crs = crs(cocoa_ugf))
names(captain_priority) <- "captain_priority"
# graph_to_grid() fills non-planning-unit cells (outside the study area) with
# 0 by default, not NA -- indistinguishable on-disk from "eligible but never
# selected". Re-mask with cocoa_ugf's boundary so the true study-area edge
# shows up correctly, matching how the prioritizr rasters are already masked.
captain_priority <- mask(captain_priority, cocoa_ugf)
captain_priority_display <- ifel(captain_priority == 0, NA, captain_priority)

output_dir <- path.expand("~/captain_test/captain2-main/outputs")

# --- conditional: protected areas overlay ---
if (add_pa) {
  protected_areas_captain <- rast(file.path(dirname(priority_csv), "..", "environmental_layers", "protected_areas.tif"))
  captain_display <- ifel(protected_areas_captain == 1, -0.2, captain_priority)
  gradient_cols <- hcl.colors(6, "YlGnBu", rev = TRUE)
  captain_cols <- c("grey40", gradient_cols)
  captain_breaks <- c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1)
  
  png(file.path(output_dir, paste0(run_label, "_priority_classified.png")), width = 2000, height = 1200, res = 150)
  plot(captain_display, col = captain_cols, breaks = captain_breaks,
       main = paste0("CAPTAIN Solution (cost ", ifelse(use_cost, "on", "off"), ", PA lock-in on)"),
       plg = list(legend = c("Existing protected", "0", "0.2", "0.4", "0.6", "0.8", "1")))
  mtext(diagnostics_text, side = 1, line = 4, cex = 0.6, adj = 0)
  dev.off()
  
} else {
  p <- ggplot() +
    geom_sf(data = st_as_sf(ugf_vect_3857), fill = "#d3d3d3", color = NA) +
    geom_spatraster(data = captain_priority_display, na.rm = TRUE) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Priority", na.value = NA) +
    coord_sf() +
    labs(title = paste0("CAPTAIN Predicted Protection Priority (cost ", ifelse(use_cost, "on", "off"), ")"), 
         caption = diagnostics_text) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
      plot.caption = element_text(size = 8, hjust = 0, family = "mono"),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_rect(fill = "white", color = NA)
    )
  p <- p + labs(caption = diagnostics_text) +
    theme(plot.caption = element_text(size = 8, hjust = 0, family = "mono"))
  ggsave(file.path(output_dir, paste0(run_label, "_priority_map.png")), p, width = 12, height = 8, dpi = 150)
}
