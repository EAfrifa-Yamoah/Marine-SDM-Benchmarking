# =====================================================================
# 09_explore_spatial.R
# Spatial exploration of the four survey datasets, for contextualisation.
# Three figures:
#   (a) sampling profiles: where hauls fall in space and how sampling
#       intensity varies within each survey extent
#   (b) occurrence in space: presence and absence of a representative
#       species per survey, showing spatial structure and prevalence
#   (c) design schematic: how the coverage levels and the conditional
#       versus marginal holdout draw training and test sets from the
#       extent, illustrated with real subsamples
# =====================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(viridis); library(patchwork)
})
source(file.path("R", "config.R"))
source(file.path("R", "02_subsample.R"))

FIG <- PATHS$figures
world <- map_data("world")

codes <- c("EBS", "NS-IBTS", "GMEX", "NEUS")
regions <- c(EBS = "Eastern Bering Sea", `NS-IBTS` = "North Sea (IBTS)",
             GMEX = "Gulf of Mexico", NEUS = "Northeast US shelf")

D <- lapply(codes, function(c) readRDS(sprintf("data_processed/%s_processed.rds", c)))
names(D) <- codes

theme_map <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_line(colour = "grey92", linewidth = 0.3),
        panel.background = element_rect(fill = "#eef3f7", colour = NA),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        legend.key.height = unit(0.9, "lines"))

# per survey bounding box with a small margin
bbox <- function(h, m = 0.5) list(
  xlim = c(min(h$lon) - m, max(h$lon) + m),
  ylim = c(min(h$lat) - m, max(h$lat) + m))

land_layer <- function(bx)
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey78", colour = "grey60", linewidth = 0.2)

# ---- (a) sampling profiles: haul point cloud with density -----------
prof_panel <- function(code) {
  h <- as.data.frame(D[[code]]$hauls); bx <- bbox(h)
  ggplot() +
    land_layer(bx) +
    geom_point(data = h, aes(lon, lat), colour = "#0b6b6b",
               size = 0.14, alpha = 0.14) +
    stat_density_2d(data = h, aes(lon, lat, colour = after_stat(level)),
                    linewidth = 0.25, bins = 6) +
    scale_colour_viridis(option = "inferno", end = 0.85, guide = "none") +
    coord_quickmap(xlim = bx$xlim, ylim = bx$ylim, expand = FALSE) +
    labs(title = regions[code],
         subtitle = sprintf("%s hauls", format(nrow(h), big.mark = ",")),
         x = NULL, y = NULL) +
    theme_map
}
p_prof <- (prof_panel("EBS") | prof_panel("NS-IBTS")) /
          (prof_panel("GMEX") | prof_panel("NEUS"))
ggsave(file.path(FIG, "fig_sampling_profiles.png"), p_prof,
       width = 9.5, height = 8.4, dpi = 200, bg = "white")
cat("wrote fig_sampling_profiles.png\n")

# ---- (b) occurrence in space: common species per survey --------------
occ_panel <- function(code) {
  r <- D[[code]]; h <- as.data.frame(r$hauls); bx <- bbox(h)
  sp <- r$species
  focal <- sp$accepted_name[sp$band == "intermediate"][1]
  y <- r$occ[[focal]][match(h$haul_id, r$occ$haul_id)]
  h$occ <- factor(ifelse(y == 1, "present", "absent"),
                  levels = c("absent", "present"))
  h <- h[order(h$occ), ]                       # draw presences on top
  ggplot() +
    land_layer(bx) +
    geom_point(data = h, aes(lon, lat, colour = occ, size = occ,
                             alpha = occ)) +
    # absent is rose rather than grey or orange: grey collided with the
    # grey landmass, and orange is reserved for the test set in the design
    # schematic. Rose and blue are mutually distinguishable under all
    # common forms of colour blindness (Okabe and Ito palette).
    scale_colour_manual(values = c(absent = "#CC79A7", present = "#0072B2"),
                        name = NULL) +
    scale_size_manual(values = c(absent = 0.26, present = 0.6), guide = "none") +
    scale_alpha_manual(values = c(absent = 0.5, present = 0.9), guide = "none") +
    coord_quickmap(xlim = bx$xlim, ylim = bx$ylim, expand = FALSE) +
    labs(title = regions[code],
         subtitle = sprintf("%s  (prevalence %.0f%%)",
                            focal, 100 * mean(y)),
         x = NULL, y = NULL) +
    theme_map +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
}
# one shared legend seated in a band between the top and bottom rows,
# rather than repeating an identical legend beside each of the four panels
occ_layout <- "
AB
CC
DE
"
p_occ <- occ_panel("EBS") + occ_panel("NS-IBTS") + guide_area() +
         occ_panel("GMEX") + occ_panel("NEUS") +
  plot_layout(design = occ_layout, heights = c(1, 0.10, 1),
              guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center",
        legend.text = element_text(size = 10))
ggsave(file.path(FIG, "fig_occurrence_space.png"), p_occ,
       width = 9.5, height = 8.6, dpi = 200, bg = "white")
cat("wrote fig_occurrence_space.png\n")

# ---- (c) design schematic: coverage and holdout on one survey --------
# GMEX is compact and reads clearly; use its common species
sc <- "GMEX"
r <- D[[sc]]; h_dt <- data.table::as.data.table(r$hauls)
h <- as.data.frame(h_dt); bx <- bbox(h)
focal <- r$species$accepted_name[r$species$band == "common"][1]
y <- r$occ[[focal]][match(h$haul_id, r$occ$haul_id)]

base_map <- function() list(
  land_layer(bx),
  geom_point(data = h, aes(lon, lat), colour = "grey75", size = 0.18,
             alpha = 0.5),
  coord_quickmap(xlim = bx$xlim, ylim = bx$ylim, expand = FALSE),
  theme_map, labs(x = NULL, y = NULL))

set.seed(GLOBAL_SEED)
cover_panel <- function(cov) {
  sp <- make_split(h_dt, y, r$env_used, n = 200, coverage = cov,
                   target = "conditional", rep_id = 1, test_size = 300)
  tr <- sp$train
  ggplot() + base_map() +
    geom_point(data = tr, aes(lon, lat), colour = "#b2182b", size = 0.7) +
    labs(title = sprintf("%s coverage", cov),
         subtitle = sprintf("%d training hauls", nrow(tr)))
}
target_panel <- function(tgt) {
  sp <- make_split(h_dt, y, r$env_used, n = 200, coverage = "distributed",
                   target = tgt, rep_id = 1, test_size = 300)
  tr <- sp$train; te <- sp$test
  g <- ggplot() + base_map()
  if (tgt == "marginal") {
    # shade the actual held out block, taken from the test point extent
    g <- g + annotate("rect",
                      xmin = min(te$lon), xmax = max(te$lon),
                      ymin = min(te$lat), ymax = max(te$lat),
                      fill = "#f1a340", alpha = 0.15)
  }
  g +
    geom_point(data = tr, aes(lon, lat, colour = "training"), size = 0.6) +
    geom_point(data = te, aes(lon, lat, colour = "test"), size = 0.6) +
    scale_colour_manual(values = c(training = "#2166ac", test = "#e08214"),
                        name = NULL) +
    labs(title = sprintf("%s target", tgt),
         subtitle = if (tgt == "conditional")
           "test within training region (interpolation)" else
           "test in held out block (extrapolation)") +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
}
row1 <- cover_panel("clustered") | cover_panel("intermediate") |
        cover_panel("distributed")
row2 <- target_panel("conditional") | target_panel("marginal")
p_design <- row1 / row2 + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = sprintf("Subsampling design in space (%s, %s)", regions[sc], focal),
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(FIG, "fig_design_schematic.png"), p_design,
       width = 11, height = 8, dpi = 200, bg = "white")
cat("wrote fig_design_schematic.png\n")

cat("spatial exploration figures complete\n")
