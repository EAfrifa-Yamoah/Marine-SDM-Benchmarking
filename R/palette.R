# =====================================================================
# palette.R
# One colour vocabulary for every figure in the manuscript and its
# supporting information, so that a colour never carries two meanings.
#
#   Hue identifies a METHOD and nothing else. The seven methods take the
#   Okabe and Ito qualitative palette, which stays distinguishable under
#   deuteranopia, protanopia and tritanopia.
#
#   Every other categorical distinction (spatial against non spatial
#   class, training against test, presence against absence, reliability
#   class, prevalence band, variance source) is carried by position,
#   facet, shape, line type or a direct label, drawn in neutral greys,
#   so it never borrows a method hue.
#
#   Continuous magnitudes (mean AUC, correlation) use the viridis scale
#   with the value printed in each cell.
# =====================================================================

METHOD_COLOURS <- c(
  RF         = "#E69F00",   # orange
  SpatialRF  = "#56B4E9",   # sky blue
  BRT        = "#009E73",   # bluish green
  MaxEntPO   = "#0072B2",   # blue
  SpatialGAM = "#D55E00",   # vermillion
  GeostatGP  = "#CC79A7",   # reddish purple
  JointSDM   = "#000000")   # black

# neutral greys for marks that do not identify a method
INK       <- "grey15"      # emphasised marks (training hauls, presences, estimates)
INK_MID   <- "grey45"      # secondary marks (reference lines, bars)
INK_LIGHT <- "grey72"      # background marks (absences, available hauls)
LAND_FILL <- "grey93"
LAND_LINE <- "grey60"
SEA_FILL  <- "#eef3f7"

# shapes for the two method classes, used wherever the class is shown
CLASS_SHAPES <- c("spatial method" = 16, "non-spatial method" = 1)
# shapes for the translation classes of the harmonisation protocol
RELIABILITY_SHAPES <- c("exact (analytic)" = 15, "reliable" = 16,
                        "context dependent" = 1, "unreliable" = 4)
