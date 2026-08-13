# create_heatmaps_en.r
# Creates Kansai municipal choropleth maps from fitted real-data results.
#
# Inputs:
#   - Fitted-value and residual CSV files (RESID_DIR/fitted_residuals_<ID>.csv)
#     columns: region, time, y1_obs, y1_pred, y1_resid, y2_obs, y2_pred, y2_resid
#   - transformed_data198_for_R.csv, region_id_mapping198_en.csv,
#     W_row_standardized_198_for_R.csv
#   - Shapefile (national N03 administrative areas)
#   - Optional Lake Biwa polygon from W09 lake data (see §1-B)


# §0. Settings

# To use pBIC results, change this to "fitted_residuals_pBIC".
RESID_DIR      <- "../Kansai_Estimation/output/fitted_residuals_pAIC"
PATH_DATA      <- "../Kansai_Estimation/data/real/transformed_data198_for_R.csv"
PATH_MAPPING   <- "region_id_mapping198_en.csv"
PATH_W         <- "../Kansai_Estimation/data/real/W_row_standardized_198_for_R.csv"
PATH_SHAPEFILE <- "N03-20210101_GML/N03-21_210101.shp"

# Optional Lake Biwa polygon. Empty or missing paths skip the lake overlay.
PATH_LAKE_BIWA <- "W09-05_GML/W09-05-g_Lake.shp"   # see §1-B

MODEL_IDS <- c("1111", "dddd")

KANSAI_PREF_CODES   <- c("25", "26", "27", "28", "29", "30")
# 25=Shiga, 26=Kyoto, 27=Osaka, 28=Hyogo, 29=Nara, 30=Wakayama

ADJACENT_PREF_CODES <- c("18", "21", "24", "23", "31", "33", "36")
# 18=Fukui, 21=Gifu, 24=Mie, 23=Aichi, 31=Tottori, 33=Okayama, 36=Tokushima

# Expand the bounding box of the six Kansai prefectures by this share.
BBOX_EXPAND <- 0.05

OUTPUT_DIR <- "heatmaps"
DPI        <- 300
MAP_WIDTH  <- 10
MAP_HEIGHT <- 12

# Water color for sea and Lake Biwa.
COLOR_WATER        <- "#cfe6f0"
COLOR_WATER_BORDER <- "#6ea3b5"


# §1. Data Acquisition Notes
# §1-A. Administrative-area N03 data:
#   https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N03-v3_1.html
#
# §1-B. Lake Biwa polygon:
#   National Land Numerical Information "Lake Data W09":
#     https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-W09-v2_2.html
#   1. Download the national ShapeFile package (for example, W09-05_GML.zip).
#   2. After extracting it, set PATH_LAKE_BIWA to the .shp path.
#   3. Attribute names vary by release. This script searches for rows whose
#      attributes contain the Japanese string for Lake Biwa, so edits are usually unnecessary.
#   Alternative: if W09 is unavailable, set PATH_LAKE_BIWA <- "" to skip the lake overlay.


# §2. Load Packages

cat("\n=== Loading packages ===\n")

required_pkgs <- c("sf", "ggplot2", "dplyr", "tidyr", "viridis",
                   "patchwork", "stringr", "spdep", "rmapshaper")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing: %s\n", pkg))
    install.packages(pkg, quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
  cat(sprintf("  ✓ %s\n", pkg))
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)


# §3. Load Data

cat("\n=== Loading data ===\n")

cat(sprintf("  Residual CSV source: %s/ (fallback: current folder)\n", RESID_DIR))

dat <- read.csv(PATH_DATA)
cat(sprintf("  Data: %d rows x %d columns (n=%d, T=%d)\n",
            nrow(dat), ncol(dat),
            length(unique(dat$region)), length(unique(dat$time))))
cat(sprintf("  Columns: %s\n", paste(colnames(dat), collapse = ", ")))

mapping <- read.csv(PATH_MAPPING)
cat(sprintf("  Mapping: %d municipalities\n", nrow(mapping)))

if (!"region_id" %in% colnames(mapping)) colnames(mapping)[1] <- "region_id"
if (!"municipality" %in% colnames(mapping) && !"municipality_name" %in% colnames(mapping)) {
  colnames(mapping)[2] <- "municipality"
}
if ("municipality" %in% colnames(mapping)) {
  mapping <- mapping %>% rename(municipality_name = municipality)
}
mapping$code_6 <- as.character(mapping$code_6)
mapping$code_5 <- substr(mapping$code_6, 1, 5)

W <- as.matrix(read.csv(PATH_W, header = TRUE))
cat(sprintf("  W matrix: %d x %d\n", nrow(W), ncol(W)))


# §4. Load Shapefile

cat("\n=== Loading shapefile ===\n")

shp_raw <- st_read(PATH_SHAPEFILE, quiet = TRUE)
shp_raw$pref_code <- substr(shp_raw$N03_007, 1, 2)
cat(sprintf("  National polygons: %d\n", nrow(shp_raw)))

# --- §4.1. Six Kansai prefectures ---
shp_kansai <- shp_raw %>% filter(pref_code %in% KANSAI_PREF_CODES)
cat(sprintf("  Six Kansai prefectures: %d polygons\n", nrow(shp_kansai)))

# Aggregate designated cities.
shp_kansai_tagged <- shp_kansai %>%
  filter(!is.na(N03_004)) %>%
  mutate(
    # N03_003 stores designated-city names in Japanese; keep the data-specific pattern.
    is_desig_ward = !is.na(N03_003) & grepl("市$", N03_003),
    city_name = ifelse(is_desig_ward, N03_003, N03_004)
  )

shp_muni <- shp_kansai_tagged %>%
  group_by(N03_001, city_name) %>%
  summarise(
    N03_003    = first(N03_003),
    N03_004    = first(N03_004),
    N03_007    = first(N03_007),
    pref_code  = first(pref_code),
    is_desig   = any(is_desig_ward),
    n_polygons = n(),
    geometry   = st_union(geometry),
    .groups    = "drop"
  ) %>%
  rename(municipality_name = city_name)

shp_muni$code_5 <- shp_muni$N03_007
desig_idx <- which(shp_muni$is_desig)
for (di in desig_idx) {
  city_nm <- shp_muni$municipality_name[di]
  m_match <- mapping %>% filter(municipality_name == city_nm)
  if (nrow(m_match) == 1) {
    shp_muni$code_5[di] <- m_match$code_5
  } else {
    ward_code_5 <- substr(as.character(shp_muni$N03_007[di]), 1, 5)
    prefix_3 <- substr(ward_code_5, 1, 3)
    code_candidates <- mapping$code_5[
      substr(mapping$code_5, 1, 3) == prefix_3 &
        suppressWarnings(as.integer(mapping$code_5) < as.integer(ward_code_5))
    ]

    if (length(code_candidates) > 0) {
      shp_muni$code_5[di] <- code_candidates[which.max(as.integer(code_candidates))]
    } else {
      warning(sprintf("Designated city %s was not found in mapping", city_nm))
    }
  }
}
cat(sprintf("  Municipalities after aggregation: %d\n", nrow(shp_muni)))


# --- §4.2. Background layer: seven adjacent prefectures ---
shp_adjacent <- shp_raw %>%
  filter(pref_code %in% ADJACENT_PREF_CODES) %>%
  ms_dissolve(field = "pref_code")
cat(sprintf("  Adjacent prefectures for background: %d\n", nrow(shp_adjacent)))

# --- §4.3. Prefecture-border lines (Kansai + adjacent prefectures) ---
shp_kansai_pref_bdy <- shp_kansai %>%
  ms_dissolve(field = "pref_code")

shp_pref_borders <- rbind(shp_kansai_pref_bdy, shp_adjacent)
cat(sprintf("  Prefecture borders (Kansai + adjacent): %d prefectures\n", nrow(shp_pref_borders)))


# --- §4.4. Bounding box ---
bb <- st_bbox(shp_kansai)
xr <- bb$xmax - bb$xmin
yr <- bb$ymax - bb$ymin
BBOX_LIMITS <- list(
  xlim = c(bb$xmin - xr * BBOX_EXPAND, bb$xmax + xr * BBOX_EXPAND),
  ylim = c(bb$ymin - yr * BBOX_EXPAND, bb$ymax + yr * BBOX_EXPAND)
)
cat(sprintf("  bbox (Kansai, +%.0f%%): x=[%.3f, %.3f], y=[%.3f, %.3f]\n",
            BBOX_EXPAND * 100,
            BBOX_LIMITS$xlim[1], BBOX_LIMITS$xlim[2],
            BBOX_LIMITS$ylim[1], BBOX_LIMITS$ylim[2]))


# --- §4.5. Lake Biwa polygon (optional) ---
shp_lake <- NULL
if (nzchar(PATH_LAKE_BIWA) && file.exists(PATH_LAKE_BIWA)) {
  lake_raw <- st_read(PATH_LAKE_BIWA, quiet = TRUE)

  # Some W09 releases omit .prj; if CRS is missing, assume JGD2000 (EPSG:4612).
  if (is.na(st_crs(lake_raw))) {
    cat("  Warning: lake data has no CRS; assuming JGD2000 (EPSG:4612)\n")
    st_crs(lake_raw) <- 4612
  }

  # Auto-detect the lake-name column (for example, W09_001).
  name_cols <- colnames(lake_raw)[sapply(lake_raw, is.character)]
  hit <- FALSE
  for (nc in name_cols) {
    v <- as.character(lake_raw[[nc]])
    idx <- grepl("琵琶", v)
    if (any(idx, na.rm = TRUE)) {
      shp_lake <- lake_raw[which(idx), ] %>% st_transform(st_crs(shp_muni))
      cat(sprintf("  Lake Biwa (column %s): %d polygons\n", nc, nrow(shp_lake)))
      hit <- TRUE
      break
    }
  }

  # Fallback when no name column identifies Lake Biwa: take the largest lake polygon in Shiga.
  if (!hit) {
    cat("  Warning: Lake Biwa not identified by name; trying area-based fallback\n")
    lake_raw <- st_transform(lake_raw, st_crs(shp_muni))
    shiga_bb <- shp_kansai %>% filter(pref_code == "25") %>% st_bbox() %>% st_as_sfc()
    st_crs(shiga_bb) <- st_crs(shp_muni)
    inter_idx <- which(lengths(st_intersects(lake_raw, shiga_bb)) > 0)
    if (length(inter_idx) > 0) {
      areas <- as.numeric(st_area(lake_raw[inter_idx, ]))
      biwa_i <- inter_idx[which.max(areas)]
      shp_lake <- lake_raw[biwa_i, ]
      cat(sprintf("  Lake Biwa (largest lake polygon in Shiga): area %.1f km²\n",
                  as.numeric(st_area(shp_lake)) / 1e6))
      hit <- TRUE
    }
  }

  if (!hit) {
    cat("  Warning: Lake Biwa was not found; please check attribute names and values\n")
  }
} else {
  cat("  Warning: Lake Biwa file is unset or missing; skipping lake overlay\n")
}


# §5. Municipality-Code Matching

cat("\n=== Matching municipality codes ===\n")

shp_merged <- shp_muni %>%
  inner_join(
    mapping %>% select(region_id, municipality_name, code_5),
    by = "code_5"
  ) %>%
  rename(
    municipality_name = municipality_name.y,
    shp_municipality  = municipality_name.x
  )
cat(sprintf("  Matched: %d / %d\n", nrow(shp_merged), nrow(mapping)))


# §6. Join Observed Values

cat("\n=== Preparing observed-value data ===\n")

dat_t1 <- dat %>% filter(time == 1) %>% arrange(region)
dat_t2 <- dat %>% filter(time == 2) %>% arrange(region)

shp_plot <- shp_merged %>%
  left_join(dat_t1 %>% select(region, y1_t1 = y1, y2_t1 = y2),
            by = c("region_id" = "region")) %>%
  left_join(dat_t2 %>% select(region, y1_t2 = y1, y2_t2 = y2),
            by = c("region_id" = "region"))

x_cols <- setdiff(colnames(dat), c("region", "time", "y1", "y2"))
# Variable definitions follow Data_I (see transformed_data198_for_R_README.md):
#   y1 = share of population aged 65+ (standardized)
#   y2 = secondary-industry employment ratio (standardized)
#   x_common1     = log taxable income
#   x_common2     = log(1 + foreign population)
#   x_specific1_1 = net in-migration rate (%)   [y1 equation]
#   x_specific2_1 = commuting inflow rate (%)    [y2 equation]
X_LABELS <- c(
  "x_common1"     = "x_common1 (Log Taxable Income)",
  "x_common2"     = "x_common2 (Log Foreign Population)",
  "x_specific1_1" = "x_specific1_1 (Net In-migration Rate %)",
  "x_specific2_1" = "x_specific2_1 (Commuting Inflow Rate %)"
)

for (xc in x_cols) {
  shp_plot <- shp_plot %>%
    left_join(dat_t1 %>% select(region, !!paste0(xc,"_t1") := !!sym(xc)),
              by = c("region_id" = "region")) %>%
    left_join(dat_t2 %>% select(region, !!paste0(xc,"_t2") := !!sym(xc)),
              by = c("region_id" = "region"))
}


# §7. Load Fitted Values and Residuals

cat("\n=== Loading fitted values and residuals (CSV) ===\n")

# Read fitted_residuals_<ID>.csv and return rows sorted by region.
#   columns: region, time, y1_obs, y1_pred, y1_resid, y2_obs, y2_pred, y2_resid
#   Prefer RESID_DIR, then the current folder. Return NULL if not found.
read_fitted_resid <- function(mid) {
  fname <- sprintf("fitted_residuals_%s.csv", mid)
  cand  <- c(file.path(RESID_DIR, fname), fname)
  fpath <- cand[file.exists(cand)][1]
  if (is.na(fpath)) return(NULL)

  df <- read.csv(fpath)
  df <- df[order(df$region), ]   # Align with W matrix and dat_t2 by ascending region.
  cat(sprintf("  ✓ %s (%d rows)\n", fpath, nrow(df)))
  list(fitted_y1 = df$y1_pred,  fitted_y2 = df$y2_pred,
       resid_y1  = df$y1_resid, resid_y2  = df$y2_resid,
       y1        = df$y1_obs,   y2        = df$y2_obs)
}

model_outputs <- list()
for (mid in MODEL_IDS) {
  out <- read_fitted_resid(mid)
  if (is.null(out)) {
    cat(sprintf("  CSV for model %s was not found; skipping\n", mid)); next
  }
  model_outputs[[mid]] <- out
}

for (mid in names(model_outputs)) {
  out <- model_outputs[[mid]]
  df_model <- data.frame(
    region_id = dat_t2$region,
    fitted_y1 = out$fitted_y1, fitted_y2 = out$fitted_y2,
    resid_y1  = out$resid_y1,  resid_y2  = out$resid_y2
  )
  colnames(df_model)[-1] <- paste0(colnames(df_model)[-1], "_", mid)
  shp_plot <- shp_plot %>% left_join(df_model, by = "region_id")
}


# §7.4. Print Summary Statistics

cat("\n=== Summary statistics (n, mean, sd, min, Q1, median, Q3, max) ===\n")

desc_stats <- function(x) {
  x <- x[!is.na(x)]
  c(n      = length(x),
    mean   = mean(x),
    sd     = stats::sd(x),
    min    = min(x),
    Q1     = as.numeric(quantile(x, 0.25)),
    median = stats::median(x),
    Q3     = as.numeric(quantile(x, 0.75)),
    max    = max(x))
}

print_desc_table <- function(stats_mat) {
  cat(sprintf("  %-32s %6s %10s %10s %10s %10s %10s %10s %10s\n",
              "variable", "n", "mean", "sd", "min", "Q1", "median", "Q3", "max"))
  for (i in seq_len(nrow(stats_mat))) {
    nm <- rownames(stats_mat)[i]
    v  <- stats_mat[i, ]
    cat(sprintf("  %-32s %6d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                nm, as.integer(v["n"]), v["mean"], v["sd"],
                v["min"], v["Q1"], v["median"], v["Q3"], v["max"]))
  }
}

cat("\n--- Response variables (y1, y2 x t1, t2) ---\n")
resp_mat <- rbind(
  "y1 (2015)" = desc_stats(dat_t1$y1),
  "y2 (2015)" = desc_stats(dat_t1$y2),
  "y1 (2020)" = desc_stats(dat_t2$y1),
  "y2 (2020)" = desc_stats(dat_t2$y2)
)
print_desc_table(resp_mat)

cat("\n--- Explanatory variables (each x x t1, t2) ---\n")
expl_rows <- list()
for (xc in x_cols) {
  expl_rows[[paste0(xc, " (2015)")]] <- desc_stats(dat_t1[[xc]])
  expl_rows[[paste0(xc, " (2020)")]] <- desc_stats(dat_t2[[xc]])
}
expl_mat <- do.call(rbind, expl_rows)
print_desc_table(expl_mat)

cat("\n--- Fitted values yhat (each model) ---\n")
fit_rows <- list()
for (mid in names(model_outputs)) {
  out <- model_outputs[[mid]]
  fit_rows[[sprintf("fit_y1 (%s)", mid)]] <- desc_stats(out$fitted_y1)
  fit_rows[[sprintf("fit_y2 (%s)", mid)]] <- desc_stats(out$fitted_y2)
}
fit_mat <- do.call(rbind, fit_rows)
print_desc_table(fit_mat)

cat("\n--- Residuals e (each model) ---\n")
res_rows <- list()
for (mid in names(model_outputs)) {
  out <- model_outputs[[mid]]
  res_rows[[sprintf("res_y1 (%s)", mid)]] <- desc_stats(out$resid_y1)
  res_rows[[sprintf("res_y2 (%s)", mid)]] <- desc_stats(out$resid_y2)
}
res_mat <- do.call(rbind, res_rows)
print_desc_table(res_mat)


# §7.5. Compute Moran's I

cat("\n=== Computing Moran's I ===\n")

# W is already row-standardized.
W_listw <- mat2listw(W, style = "W", zero.policy = TRUE)

compute_moran <- function(x, listw, label = "") {
  mt <- moran.test(x, listw, zero.policy = TRUE, randomisation = TRUE)
  variance_I <- as.numeric(mt$estimate["Variance"])
  sd_I       <- sqrt(variance_I)
  res <- list(I  = as.numeric(mt$estimate["Moran I statistic"]),
              sd = sd_I,
              p  = mt$p.value,
              z  = as.numeric(mt$statistic))
  if (nchar(label) > 0) {
    cat(sprintf("  %-32s  I = %+.4f   SD = %.4f   z = %+.3f   p = %.3g\n",
                label, res$I, res$sd, res$z, res$p))
  }
  res
}

moran_results <- list()

moran_results$y1_t1 <- compute_moran(dat_t1$y1, W_listw, "y1 2015 (std_aging_rate_65plus)")
moran_results$y2_t1 <- compute_moran(dat_t1$y2, W_listw, "y2 2015 (std_secondary_industry_ratio)")
moran_results$y1_t2 <- compute_moran(dat_t2$y1, W_listw, "y1 2020 (std_aging_rate_65plus)")
moran_results$y2_t2 <- compute_moran(dat_t2$y2, W_listw, "y2 2020 (std_secondary_industry_ratio)")

for (mid in names(model_outputs)) {
  out <- model_outputs[[mid]]
  moran_results[[paste0("fitted_y1_", mid)]] <-
    compute_moran(out$fitted_y1, W_listw, sprintf("ŷ1 (%s)", mid))
  moran_results[[paste0("fitted_y2_", mid)]] <-
    compute_moran(out$fitted_y2, W_listw, sprintf("ŷ2 (%s)", mid))
  moran_results[[paste0("resid_y1_", mid)]] <-
    compute_moran(out$resid_y1, W_listw, sprintf("e1 (%s)", mid))
  moran_results[[paste0("resid_y2_", mid)]] <-
    compute_moran(out$resid_y2, W_listw, sprintf("e2 (%s)", mid))
}

for (xc in x_cols) {
  moran_results[[paste0(xc, "_t1")]] <-
    compute_moran(dat_t1[[xc]], W_listw, sprintf("%s 2015", xc))
  moran_results[[paste0(xc, "_t2")]] <-
    compute_moran(dat_t2[[xc]], W_listw, sprintf("%s 2020", xc))
}

# §8. Heatmap Plotting Function

plot_heatmap <- function(shp, var_col, title, subtitle = NULL,
                         palette = "viridis", midpoint = NULL,
                         legend_title = "value", limits = NULL) {

  p <- ggplot()

  p <- p + geom_sf(data = shp_adjacent,
                   fill = "grey93", color = "grey70", linewidth = 0.2)

  p <- p + geom_sf(data = shp, aes(fill = .data[[var_col]]),
                   color = "grey55", linewidth = 0.1)

  if (!is.null(shp_lake)) {
    p <- p + geom_sf(data = shp_lake, fill = COLOR_WATER,
                     color = COLOR_WATER_BORDER, linewidth = 0.3)
  }

  p <- p + geom_sf(data = shp_pref_borders, fill = NA,
                   color = "black", linewidth = 0.45)

  p <- p + coord_sf(xlim = BBOX_LIMITS$xlim,
                    ylim = BBOX_LIMITS$ylim, expand = FALSE)

  if (!is.null(midpoint)) {
    # Use explicit limits when supplied, for shared residual and difference scales.
    if (!is.null(limits)) {
      use_limits <- limits
    } else {
      max_abs <- max(abs(shp[[var_col]]), na.rm = TRUE)
      use_limits <- c(-max_abs, max_abs)
    }
    p <- p + scale_fill_distiller(
      palette = "RdBu", direction = 1,
      limits = use_limits,
      name = legend_title, na.value = "grey85"
    )
  } else {
    p <- p + scale_fill_viridis_c(
      option = palette, name = legend_title,
      limits = limits, na.value = "grey85"
    )
  }

  cap <- title
  if (!is.null(subtitle)) cap <- paste0(title, "\n", subtitle)

  p + theme_minimal() +
    theme(
      plot.caption  = element_text(size = 11, face = "bold", hjust = 0.5,
                                   color = "grey20", margin = margin(t = 6)),
      legend.position = "right",
      legend.title  = element_text(size = 9),
      legend.text   = element_text(size = 8),
      axis.text     = element_blank(),
      axis.ticks    = element_blank(),
      panel.grid    = element_blank(),
      panel.background = element_rect(fill = COLOR_WATER, color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin   = margin(5, 5, 5, 5)
    ) +
    labs(caption = cap)
}


# §8.5. Compute Shared Color Scales

cat("\n=== Computing shared color scales ===\n")

y1_all <- c(shp_plot$y1_t1, shp_plot$y1_t2)
y2_all <- c(shp_plot$y2_t1, shp_plot$y2_t2)
for (mid in names(model_outputs)) {
  y1_all <- c(y1_all, shp_plot[[paste0("fitted_y1_", mid)]])
  y2_all <- c(y2_all, shp_plot[[paste0("fitted_y2_", mid)]])
}
y1_range <- range(y1_all, na.rm = TRUE)
y2_range <- range(y2_all, na.rm = TRUE)
cat(sprintf("  y1 shared range (obs+fitted): [%.3f, %.3f]\n", y1_range[1], y1_range[2]))
cat(sprintf("  y2 shared range (obs+fitted): [%.3f, %.3f]\n", y2_range[1], y2_range[2]))

resid_y1_all <- numeric(0); resid_y2_all <- numeric(0)
for (mid in names(model_outputs)) {
  resid_y1_all <- c(resid_y1_all, shp_plot[[paste0("resid_y1_", mid)]])
  resid_y2_all <- c(resid_y2_all, shp_plot[[paste0("resid_y2_", mid)]])
}
# Keep resid_y1_limits and resid_y2_limits for downstream compatibility, but set them equal.
resid_max <- max(abs(c(resid_y1_all, resid_y2_all)), na.rm = TRUE)
resid_y1_limits <- c(-resid_max, resid_max)
resid_y2_limits <- c(-resid_max, resid_max)
cat(sprintf("  e1/e2 shared symmetric range (all outcomes/models): +/- %.3f\n", resid_max))


# §9. Generate Choropleth Maps

cat("\n=== Generating choropleth maps ===\n")

# --- §9a. Observed-value maps (y1, y2 x t1, t2) ---
p_y1_t1 <- plot_heatmap(shp_plot, "y1_t1",
                        title = "Aging Rate 65+ (standardized)",
                        subtitle = "2015", palette = "plasma",
                        legend_title = "std_y1",
                        limits = y1_range)
p_y1_t2 <- plot_heatmap(shp_plot, "y1_t2",
                        title = "Aging Rate 65+ (standardized)",
                        subtitle = "2020", palette = "plasma",
                        legend_title = "std_y1",
                        limits = y1_range)
p_y2_t1 <- plot_heatmap(shp_plot, "y2_t1",
                        title = "Secondary Industry Ratio (standardized)",
                        subtitle = "2015", palette = "mako",
                        legend_title = "std_y2",
                        limits = y2_range)
p_y2_t2 <- plot_heatmap(shp_plot, "y2_t2",
                        title = "Secondary Industry Ratio (standardized)",
                        subtitle = "2020", palette = "mako",
                        legend_title = "std_y2",
                        limits = y2_range)

p_observed <- (p_y1_t1 | p_y1_t2) / (p_y2_t1 | p_y2_t2)

ggsave(file.path(OUTPUT_DIR, "observed_values.png"),
       plot = p_observed, width = MAP_WIDTH, height = MAP_HEIGHT, dpi = DPI)
cat("  ✓ observed_values.png\n")


# --- §9a-1. 2020 only (y1_2020 and y2_2020 side by side) ---
p_observed_2020 <- (p_y1_t2 | p_y2_t2)
ggsave(file.path(OUTPUT_DIR, "observed_values_2020.png"),
       plot = p_observed_2020, width = 12, height = 7, dpi = DPI)
cat("  ✓ observed_values_2020.png\n")


# --- §9a-2. Explanatory-variable maps (one 2020 map per x variable) ---
cat("\n  Explanatory variable maps...\n")

for (xc in x_cols) {

  label  <- if (xc %in% names(X_LABELS)) X_LABELS[[xc]] else xc
  col_t2 <- paste0(xc, "_t2")

  val_range <- range(shp_plot[[col_t2]], na.rm = TRUE)

  p_x_t2 <- plot_heatmap(shp_plot, col_t2,
                         title = label, subtitle = "2020",
                         palette = "viridis", legend_title = xc,
                         limits = val_range)

  fname_x <- sprintf("xvar_%s.png", xc)
  ggsave(file.path(OUTPUT_DIR, fname_x),
         plot = p_x_t2, width = 7, height = 8, dpi = DPI)
  cat(sprintf("  ✓ %s  (%s)\n", fname_x, label))
}


# --- §9b. By model: fitted values + residuals ---
for (mid in names(model_outputs)) {
  model_label <- mid

  p_fit_y1 <- plot_heatmap(shp_plot, paste0("fitted_y1_", mid),
                           title = "Fitted y1: Aging Rate 65+",
                           subtitle = sprintf("Model %s, 2020", model_label),
                           palette = "plasma",
                           legend_title = expression(hat(y)[1]),
                           limits = y1_range)
  p_fit_y2 <- plot_heatmap(shp_plot, paste0("fitted_y2_", mid),
                           title = "Fitted y2: Secondary Industry Ratio",
                           subtitle = sprintf("Model %s, 2020", model_label),
                           palette = "mako",
                           legend_title = expression(hat(y)[2]),
                           limits = y2_range)
  p_res_y1 <- plot_heatmap(shp_plot, paste0("resid_y1_", mid),
                           title = "Residual: Aging Rate 65+",
                           subtitle = sprintf("Model %s, 2020", model_label),
                           midpoint = 0,
                           legend_title = expression(y[1]-hat(y)[1]),
                           limits = resid_y1_limits)
  p_res_y2 <- plot_heatmap(shp_plot, paste0("resid_y2_", mid),
                           title = "Residual: Secondary Industry Ratio",
                           subtitle = sprintf("Model %s, 2020", model_label),
                           midpoint = 0,
                           legend_title = expression(y[2]-hat(y)[2]),
                           limits = resid_y2_limits)

  p_model <- (p_fit_y1 | p_fit_y2) / (p_res_y1 | p_res_y2)

  ggsave(file.path(OUTPUT_DIR, sprintf("model_%s_fitted_resid.png", mid)),
         plot = p_model, width = MAP_WIDTH, height = MAP_HEIGHT, dpi = DPI)
  cat(sprintf("  ✓ model_%s_fitted_resid.png\n", mid))
}


# §10. Residual-Difference Map Between Two Models (Signed)

if (length(names(model_outputs)) >= 2) {

  cat("\n  Residual difference map (signed) ...\n")

  mid1 <- names(model_outputs)[1]
  mid2 <- names(model_outputs)[2]

  shp_plot[["resid_diff_y1"]] <-
    shp_plot[[paste0("resid_y1_", mid2)]] - shp_plot[[paste0("resid_y1_", mid1)]]
  shp_plot[["resid_diff_y2"]] <-
    shp_plot[[paste0("resid_y2_", mid2)]] - shp_plot[[paste0("resid_y2_", mid1)]]

  diff_y1_vec <- model_outputs[[mid2]]$resid_y1 - model_outputs[[mid1]]$resid_y1
  diff_y2_vec <- model_outputs[[mid2]]$resid_y2 - model_outputs[[mid1]]$resid_y2
  compute_moran(diff_y1_vec, W_listw, sprintf("e_%s - e_%s (y1)", mid2, mid1))
  compute_moran(diff_y2_vec, W_listw, sprintf("e_%s - e_%s (y2)", mid2, mid1))

  resid_diff_max <- max(abs(c(shp_plot[["resid_diff_y1"]],
                              shp_plot[["resid_diff_y2"]])), na.rm = TRUE)
  resid_diff_limits <- c(-resid_diff_max, resid_diff_max)
  cat(sprintf("  Residual-difference shared symmetric range (y1/y2): +/- %.3f\n", resid_diff_max))

  p_diff_y1 <- plot_heatmap(
    shp_plot, "resid_diff_y1",
    title = "Residual Difference (y1, signed)",
    subtitle = sprintf("e_%s - e_%s", mid2, mid1),
    midpoint = 0,
    legend_title = sprintf("e_%s - e_%s", mid2, mid1),
    limits = resid_diff_limits
  )
  p_diff_y2 <- plot_heatmap(
    shp_plot, "resid_diff_y2",
    title = "Residual Difference (y2, signed)",
    subtitle = sprintf("e_%s - e_%s", mid2, mid1),
    midpoint = 0,
    legend_title = sprintf("e_%s - e_%s", mid2, mid1),
    limits = resid_diff_limits
  )

  p_diff <- (p_diff_y1 | p_diff_y2)

  ggsave(file.path(OUTPUT_DIR, "model_comparison_residuals_signed.png"),
         plot = p_diff, width = MAP_WIDTH, height = MAP_HEIGHT * 0.55, dpi = DPI)
  cat("  ✓ model_comparison_residuals_signed.png\n")
}


# §11. Generate Scatter Plots

cat("\n=== Generating scatter plots ===\n")

scatter_theme <- theme_minimal() +
  theme(plot.caption = element_text(size = 11, face = "bold", hjust = 0.5,
                                    color = "grey20", margin = margin(t = 6)),
        panel.grid.minor = element_blank())

df_sc <- data.frame(
  region_id = dat_t2$region,
  y1 = dat_t2$y1, y2 = dat_t2$y2
)
for (mid in names(model_outputs)) {
  out <- model_outputs[[mid]]
  df_sc[[paste0("fit_y1_", mid)]] <- out$fitted_y1
  df_sc[[paste0("fit_y2_", mid)]] <- out$fitted_y2
  df_sc[[paste0("res_y1_", mid)]] <- out$resid_y1
  df_sc[[paste0("res_y2_", mid)]] <- out$resid_y2
}

lag_vec <- function(x) as.numeric(W %*% x)
df_sc$Wy1 <- lag_vec(df_sc$y1)
df_sc$Wy2 <- lag_vec(df_sc$y2)
for (mid in names(model_outputs)) {
  df_sc[[paste0("Wres_y1_", mid)]] <- lag_vec(df_sc[[paste0("res_y1_", mid)]])
  df_sc[[paste0("Wres_y2_", mid)]] <- lag_vec(df_sc[[paste0("res_y2_", mid)]])
}


# --- #1 Observed vs Fitted ---
cat("  [#1] Observed vs Fitted (colored by residual) ...\n")

make_obs_fit <- function(df, fit_col, obs_col, title,
                         color_col = NULL, color_limits = NULL,
                         color_label = "residual",
                         console_label = NULL) {
  r <- cor(df[[obs_col]], df[[fit_col]], use = "complete.obs")
  rng <- range(c(df[[obs_col]], df[[fit_col]]), na.rm = TRUE)
  sub <- sprintf("r = %.3f", r)

  if (!is.null(console_label)) {
    cat(sprintf("    %-32s r = %+.4f   R² = %.4f\n",
                console_label, r, r^2))
  }

  if (!is.null(color_col)) {
    p <- ggplot(df, aes(x = .data[[obs_col]], y = .data[[fit_col]],
                        fill = .data[[color_col]])) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(shape = 21, color = "grey25", stroke = 0.35,
                 size = 2.0, alpha = 0.9) +
      scale_fill_distiller(palette = "RdBu", direction = 1,
                           limits = color_limits, name = color_label)
  } else {
    p <- ggplot(df, aes(x = .data[[obs_col]], y = .data[[fit_col]])) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(shape = 21, fill = "steelblue", color = "grey25",
                 stroke = 0.3, size = 1.7, alpha = 0.85)
  }

  p + coord_fixed(xlim = rng, ylim = rng) +
    labs(caption = paste0(title, "\n", sub),
         x = "Observed y", y = "Fitted ŷ") +
    scatter_theme
}

cat("\n  --- Correlation coefficients (#1 Obs vs Fit) ---\n")

panels_1 <- list()
for (mid in names(model_outputs)) {
  panels_1[[paste0("y1_", mid)]] <-
    make_obs_fit(df_sc, paste0("fit_y1_", mid), "y1",
                 sprintf("y1: Fit vs Obs (Model %s)", mid),
                 color_col = paste0("res_y1_", mid),
                 color_limits = resid_y1_limits,
                 color_label = "e1",
                 console_label = sprintf("y1 vs ŷ1 [Model %s]", mid))
  panels_1[[paste0("y2_", mid)]] <-
    make_obs_fit(df_sc, paste0("fit_y2_", mid), "y2",
                 sprintf("y2: Fit vs Obs (Model %s)", mid),
                 color_col = paste0("res_y2_", mid),
                 color_limits = resid_y2_limits,
                 color_label = "e2",
                 console_label = sprintf("y2 vs ŷ2 [Model %s]", mid))
}
p_scatter_1 <- wrap_plots(panels_1, ncol = 2)
ggsave(file.path(OUTPUT_DIR, "scatter_1_obs_vs_fitted.png"),
       plot = p_scatter_1, width = 10, height = 10, dpi = DPI)
cat("  ✓ scatter_1_obs_vs_fitted.png\n")


# --- #2 Residuals vs Fitted ---
cat("  [#2] Residuals vs Fitted (colored by observed y) ...\n")
cat("\n  --- Correlation coefficients (#2 Resid vs Fit) ---\n")
cat("  (r should be close to 0 for a well-specified model)\n")

make_res_fit <- function(df, fit_col, res_col, title,
                         color_col = NULL, color_limits = NULL,
                         color_palette = "plasma", color_label = "y",
                         console_label = NULL) {
  r <- cor(df[[fit_col]], df[[res_col]], use = "complete.obs")
  if (!is.null(console_label)) {
    cat(sprintf("    %-32s r = %+.4f   R² = %.4f\n",
                console_label, r, r^2))
  }
  sub <- sprintf("r(ŷ, e) = %+.3f", r)

  if (!is.null(color_col)) {
    # Put fill inside geom_point() so geom_smooth() does not inherit it.
    p <- ggplot(df, aes(x = .data[[fit_col]], y = .data[[res_col]])) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(fill = .data[[color_col]]),
                 shape = 21, color = "grey25", stroke = 0.35,
                 size = 2.0, alpha = 0.9) +
      geom_smooth(method = "loess", se = FALSE, color = "darkred",
                  linewidth = 0.6, formula = y ~ x) +
      scale_fill_viridis_c(option = color_palette,
                           limits = color_limits, name = color_label)
  } else {
    p <- ggplot(df, aes(x = .data[[fit_col]], y = .data[[res_col]])) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(shape = 21, fill = "tomato", color = "grey25",
                 stroke = 0.3, size = 1.7, alpha = 0.85) +
      geom_smooth(method = "loess", se = FALSE, color = "darkred",
                  linewidth = 0.6, formula = y ~ x)
  }
  p + labs(caption = paste0(title, "\n", sub),
           x = "Fitted ŷ", y = "Residual e") + scatter_theme
}

panels_2 <- list()
for (mid in names(model_outputs)) {
  panels_2[[paste0("y1_", mid)]] <-
    make_res_fit(df_sc, paste0("fit_y1_", mid), paste0("res_y1_", mid),
                 sprintf("y1: Resid vs Fit (Model %s)", mid),
                 color_col = "y1", color_limits = y1_range,
                 color_palette = "plasma", color_label = "y1",
                 console_label = sprintf("ŷ1 vs e1 [Model %s]", mid))
  panels_2[[paste0("y2_", mid)]] <-
    make_res_fit(df_sc, paste0("fit_y2_", mid), paste0("res_y2_", mid),
                 sprintf("y2: Resid vs Fit (Model %s)", mid),
                 color_col = "y2", color_limits = y2_range,
                 color_palette = "mako", color_label = "y2",
                 console_label = sprintf("ŷ2 vs e2 [Model %s]", mid))
}
p_scatter_2 <- wrap_plots(panels_2, ncol = 2)
ggsave(file.path(OUTPUT_DIR, "scatter_2_resid_vs_fitted.png"),
       plot = p_scatter_2, width = 10, height = 10, dpi = DPI)
cat("  ✓ scatter_2_resid_vs_fitted.png\n")


# --- #3 Moran Scatter ---
cat("  [#3] Moran scatter ...\n")

make_moran_scatter <- function(df, xcol, wcol, title, I) {
  mx <- mean(df[[xcol]], na.rm = TRUE); mw <- mean(df[[wcol]], na.rm = TRUE)
  ggplot(df, aes(x = .data[[xcol]], y = .data[[wcol]])) +
    geom_hline(yintercept = mw, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = mx, linetype = "dotted", color = "grey60") +
    geom_point(shape = 21, fill = "darkgreen", color = "grey25",
               stroke = 0.3, size = 1.7, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "red",
                linewidth = 0.7, formula = y ~ x) +
    labs(caption = paste0(title, "\n", sprintf("Moran's I (slope) = %+.3f", I)),
         x = "x", y = "W · x (spatial lag)") +
    scatter_theme
}

panels_3 <- list()
panels_3$y1 <- make_moran_scatter(df_sc, "y1", "Wy1",
                                  "y1 (obs, 2020)",
                                  moran_results$y1_t2$I)
panels_3$y2 <- make_moran_scatter(df_sc, "y2", "Wy2",
                                  "y2 (obs, 2020)",
                                  moran_results$y2_t2$I)
for (mid in names(model_outputs)) {
  panels_3[[paste0("e1_", mid)]] <- make_moran_scatter(
    df_sc, paste0("res_y1_", mid), paste0("Wres_y1_", mid),
    sprintf("e1 (%s)", mid),
    moran_results[[paste0("resid_y1_", mid)]]$I)
  panels_3[[paste0("e2_", mid)]] <- make_moran_scatter(
    df_sc, paste0("res_y2_", mid), paste0("Wres_y2_", mid),
    sprintf("e2 (%s)", mid),
    moran_results[[paste0("resid_y2_", mid)]]$I)
}
p_scatter_3 <- wrap_plots(panels_3, ncol = 2)
ggsave(file.path(OUTPUT_DIR, "scatter_3_moran.png"),
       plot = p_scatter_3, width = 10, height = 13, dpi = DPI)
cat("  ✓ scatter_3_moran.png\n")


# --- #4 Residual vs Residual (cross-model) ---
cat("  [#4] Residual vs Residual (cross-model) ...\n")

if (length(names(model_outputs)) >= 2) {
  mid1 <- names(model_outputs)[1]; mid2 <- names(model_outputs)[2]

  cat("\n  --- Correlation coefficients (#4 Resid cross-model) ---\n")

  make_rr <- function(df, col1, col2, title, console_label = NULL) {
    r <- cor(df[[col1]], df[[col2]], use = "complete.obs")
    if (!is.null(console_label)) {
      cat(sprintf("    %-32s r = %+.4f   R² = %.4f\n",
                  console_label, r, r^2))
    }
    rng <- range(c(df[[col1]], df[[col2]]), na.rm = TRUE)
    ggplot(df, aes(x = .data[[col1]], y = .data[[col2]])) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      geom_point(shape = 21, fill = "purple", color = "grey25",
                 stroke = 0.3, size = 1.7, alpha = 0.85) +
      coord_fixed(xlim = rng, ylim = rng) +
      labs(caption = paste0(title, "\n", sprintf("r = %.3f", r)),
           x = sprintf("Residual (%s)", mid1),
           y = sprintf("Residual (%s)", mid2)) +
      scatter_theme
  }

  p4_y1 <- make_rr(df_sc, paste0("res_y1_", mid1), paste0("res_y1_", mid2),
                   "y1: Residual comparison",
                   console_label = sprintf("e1 [%s] vs e1 [%s]", mid1, mid2))
  p4_y2 <- make_rr(df_sc, paste0("res_y2_", mid1), paste0("res_y2_", mid2),
                   "y2: Residual comparison",
                   console_label = sprintf("e2 [%s] vs e2 [%s]", mid1, mid2))

  p_scatter_4 <- (p4_y1 | p4_y2)

  ggsave(file.path(OUTPUT_DIR, "scatter_4_resid_vs_resid.png"),
         plot = p_scatter_4, width = 11, height = 5.5, dpi = DPI)
  cat("  ✓ scatter_4_resid_vs_resid.png\n")
}


# --- #6 Between-equation residuals: e1 vs e2 ---
cat("  [#6] e1 vs e2 (between-equation) ...\n")
cat("\n  --- Correlation coefficients (#6 e1 vs e2 between-equation) ---\n")

make_e12 <- function(df, c1, c2, title, console_label = NULL) {
  r <- cor(df[[c1]], df[[c2]], use = "complete.obs")
  if (!is.null(console_label)) {
    cat(sprintf("    %-32s r = %+.4f   R² = %.4f\n",
                console_label, r, r^2))
  }
  ggplot(df, aes(x = .data[[c1]], y = .data[[c2]])) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    geom_point(shape = 21, fill = "orange3", color = "grey25",
               stroke = 0.3, size = 1.7, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, color = "red",
                linewidth = 0.7, formula = y ~ x) +
    labs(caption = paste0(title, "\n", sprintf("r(e1, e2) = %.3f", r)),
         x = "e1 (Aging Rate 65+)", y = "e2 (Secondary Industry)") +
    scatter_theme
}

panels_6 <- list()
for (mid in names(model_outputs)) {
  panels_6[[mid]] <- make_e12(df_sc, paste0("res_y1_", mid), paste0("res_y2_", mid),
                              sprintf("Model %s", mid),
                              console_label = sprintf("e1 vs e2 [Model %s]", mid))
}
p_scatter_6 <- wrap_plots(panels_6, nrow = 1)

ggsave(file.path(OUTPUT_DIR, "scatter_6_e1_vs_e2.png"),
       plot = p_scatter_6, width = 11, height = 5.5, dpi = DPI)
cat("  ✓ scatter_6_e1_vs_e2.png\n")


# §12. Completion Report

cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("Done: ", OUTPUT_DIR, "/\n", sep = "")
cat(paste(rep("=", 60), collapse = ""), "\n")

for (f in list.files(OUTPUT_DIR, full.names = TRUE)) {
  cat(sprintf("  %-40s  %.1f KB\n", basename(f), file.info(f)$size / 1024))
}

cat("\nMoran's I summary (reprinted):\n")
for (nm in names(moran_results)) {
  mi <- moran_results[[nm]]
  cat(sprintf("  %-25s  I = %+.4f   SD = %.4f   z = %+.3f   p = %.3g\n",
              nm, mi$I, mi$sd, mi$z, mi$p))
}
