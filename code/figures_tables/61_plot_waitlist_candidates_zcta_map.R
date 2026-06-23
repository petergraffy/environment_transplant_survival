#!/usr/bin/env Rscript

runtime_source <- file.path("code", "r_runtime.R")
if (!file.exists(runtime_source)) runtime_source <- "r_runtime.R"
if (file.exists(runtime_source)) {
  source(runtime_source)
  ensure_user_library()
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(sf)
  library(tibble)
  library(viridis)
})

sf_use_s2(FALSE)

analysis_path <- file.path(
  "output",
  "primary_waitlist_period_pollution_cox",
  "primary_waitlist_period_pollution_analysis_dataset.csv.gz"
)
zcta_path <- file.path("data", "cache", "cb_2020_us_zcta520_500k", "cb_2020_us_zcta520_500k.shp")
fig_dir <- file.path("output", "figures", "waitlist_candidates_zcta_map")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_path)) stop("Missing analysis dataset: ", analysis_path, call. = FALSE)
if (!file.exists(zcta_path)) stop("Missing cached ZCTA shapefile: ", zcta_path, call. = FALSE)

clean_zip <- function(x) sprintf("%05s", as.character(x))

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

log_msg("Reading candidate ZCTA counts")
candidate_counts <- read_csv(
  analysis_path,
  col_select = c(PERS_ID, candidate_zip),
  col_types = cols(PERS_ID = col_character(), candidate_zip = col_character())
) %>%
  filter(!is.na(candidate_zip), nzchar(candidate_zip)) %>%
  mutate(zip = clean_zip(candidate_zip)) %>%
  group_by(zip) %>%
  summarise(
    waitlist_candidates = n_distinct(PERS_ID),
    waitlist_episodes = n(),
    .groups = "drop"
  )

count_breaks <- c(-Inf, 49, 99, 249, 499, 999, 2499, Inf)
count_labels <- c("<50", "50-99", "100-249", "250-499", "500-999", "1,000-2,499", "≥2,500")

candidate_counts <- candidate_counts %>%
  mutate(
    count_bin = cut(
      waitlist_candidates,
      breaks = count_breaks,
      labels = count_labels,
      include.lowest = TRUE,
      right = TRUE
    )
  )

write_csv(candidate_counts, file.path(fig_dir, "waitlist_candidates_by_zcta_counts.csv"))

summary <- candidate_counts %>%
  summarise(
    zctas_with_candidates = n(),
    total_distinct_candidates_across_zctas = sum(waitlist_candidates),
    total_waitlist_episodes = sum(waitlist_episodes),
    min_candidates = min(waitlist_candidates),
    median_candidates = median(waitlist_candidates),
    max_candidates = max(waitlist_candidates)
  )
write_csv(summary, file.path(fig_dir, "waitlist_candidates_by_zcta_summary.csv"))

log_msg("Reading CONUS ZCTA boundaries")
zcta <- st_read(zcta_path, quiet = TRUE)[, c("ZCTA5CE20", "geometry")]
names(zcta)[1] <- "zip"
zcta$zip <- clean_zip(zcta$zip)
conus_bbox <- st_as_sfc(st_bbox(c(xmin = -125, ymin = 24, xmax = -66, ymax = 50), crs = 4326))
zcta <- suppressWarnings(st_crop(zcta, st_transform(conus_bbox, st_crs(zcta))))
zcta <- st_transform(zcta, 5070)

map_data <- merge(zcta, candidate_counts, by = "zip", all.x = TRUE)
map_data$count_bin <- factor(map_data$count_bin, levels = count_labels)
legend_breaks <- count_labels[count_labels %in% as.character(stats::na.omit(unique(candidate_counts$count_bin)))]

palette <- c(
  "<50" = "#d7e8f4",
  "50-99" = "#a6bddb",
  "100-249" = "#67a9cf",
  "250-499" = "#2166ac",
  "500-999" = "#762a83",
  "1,000-2,499" = "#4d004b",
  "≥2,500" = "#250026"
)

p <- ggplot(map_data) +
  geom_sf(aes(fill = count_bin), color = NA) +
  scale_fill_manual(
    values = palette,
    breaks = legend_breaks,
    na.value = "#f8f8f8",
    drop = TRUE,
    name = "Candidates",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(20, "pt"),
      keyheight = unit(9, "pt")
    )
  ) +
  coord_sf(datum = NA) +
  labs(title = "Waitlist Candidates by ZCTA") +
  theme_void(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

png_path <- file.path(fig_dir, "waitlist_candidates_by_zcta_map.png")
pdf_path <- file.path(fig_dir, "waitlist_candidates_by_zcta_map.pdf")
ggsave(png_path, p, width = 9.5, height = 6.2, dpi = 600, bg = "white")
ggsave(pdf_path, p, width = 9.5, height = 6.2, device = grDevices::pdf, bg = "white")

write_csv(
  tibble(figure_path = c(png_path, pdf_path)),
  file.path(fig_dir, "waitlist_candidates_by_zcta_map_manifest.csv")
)

log_msg("Wrote waitlist candidate ZCTA map to ", normalizePath(fig_dir, winslash = "/"))
