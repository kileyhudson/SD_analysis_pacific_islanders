# violin_SD_pair_counts.R
# Plot number of SD pairs for HPRC vs Pacific Islanders

# auto-install minimal deps -------------------------------------------
pkgs <- c("ggplot2", "dplyr", "readr", "scales", "viridisLite")
to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(viridisLite)
})

# read data ------------------------------------------------------------
cols_needed <- c(1, 3, 4)          # Sample | SD pair counts | Superpopulation/subpopulation
hprc <- read_tsv("HPRC_SD_counts.tsv",
                 col_select = all_of(cols_needed),
                 show_col_types = FALSE)
pi   <- read_tsv("PI_SD_counts.tsv",
                 col_select = all_of(cols_needed),
                 show_col_types = FALSE)
colnames(hprc) <- colnames(pi) <- c("Sample", "SD_pairs", "Superpopulation")

#hprc <- filter(hprc, !Sample %in% c("hg38", "T2T"))
#pi   <- filter(pi,   !Sample %in% c("hg38", "T2T"))

hprc <- mutate(hprc, Group = Superpopulation, Island = "-")
pi   <- mutate(pi,   Group = "PI",              Island = Superpopulation)

df <- bind_rows(hprc, pi) %>%
      mutate(SD_pairs = as.numeric(SD_pairs))

# colour palettes ------------------------------------------------------
hprc_groups  <- sort(unique(hprc$Group))
group_levels <- c("PI", hprc_groups)           # PI violin on the left

violin_cols  <- c("PI" = "aliceblue",
                  setNames(c("#FFB366", "#FFE066", "#B3FF66",
                             "#66FFB3", "#66B3FF")[seq_along(hprc_groups)],
                           hprc_groups))

pi_islands   <- sort(unique(pi$Island))
island_cols <- c(
  "Samoa"       = "#E41A1C",  # bright red
  "Fiji"        = "#377EB8",  # clear blue
  "Philippines" = "#4DAF4A",  # green
  "Guam"        = "#FF7F00",  # orange
  "Tonga"       = "#984EA3",  # purple
  "Marshall"    = "#FFD92F",  # yellow
  "Pohnpei"     = "#A65628",  # brown
  "Tahiti"      = "#999999"   # grey
)
# plot -----------------------------------------------------------------
p <- ggplot(df, aes(x = factor(Group, levels = group_levels),
                    y = SD_pairs)) +
  geom_violin(aes(fill = Group),
              width = 0.9, alpha = 0.6, colour = "grey40") +
  geom_jitter(data = hprc,
              width = 0.1, size = 1.2, colour = "black", alpha = 0.8) +
  geom_point(
  data     = pi,
  aes(colour = Island, group = Island),
  position = position_dodge(width = 0.5),  # small parallel lanes
  size     = 1.5,
  alpha    = 0.7
  )+
  scale_fill_manual(values = violin_cols, guide = "none") +
  scale_colour_manual(values = island_cols, name = "PI population") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +  # auto-scales
  labs(title = "Number of Segmental Duplication Pairs",
       subtitle = "Pacific Ancestry assemblies vs. HPRC",
       x = "Population group", y = "Number of SD pairs") +
  theme_classic(base_size = 14)

# save -----------------------------------------------------------------
ggsave("violin_SD_pair_counts.png", p, width = 8, height = 5, dpi = 300)
