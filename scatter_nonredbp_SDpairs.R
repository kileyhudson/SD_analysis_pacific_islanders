# scatter_nonredbp_SD pairs.R
# Creates a scatter plot of the number of nonredundant basepairs and the number of SD pairs for Pacific islander and HPRC assemblies
# HPRC in black, PI stratified by island subpopulations

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
cols_needed <- c(1, 3, 4, 5)          # Sample | SD pair counts | Nonredundant bp | Superpopulation/subpopulation
hprc <- read_tsv("HPRC_SD_info.tsv",
                 col_select = all_of(cols_needed),
                 show_col_types = FALSE)
pi   <- read_tsv("PI_SD_info.tsv",
                 col_select = all_of(cols_needed),
                 show_col_types = FALSE)
colnames(hprc) <- colnames(pi) <- c("Sample", "SD_pairs", "Nonredundant_bp", "Superpopulation")

# Process data
hprc <- mutate(hprc, Group = "HPRC", Island = "-")
pi   <- mutate(pi,   Group = "PI",   Island = Superpopulation)

df <- bind_rows(hprc, pi) %>%
  mutate(
    SD_pairs = as.numeric(SD_pairs),
    Nonredundant_bp = as.numeric(Nonredundant_bp),
    Nonredundant_Mb = Nonredundant_bp / 1e6
  )

# colour palettes ------------------------------------------------------
# Same island colors as in the previous script
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

# Only keep colors for islands that exist in the data
island_cols <- island_cols[names(island_cols) %in% pi_islands]

# Assign color mapping by group
# HPRC gets "HPRC" label, PI gets island names
df <- df %>%
  mutate(ColorGroup = ifelse(Group == "PI", Island, Group))

# Create color map with black for HPRC and island colors for PI
color_map <- c("HPRC" = "black", island_cols)

# Calculate correlations separately for PI and HPRC
pi_only <- filter(df, Group == "PI")
pi_correlation <- cor(pi_only$Nonredundant_Mb, pi_only$SD_pairs, use = "complete.obs")

hprc_only <- filter(df, Group != "PI")
hprc_correlation <- cor(hprc_only$Nonredundant_Mb, hprc_only$SD_pairs, use = "complete.obs")

# Create scatter plot
p <- ggplot(df, aes(x = Nonredundant_Mb, y = SD_pairs, color = ColorGroup)) +
  geom_point(alpha = 0.85, size = 2.6) +
  scale_color_manual(values = color_map, 
                     name = "Population",
                     breaks = c("HPRC", names(island_cols))) +  # Order legend items
  labs(
    x = expression(paste("Nonredundant SD bases (×10"^6, " bp)")),
    y = "Number of SD pairs",
    title = "Nonredundant SD bases vs. SD pairs\nin Pacific Islander and HPRC assemblies"
  ) +
  annotate("text",
           x = min(df$Nonredundant_Mb, na.rm = TRUE) + 2,
           y = max(df$SD_pairs, na.rm = TRUE) - 500,
           label = paste0("Pacific Islander r = ", round(pi_correlation, 3)),
           size = 5,
           hjust = 0,
           color = "black") +
  annotate("text",
           x = min(df$Nonredundant_Mb, na.rm = TRUE) + 2,
           y = max(df$SD_pairs, na.rm = TRUE) +250,
           label = paste0("HPRC r = ", round(hprc_correlation, 3)),
           size = 5,
           hjust = 0,
           color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))  # Make legend points larger

# Save plot
ggsave("scatter_nonredbp_SDpairs_withsubpop1.png", plot = p, width = 9, height = 6, dpi = 300)
