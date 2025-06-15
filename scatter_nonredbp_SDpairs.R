# scatter_nonredbp_SD pairs.R
# Creates a scatter plot of the number of nonredundant basepairs and the number of SD pairs for Pacific islander and HPRC assemblies
# option to stratify by pacific islander subpopulations-- 

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

#hprc <- filter(hprc, !Sample %in% c("hg38", "T2T"))
#pi   <- filter(pi,   !Sample %in% c("hg38", "T2T"))

hprc <- mutate(hprc, Group = Superpopulation, Island = "-")
pi   <- mutate(pi,   Group = "PI",              Island = Superpopulation)

df <- bind_rows(hprc, pi) %>%
  mutate(
    SD_pairs = as.numeric(SD_pairs),
    Nonredundant_bp = as.numeric(Nonredundant_bp),
    Nonredundant_Mb = Nonredundant_bp / 1e6
  )	
# colour palettes ------------------------------------------------------
hprc_groups  <- sort(unique(hprc$Group))
superpopulation_cols <- c(
  #"PI"		= "gray25",	
  "AFR"		= "#FFB366",  #orange
  "AMR"		= "red",			# "#FFE066",  #yellow
  "EAS" 	= "green",				#"#B3FF66",  #green
 # "EUR"	= "#66FFB3",  #aquamarine
  "EUR"         = "purple",
  "SAS"		= "#66B3FF"   #blue
)

pi_islands   <- sort(unique(pi$Island))
island_cols <- c(
  "Samoa"       = "black",  # bright red
  "Fiji"        = "cadetblue4",  # clear blue
  "Philippines" = "lightpink",  # green
  "Guam"        = "violetred4",  # orange
  "Tonga"       = "yellowgreen",  # purple
  "Marshall"    = "deeppink",  # yellow
  "Pohnpei"     = "darkgoldenrod2",  # brown
  "Tahiti"      = "cyan1"   # grey
)
# Assign color mapping by group
# for colorful PI subpopulations
df <- df %>%
  mutate(ColorGroup = ifelse(Group == "PI", Island, Group))
color_map <- c(superpopulation_cols, island_cols)

#for all gray PI superpopulation
#df <- df %>%
#  mutate(ColorGroup = ifelse(Group == "PI", "PI", Group))
#color_map <- c(superpopulation_cols)

# Calculate correlation
# correlation <- cor(df$Nonredundant_Mb, df$SD_pairs)

#calculate correlation for PI and hprc separately
pi_only <- filter(df, Group == "PI")
correlation <- cor(pi_only$Nonredundant_Mb, pi_only$SD_pairs, use = "complete.obs")

hprc_only <-filter(df, Group != "PI")
hprc_correlation <- cor(hprc_only$Nonredundant_Mb, hprc_only$SD_pairs, use = "complete.obs")


# Create scatter plot
p <- ggplot(df, aes(x = Nonredundant_Mb, y = SD_pairs, color = ColorGroup)) +
  geom_point(alpha = 0.85, size = 2.6) +
  scale_color_manual(values = color_map, name = "Superpopulation") +
  labs(
    x = expression(paste("Nonredundant SD bases (×10"^6, " bp)")),
    y = "Number of SD pairs",
    title = "Nonredundant SD bases vs. SD pairs\nin Pacific Islander and HPRC assemblies"
  ) +
  annotate("text",
           x = min(df$Nonredundant_Mb, na.rm = TRUE) + 2,
           y = max(df$SD_pairs, na.rm = TRUE) - 500,
           label = paste0("Pacific Islander r = ", round(correlation, 3)),
           size = 4,
           hjust = 0,
           color = "black") +
  annotate("text",
         x = min(df$Nonredundant_Mb, na.rm = TRUE) + 2,
         y = max(df$SD_pairs, na.rm = TRUE) +100,
         label = paste0("HPRC r = ", round(hprc_correlation, 3)),
         size = 4,
         hjust = 0,
         color = "black") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 14)
  )
#colorul subpopulations:
ggsave("scatter_nonredbp_SDpairs_withsubpop.png", plot = p, width = 8, height = 6, dpi = 300)

#only see superpopulation (no subpopulation stratification)
#ggsave("scatter_nonredbp_SDpairs.png", plot = p, width = 8, height = 6, dpi = 300)

