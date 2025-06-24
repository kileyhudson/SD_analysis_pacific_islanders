# violin_pi_filtering_comparison.R
# Plot PI nonredundant bps across different filtering conditions (unfiltered, filter SDs overlapping structural errors, structural+small (all errors)
# doesn't include HPRC

pkgs <- c("ggplot2", "dplyr", "readr", "scales")
to_get <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

# read data ------------------------------------------------------------
structural <- read_tsv("filter_SDs_by_structural_errors.tsv", show_col_types = FALSE)
all_errors <- read_tsv("filter_SDs_by_all_errors.tsv", show_col_types = FALSE)
ancestry_map <- read_tsv("PI_ancestry_map.tsv", col_names = c("Sample", "Subpopulation"), show_col_types = FALSE)

# prepare data for each filtering condition ----------------------------
before <- structural %>%
  select(Sample, Haplotype, Nonredundant_bp_before_filtering) %>%
  mutate(Condition = "unfiltered",
         Nonredundant_bp = Nonredundant_bp_before_filtering) %>%
  select(-Nonredundant_bp_before_filtering)

struct_only <- structural %>%
  select(Sample, Haplotype, nonredundat_bp_after_filtering) %>%
  mutate(Condition = "structural error filter",
         Nonredundant_bp = nonredundat_bp_after_filtering) %>%
  select(-nonredundat_bp_after_filtering)

both <- all_errors %>%
  select(Sample, Haplotype, nonredundat_bp_after_filtering) %>%
  mutate(Condition = "all errors filter",
         Nonredundant_bp = nonredundat_bp_after_filtering) %>%
  select(-nonredundat_bp_after_filtering)

# combine and add subpopulation info -----------------------------------
df <- bind_rows(before, struct_only, both) %>%
  left_join(ancestry_map, by = "Sample") %>%
  mutate(Nonredundant_Mb = Nonredundant_bp / 1e6,
         Condition = factor(Condition, levels = c("unfiltered", "structural error filter", "all errors filter")))

# color palette -------------------------------------------------------
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
p <- ggplot(df, aes(x = Condition, y = Nonredundant_Mb)) +
  geom_violin(fill = "aliceblue", width = 0.9, alpha = 0.6, scale = "width", colour = "grey40", adjust = 2) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(colour = Subpopulation), 
             position = position_dodge(width = 0.5), 
             size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = island_cols, name = "Subpopulation") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(title = "Nonredundant bps across different filtering conditions",
       x=NULL, y = "Number of SD nonredundant bases (Mb)") +
  theme_classic(base_size = 15)

# save -----------------------------------------------------------------
ggsave("violin_error_filtering_comparison.png", p, width = 7.5, height = 5, dpi = 300)
