# scatter_total_vs_filtered_bp.R
# Scatter plot: total bp vs. bp after filtering
# HPRC points coloured by super-population; PI points in grey.

# auto-install minimal deps --------------------------------------------
pkgs <- c("ggplot2", "dplyr", "readr", "scales")
to_add <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_add)) install.packages(to_add, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

# read data -------------------------------------------------------------
cols <- c("Sample", "Haplotype", "Total_Contigs", "Pct_Contigs_Under1MB",
          "N50", "Superpopulation", "total_bp", "bp_after_filtering")

hprc <- read_tsv("HPRC_contig_length_info1.tsv",
                 col_names = cols,
                 col_types = cols(
                   .default = "c",
                   Total_Contigs = "d",
                   Pct_Contigs_Under1MB = "d",
                   N50 = "d",
                   total_bp = "c",
                   bp_after_filtering = "c"
                 ),
                 show_col_types = FALSE)

pi <- read_tsv("PI_contig_length_info1.tsv",
               col_names = cols,
               col_types = cols(
                 .default = "c",
                 Total_Contigs = "d",
                 Pct_Contigs_Under1MB = "d",
                 N50 = "d",
                 total_bp = "c",
                 bp_after_filtering = "c"
               ),
               show_col_types = FALSE)

# drop hg38 / T2T -------------------------------------------------------
hprc <- filter(hprc, !Sample %in% c("hg38", "T2T"))
pi   <- filter(pi,   !Sample %in% c("hg38", "T2T"))

# parse numbers (strip commas) -----------------------------------------
parse_bp <- function(x) {
  # Handle NA and non-numeric values
  result <- suppressWarnings(parse_number(x, locale = locale(grouping_mark = ",")))
  # Replace any remaining parsing failures with NA
  result[is.na(result) & !is.na(x) & x != ""] <- NA
  return(result)
}

hprc <- mutate(hprc,
               total_bp = parse_bp(total_bp),
               bp_after_filtering = parse_bp(bp_after_filtering),
               Group = "HPRC",  # All HPRC samples get the same group
               Category = "HPRC")

pi <- mutate(pi,
             total_bp = parse_bp(total_bp),
             bp_after_filtering = parse_bp(bp_after_filtering),
             Group = Superpopulation,  # PI samples grouped by island
             Category = "PI")

# Remove rows with NA values in bp columns
hprc <- filter(hprc, !is.na(total_bp) & !is.na(bp_after_filtering))
pi <- filter(pi, !is.na(total_bp) & !is.na(bp_after_filtering))

df <- bind_rows(hprc, pi)

# colour palettes -------------------------------------------------------
# Get unique islands from PI data
pi_islands <- sort(unique(pi$Group))

# Define colors for each island
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

# Combine HPRC (gray) with island colors
all_cols <- c("HPRC" = "black", island_cols)  # Changed to a more visible gray

# Define what should appear in the legend
legend_order <- c("HPRC", names(island_cols))

# plot -----------------------------------------------------------------
p <- ggplot(df,
            aes(x = total_bp, y = bp_after_filtering,
                colour = Group)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey70") +
  geom_point(alpha = 0.7, size = 2.5) +  # Slightly larger and more opaque
  scale_colour_manual(values = all_cols, 
                      breaks = legend_order,  # Show HPRC and islands in legend
                      name = "Population") +
  scale_x_continuous(labels = label_number(scale = 1 / 1e9,
                                           accuracy = 0.01,
                                           suffix = " Gb")) +
  scale_y_continuous(labels = label_number(scale = 1 / 1e9,
                                           accuracy = 0.01,
                                           suffix = " Gb")) +
  labs(title = "Total base pairs vs. base pairs after filtering",
       subtitle = "HPRC and Pacific Ancestry Assemblies",
       x = "Total bp (Gb)",
       y = "After filtering bp (Gb)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

# save -----------------------------------------------------------------
ggsave("scatter_total_vs_filtered_bp1.png", p, width = 9, height = 6, dpi = 300)
