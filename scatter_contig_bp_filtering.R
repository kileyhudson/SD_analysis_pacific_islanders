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
parse_bp <- function(x) parse_number(x, locale = locale(grouping_mark = ","))
hprc <- mutate(hprc,
               total_bp = parse_bp(total_bp),
               bp_after_filtering = parse_bp(bp_after_filtering),
               Group = Superpopulation,
               Category = "HPRC")

pi <- mutate(pi,
             total_bp = parse_bp(total_bp),
             bp_after_filtering = parse_bp(bp_after_filtering),
             Group = "PI",
             Category = "PI")

df <- bind_rows(hprc, pi)

# colour palettes -------------------------------------------------------
hprc_groups <- sort(unique(hprc$Group))
hprc_cols <- setNames(c("red", "orange", "green",
                        "#66FFB3", "purple")[seq_along(hprc_groups)],
                      hprc_groups)
all_cols <- c(hprc_cols, PI = "#555555")

# plot -----------------------------------------------------------------
wanted <- c("AFR", "AMR", "EAS", "EUR", "SAS", "PI")
p <- ggplot(df,
            aes(x = total_bp, y = bp_after_filtering,
                colour = Group)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey70") +
  geom_point(alpha = 0.6, size = 2) +
  scale_colour_manual(values = all_cols, breaks= wanted,
                     name = "Populations") +
  scale_x_continuous(labels = label_number(scale = 1 / 1e9,
                                           accuracy = 0.01,
                                           suffix = " Gb")) +
  scale_y_continuous(labels = label_number(scale = 1 / 1e9,
                                           accuracy = 0.01,
                                           suffix = " Gb")) +
  labs(title = "Total base pairs vs. base pairs after filtering",
       subtitle= "HPRC and Pacific Ancestry Assemblies",
	x = "Total bp (Gb)",
       y = "After filtering bp (Gb)") +
  theme_classic(base_size = 14)

# save -----------------------------------------------------------------
ggsave("scatter_total_vs_filtered_bp.png", p, width = 8, height = 5, dpi = 300)
