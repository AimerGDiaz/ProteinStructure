#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
data_dir <- file.path(root, "Data")
figure_dir <- file.path(root, "Figures", "Reproduced")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_rds <- function(name) readRDS(file.path(data_dir, name))

pae <- bind_rows(
  read_rds("P6_dsRNA_DF_ALL.rds") %>% mutate(Complex = "P6-dsRNA"),
  read_rds("P6_TOR_DF_ALL.rds") %>% mutate(Complex = "P6-TOR")
) %>%
  mutate(Model = factor(Model))

pae_summary <- pae %>%
  group_by(Complex, Interaction, Residue, Type) %>%
  summarise(
    mean_value = mean(MinPAE, na.rm = TRUE),
    se_value = sd(MinPAE, na.rm = TRUE) / sqrt(sum(!is.na(MinPAE))),
    .groups = "drop"
  )

pae_plot <- ggplot(
  pae_summary,
  aes(Residue, mean_value, colour = Interaction, group = Interaction)
) +
  geom_ribbon(
    aes(ymin = mean_value - se_value, ymax = mean_value + se_value,
        fill = Interaction),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.6) +
  facet_grid(Complex + Type ~ ., scales = "free_y") +
  labs(x = "P6 residue", y = "Minimum predicted aligned error (A)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "P6_interface_PAE.svg"),
  pae_plot,
  width = 11,
  height = 9
)

contact <- read_tsv(
  file.path(root, "Results", "ContactProb_Integration_Summary.txt"),
  show_col_types = FALSE
)

contact_plot <- contact %>%
  filter(Contactprobability > 0) %>%
  ggplot(aes(Residue, Contactprobability, fill = Interaction)) +
  geom_col(width = 1) +
  facet_grid(Type ~ ., scales = "free_y") +
  labs(x = "P6 residue", y = "Maximum contact probability") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "P6_interface_contact_probability.svg"),
  contact_plot,
  width = 11,
  height = 9
)

write_tsv(pae_summary, file.path(root, "Results", "PAE_summary_reproduced.tsv"))

