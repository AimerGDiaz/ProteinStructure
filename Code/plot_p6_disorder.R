#!/usr/bin/env Rscript
# Run from the repository root, or source from Readme.Rmd.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(dplyr)
  library(tidyr)
  library(readr)
})
dir.create("Results/Disorder", recursive = TRUE, showWarnings = FALSE)
dir.create("Figures/Disorder", recursive = TRUE, showWarnings = FALSE)
idr <- read_tsv("Data/af_disorder_AllP6s_pred.tsv", show_col_types = FALSE)
required <- c("name", "pos", "aa", "lddt", "disorder", "rsa",
              "disorder-25", "binding-25-0.581")
stopifnot(all(required %in% names(idr)))
idr <- idr %>% extract(name, c("condition", "seed", "model", "chain"),
  "^(.*)_seed_([0-9]+)_sample_([0-4])_chain_([AB])$", remove = FALSE,
  convert = TRUE)
labels <- c(P6_Monomer = "P6 monomer", P6_Dimer = "P6 dimer",
            P6_dimer_dsRNA = "P6-dsRNA", P6_dimer_TOR = "P6-TOR",
            P6_dimer_ivIDR_dsRNA = "P6IDR-dsRNA")
stopifnot(!anyNA(idr), all(idr$condition %in% names(labels)),
          !anyDuplicated(idr[c("name", "pos")]))
idr$condition <- factor(unname(labels[idr$condition]), levels = unname(labels))
coverage <- idr %>% group_by(condition, seed, model, chain, name) %>%
  summarise(n_residues = n(), first = min(pos), last = max(pos), .groups = "drop")
stopifnot(nrow(coverage) == 45L, all(coverage$first == 1),
          all(coverage$n_residues == coverage$last),
          all(coverage$n_residues == ifelse(coverage$condition == "P6IDR-dsRNA", 509, 520)))
for (condition_label in levels(idr$condition)) {
  expected <- expand.grid(model = 0:4,
    chain = if (condition_label == "P6 monomer") "A" else c("A", "B"))
  observed <- filter(coverage, .data$condition == condition_label)
  stopifnot(setequal(paste(expected$model, expected$chain),
                     paste(observed$model, observed$chain)))
}
stopifnot(all(vapply(idr[required[4:8]], function(x)
  all(is.finite(x) & x >= 0 & x <= 1), logical(1))),
  max(abs(idr$disorder + idr$lddt - 1)) < 0.002)
long <- idr %>% pivot_longer(c(`disorder-25`, `binding-25-0.581`, disorder),
  names_to = "metric", values_to = "score") %>%
  mutate(metric = factor(metric, levels = c("disorder-25", "binding-25-0.581", "disorder"),
    labels = c("Disorder score (RSA-25)", "IDR-binding score", "1 - pLDDT")))
# Display the 13-aa replacement at 224-236 and pad the remaining WT span.
# Preserve original positions and flag synthetic zeros; never average them
# into whole-chain scores.
del_start <- 237L
del_end <- 247L
del_len <- del_end - del_start + 1L  # 11 residues
long <- long %>% mutate(construct_pos = pos, is_placeholder = FALSE,
  pos = if_else(condition == "P6IDR-dsRNA" & pos >= del_start,
                pos + del_len, pos))
deleted_rows <- long %>% filter(condition == "P6IDR-dsRNA") %>%
  distinct(name, condition, seed, model, chain, metric) %>%
  crossing(pos = del_start:del_end) %>%
  mutate(construct_pos = NA_integer_, is_placeholder = TRUE,
         aa = "", lddt = 0, disorder = 0, rsa = 0, ss = "", score = 0)
long <- bind_rows(long, deleted_rows) %>% arrange(condition, model, chain, metric, pos)
stopifnot(nrow(deleted_rows) == 5L * 2L * 3L * del_len,
          !anyDuplicated(long[c("name", "metric", "pos")]),
          all((long %>% count(name, metric))$n == 520L),
          all(long$pos[long$is_placeholder] %in% del_start:del_end))
profiles <- long %>% group_by(condition, chain, pos, metric, is_placeholder) %>%
  summarise(n_models = n(), mean = mean(score), sd = sd(score),
            minimum = min(score), maximum = max(score), .groups = "drop")
stopifnot(all(profiles$n_models == 5))
model_summary <- long %>% filter(!is_placeholder) %>%
  group_by(condition, model, chain, metric) %>%
  summarise(mean_score = mean(score), .groups = "drop")
# Pair A/B within each model before summarising the five model predictions.
ensemble <- model_summary %>% group_by(condition, model, metric) %>%
  summarise(mean_score = mean(mean_score), .groups = "drop") %>%
  group_by(condition, metric) %>% summarise(n_models = n(),
    mean = mean(mean_score), sd = sd(mean_score), .groups = "drop")
write_csv(coverage, "Results/Disorder/coverage.csv")
write_csv(profiles, "Results/Disorder/residue_summary.csv")
write_csv(long, "Results/Disorder/aligned_model_chain_scores.csv")
write_csv(model_summary, "Results/Disorder/model_chain_summary.csv")
write_csv(ensemble, "Results/Disorder/ensemble_summary.csv")
prism_style <- function() theme_prism(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
for (metric_label in levels(long$metric)) {
  p <- ggplot(filter(profiles, metric == metric_label),
              aes(pos, mean, colour = chain, fill = chain)) +
    geom_ribbon(aes(ymin = minimum, ymax = maximum), alpha = 0.14, colour = NA) +
    geom_line(data = filter(long, metric == metric_label),
      aes(y = score, group = interaction(model, chain)), linewidth = 0.25, alpha = 0.22) +
    geom_line(linewidth = 0.75) + facet_wrap(~condition, ncol = 2) +
    scale_colour_manual(values = c(A = "#0072B2", B = "#D55E00")) +
    scale_fill_manual(values = c(A = "#0072B2", B = "#D55E00")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = metric_label, x = "Residue (WT-aligned numbering)", y = "Score",
      colour = "P6 chain", fill = "P6 chain",
      caption = paste("Thin lines: models 0-4; thick lines: mean; shading: model range.",
        "P6IDR positions 237-247: inserted zero placeholders.", sep = "\n")) + prism_style()
  stem <- c("disorder_score", "idr_binding_score", "one_minus_plddt")[match(metric_label, levels(long$metric))]
  ggsave(paste0("Figures/Disorder/", stem, ".svg"), p, width = 13, height = 11)
  ggsave(paste0("Figures/Disorder/", stem, ".png"), p, width = 13, height = 11, dpi = 180)
}
# Original IntrinsicallyDisorderRegions/Readme.Rmd condition colours.
condition_colours <- c("P6 monomer" = "#E67300", "P6 dimer" = "#E64000",
  "P6-dsRNA" = "#28CC48", "P6-TOR" = "#005770", "P6IDR-dsRNA" = "blue")
comparison_long <- long %>% filter(metric != "1 - pLDDT")
comparison_profiles <- profiles %>% filter(metric != "1 - pLDDT")
comparison <- ggplot(comparison_profiles, aes(pos, mean, colour = condition)) +
  annotate("rect", xmin = 224, xmax = 247, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "#63666A") +
  annotate("rect", xmin = 486, xmax = 520, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "#63666A") +
  geom_line(data = comparison_long,
    aes(y = score, group = interaction(condition, model, chain)),
    linewidth = 0.25, alpha = 0.18) +
  geom_line(linewidth = 0.85) +
  facet_grid(metric ~ chain, labeller = labeller(chain = c(A = "Monomer A", B = "Monomer B"))) +
  scale_colour_manual(values = condition_colours, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "P6 disorder and IDR-binding comparison",
       x = "P6 residues (WT-aligned numbering)", y = "Score", colour = "P6 condition",
       caption = paste("Thin lines: all five models; thick lines: five-model mean.",
         "P6IDR 237-247: zero placeholders. Monomer-only condition has chain A only.", sep = "\n")) +
  prism_style() + guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1)))
ggsave("Figures/Disorder/all_p6_by_monomer.svg", comparison, width = 14, height = 10)
ggsave("Figures/Disorder/all_p6_by_monomer.png", comparison, width = 14, height = 10, dpi = 180)
writeLines(capture.output(sessionInfo()), "Results/Disorder/sessionInfo.txt")
