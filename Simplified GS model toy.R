# ============================================================
# Gulf Sturgeon – Simple Age-Structured Model (SturMod-lite from FloPop/Brett)
# Purpose:
#   Small population (~ABCDE total fish) with SturMod life history.
#   Compare benefits to adult abundance (age ≥ 10) from:
#     (1) Increasing recruitment (R0 up)
#     (2) Decreasing adult mortality (Z down)
# Years: 2000–2075
# 100 simulations per scenario; plots show mean + 95% CI across sims.
# ============================================================

library(tidyverse)
library(cowplot)

set.seed(1)

# ----------------------------
# GLOBAL INPUTS
# ----------------------------
years <- 2000:2075
ny    <- length(years)

ages <- 1:50        # SturMod oldest modeled age-class A = 50
A    <- length(ages)

adult_age <- 10     # define "adult" as age >= 10

# Year where management change occurs
change_year <- 2020

# Number of simulations
n_sims <- 100

# ----------------------------
# GULF STURGEON LIFE HISTORY (from SturMod)
# ----------------------------

# Growth / length-at-age
Linf <- 220        # von Bertalanffy asymptotic length
K    <- 0.13       # von Bertalanffy K
t0   <- 0

Length <- Linf * (1 - exp(-K * (ages - t0)))   # length-at-age
aLW   <- 6.11e-6
bLW   <- 3
Wt    <- aLW * Length^bLW                     # weight-at-age

# Weight at maturity (mai = 8), as in SturMod
mai   <- 8
Wmat  <- aLW * (Linf * (1 - exp(-K * (mai - 1))))^bLW

# Fecundity-at-age
fec <- pmax(0, Wt - Wmat)

# Recruitment (Beverton–Holt)
recK   <- 2.9     # recruitment compensation ratio (SturMod-ish)
RecSD  <- 0.2     # lognormal recruitment SD

# Mortality (Z) – adult mortality at Linf
M_base <- 0.095   # Mad from SturMod (minimum adult mortality)

# Flatten Z for ages >= flat_age (Z similar beyond ~age 3–4)
flat_age   <- 4
Length_Z   <- Length
Length_Z[ages >= flat_age] <- Length[flat_age]

# ------------------------------------------------------------
# CALCULATE BASELINE R0 TO GET ~5000 TOTAL FISH AT EQUILIBRIUM
# ------------------------------------------------------------

Mu0  <- M_base * (Linf / Length_Z)
Sa0  <- exp(-Mu0)
la0  <- cumprod(c(1, Sa0[1:(A - 1)]))  # survivorship-at-age (age 1 starts at 1)

target_total <- 5000
Ro_base      <- target_total / sum(la0)

# Management options:
Ro_factor_up  <- 1.05   # 5% multiply R0 by this after change_year (recruitment better)
M_factor_down <- 0.95   # 5% improvement multiply M_base by this after change_year (Z lower)
M_low         <- M_base * M_factor_down

# ----------------------------
# CORE GULF STURGEON MODEL 
# ----------------------------
run_model_GS <- function(M_vec, Ro_vec, label) {
  
  N <- matrix(
    0,
    nrow = ny,
    ncol = A,
    dimnames = list(years, paste0("age", ages))
  )
  
  # Initial equilibrium-ish age structure under baseline M_base and Ro_vec[1]
  N[1, ] <- Ro_vec[1] * la0
  
  for (t in 2:ny) {
    
    # Age-specific mortality and survival for this year (flattened above flat_age)
    Mu_t <- M_vec[t] * (Linf / Length_Z)
    Sa_t <- exp(-Mu_t)
    
    # Eggs-per-recruit and BH parameters for this year's M and Ro
    la_t <- cumprod(c(1, Sa_t[1:(A - 1)]))
    epro <- sum(fec * la_t)                # eggs per recruit
    
    bha <- recK / epro                     # BH alpha
    bhb <- (recK - 1) / (Ro_vec[t] * epro) # BH beta (depends on Ro_t)
    
    # Eggs from previous year
    eggs_prev <- sum(fec * N[t - 1, ])
    
    # Lognormal recruitment anomaly
    rec_anom <- exp(rnorm(1, 0, RecSD))
    
    # Recruits to age 1
    N[t, 1] <- (bha * eggs_prev) / (1 + bhb * eggs_prev) * rec_anom
    
    # Survival transitions age a -> a+1
    for (a in 2:A) {
      N[t, a] <- N[t - 1, a - 1] * Sa_t[a - 1]
    }
  }
  
  as.data.frame(N) %>%
    rownames_to_column("year") %>%
    mutate(year = as.integer(year)) %>%
    pivot_longer(
      starts_with("age"),
      names_to  = "age",
      values_to = "N"
    ) %>%
    mutate(
      age_num  = as.numeric(gsub("age", "", age)),
      scenario = label
    )
}

# ----------------------------
# SIMULATE RUNS FOR ONE SCENARIO
# ----------------------------
simulate_scenario <- function(M_vec, Ro_vec, label) {
  out <- vector("list", n_sims)
  for (s in 1:n_sims) {
    out[[s]] <- run_model_GS(M_vec, Ro_vec, label) %>%
      mutate(sim = s)
  }
  bind_rows(out)
}

# ----------------------------
# SCENARIOS (focus on adult abundance, age >= adult_age)
# ----------------------------

# 1) Baseline – constant R0 and Z
M_1  <- rep(M_base, ny)
Ro_1 <- rep(Ro_base, ny)
df_1_sims <- simulate_scenario(M_1, Ro_1, "1) Baseline")

# 2) R0 ↑ (persistent) – R0 increases at change_year and stays high
Ro_2 <- ifelse(years >= change_year, Ro_base * Ro_factor_up, Ro_base)
M_2  <- M_1
df_2_sims <- simulate_scenario(M_2, Ro_2, "2) R0 ↑ (persistent)")

# 3) Z ↓ (persistent) – Z decreases at change_year and stays low
M_3  <- ifelse(years >= change_year, M_low, M_base)
Ro_3 <- Ro_1
df_3_sims <- simulate_scenario(M_3, Ro_3, "3) Z ↓ (persistent)")

# ----------------------------
# SUMMARISE ADULT / JUVENILE / TOTAL ABUNDANCE (BY SIM)
# ----------------------------

df_all_sims <- bind_rows(df_1_sims, df_2_sims, df_3_sims)

df_all_sims <- df_all_sims %>%
  mutate(is_adult = age_num >= adult_age)

summary_by_sim <- df_all_sims %>%
  group_by(scenario, sim, year) %>%
  summarise(
    N_total    = sum(N),
    N_adults   = sum(N[is_adult]),
    N_juvenile = sum(N[!is_adult]),
    .groups    = "drop"
  )

# Ensure scenarios are ordered and baseline is first
scenario_levels <- c(
  "1) Baseline",
  "2) R0 ↑ (persistent)",
  "3) Z ↓ (persistent)"
)

summary_by_sim <- summary_by_sim %>%
  mutate(scenario = factor(scenario, levels = scenario_levels, ordered = TRUE))

# ----------------------------
# MEAN TRAJECTORIES + 95% INTERVALS ACROSS SIMS
# ----------------------------

summary_stats <- summary_by_sim %>%
  group_by(scenario, year) %>%
  summarise(
    N_total_mean    = mean(N_total),
    N_total_lwr     = quantile(N_total, 0.025),
    N_total_upr     = quantile(N_total, 0.975),
    N_adults_mean   = mean(N_adults),
    N_adults_lwr    = quantile(N_adults, 0.025),
    N_adults_upr    = quantile(N_adults, 0.975),
    N_juvenile_mean = mean(N_juvenile),
    N_juvenile_lwr  = quantile(N_juvenile, 0.025),
    N_juvenile_upr  = quantile(N_juvenile, 0.975),
    .groups         = "drop"
  ) %>%
  mutate(scenario = factor(scenario, levels = scenario_levels, ordered = TRUE))

# Adults relative to baseline, per sim, then summarise
rel_by_sim <- summary_by_sim %>%
  group_by(year, sim) %>%
  mutate(
    rel_adults = N_adults / N_adults[scenario == "1) Baseline"]
  ) %>%
  ungroup()

summary_rel_stats <- rel_by_sim %>%
  group_by(scenario, year) %>%
  summarise(
    rel_adults_mean = mean(rel_adults),
    rel_adults_lwr  = quantile(rel_adults, 0.025),
    rel_adults_upr  = quantile(rel_adults, 0.975),
    .groups         = "drop"
  ) %>%
  mutate(scenario = factor(scenario, levels = scenario_levels, ordered = TRUE))

# ----------------------------
# COLOR PALETTE (colourblind-friendly; Baseline = black)
# ----------------------------

cb_palette <- c(
  "1) Baseline"           = "black",
  "2) R0 ↑ (persistent)"  = "#0072B2",  # blue
  "3) Z ↓ (persistent)"   = "#D55E00"   # orange
)

# ----------------------------
# Y-LIMITS FOR CONSISTENT AXES WITHIN EACH GRID
# ----------------------------

adults_ylim <- c(
  0,
  max(summary_stats$N_adults_upr, na.rm = TRUE)
)

juveniles_ylim <- c(
  0,
  max(summary_stats$N_juvenile_upr, na.rm = TRUE)
)

total_ylim <- c(
  0,
  max(summary_stats$N_total_upr, na.rm = TRUE)
)

rel_ylim <- c(
  min(summary_rel_stats$rel_adults_lwr, na.rm = TRUE),
  max(summary_rel_stats$rel_adults_upr, na.rm = TRUE)
)

# ============================================================
# PANEL HELPERS (3-panel style with dotted change_year line)
# ============================================================

# Adults: one-panel helper
make_adults_panel <- function(scen_label, panel_title) {
  df_plot <- summary_stats %>%
    filter(scenario == scen_label)
  
  this_col <- cb_palette[scen_label]
  
  ggplot(df_plot,
         aes(x = year, y = N_adults_mean)) +
    geom_ribbon(aes(ymin = N_adults_lwr, ymax = N_adults_upr),
                fill = this_col, alpha = 0.2, colour = NA) +
    geom_line(colour = this_col, linewidth = 1) +
    geom_vline(xintercept = change_year,
               linetype = "dotted",
               colour   = "black",
               linewidth = 0.7) +
    coord_cartesian(ylim = adults_ylim) +
    theme_minimal() +
    labs(
      title    = panel_title,
      subtitle = paste0("Adults (age ≥ ", adult_age, "); ", n_sims, " sims"),
      x        = "Year",
      y        = "Number of adults"
    )
}

# Juveniles: one-panel helper
make_juveniles_panel <- function(scen_label, panel_title) {
  df_plot <- summary_stats %>%
    filter(scenario == scen_label)
  
  this_col <- cb_palette[scen_label]
  
  ggplot(df_plot,
         aes(x = year, y = N_juvenile_mean)) +
    geom_ribbon(aes(ymin = N_juvenile_lwr, ymax = N_juvenile_upr),
                fill = this_col, alpha = 0.2, colour = NA) +
    geom_line(colour = this_col, linewidth = 1) +
    geom_vline(xintercept = change_year,
               linetype = "dotted",
               colour   = "black",
               linewidth = 0.7) +
    coord_cartesian(ylim = juveniles_ylim) +
    theme_minimal() +
    labs(
      title    = panel_title,
      subtitle = paste0("Juveniles (age < ", adult_age, "); ", n_sims, " sims"),
      x        = "Year",
      y        = "Number of juveniles"
    )
}

# Total: one-panel helper
make_total_panel <- function(scen_label, panel_title) {
  df_plot <- summary_stats %>%
    filter(scenario == scen_label)
  
  this_col <- cb_palette[scen_label]
  
  ggplot(df_plot,
         aes(x = year, y = N_total_mean)) +
    geom_ribbon(aes(ymin = N_total_lwr, ymax = N_total_upr),
                fill = this_col, alpha = 0.2, colour = NA) +
    geom_line(colour = this_col, linewidth = 1) +
    geom_vline(xintercept = change_year,
               linetype = "dotted",
               colour   = "black",
               linewidth = 0.7) +
    coord_cartesian(ylim = total_ylim) +
    theme_minimal() +
    labs(
      title    = panel_title,
      subtitle = paste0("Total abundance; ", n_sims, " sims"),
      x        = "Year",
      y        = "Number of fish"
    )
}

# Adults relative: one-panel helper
make_rel_panel <- function(scen_label, panel_title) {
  df_plot <- summary_rel_stats %>%
    filter(scenario == scen_label)
  
  this_col <- cb_palette[scen_label]
  
  ggplot(df_plot,
         aes(x = year, y = rel_adults_mean)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(aes(ymin = rel_adults_lwr, ymax = rel_adults_upr),
                fill = this_col, alpha = 0.2, colour = NA) +
    geom_line(colour = this_col, linewidth = 1) +
    geom_vline(xintercept = change_year,
               linetype = "dotted",
               colour   = "black",
               linewidth = 0.7) +
    coord_cartesian(ylim = rel_ylim) +
    theme_minimal() +
    labs(
      title    = panel_title,
      subtitle = paste0("Adults (age ≥ ", adult_age,
                        ") / Baseline; ", n_sims, " sims"),
      x        = "Year",
      y        = "Adult abundance / Baseline"
    )
}

# ============================================================
# 3-PANEL GRIDS (2x2 with bottom-left blank)
# ------------------------------------------------------------

# Adults grid
p_base_adults <- make_adults_panel("1) Baseline",          "Baseline")
p_R0_adults   <- make_adults_panel("2) R0 ↑ (persistent)", "R0 ↑ (persistent)")
p_Z_adults    <- make_adults_panel("3) Z ↓ (persistent)",  "Z ↓ (persistent)")
p_blank       <- ggplot() + theme_void()

p_adults_mean <- plot_grid(
  p_base_adults, p_R0_adults,
  p_blank,       p_Z_adults,
  nrow  = 2,
  ncol  = 2,
  labels = c("A", "B", "", "C"),
  label_size = 12
)

# Juveniles grid
p_base_juv <- make_juveniles_panel("1) Baseline",          "Baseline")
p_R0_juv   <- make_juveniles_panel("2) R0 ↑ (persistent)", "R0 ↑ (persistent)")
p_Z_juv    <- make_juveniles_panel("3) Z ↓ (persistent)",  "Z ↓ (persistent)")

p_juveniles_mean <- plot_grid(
  p_base_juv, p_R0_juv,
  p_blank,    p_Z_juv,
  nrow  = 2,
  ncol  = 2,
  labels = c("A", "B", "", "C"),
  label_size = 12
)

# Total grid
p_base_tot <- make_total_panel("1) Baseline",          "Baseline")
p_R0_tot   <- make_total_panel("2) R0 ↑ (persistent)", "R0 ↑ (persistent)")
p_Z_tot    <- make_total_panel("3) Z ↓ (persistent)",  "Z ↓ (persistent)")

p_total_mean <- plot_grid(
  p_base_tot, p_R0_tot,
  p_blank,    p_Z_tot,
  nrow  = 2,
  ncol  = 2,
  labels = c("A", "B", "", "C"),
  label_size = 12
)

# Adults relative grid
p_base_rel <- make_rel_panel("1) Baseline",          "Baseline")
p_R0_rel   <- make_rel_panel("2) R0 ↑ (persistent)", "R0 ↑ (persistent)")
p_Z_rel    <- make_rel_panel("3) Z ↓ (persistent)",  "Z ↓ (persistent)")

p_adults_rel_mean <- plot_grid(
  p_base_rel, p_R0_rel,
  p_blank,    p_Z_rel,
  nrow  = 2,
  ncol  = 2,
  labels = c("A", "B", "", "C"),
  label_size = 12
)

# ------------------------------------------------------------
# Print the four grids
# ------------------------------------------------------------
p_adults_mean
p_juveniles_mean
p_total_mean
p_adults_rel_mean
