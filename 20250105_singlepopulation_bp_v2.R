############################################################
##GULF STURGEON PVA CALCULATOR 
##Single-population PVA (self-contained within inputs in the code here)
##
## Endpoint:
##   Extirpation = min(adults) <= threshold at any point in first 100 years
##
## Plots (user-chosen in INPUTS):
##   (1) Age plot_ageA+ abundance over years (light traces + mean line + 95% interval ribbon)
##   (2) Age plot_ageB+ abundance over years (light traces + mean line + 95% interval ribbon)
##
## Clarifications included:
##   1) Explicit "vulnerable age" definition (vuln_age)
##   2) Optional maturity age override (amat_override)
##   3) Initialization bookkeeping (vulnerable vs adult, per-sim + summaries)
##   4) Scales initial total abundance so mean vulnerable abundance ~= No
##   5) Saves FULL per-simulation initialization table (start_df_sims) to CSV (optional)
############################################################

rm(list = ls())
library(ggplot2)

# ============================================================
# INPUTS  (EDIT ONLY THIS SECTION)
# ============================================================

# Reproducibility
SEED <- 1

# Simulation controls
nT    <- 200    # number of years
A     <- 60      # max age
AR    <- 1       # recruitment age (age class recruits enter simulation)
n_sim <- 2000     # number of simulations (>=500 recommended for probability work)

# --- NEW: Choose TWO ages to plot (Age X+ and Age Y+) ---
plot_ageA <- 1    # e.g., 1 means Age 1+
plot_ageB <- 10   # e.g., 10 means Age 10+

# --- Clarification (1): definition of "vulnerable" fish (for reporting + init scaling)
vuln_age <- 4    # "vulnerable abundance" refers to females age vuln_age+
#vuln_age is a reporting and initialization convention, not related to gear
#we interpret the input No as the target female abundance age vuln_age+ 
#(the “vulnerable” pool). We then scale the starting age structure so the average 
#initial number of females age vuln_age+ across simulations is approximately No, 
#these are the starting conditions.


# Episodic event settings
epfreq <- 10     # years between episodic events (smaller => more frequent things happen)
epAM   <- 0.1   # adult episodic mortality severity (0 to 1) in episodic years, adult survival is multiplied by (1-epAM)
epJM   <- 0.1   # subadult episodic mortality severity (0 to 1) in episodic years, subadult survival is multiplied by (1-epJM)

# Extirpation definition (adult threshold)
extir_threshold <- 50

# Global biological scalars
reck <- 2.8
K    <- 0.13
afec <- 100000
Wmat <- 0.15
sd_S <- 0.6

# --- Clarification (2): optional maturity age override
amat_override <- NA
# e.g., set to 12 for "adults are age 12+"; leave NA to use the Wmat-based formula

# --- Clarification (3): print initialization bookkeeping
report_initial_composition <- TRUE

# --- Clarification (5): save the full per-simulation initialization table
save_start_table <- TRUE
start_table_file <- NA_character_   # will be auto-set after pop_name is chosen

# Plotting controls
plot_trace_sims <- 250        # number of simulations to draw as light traces (set <= n_sim)
save_plots <- FALSE

# Output plot filenames (auto-set after pop_name is chosen unless you override)
plot_A_file <- NA_character_
plot_B_file <- NA_character_

# ------------------------------------------------------------
# Population-specific scalars (UNCOMMENT EXACTLY ONE BLOCK)
# ------------------------------------------------------------

# # Example
# pop_name <- "Example"
# R0       <- 10000
# Madult   <- 0.15
# No       <- 5000

# # Pearl
# pop_name <- "Pearl"
# R0       <- 1832
# Madult   <- 0.274
# No       <- 797

# # Pascagoula
#pop_name <- "Pascagoula"
R0       <- 1869
Madult   <- 0.139
No       <- 161

# # Escambia
# pop_name <- "Escambia"
# R0       <- 498
# Madult   <- 0.139
# No       <- 339

# # Yellow
# pop_name <- "Yellow"
# R0       <- 327
# Madult   <- 0.139
# No       <- 2962

# # Choctawhatchee
# pop_name <- "Choctawhatchee"
# R0       <- 461
# Madult   <- 0.117
# No       <- 2183

# # Apalachicola
# pop_name <- "Apalachicola"
# R0       <- 537
# Madult   <- 0.117
# No       <- 1359

# # Suwannee
 pop_name <- "Suwannee"
 R0       <- 655
 Madult   <- 0.117
 No       <- 6295


# ============================================================
# SAFETY CHECKS
# ============================================================

if (!exists("pop_name") || !exists("R0") || !exists("Madult") || !exists("No")) {
  stop("You must uncomment EXACTLY ONE population block in the INPUTS section (pop_name, R0, Madult, No).")
}
if (vuln_age < 1 || vuln_age > A) stop("vuln_age must be between 1 and A.")
if (AR < 1 || AR > A) stop("AR must be between 1 and A.")
if (n_sim < 1) stop("n_sim must be >= 1.")
if (plot_ageA < 1 || plot_ageA > A) stop("plot_ageA must be between 1 and A.")
if (plot_ageB < 1 || plot_ageB > A) stop("plot_ageB must be between 1 and A.")

plot_trace_sims <- min(plot_trace_sims, n_sim)

if (is.na(start_table_file)) {
  start_table_file <- paste0("start_composition_sims_", pop_name, ".csv")
}
if (is.na(plot_A_file)) {
  plot_A_file <- paste0("timeseries_age", plot_ageA, "plus_", pop_name, ".png")
}
if (is.na(plot_B_file)) {
  plot_B_file <- paste0("timeseries_age", plot_ageB, "plus_", pop_name, ".png")
}

# ============================================================
# MODEL DO NOT MODIFY BELOW
# ============================================================

set.seed(SEED)

# Ages
ages <- 1:A
la   <- (1 - exp(-K * ages))
wa   <- la^3
fec  <- pmax(0, wa - Wmat) * afec

# Adult cutoff age (MetaPVA formula unless overridden)
amat <- as.integer(-log(1 - Wmat^(1/3)) / K)
amat <- max(1, min(A, amat))
if (!is.na(amat_override)) {
  amat <- as.integer(amat_override)
  amat <- max(1, min(A, amat))
}

# Length-based survival-at-age (MetaPVA uses *0.66)
Sa <- exp(-Madult / la * 0.66)   # transition a->a+1 uses Sa[a]

# Survivorship lx
lx <- numeric(A)
lx[1] <- 1
if (A > 1) {
  for (a in 2:A) lx[a] <- lx[a - 1] * Sa[a - 1]
}

# Eggs per recruit + Beverton-Holt params
phie <- sum(lx * fec)
R_A  <- reck / phie
R_B  <- (reck - 1) / (R0 * phie)

# ------------------------------------------------------------
# Clarification (4): scale initial total abundance so mean vulnerable ~= No
# ------------------------------------------------------------
# Expected vulnerable fraction among ages AR..A (without deviates):
vuln_frac <- sum(lx[vuln_age:A]) / sum(lx[AR:A])
vuln_frac <- max(vuln_frac, 1e-12)
N1 <- as.integer(round(No / vuln_frac))

# Episodic years per simulation
ep_flag <- matrix(FALSE, nrow = nT, ncol = n_sim)
if (epfreq <= nT) {
  n_samp <- as.integer(nT / epfreq)
  for (j in 1:n_sim) {
    yrs <- sample.int(nT, n_samp, replace = FALSE)
    ep_flag[yrs, j] <- TRUE
  }
}

# Episodic multipliers for survival transitions (ages AR..A-1)
# Subadults: AR < age < amat
# Adults:    age >= amat
ages_surv <- AR:(A - 1)
mult_episode <- rep(1, length(ages_surv))
if (length(ages_surv) > 0) {
  for (k in seq_along(ages_surv)) {
    a <- ages_surv[k]
    if (a == AR) {
      mult_episode[k] <- 1
    } else if (a < amat) {
      mult_episode[k] <- 1 - epJM
    } else {
      mult_episode[k] <- 1 - epAM
    }
  }
}

# Environmental effects
rec_dev <- matrix(
  exp(rnorm((A - AR + 1) * n_sim, 0, sd_S) + 0.5 * sd_S^2),
  nrow = (A - AR + 1), ncol = n_sim
)
anom <- matrix(
  exp(rnorm((nT - 1) * n_sim, 0, sd_S)),
  nrow = (nT - 1), ncol = n_sim
)

# State: numbers-at-age each year (A x n_sim)
N <- matrix(0L, nrow = A, ncol = n_sim)

# Initialize populations across ages AR..A
base_w <- lx[AR:A]
for (j in 1:n_sim) {
  w <- base_w * rec_dev[, j]
  w <- w / sum(w)
  N[AR:A, j] <- as.vector(rmultinom(1, N1, w))
}

# ============================================================
# Clarification (3 + 5): FULL per-simulation initialization table
# ============================================================

init_total <- colSums(N[1:A, , drop = FALSE])
init_vuln  <- colSums(N[vuln_age:A, , drop = FALSE])
init_adult <- colSums(N[amat:A,    , drop = FALSE])

start_df_sims <- data.frame(
  pop_name = pop_name,
  sim = 1:n_sim,
  init_total_age1plus = init_total,
  init_vulnerable_ageplus = init_vuln,
  init_adults_ageplus = init_adult,
  adult_frac_within_vulnerable = init_adult / pmax(1, init_vuln)
)

if (report_initial_composition) {
  cat("\n--- Initialization bookkeeping (Year 1) ---\n")
  cat("Population:", pop_name, "\n")
  cat("Interpretation: No is target vulnerable abundance (age", vuln_age, "+ females) =", No, "\n")
  cat("Derived total initial abundance (N1) =", N1, "\n")
  cat("Adult cutoff age (amat) =", amat, "\n")
  cat("Observed mean vulnerable (age", vuln_age, "+) =", round(mean(init_vuln), 1), "\n")
  cat("Observed mean adults (age", amat, "+) =", round(mean(init_adult), 1), "\n")
  cat("Adults as % of vulnerable (mean across sims) =",
      round(100 * mean(start_df_sims$adult_frac_within_vulnerable), 1), "%\n")
  cat("-----------------------------------------\n\n")
  print(utils::head(start_df_sims, 10))
}

if (isTRUE(save_start_table)) {
  utils::write.csv(start_df_sims, file = start_table_file, row.names = FALSE)
  cat("Saved per-simulation start table to:", start_table_file, "\n\n")
}

# ============================================================
# TIME SERIES STORAGE
# ============================================================

# two user-chosen age thresholds:
N_ageAplus <- matrix(0, nrow = nT, ncol = n_sim)
N_ageBplus <- matrix(0, nrow = nT, ncol = n_sim)

# kept for extirpation endpoint + bookkeeping
N_adults   <- matrix(0, nrow = nT, ncol = n_sim)  # adults = amat+
N_vuln     <- matrix(0, nrow = nT, ncol = n_sim)  # vulnerable = vuln_age+

record_metrics <- function(t, Nmat) {
  N_ageAplus[t, ] <<- colSums(Nmat[plot_ageA:A, , drop = FALSE])
  N_ageBplus[t, ] <<- colSums(Nmat[plot_ageB:A, , drop = FALSE])
  N_adults[t, ]   <<- colSums(Nmat[amat:A, , drop = FALSE])
  N_vuln[t, ]     <<- colSums(Nmat[vuln_age:A, , drop = FALSE])
}

# Record year 1
record_metrics(1, N)

# ============================================================
# DYNAMICS OVER YEARS
# ============================================================

for (t in 2:nT) {
  
  # Eggs from prior year
  Et <- colSums(N[AR:A, , drop = FALSE] * fec[AR:A])
  
  # Recruitment to age AR
  recruits <- as.integer((R_A * Et * anom[t - 1, ]) / (1 + R_B * Et))
  
  # Next year state
  N_next <- matrix(0L, nrow = A, ncol = n_sim)
  N_next[AR, ] <- recruits
  
  # Survival transitions a -> a+1 for a = AR..A-1
  if (A > AR) {
    surv_base <- Sa[AR:(A - 1)]
    for (k in seq_along(ages_surv)) {
      mult_sim <- ifelse(ep_flag[t - 1, ], mult_episode[k], 1)
      pvec <- pmin(pmax(surv_base[k] * mult_sim, 0), 1)
      
      N_next[ages_surv[k] + 1, ] <- rbinom(
        n    = n_sim,
        size = pmax(0L, N[ages_surv[k], ]),
        prob = pvec
      )
    }
  }
  
  N <- N_next
  record_metrics(t, N)
}

# ============================================================
# EXTIRPATION PROBABILITY
# ============================================================

ext_window <- min(100, nT)
min_adults <- apply(N_adults[1:ext_window, , drop = FALSE], 2, min)
extirp_ind <- (min_adults <= extir_threshold)
p_extirp   <- mean(extirp_ind)

cat("\nPopulation:", pop_name, "\n")
cat("Years:", nT, " | Sims:", n_sim, " | A:", A, " | AR:", AR, "\n")
cat("Vulnerable age (vuln_age) =", vuln_age, " | Adult age cutoff (amat) =", amat, "\n")
cat("Plotted ages: Age", plot_ageA, "+ and Age", plot_ageB, "+\n")
cat("epfreq =", epfreq, " | epAM =", epAM, " | epJM =", epJM, "\n")
cat("Extirpation threshold (adults) =", extir_threshold, "\n")
cat("P(extirpation) = P(min adults <= threshold in years 1..", ext_window, ") = ",
    round(p_extirp, 4), "\n\n", sep = "")

# ============================================================
# PLOTTING (traces + mean + 95% interval ribbon)
# ============================================================

years <- 1:nT
trace_ids <- if (plot_trace_sims < n_sim) sample.int(n_sim, plot_trace_sims) else 1:n_sim

summ_by_year <- function(mat) {
  data.frame(
    year = years,
    mean = rowMeans(mat),
    lwr  = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    upr  = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# Helper to build the plot for an arbitrary "Age X+ matrix"
make_plot <- function(mat, age_label, trace_ids, title_prefix, p_extirp, epfreq, epAM, extir_threshold) {
  
  df_trace <- data.frame(
    year = rep(years, times = length(trace_ids)),
    sim  = rep(trace_ids, each = nT),
    N    = as.vector(mat[, trace_ids, drop = FALSE])
  )
  
  sum_tbl <- summ_by_year(mat)
  
  ggplot() +
    geom_line(data = df_trace, aes(x = year, y = N, group = sim),
              alpha = 0.12, linewidth = 0.3, color = "grey60") +
    geom_ribbon(data = sum_tbl, aes(x = year, ymin = lwr, ymax = upr),
                alpha = 0.20) +
    geom_line(data = sum_tbl, aes(x = year, y = mean),
              linewidth = 1.0, color = "black") +
    labs(
      title = paste0(title_prefix, ": ", age_label, " abundance"),
      subtitle = paste0("Mean (line) + 95% interval (ribbon); traces show ", length(trace_ids),
                        " sims.  P(extirp)=", round(p_extirp, 3),
                        " | epfreq=", epfreq, ", epAM=", epAM, ", threshold=", extir_threshold),
      x = "Year",
      y = paste0("Abundance (", age_label, ")")
    ) +
    theme_minimal()
}

pA <- make_plot(
  mat = N_ageAplus,
  age_label = paste0("Age ", plot_ageA, "+"),
  trace_ids = trace_ids,
  title_prefix = pop_name,
  p_extirp = p_extirp,
  epfreq = epfreq,
  epAM = epAM,
  extir_threshold = extir_threshold
)

pB <- make_plot(
  mat = N_ageBplus,
  age_label = paste0("Age ", plot_ageB, "+"),
  trace_ids = trace_ids,
  title_prefix = pop_name,
  p_extirp = p_extirp,
  epfreq = epfreq,
  epAM = epAM,
  extir_threshold = extir_threshold
)

print(pA)
print(pB)

if (isTRUE(save_plots)) {
  ggsave(plot_A_file, pA, width = 9, height = 4.8, dpi = 300)
  ggsave(plot_B_file, pB, width = 9, height = 4.8, dpi = 300)
  cat("Saved plots:\n  ", plot_A_file, "\n  ", plot_B_file, "\n", sep = "")
}

# NOTE:
# - Keep n_sim large for stable extirpation probabilities.
# - For plotting, use plot_trace_sims << n_sim to keep figures fast and readable.
# - start_df_sims is the full per-simulation initialization table (Year 1).
############################################################
