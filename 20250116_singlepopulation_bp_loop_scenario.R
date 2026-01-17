############################################################
## GULF STURGEON PVA CALCULATOR — 2-scenario comparison
## Scenario 1: baseline (epJM as set in INPUTS)
## Scenario 2: epJM = 0.5
## Output: 2x2 plot:
##   [pA_baseline | pA_epJM0.5]
##   [pB_baseline | pB_epJM0.5]
## With matched y-axes across columns within each row.
############################################################

rm(list = ls())

library(ggplot2)
library(grid)

# ============================================================
# INPUTS  (EDIT ONLY THIS SECTION)
# ============================================================

# Reproducibility
SEED <- 1

# Simulation controls
nT    <- 50
A     <- 60
AR    <- 1
n_sim <- 2000

# Choose TWO ages to plot (Age X+ and Age Y+)
plot_ageA <- 1
plot_ageB <- 10

# "Vulnerable" fish definition (reporting + init scaling only)
vuln_age <- 4

# Episodic event settings
epfreq <- 5
epAM   <- 0.01
epJM   <- 0.01   # <-- baseline value (Scenario 1 uses this)

# Extirpation definition (adult threshold)
extir_threshold <- 50

# Global biological scalars
reck <- 2.8
K    <- 0.13
afec <- 100000
Wmat <- 0.15
sd_S <- 0.6

# Optional maturity age override
amat_override <- NA

# Initialization bookkeeping
report_initial_composition <- TRUE

# Save per-simulation initialization table
save_start_table <- FALSE   # recommend FALSE for quick comparisons

# Plotting controls
plot_trace_sims <- 250
save_plots <- FALSE

# ------------------------------------------------------------
# Population-specific scalars (UNCOMMENT EXACTLY ONE BLOCK)
# ------------------------------------------------------------

pop_name <- "Example"
R0       <- 1000
Madult   <- 0.10
No       <- 5000

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

# ============================================================
# HELPER: run one scenario and return matrices + summaries
# ============================================================

summ_by_year <- function(mat, years) {
  data.frame(
    year = years,
    mean = rowMeans(mat),
    lwr  = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    upr  = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_plot <- function(mat, years, trace_ids, age_label, title_text, subtitle_text, y_max) {
  df_trace <- data.frame(
    year = rep(years, times = length(trace_ids)),
    sim  = rep(trace_ids, each = length(years)),
    N    = as.vector(mat[, trace_ids, drop = FALSE])
  )
  
  sum_tbl <- summ_by_year(mat, years)
  
  ggplot() +
    geom_line(data = df_trace, aes(x = year, y = N, group = sim),
              alpha = 0.12, linewidth = 0.3, color = "grey60") +
    geom_ribbon(data = sum_tbl, aes(x = year, ymin = lwr, ymax = upr),
                alpha = 0.20) +
    geom_line(data = sum_tbl, aes(x = year, y = mean),
              linewidth = 1.0, color = "black") +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      x = "Year",
      y = paste0("Abundance (", age_label, ")")
    ) +
    theme_minimal()
}

run_one <- function(epJM_value, scenario_tag) {
  
  # Reset seed so both scenarios share the same random draws (fair compare)
  set.seed(SEED)
  
  # Local scenario value
  epJM_local <- epJM_value
  
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
  Sa <- exp(-Madult / la * 0.66)
  
  # Survivorship lx
  lx <- numeric(A)
  lx[1] <- 1
  if (A > 1) {
    for (a in 2:A) lx[a] <- lx[a - 1] * Sa[a - 1]
  }
  
  # Beverton-Holt params
  phie <- sum(lx * fec)
  R_A  <- reck / phie
  R_B  <- (reck - 1) / (R0 * phie)
  
  # Scale initial total abundance so mean vulnerable ~= No
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
  ages_surv <- AR:(A - 1)
  mult_episode <- rep(1, length(ages_surv))
  if (length(ages_surv) > 0) {
    for (k in seq_along(ages_surv)) {
      a <- ages_surv[k]
      if (a == AR) {
        mult_episode[k] <- 1
      } else if (a < amat) {
        mult_episode[k] <- 1 - epJM_local
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
  
  # Init composition table (optional)
  init_vuln  <- colSums(N[vuln_age:A, , drop = FALSE])
  init_adult <- colSums(N[amat:A,    , drop = FALSE])
  
  if (report_initial_composition) {
    cat("\n--- Initialization bookkeeping (Year 1) ---\n")
    cat("Scenario:", scenario_tag, "\n")
    cat("Population:", pop_name, "\n")
    cat("No target vulnerable (age", vuln_age, "+) =", No, "\n")
    cat("Derived total initial abundance (N1) =", N1, "\n")
    cat("Adult cutoff age (amat) =", amat, "\n")
    cat("Observed mean vulnerable =", round(mean(init_vuln), 1), "\n")
    cat("Observed mean adults     =", round(mean(init_adult), 1), "\n")
    cat("-----------------------------------------\n\n")
  }
  
  # Time series storage
  N_ageAplus <- matrix(0, nrow = nT, ncol = n_sim)
  N_ageBplus <- matrix(0, nrow = nT, ncol = n_sim)
  N_adults   <- matrix(0, nrow = nT, ncol = n_sim)
  
  record_metrics <- function(t, Nmat) {
    N_ageAplus[t, ] <<- colSums(Nmat[plot_ageA:A, , drop = FALSE])
    N_ageBplus[t, ] <<- colSums(Nmat[plot_ageB:A, , drop = FALSE])
    N_adults[t, ]   <<- colSums(Nmat[amat:A,    , drop = FALSE])
  }
  
  years_vec <- 1:nT
  
  # Record year 1
  record_metrics(1, N)
  
  # Dynamics
  for (t in 2:nT) {
    
    Et <- colSums(N[AR:A, , drop = FALSE] * fec[AR:A])
    
    recruits <- as.integer((R_A * Et * anom[t - 1, ]) / (1 + R_B * Et))
    
    N_next <- matrix(0L, nrow = A, ncol = n_sim)
    N_next[AR, ] <- recruits
    
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
  
  # Extirpation probability
  ext_window <- min(100, nT)
  min_adults <- apply(N_adults[1:ext_window, , drop = FALSE], 2, min)
  p_extirp   <- mean(min_adults <= extir_threshold)
  
  list(
    scenario_tag = scenario_tag,
    epJM = epJM_local,
    amat = amat,
    years = years_vec,
    N_ageAplus = N_ageAplus,
    N_ageBplus = N_ageBplus,
    p_extirp = p_extirp
  )
}

# ============================================================
# RUN BOTH SCENARIOS
# ============================================================

res_base <- run_one(epJM_value = epJM,  scenario_tag = "Baseline")
res_hi   <- run_one(epJM_value = 0.5,   scenario_tag = "epJM = 0.5")

# Use identical trace IDs in both panels for visual comparability
trace_ids <- 1:plot_trace_sims

# Compute shared y-axis maxima (per row) using the 97.5% ribbon bound across BOTH scenarios
sumA_base <- summ_by_year(res_base$N_ageAplus, res_base$years)
sumA_hi   <- summ_by_year(res_hi$N_ageAplus,   res_hi$years)
ymax_A    <- max(sumA_base$upr, sumA_hi$upr)

sumB_base <- summ_by_year(res_base$N_ageBplus, res_base$years)
sumB_hi   <- summ_by_year(res_hi$N_ageBplus,   res_hi$years)
ymax_B    <- max(sumB_base$upr, sumB_hi$upr)

# Build plots (matched y within each row)
pA_base <- make_plot(
  mat = res_base$N_ageAplus, years = res_base$years, trace_ids = trace_ids,
  age_label = paste0("Age ", plot_ageA, "+"),
  title_text = paste0(pop_name, " — Baseline"),
  subtitle_text = paste0("P(extirp)=", round(res_base$p_extirp, 3),
                         " | epJM=", res_base$epJM, ", epAM=", epAM,
                         " | epfreq=", epfreq, ", threshold=", extir_threshold),
  y_max = ymax_A
)

pA_hi <- make_plot(
  mat = res_hi$N_ageAplus, years = res_hi$years, trace_ids = trace_ids,
  age_label = paste0("Age ", plot_ageA, "+"),
  title_text = paste0(pop_name, " — epJM = 0.5"),
  subtitle_text = paste0("P(extirp)=", round(res_hi$p_extirp, 3),
                         " | epJM=", res_hi$epJM, ", epAM=", epAM,
                         " | epfreq=", epfreq, ", threshold=", extir_threshold),
  y_max = ymax_A
)

pB_base <- make_plot(
  mat = res_base$N_ageBplus, years = res_base$years, trace_ids = trace_ids,
  age_label = paste0("Age ", plot_ageB, "+"),
  title_text = paste0(pop_name, " — Baseline"),
  subtitle_text = paste0("P(extirp)=", round(res_base$p_extirp, 3),
                         " | epJM=", res_base$epJM, ", epAM=", epAM,
                         " | epfreq=", epfreq, ", threshold=", extir_threshold),
  y_max = ymax_B
)

pB_hi <- make_plot(
  mat = res_hi$N_ageBplus, years = res_hi$years, trace_ids = trace_ids,
  age_label = paste0("Age ", plot_ageB, "+"),
  title_text = paste0(pop_name, " — epJM = 0.5"),
  subtitle_text = paste0("P(extirp)=", round(res_hi$p_extirp, 3),
                         " | epJM=", res_hi$epJM, ", epAM=", epAM,
                         " | epfreq=", epfreq, ", threshold=", extir_threshold),
  y_max = ymax_B
)

# ============================================================
# 2x2 LAYOUT (baseline left column, epJM=0.5 right column)
# ============================================================

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

print(pA_base, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pA_hi,   vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pB_base, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pB_hi,   vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

############################################################
# ---- SAVE the 2x2 grid plot ----
out_dir  <- "C:/Users/flori/Documents/GS/20251210_GS_Connor_SDM"
out_file <- file.path(out_dir, "pva_compare_baseline_vs_epJM0p5.png")

png(filename = out_file, width = 14, height = 8, units = "in", res = 300)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

print(pA_base, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pA_hi,   vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pB_base, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pB_hi,   vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

dev.off()

cat("Saved:", out_file, "\n")
