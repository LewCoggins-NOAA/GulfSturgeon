############################################################
## Single-population PVA (self-contained)
## Endpoint: min(adults) <= threshold at any point in first 100 years
## Plots:
##   (1) Age 1+ abundance over years (traces + mean)
##   (2) Age 4+ abundance over years (traces + mean)
############################################################

rm(list = ls())

library(ggplot2)

# ============================================================
# INPUTS  (EDIT ONLY THIS SECTION)
# ============================================================

# Reproducibility
SEED <- 1

# Simulation controls
nT    <- 100     # number of years
A     <- 60      # max age
AR    <- 1       # recruitment age (age class recruits enter)
n_sim <- 2000    # number of simulations (>=500 recommended for probability work)

# Episodic event settings
epfreq <- 10     # years between episodic events (smaller => more frequent)
epAM   <- 0.50   # adult episodic mortality severity (0 to 1), adult survival reduced by X% in episodic years
epJM   <- 0.20   # subadult episodic mortality severity (0 to 1)

# Extirpation definition (adult threshold)
extir_threshold <- 50

# Global biological scalars (shared across pops in your meta inputs)
reck <- 2.8
K    <- 0.13
afec <- 100000
Wmat <- 0.15
sd_S <- 0.6

# ------------------------------------------------------------
# Population-specific scalars (UNCOMMENT EXACTLY ONE BLOCK)
# ------------------------------------------------------------

# # Example
 pop_name <- "Example"
 R0       <- 10000
 Madult   <- 0.15
 No       <- 5000

# # Pearl
# pop_name <- "Pearl"
# R0       <- 1832
# Madult   <- 0.274
# No       <- 797

# # Pascagoula
# pop_name <- "Pascagoula"
# R0       <- 1869
# Madult   <- 0.139
# No       <- 161

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
# pop_name <- "Suwannee"
# R0       <- 655
# Madult   <- 0.117
# No       <- 6295

# IMPORTANT: Uncomment one block above. If you forget, you'll get an error below.

# Plot options
save_plots <- FALSE
plot_age1_file <- "timeseries_age1plus.png"
plot_age4_file <- "timeseries_age4plus.png"


# ============================================================
# Make sure to have a population above or errors occur
# ============================================================
if (!exists("pop_name") || !exists("R0") || !exists("Madult") || !exists("No")) {
  stop("You must uncomment EXACTLY ONE population block in the INPUTS section (pop_name, R0, Madult, No).")
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

# Adult cutoff age (matches MetaPVA formula)
amat <- as.integer(-log(1 - Wmat^(1/3)) / K)
amat <- max(1, min(A, amat))

# Length-based survival-at-age (MetaPVA uses *0.66)
Sa <- exp(-Madult / la * 0.66)   # length A; transition a->a+1 uses Sa[a]

# Survivorship lx (length A)
lx <- numeric(A)
lx[1] <- 1
if (A > 1) {
  for (a in 2:A) lx[a] <- lx[a - 1] * Sa[a - 1]
}

# Eggs per recruit + Beverton-Holt params
phie <- sum(lx * fec)
R_A  <- reck / phie
R_B  <- (reck - 1) / (R0 * phie)

# Initial population size
Rinit <- No / sum(lx)
N1    <- as.integer(round(sum(Rinit * lx)))

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

# Initialize populations
base_w <- lx[AR:A]
for (j in 1:n_sim) {
  w <- base_w * rec_dev[, j]
  w <- w / sum(w)
  N[AR:A, j] <- as.vector(rmultinom(1, N1, w))
}

# Time series storage
N_age1plus <- matrix(0, nrow = nT, ncol = n_sim)
N_age4plus <- matrix(0, nrow = nT, ncol = n_sim)
N_adults   <- matrix(0, nrow = nT, ncol = n_sim)  # adults = amat+

record_metrics <- function(t, Nmat) {
  N_age1plus[t, ] <<- colSums(Nmat)
  if (A >= 4) N_age4plus[t, ] <<- colSums(Nmat[4:A, , drop = FALSE]) else N_age4plus[t, ] <<- 0
  N_adults[t, ] <<- colSums(Nmat[amat:A, , drop = FALSE])
}

# Record year 1
record_metrics(1, N)

# Dynamics over years
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
    surv_base <- Sa[AR:(A - 1)]  # baseline per-age transition survival
    
    for (k in seq_along(ages_surv)) {
      a <- ages_surv[k]
      
      mult_sim <- ifelse(ep_flag[t - 1, ], mult_episode[k], 1)
      pvec <- pmin(pmax(surv_base[k] * mult_sim, 0), 1)
      
      N_next[a + 1, ] <- rbinom(
        n    = n_sim,
        size = pmax(0L, N[a, ]),
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
cat("epfreq =", epfreq, " | epAM =", epAM, " | epJM =", epJM, "\n")
cat("Adult cutoff amat =", amat, " | threshold =", extir_threshold, "\n")
cat("P(extirpation) = P(min adults <= threshold in years 1..", ext_window, ") = ",
    round(p_extirp, 4), "\n\n", sep = "")

# ============================================================
# PLOTTING (traces + mean)
# ============================================================

years <- 1:nT

df1 <- data.frame(
  year = rep(years, times = n_sim),
  sim  = rep(1:n_sim, each = nT),
  N    = as.vector(N_age1plus)
)

df4 <- data.frame(
  year = rep(years, times = n_sim),
  sim  = rep(1:n_sim, each = nT),
  N    = as.vector(N_age4plus)
)

mean1 <- data.frame(year = years, N = rowMeans(N_age1plus))
mean4 <- data.frame(year = years, N = rowMeans(N_age4plus))

p1 <- ggplot(df1, aes(x = year, y = N, group = sim)) +
  geom_line(alpha = 0.12, linewidth = 0.3, color = "grey60") +
  geom_line(data = mean1, aes(x = year, y = N), inherit.aes = FALSE,
            linewidth = 1.1, color = "black") +
  labs(
    title = paste0(pop_name, ": Age 1+ abundance (traces + mean)"),
    subtitle = paste0("P(extirp) = ", round(p_extirp, 3),
                      " | epfreq=", epfreq, ", epAM=", epAM, ", threshold=", extir_threshold),
    x = "Year",
    y = "Abundance (Age 1+)"
  ) +
  theme_minimal()

p4 <- ggplot(df4, aes(x = year, y = N, group = sim)) +
  geom_line(alpha = 0.12, linewidth = 0.3, color = "grey60") +
  geom_line(data = mean4, aes(x = year, y = N), inherit.aes = FALSE,
            linewidth = 1.1, color = "black") +
  labs(
    title = paste0(pop_name, ": Age 4+ abundance (traces + mean)"),
    subtitle = paste0("P(extirp) = ", round(p_extirp, 3),
                      " | epfreq=", epfreq, ", epAM=", epAM, ", threshold=", extir_threshold),
    x = "Year",
    y = "Abundance (Age 4+)"
  ) +
  theme_minimal()

print(p1)
print(p4)

if (save_plots) {
  ggsave(plot_age1_file, p1, width = 9, height = 4.8, dpi = 300)
  ggsave(plot_age4_file, p4, width = 9, height = 4.8, dpi = 300)
  cat("Saved plots:\n  ", plot_age1_file, "\n  ", plot_age4_file, "\n", sep = "")
}

# NOTE:
# plotting becomes slow for very large n_sim (e.g., 20000),
# keep n_sim large for probability estimation but plot traces using a smaller
# subset of simulations (e.g., sample 200 sims for plotting).
