###############################################################################
# Compound Poisson-Lognormal P&C loss model with reinsurance
#
# Simulates a property & casualty portfolio's annual aggregate loss, applies
# quota-share and excess-of-loss treaties, prices them by Monte Carlo, and
# compares expected loss, VaR and TVaR gross and net. Also runs a retention
# sensitivity analysis for the excess-of-loss layer.
#
# Base R only -- no external packages -- so it runs anywhere with just R.
# Running it writes two CSV summaries to outputs/ and four PNG charts to
# figures/. Author: Muhammad Muavia.
###############################################################################

set.seed(42)

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------
n_years   <- 100000       # Monte Carlo years
lambda    <- 500          # expected number of claims per year (Poisson)
meanlog   <- 8.0          # lognormal mu   (severity)
sdlog     <- 1.6          # lognormal sigma (heavy-tailed severity)
qs_cede   <- 0.30         # quota share: fraction of every claim ceded (retain 70%)
xol_ret   <- 250000       # excess-of-loss per-claim retention (attachment point)
xol_lim   <- 1000000      # excess-of-loss per-claim limit
alpha     <- 0.995        # VaR / TVaR confidence level (1-in-200)
loading   <- 0.20         # reinsurance premium loading over expected ceded loss
ret_grid  <- seq(100000, 1000000, by = 50000)   # retention sensitivity grid

dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Risk measures
## ---------------------------------------------------------------------------
# Value-at-Risk: the alpha-quantile of the annual loss.
var_fun  <- function(x, a = alpha) as.numeric(quantile(x, a))
# Tail Value-at-Risk (expected shortfall): the mean loss beyond the VaR.
tvar_fun <- function(x, a = alpha) { v <- var_fun(x, a); mean(x[x >= v]) }

## ---------------------------------------------------------------------------
## Monte Carlo simulation of one year of aggregate losses
## ---------------------------------------------------------------------------
# For each simulated year we draw a Poisson number of claims, then that many
# lognormal severities, and sum them. The excess-of-loss treaty is applied per
# claim (reinsurer pays min(max(X - d, 0), limit)); quota share is proportional.
#
# We keep per-year aggregates only, and (for the sensitivity) the per-year ceded
# and net amounts at every retention on the grid. Severities are drawn once and
# reused across retentions so the sensitivity curves are smooth.

simulate <- function() {
  gross   <- numeric(n_years)
  net_qs  <- numeric(n_years)
  ceded_grid <- matrix(0, nrow = n_years, ncol = length(ret_grid))
  net_grid   <- matrix(0, nrow = n_years, ncol = length(ret_grid))

  chunk <- 10000
  start <- 1
  while (start <= n_years) {
    idx <- start:min(start + chunk - 1, n_years)
    c_n <- length(idx)

    n_claims <- rpois(c_n, lambda)                 # claims per year in this chunk
    total    <- sum(n_claims)
    sev      <- rlnorm(total, meanlog, sdlog)      # individual claim severities
    yr       <- rep(seq_len(c_n), n_claims)        # year index for each claim

    g <- tapply(sev, factor(yr, levels = seq_len(c_n)), sum)
    g <- as.numeric(g); g[is.na(g)] <- 0           # years with zero claims -> 0
    gross[idx]  <- g
    net_qs[idx] <- (1 - qs_cede) * g

    for (k in seq_along(ret_grid)) {
      d      <- ret_grid[k]
      ceded  <- pmin(pmax(sev - d, 0), xol_lim)     # reinsurer's per-claim payment
      cy     <- tapply(ceded, factor(yr, levels = seq_len(c_n)), sum)
      cy     <- as.numeric(cy); cy[is.na(cy)] <- 0
      ceded_grid[idx, k] <- cy
      net_grid[idx, k]   <- g - cy
    }
    start <- start + chunk
  }
  list(gross = gross, net_qs = net_qs,
       ceded_grid = ceded_grid, net_grid = net_grid)
}

cat("Simulating", n_years, "years ...\n")
sim <- simulate()
gross <- sim$gross
xol_col <- which(ret_grid == xol_ret)              # headline XoL column
net_xol <- sim$net_grid[, xol_col]

g_var  <- var_fun(gross)
g_tvar <- tvar_fun(gross)

## ---------------------------------------------------------------------------
## Headline comparison: gross vs net of each treaty
## ---------------------------------------------------------------------------
ceded_qs  <- qs_cede * gross
ceded_xol <- sim$ceded_grid[, xol_col]

make_row <- function(name, net, ceded) {
  v <- var_fun(net); t <- tvar_fun(net)
  data.frame(
    Programme                     = name,
    Expected_retained_loss        = mean(net),
    Std_dev                       = sd(net),
    VaR_99_5                       = v,
    TVaR_99_5                      = t,
    VaR_reduction_vs_gross_pct     = 100 * (g_var  - v) / g_var,
    TVaR_reduction_vs_gross_pct    = 100 * (g_tvar - t) / g_tvar,
    Expected_ceded_pure_premium    = mean(ceded),
    Reinsurance_premium_loaded     = mean(ceded) * (1 + loading),
    stringsAsFactors = FALSE)
}

# Human-readable CSV headers, identical to the committed CSV/Excel, so a re-run
# reproduces the same schema. Internal data frames keep simple names for plotting.
risk_names <- c("Programme", "Expected retained loss", "Std dev", "VaR 99.5%",
                "TVaR 99.5%", "VaR reduction vs gross %", "TVaR reduction vs gross %",
                "Expected ceded (pure premium)", "Reinsurance premium (loaded)")
sens_names <- c("Retention (EUR)", "Expected retained loss", "Net VaR 99.5%",
                "Net TVaR 99.5%", "Expected ceded (pure premium)",
                "Reinsurance premium (loaded)", "TVaR reduction vs gross %")
write_csv_named <- function(df, names_vec, path) {
  out <- df; names(out) <- names_vec
  write.csv(out, path, row.names = FALSE)
}

risk_measures <- rbind(
  make_row("Gross (no reinsurance)",           gross,   rep(0, n_years)),
  make_row("Net of quota share (30% ceded)",   sim$net_qs, ceded_qs),
  make_row(sprintf("Net of XoL (%s xs, limit %s)",
                   format(xol_ret, big.mark = ","), format(xol_lim, big.mark = ",")),
           net_xol, ceded_xol))

write_csv_named(risk_measures, risk_names, "outputs/risk_measures.csv")

## ---------------------------------------------------------------------------
## Retention sensitivity (excess-of-loss)
## ---------------------------------------------------------------------------
sens <- data.frame(
  Retention_EUR              = ret_grid,
  Expected_retained_loss     = apply(sim$net_grid,   2, mean),
  Net_VaR_99_5               = apply(sim$net_grid,   2, var_fun),
  Net_TVaR_99_5              = apply(sim$net_grid,   2, tvar_fun),
  Expected_ceded_pure_prem   = apply(sim$ceded_grid, 2, mean))
sens$Reinsurance_premium_loaded <- sens$Expected_ceded_pure_prem * (1 + loading)
sens$TVaR_reduction_vs_gross_pct <- 100 * (g_tvar - sens$Net_TVaR_99_5) / g_tvar
write_csv_named(sens, sens_names, "outputs/retention_sensitivity.csv")

cat("\nExpected annual gross loss: ", format(mean(gross), big.mark = ",",
    nsmall = 0, scientific = FALSE), "\n")
cat("Gross VaR 99.5%:  ", format(round(g_var),  big.mark = ","), "\n")
cat("Gross TVaR 99.5%: ", format(round(g_tvar), big.mark = ","), "\n\n")
print(risk_measures[, c("Programme", "Expected_retained_loss", "VaR_99_5",
                        "TVaR_99_5", "Expected_ceded_pure_premium")], row.names = FALSE)

## ---------------------------------------------------------------------------
## Charts (base graphics)
## ---------------------------------------------------------------------------
col_gross <- "#2b6cb0"; col_qs <- "#dd6b20"; col_xol <- "#2f855a"; col_mut <- "#4a5568"
mm <- function(x) x / 1e6

# Figure 1 -- aggregate loss distribution, gross vs net
png("figures/01_loss_distribution.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 3, 1), family = "sans")
hi   <- var_fun(gross, 0.995)
brk  <- seq(0, hi, length.out = 70)
h_g  <- hist(pmin(gross,        hi), breaks = brk, plot = FALSE)
h_q  <- hist(pmin(sim$net_qs,   hi), breaks = brk, plot = FALSE)
h_x  <- hist(pmin(net_xol,      hi), breaks = brk, plot = FALSE)
ymax <- max(h_g$density, h_q$density, h_x$density)
plot(0, 0, type = "n", xlim = c(0, mm(hi)), ylim = c(0, ymax),
     xlab = "Annual aggregate loss (EUR m)", ylab = "Density", yaxt = "n",
     main = "Reinsurance compresses the annual loss distribution")
lines(mm(h_g$mids), h_g$density, type = "s", col = col_gross, lwd = 2)
lines(mm(h_q$mids), h_q$density, type = "s", col = col_qs,    lwd = 2)
lines(mm(h_x$mids), h_x$density, type = "s", col = col_xol,   lwd = 2)
abline(v = mm(g_var), col = col_gross, lty = 3)
legend("topright", bty = "n",
       legend = c("Gross", "Net of quota share", "Net of excess-of-loss"),
       col = c(col_gross, col_qs, col_xol), lwd = 2)
dev.off()

# Figure 2 -- risk measures gross vs net of each treaty
png("figures/02_risk_measures.png", width = 1600, height = 1000, res = 200)
par(mar = c(4, 4.5, 3, 1), family = "sans")
M <- mm(rbind(risk_measures$Expected_retained_loss,
              risk_measures$VaR_99_5, risk_measures$TVaR_99_5))
colnames(M) <- c("Gross", "Net QS", "Net XoL")
bp <- barplot(t(M), beside = TRUE, col = c(col_gross, col_qs, col_xol),
              names.arg = c("Expected loss", "VaR 99.5%", "TVaR 99.5%"),
              ylab = "EUR m", ylim = c(0, max(M) * 1.15),
              main = "Risk measures: gross vs net of each treaty")
text(bp, t(M), labels = sprintf("%.1f", t(M)), pos = 3, cex = 0.7)
legend("topleft", bty = "n", fill = c(col_gross, col_qs, col_xol),
       legend = c("Gross", "Net quota share", "Net excess-of-loss"))
dev.off()

# Figure 3 -- retention sensitivity
png("figures/03_retention_sensitivity.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 3, 1), family = "sans")
plot(sens$Retention_EUR, mm(sens$Net_TVaR_99_5), type = "b", pch = 19, col = col_xol,
     lwd = 2, ylim = c(min(mm(sens$Expected_retained_loss)) * 0.95, mm(g_tvar) * 1.03),
     xlab = "Excess-of-loss retention / attachment (EUR)", ylab = "EUR m",
     main = "Retention sensitivity: lower retention buys more tail relief")
lines(sens$Retention_EUR, mm(sens$Net_VaR_99_5), type = "b", pch = 15, col = col_gross, lwd = 2)
lines(sens$Retention_EUR, mm(sens$Expected_retained_loss), type = "b", pch = 17, col = col_mut, lwd = 1.6)
abline(h = mm(g_tvar), col = col_xol, lty = 3)
legend("right", bty = "n", col = c(col_xol, col_gross, col_mut), lwd = 2, pch = c(19, 15, 17),
       legend = c("Net TVaR 99.5%", "Net VaR 99.5%", "Expected retained loss"))
dev.off()

# Figure 4 -- cost vs tail relief (efficiency frontier)
png("figures/04_efficiency_frontier.png", width = 1600, height = 1000, res = 200)
par(mar = c(4.5, 4.8, 3, 1), family = "sans")
qs_cost <- mean(ceded_qs)
plot(mm(sens$Expected_ceded_pure_prem), sens$TVaR_reduction_vs_gross_pct,
     type = "b", pch = 19, col = col_xol, lwd = 2,
     xlim = c(0, mm(qs_cost) * 1.10),
     ylim = c(0, 33),
     xlab = "Reinsurance pure premium ceded (EUR m)",
     ylab = "TVaR 99.5% reduction vs gross (%)",
     main = "Cost vs tail relief: XoL is more capital-efficient at low cost")
points(mm(qs_cost), risk_measures$TVaR_reduction_vs_gross_pct[2],
       pch = 19, col = col_qs, cex = 1.6)
text(mm(qs_cost), risk_measures$TVaR_reduction_vs_gross_pct[2],
     "30% QS", pos = 2, col = col_qs, cex = 0.9)
legend("bottomright", bty = "n", col = c(col_xol, col_qs), pch = 19,
       legend = c("Excess-of-loss (varying retention)", "Quota share (30% ceded)"))
dev.off()

cat("\nDone. CSV summaries in outputs/, charts in figures/.\n")
