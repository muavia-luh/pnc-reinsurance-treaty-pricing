# Compound Poisson–Lognormal P&C Loss Model with Reinsurance

An actuarial loss-and-reinsurance model in R. It simulates a property-and-casualty
portfolio's annual claims, applies quota-share and excess-of-loss treaties, prices them
by Monte Carlo, and compares the expected loss, Value-at-Risk and Tail Value-at-Risk the
cedent carries before and after each treaty. A retention sensitivity analysis then shows
how the excess-of-loss layer trades premium for tail protection.

The theme of the project is a practical one: **quota share and excess-of-loss protect
against very different risks, and buy very different amounts of tail relief for the
premium they cost.**

## Model

Annual aggregate loss follows a compound Poisson–lognormal:

- number of claims per year `N ~ Poisson(λ)`,
- individual claim size `X ~ Lognormal(μ, σ)` (heavy-tailed),
- annual loss `S = X₁ + … + X_N`.

Two reinsurance structures are layered on top:

- **Quota share** — the cedent retains a fixed proportion of every claim (70% here) and
  cedes the rest. Being proportional, it scales the entire loss distribution.
- **Excess of loss** — the reinsurer pays the part of each claim above a retention, up to
  a limit (€1,000,000 excess of €250,000 here). It responds only to large claims, so it
  reshapes the tail.

Risk is measured by **VaR and TVaR at the 99.5% level** (a one-in-200-year outcome), and
each treaty is priced as its expected ceded loss plus a 20% loading.

## Results

100,000 simulated years on the default parameters:

| Programme | Expected loss | VaR 99.5% | TVaR 99.5% | Ceded premium |
|---|--:|--:|--:|--:|
| Gross | €5.36m | €8.66m | €9.94m | – |
| Net of quota share (30% ceded) | €3.76m | €6.06m | €6.96m | €1.61m |
| Net of excess-of-loss (250k xs, 1m limit) | €5.09m | €7.21m | €8.46m | €0.27m |

Quota share reduces expected loss and capital by the same 30%, since it scales every
claim, and it cedes €1.61m of premium to achieve that. Excess of loss leaves expected
loss almost unchanged — large claims are rare — yet removes about 15% of the 99.5% tail
for only €0.27m of ceded premium. Measured per euro of premium, the excess-of-loss layer
delivers far more tail protection, while quota share provides broad proportional relief
at a higher cost. The efficiency chart draws that contrast directly, and the retention
analysis shows how far the excess-of-loss trade-off can be pushed.

## Figures

| File | Shows |
|---|---|
| `figures/01_loss_distribution.png` | Annual aggregate loss distribution, gross vs net of each treaty |
| `figures/02_risk_measures.png` | Expected loss, VaR and TVaR side by side |
| `figures/03_retention_sensitivity.png` | Net VaR and TVaR as the excess-of-loss retention moves |
| `figures/04_efficiency_frontier.png` | Ceded premium against tail relief — the cost/benefit frontier |

## Running the model

```r
Rscript reinsurance_model.R
```

Base R only, with no packages to install. The script writes the two CSV summaries to
`outputs/` and the four charts to `figures/`. The random seed is fixed, so results
reproduce exactly.

## Repository

```
reinsurance_model.R                  the model: simulation, treaties, pricing, charts
outputs/risk_measures.csv            gross vs net comparison
outputs/retention_sensitivity.csv    the excess-of-loss retention grid
outputs/reinsurance_summary.xlsx     formatted summary workbook (both tables)
figures/                             the four charts
```

## Scope and assumptions

The frequency and severity parameters are set to be representative of a mid-sized
non-life book rather than fitted to a particular portfolio; adjust `λ`, `μ`, `σ` and the
treaty terms at the top of the script for a specific case. The model covers one line of
business over a one-year horizon and does not include reinstatements, aggregate limits,
reinsurer default, expenses, or dependence between claims. TVaR at 99.5% is estimated
from the 0.5% tail of the simulation and carries some Monte Carlo noise; increasing
`n_years` tightens it.

*Author: Muhammad Muavia · MIT licensed.*
