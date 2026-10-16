# Mathematical Modelling of Drug Interaction on Evolution of Antibiotic Resistance

MATLAB code for an ODE model of antibiotic-sensitive and antibiotic-resistant bacteria under combination therapy with two drugs. Interaction parameters λ₁ and λ₂ control how strongly the two drugs interact on the sensitive and the resistant population.

## Model

The model is implemented in [`src/Evolution_of_antibiotic_resistance.m`](src/Evolution_of_antibiotic_resistance.m). State vector `y = [S, R, A, B]`:

$$
\begin{aligned}
\frac{dS}{dt} &= \beta_s S\left(1-\frac{S+R}{K}\right) - mS - \left(\lambda_1\alpha_{11}\alpha_{12}AB + \alpha_{11}A + \alpha_{12}B\right)S - \mu_s S \\
\frac{dR}{dt} &= \beta_r R\left(1-\frac{S+R}{K}\right) + mS - \left(\lambda_2\alpha_{21}\alpha_{22}AB + \alpha_{21}A + \alpha_{22}B\right)R - \mu_r R \\
\frac{dA}{dt} &= \theta_1 - \mu_1 A \\
\frac{dB}{dt} &= \theta_2 - \mu_2 B
\end{aligned}
$$

| Symbol | Code name | Meaning |
|--------|-----------|---------|
| S, R | `y(1)`, `y(2)` | Sensitive and resistant bacteria |
| A, B | `y(3)`, `y(4)` | Concentrations of antibiotics A and B |
| β_s, β_r | `beta_s`, `beta_r` | Growth rates |
| μ_s, μ_r | `moo_s`, `moo_r` | Natural death rates |
| K | `K` | Carrying capacity |
| m | `m` (= `q1 + q2`) | Mutation rate from S to R |
| α_ij | `alpha11` … `alpha22` | Effect of antibiotic j on population i (1 = sensitive, 2 = resistant) |
| λ₁, λ₂ | `landa1`, `landa2` | Drug interaction strength on S and on R (0 = no interaction term) |
| θ₁, θ₂ | `teta_1`, `teta_2` | Antibiotic supply rates |
| μ₁, μ₂ | `moo_1`, `moo_2` | Antibiotic elimination rates |

## Repository structure

```
.
├── src/
│   └── Evolution_of_antibiotic_resistance.m    ODE right-hand side
├── scripts/
│   ├── simulate_combination_therapy.m          Solve the ODE model and plot S, R, A, B
│   └── plot_resistant_growth_rate_surface.m    Resistant growth rate over (λ₁, λ₂)
└── results/                                    Saved simulation runs (.mat)
    ├── fig6_combination_therapy/
    ├── interaction_scenarios/
    └── lambda2_MIC.mat
```

## Requirements

MATLAB only, with no toolboxes required. Tested with R2026a.

## How to run

Run either script from any folder, e.g. `run('scripts/simulate_combination_therapy.m')`, or open it and press **Run**.

- **`simulate_combination_therapy.m`**: set the parameters at the top (in particular `landa1` and `landa2`), then run. It solves the model with `ode45` from `y0 = [0.9e9, 0.1e9, 0, 0]` and plots S and R (scaled by 10⁹) with the scaled drug concentrations on a second axis. It also computes `RS` and `RR`, the ratio of growth rate to total loss rate for each population.
- **`plot_resistant_growth_rate_surface.m`**: plots the resistant-strain growth rate over a grid of λ₁, λ₂ ∈ [−1.66, 1.66].

## Saved results

Each `.mat` file in `results/` is a complete MATLAB workspace saved after a simulation run. It contains every parameter plus:

- `t`, `y`: the raw ODE solution, with `y` columns `[S, R, A, B]`
- `ynew1`: the plotted quantities, `[S/1e9, R/1e9, A·μ₁/θ₁, B·μ₂/θ₂]`
- `RS`, `RR`: the ratios printed in the tables below

These runs were produced by earlier versions of the simulation script. Load a file to see its exact settings. In particular, all runs in `fig6_combination_therapy/` used `alpha11 = alpha12 = Emax_s/IC50_s`, while the current script multiplies this by `moo_1/teta_1`. See the note at the top of the script.

### `fig6_combination_therapy/`

| File | λ₁ | λ₂ | IC50_s | IC50_r | RS | RR | t end |
|------|----|----|--------|--------|----|----|-------|
| `FIG6_1_com_therapy.mat` | -1 | 1 | 1.25 | 25 | 0.743 | 1.19 | 100 |
| `FIG6_2_com_therapy.mat` | 1 | -1 | 1.25 | 25 | 0.479 | 1.19 | 100 |
| `FIG6_3_com_therapy.mat` | -1 | -1 | 1.25 | 25 | 0.743 | 1.19 | 100 |
| `FIG6_4_com_therapy.mat` | 1 | 1 | 1.25 | 25 | 0.479 | 1.19 | 100 |

### `interaction_scenarios/`

Three (λ₁, λ₂) pairs, each run in three regimes. The file names describe the regime.

| File | λ₁ | λ₂ | IC50_s | IC50_r | RS | RR | t end |
|------|----|----|--------|--------|----|----|-------|
| **RS < 1, RR < 1** | | | | | | | |
| `lambda1_-1_lambda2_1_Rs_Rr_less_1_increase_IC50_s.mat` | -1 | 1 | 1.25 | 5 | 0.743 | 0.651 | 20 |
| `lambda1_0_lambda2_0_Rs_Rr_less_1_increase_IC50_s.mat` | 0 | 0 | 1.25 | 5 | 0.583 | 0.686 | 20 |
| `lambda1_1_lambda2_-1_Rs_Rr_less_1_increase_IC50_s.mat` | 1 | -1 | 1.25 | 5 | 0.479 | 0.724 | 20 |
| **RS < 1, RR > 1** | | | | | | | |
| `lambda1_-1_lambda2_1_Rs_less_1_Rr_more_1_inc_IC50r.mat` | -1 | 1 | 1.25 | 25 | 0.743 | 1.19 | 20 |
| `lambda1_0_lambda2_0_Rs_less_1_Rr_more_1_inc_IC50r.mat` | 0 | 0 | 1.25 | 25 | 0.583 | 1.19 | 20 |
| `lambda1_1_lambda2_-1_Rs_less_1_Rr_more_1_inc_IC50r.mat` | 1 | -1 | 1.25 | 25 | 0.479 | 1.19 | 20 |
| **RS > 1, RR > 1** | | | | | | | |
| `lambda1_-1_lambda2_1_Rs_bigger_Rr_more_1_inc_IC50r.mat` | -1 | 1 | 2.5 | 25 | 1.62 | 1.28 | 1000 |
| `lambda1_0_lambda2_0_Rs_bigger_Rr_more_1_inc_IC50r.mat` | 0 | 0 | 2.5 | 25 | 1.61 | 1.28 | 1000 |
| `lambda1_1_lambda2_-1_Rs_bigger_Rr_more_1_inc_IC50r.mat` | 1 | -1 | 2.5 | 25 | 1.6 | 1.28 | 1000 |

### `lambda2_MIC.mat`

A single vector `lambda2_MIC` of 15 λ₂ values from −1.667 to 1.6. No script in the repository loads it.

## File name changes

| Previous name | New path |
|---------------|----------|
| `ODE45_solver_EV_AB_R.m` | `src/Evolution_of_antibiotic_resistance.m` (renamed to match the function it defines, so MATLAB can find it) |
| `trial.m` | `scripts/simulate_combination_therapy.m` |
| `Gr_lambda1_lambda2_2.m` | `scripts/plot_resistant_growth_rate_surface.m` |
| `FIG6_*_com_therapy.mat` | `results/fig6_combination_therapy/` |
| `lambda1_*.mat` | `results/interaction_scenarios/` |
| `lambda2_MIC.mat` | `results/lambda2_MIC.mat` |
