# A Practical Introduction to Bayesian Estimation of Causal Effects: Parametric and Nonparametric Approaches
## Companion GitHub Repository

The working paper (currently under review) can be found on arXiv: https://arxiv.org/abs/2004.07375 .

### Topics

---

The paper touches on the following topics
- Standardization (i.e. g-computation in the point-treatment setting).
- G-computation for time-varying treatments.
  - Estimating effects of static regimes.
  - Estimating effects of dynamic regimes.
- Performing sensitivity analyses around causal assumptions via priors.

In terms of Bayesian models we touch upon
- Bayesian bootstrapping. 
- Ridge-like and horseshoe priors for sparsity in high-dimensional regressions. 
- Hierarchical priors that induce partial pooling of conditional causal effects.
- Dirichlet Process (DP) priors.
- Bayesian Additive Regression Trees (BART).
- Gaussian process (GP) priors.

### Software Dependencies.

---

All code and analyses generated in R version 3.6.3. Particularly, we rely on 
- `Stan` to fit GPs and parametric causal models. 
- `ChiRP` (https://stablemarkets.github.io/ChiRPsite/index.html) to fit DP models. 
- `BayesTree` to fit BART models.

### Directory

---

- `dose_response`: contains code implementing model discussed in Section 3.1 of the paper. Code generates Figure 1a.
- `partial_pool`: contains code implementing model discussed in Section 3.2. Code generates Figure 1b.
- `partial_pool`: contains code implementing model discussed in Section 3.2. Code generates Figure 1b.
- `g_comp`: contains code estimating model discussed in Section 4.1 using Ridge prior in Equation (10). Generates Figure 2a.
- `sensitivity`: contains code implementing sensitivity analysis for ignorability violations (Section 5). Generates Figure 2b.
- `Nonparametrics`: contains code implementing DP, GP, and BART models. The file `npbayes.R` generates Figure 3a-c. The filte `npbayes_ATE.R` uses specified models to estimate average treatment effects (ATEs) and generates Figure 3d.


