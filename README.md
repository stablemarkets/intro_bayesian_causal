<p float="center">
  <img src="/Nonparametrics/dp_oscil.png" width="250" />
  <img src="/Nonparametrics/gp_oscil.png" width="250" />
  <img src="/Nonparametrics/bart_oscil.png" width="250" /> 
</p>

# A Practical Introduction to Bayesian Estimation of Causal Effects: Parametric and Nonparametric Approaches

This is the companion GitHub repository for the paper here: https://onlinelibrary.wiley.com/doi/10.1002/sim.8761.

Please cite the code examples here and discussed in the paper by citing the paper. The BibTex is

```
@article{doi:10.1002/sim.8761,
author = {Oganisian, Arman and Roy, Jason A.},
title = {A practical introduction to Bayesian estimation of causal effects: Parametric and nonparametric approaches},
journal = {Statistics in Medicine},
volume = {n/a},
number = {n/a},
pages = {1-34},
keywords = {BART, Bayesian, Bayesian nonparametric, causal inference, confounding, Dirichlet process, g-computation, Gaussian process},
doi = {10.1002/sim.8761},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.8761},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/sim.8761}}
```

### Software Dependencies

---

All code and analyses generated in R version 3.6.3. We particularly rely on 
- `Stan` to fit GPs and parametric causal models. 
- `ChiRP` (https://stablemarkets.github.io/ChiRPsite/index.html) to fit DP models. 
- `BayesTree` to fit BART models.

### Directory

---

- `dose_response`: contains code implementing model discussed in Section 3.1 of the paper. Code generates Figure 1a.
- `partial_pool`: contains code implementing model discussed in Section 3.2. Code generates Figure 1b.
- `partial_pool`: contains code implementing model discussed in Section 3.2. Code generates Figure 1b.
- `g_comp`: contains code estimating model discussed in Section 4.1 using Ridge prior in Equation (10). Generates Figure 3a.
- `sensitivity`: contains code implementing sensitivity analysis for ignorability violations (Section 5). Generates Figure 3b.
- `Nonparametrics`: contains code implementing DP, GP, and BART models. The file `npbayes.R` generates Figure 4a-c. The filte `npbayes_ATE.R` uses specified models to estimate average treatment effects (ATEs) and generates Figure 4d.


### Causal and Bayesian Topics

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
