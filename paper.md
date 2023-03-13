---
title: 'CRE: A R package for interpretable discovery and inference of heterogeneous treatment effects'
tags:
  - R
  - causal inference
  - heterogeneous effect 
  - interpretability
  - machine learning
authors:
  - name: Riccardo Cadei
    orcid: 0000-0003-2416-8943
    equal-contrib: true
    affiliation: "1, 2"
  - name: Naeem Khoshnevis
    orcid: 0000-0003-4315-1426
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1"
  - name: Falco J. Bargagli Stoffi
    orcid: 0000-0002-6131-8165
    affiliation: "1"
  - name: Kwonsang Lee
    orcid: 0000-0002-5823-4331
    affiliation: "1"
  - name: Daniela Maria Garcia
    orcid: 0000-0003-3226-3561
    affiliation: "1"
affiliations:
 - name: Department of Biostatistics, Harvard School of Public Health
   index: 1
 - name: School of Computer and Communication Sciences, EPFL
   index: 2
date: 2 November 2022
bibliography: paper.bib
header-includes: 
 - \usepackage{algorithm}
 - \usepackage{bm}
 - \usepackage{amsmath}

---

# Summary

In health and social sciences, it is critically important to identify subgroups of the study population where a treatment has notable heterogeneity in the causal effects with respect to the average treatment effect. The bulk of heterogeneous treatment effect (HTE) literature focuses on two major tasks [@dwivedi2020stable]: (i) estimating HTEs by examining the conditional average treatment effect (CATE); (ii) discovering subgroups of a population characterized by HTE. 
Several methodologies have been proposed for both tasks, but providing interpretability in the results is still an open challenge. Interpretability is a non-mathematical concept, yet it is often defined as the degree to which a human can understand the cause of a decision [@kim2016examples], [@miller2018explanation], [@lakkaraju2016interpretable], [@wang2022causal]. Honest Causal Tree [@athey2016] fits this definition perfectly, but despite its high interpretability, it tends to be highly unstable and to find an oversimplified representation of treatment heterogeneity. To accommodate these shortcomings, @bargagli2023causal proposed Causal Rule Ensemble, a new method for HTE characterization in terms of decision rules, via an extensive exploration of heterogeneity patterns by an ensemble-of-trees approach, enforcing high stability in the discovery. `CRE` is an R Package providing a flexible implementation of the Causal Rule Ensemble algorithm.


# Algorithm

Causal Rule Ensemble relies on the Treatment Effect linear decomposition assumption, characterizing the Conditional Average Treatment Effect (CATE) by $M+1$ distinct contributions:
$$\tau(\boldsymbol{x}) = \mathbb{E}[\tau_i | X_i=\boldsymbol{x}] = \bar{\tau} + \sum_{m=1}^M \alpha_m \cdot r_m(\boldsymbol{x})$$
where $\bar{\tau}$ is the Average Treatment Effect (ATE), and for each $m$ in $\{1,..., M\}$, $r_m$ is an interpretable decision rule characterizing a specific subset of the covariate space, and $\alpha_m$ is the corresponding Additive Average Treatment Effect.
`CRE` procedure is divided into two steps, discovery and estimation, and each observation is used for only one of the two steps (honest splitting).
During the discovery step, `CRE` retrieves the $M$ decision rules characterizing the heterogeneity in the treatment effect. A set of candidate decision rules is extracted by an ensemble of trees trained by a fit-the-fit procedure to model some Individual Treatment Effect (ITE) estimates, and among these, only a simple and robust subset of rules is selected for the linear decomposition by the Stability Selection algorithm via LASSO.
During the estimation step, `CRE` estimates the ATE and AATEs, by the normal equations to model some ITE estimates.
In both steps, `CRE` is agnostic concerning the method used for ITE estimation.


# Usage
`CRE` is available both on [CRAN](https://cran.r-project.org/web/packages/CRE/index.html) and [GitHub](https://github.com/NSAPH-Software/CRE) and can be installed and loaded into the R session
using:
```R
install.packages("CRE")
library("CRE")
```

`generate_cre_dataset()` is a flexible synthetic dataset generator, which can be used for simulations before applying CRE to real-world observational data sets. 
```R
dataset <- generate_cre_dataset(n = 2000, 
                                rho = 0, 
                                n_rules = 4, 
                                p = 10,
                                effect_size = 2, 
                                binary_covariates = TRUE,
                                binary_outcome = FALSE,
                                confounding = "no")
y <- dataset$y
z <- dataset$z
X <- dataset$X
```

We propose here three examples of how to run the Causal Rule Esemble algorithm by the `CRE` package.

**Example 1.** Running Causal Rule Ensemble with default parameters described in @bargagli2023causal.
```R
cre_results <- cre(y, z, X)
```

**Example 2.** Running Causal Rule Ensemble with customized ITE estimator.
```R
ite_pred <- ... # personalized ite estimation
cre_results <- cre(y, z, X, ite = ite_pred)
```

**Example 3.** Running Causal Rule Ensemble with customized parameters.
```R
method_params <- list(ratio_dis = 0.25,
                      ite_method_dis="aipw",
                      ps_method_dis = "SL.xgboost",
                      oreg_method_dis = "SL.xgboost",
                      ite_method_inf = "aipw",
                      ps_method_inf = "SL.xgboost",
                      oreg_method_inf = "SL.xgboost")

hyper_params <- list(intervention_vars = c("x1","x2","x3","x4"),
                     offset = NULL,
                     ntrees_rf = 20,
                     ntrees_gbm = 20,
                     node_size = 20,
                     max_nodes = 5,
                     max_depth = 3,
                     t_decay = 0.025,
                     t_ext = 0.025,
                     t_corr = 1,
                     t_pvalue = 0.05,
                     replace = FALSE,
                     stability_selection = TRUE,
                     cutoff = 0.8,
                     pfer = 0.1,
                     penalty_rl = 1)

cre_results <- cre(y, z, X, method_params, hyper_params)
```

The results are collected in an S3 object containing: the number of decision rules extracted at each step (`M`), the `data.frame` of the CATE decomposition estimates with corresponding uncertainty quantification (`CATE`), the list of selected parameters (`method_params` and `hyper_params`), and the predicted ITEs (`ite_pred`). 

`summarize()` and `print()` display a summary of these results, and `plot()` visualizes the CATE decomposition estimates in a range bar plot. Figure 1 reports an example of the proposed results visualization for Example 1. 

![Visualization of Causal Rule Ensemble HTE linear decomposition for Example 1. For each decision rule discovered, the corresponding AATE estimate with 95% confidence interval is reported in a range bar plot. The decision rules are ordered from the most vulnerable (high AATE) to the least, and the ATE is reported on top of the plot.](images/example.pdf)

# Acknowledgements

This work was partially funded by the following grants: NIH: R01ES026217, R01MD012769, R01ES028033, 1R01ES030616, 1R01AG066793, 1R01MD016054-01A1; Sloan Foundation: G-2020-13946.

# References
