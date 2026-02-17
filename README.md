# ROC-Cutoff-Simulation-Study

This repository contains R code for a Monte Carlo simulation study comparing commonly used ROC-based cut-off selection methods under varying distributional assumptions and realistic diagnostic conditions.

The objective is to identify cut-off criteria that provide **stable, balanced, and reliable diagnostic performance** across effect sizes, variance ratios, sample sizes, and distribution types.

---

## Study Overview

Let X₁ and X₂ denote predictor variables from control and case populations with parameters:

- Control group: (μ₁, σ₁)
- Case group: (μ₂, σ₂)

The goal is to determine an optimal cut-off `c` such that:

- x ≥ c → classified as case  
- x < c → classified as control  

Diagnostic performance is evaluated primarily using sensitivity and specificity.

---

## Simulation Design

A factorial design was used:

- **Effect size (d):** {0.2, 0.5, 0.8}  
- **Variance ratio (k = σ₂ / σ₁):** {1, 1.5, 2}  
- **Sample size per group (n):** {30, 50, 100, ..., 1000}  
- **Monte Carlo replications per scenario:** 1000  

Control group parameters were fixed at:

- μ₁ = 2.5  
- σ₁ = 0.6  

---

## Simulation Experiments

### Experiment I: Normal – Normal (Symmetric Case)

Both groups follow Normal distributions:

X₁ ~ N(μ₁, σ₁²)  
X₂ ~ N(μ₂, σ₂²)

Effect size defined using Cohen’s d.

---

### Experiment II: Gamma – Gamma (Skewed Case)

Both groups follow Gamma distributions:

Xⱼ ~ Gamma(αⱼ, θⱼ)

Effect size defined using a median–IQR standardized difference.  
Case parameters are obtained using numerical optimization (Brent’s method).

---

### Experiment III: Normal – Gamma (Mixed Case)

Control group: Normal distribution  
Case group: Gamma distribution  

Effect size calibrated using the same median–IQR approach.

---

## Cut-off Selection Methods Compared

The following ROC-based criteria were evaluated:

- Youden’s Index  
- Accuracy  
- Sum of Sensitivity and Specificity  
- Product of Sensitivity and Specificity  
- Cohen’s Kappa  
- F1 Score  
- ROC01  
- Minimum p-value  
- Odds Ratio  
- Risk Ratio  

---

## Evaluation Criteria

Performance metrics computed:

- Sensitivity  
- Specificity  
- Accuracy  
- AUC  

Classification of performance:

- **Top-performing:** Sensitivity and Specificity ≥ 70%  
- **Moderate:** Between 40% and 70%  
- **Unstable:** Substantial imbalance between sensitivity and specificity  

Stability was assessed qualitatively across varying sample sizes.

---

## Software and Reproducibility

- R version 4.5.1  
- Package: cutpointr  

A fixed random seed was used:

```r
set.seed(1234)
