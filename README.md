# GBM Histomolecular Survival Analysis

This repository contains the R analysis script for a retrospective cohort study examining the integrated histomolecular diagnosis of glioblastoma (GBM). The analysis characterises patient demographics, surgical and treatment patterns, radiological follow-up, and the prognostic impact of molecular biomarkers.

Overall survival (OS) was measured from the date of surgery to death or censorship at the date of last follow-up. Survival time was expressed in months (days ÷ 30.4167). The event indicator was coded 1 (death) or 0 (censored/alive).

Unadjusted OS was estimated using the Kaplan-Meier method (survival::survfit). Between-group differences were tested with the log-rank test (survdiff). KM curves were stratified by:
* Surgical procedure (debulking vs. biopsy)
* First-line treatment (concomitant chemoradiation vs. other)
* Age category (<50, 50–59, 60–69, ≥70 years)

Cox Proportional Hazards Regression
A multivariable Cox proportional hazards model was fitted with age (continuous), sex, and surgical procedure as covariates to produce adjusted hazard ratios (HR) with 95% confidence intervals. Regression tables were formatted using gtsummary::tbl_regression.

Covariate-Adjusted Survival Curves
Covariate-adjusted survival curves were generated using the direct standardisation method (adjustedCurves::adjustedsurv), with 1,000 bootstrap replicates for confidence interval estimation. Adjusted median survival and between-group differences at 12 months were quantified.
