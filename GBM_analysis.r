# =============================================================================
# Integrated Histomolecular Diagnosis of Glioblastoma Multiforme (GBM)
# Survival and Molecular Biomarker Analysis
# =============================================================================
# Description:
#   This script performs comprehensive survival analysis and molecular
#   biomarker profiling of a retrospective GBM cohort. Analyses include:
#     - Patient demographics and clinical characteristics
#     - Kaplan-Meier overall survival (unadjusted and stratified)
#     - Cox proportional hazards regression (age/gender-adjusted)
#     - Covariate-adjusted survival curves
#     - Radiological MRI follow-up analysis
#     - Molecular panel reporting completeness
#     - MGMT promoter methylation and TERT promoter mutation survival analysis
#
# Input:
#   A CSV file conforming to the required variable structure (see README.md).
#   Update DATA_PATH below to point to your local data file.
#
# Output:
#   Kaplan-Meier plots, density plots, stacked bar charts, forest tables,
#   and adjusted survival curves rendered to the active graphics device.
# 


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# UPDATE THIS PATH to point to your local copy of the anonymised data file.
DATA_PATH <- "path/to/your/GBM_data.csv"


# =============================================================================
# 1. PACKAGE INSTALLATION (run once)
# =============================================================================

required_packages <- c(
  "survival", "ranger", "ggfortify", "ggsurvfit", "adjustedCurves",
  "survminer", "gtsummary", "dplyr", "riskRegression", "pammtools",
  "gtable", "wesanderson", "ggsci", "scales", "forcats", "ggplot2",
  "tidyr", "ggpubr"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# survminer: install from GitHub for latest version
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!"survminer" %in% installed.packages()[, "Package"]) {
  devtools::install_github("kassambara/survminer", build_vignettes = FALSE)
}

# ComplexHeatmap via Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!"ComplexHeatmap" %in% installed.packages()[, "Package"]) {
  BiocManager::install("ComplexHeatmap")
}


# =============================================================================
# 2. LOAD LIBRARIES
# =============================================================================

library(survival)
library(dplyr)
library(ranger)
library(ggplot2)
library(tidyr)
library(ggfortify)
library(ggsurvfit)
library(adjustedCurves)
library(survminer)
library(gtsummary)
library(riskRegression)
library(pammtools)
library(gtable)
library(wesanderson)
library(ggsci)
library(scales)
library(forcats)
library(ggpubr)
library(ComplexHeatmap)


# =============================================================================
# 3. COLOUR PALETTES
# =============================================================================

# NPG (Nature Publishing Group) palette — used throughout
show_col(pal_npg("nrc", alpha = 0.7)(8))

# Named colour constants for consistency
COL_MALE       <- "#4DBBD5B2"
COL_FEMALE     <- "#E64B35B2"
COL_BIOPSY     <- "#3C5488B2"
COL_DEBULKING  <- "#00A087B2"
COL_CRT        <- "#91D1C2B2"
COL_NCRT       <- "#F39B7FB2"


# =============================================================================
# 4. DATA IMPORT
# =============================================================================

GBM <- read.csv(DATA_PATH)


# =============================================================================
# 5. PATIENT DEMOGRAPHICS
# =============================================================================

# --- 5.1 Vital status at censorship ---
cat("Patients alive at censorship (n):", sum(GBM$Alive. == "Y"), "\n")
cat("Patients alive at censorship (%):", 100 * sum(GBM$Alive. == "Y") / nrow(GBM), "\n")

# --- 5.2 Sex distribution ---
cat("Male patients (n):", sum(GBM$Gender == "M"), "\n")
cat("Female patients (n):", sum(GBM$Gender == "F"), "\n")
chisq.test(table(GBM$Gender))

# --- 5.3 Age: summary and sex comparison ---
summary(GBM$Age)

summary(filter(GBM, Gender == "M")$Age)
summary(filter(GBM, Gender == "F")$Age)

shapiro.test(filter(GBM, Gender == "M")$Age)
shapiro.test(filter(GBM, Gender == "F")$Age)

wilcox.test(
  x = filter(GBM, Gender == "M")$Age,
  y = filter(GBM, Gender == "F")$Age,
  alternative = "two.sided"
)

# --- 5.4 Age visualisations ---

# Overall age density
ggDensity.Age.All <- ggplot(data = GBM, aes(x = Age)) +
  geom_histogram(
    aes(y = after_stat(density)), binwidth = 2,
    position = "identity", alpha = 0.7, colour = "black", fill = "grey"
  ) +
  labs(x = "Age (Years)", y = "Density") +
  theme_pubr()

# Age by sex
ggDensity.Age.Gender <- ggplot(data = GBM, aes(x = Age, colour = Gender, fill = Gender)) +
  geom_histogram(
    aes(y = after_stat(density)), binwidth = 2,
    position = "identity", alpha = 0.7
  ) +
  scale_colour_manual(values = c(COL_FEMALE, COL_MALE)) +
  scale_fill_manual(values = c(COL_FEMALE, COL_MALE)) +
  labs(x = "Age (Years)", y = "Density", colour = "Sex", fill = "Sex") +
  annotate("text", label = "p = 0.43", x = 33, y = 0.06) +
  geom_vline(
    xintercept = median(filter(GBM, Gender == "M")$Age),
    linetype = "dashed", colour = COL_MALE, linewidth = 1
  ) +
  geom_vline(
    xintercept = median(filter(GBM, Gender == "F")$Age),
    linetype = "dashed", colour = COL_FEMALE, linewidth = 1
  ) +
  theme_pubr(legend = "bottom")

# --- 5.5 Surgical procedure: summary and age comparison ---

cat("Debulking (n):", sum(GBM$Extent.of.Surgery == "Debulking"), "\n")
cat("Debulking (%):", sum(GBM$Extent.of.Surgery == "Debulking") / nrow(GBM), "\n")
cat("Biopsy (n):",    sum(GBM$Extent.of.Surgery == "Biopsy"), "\n")
cat("Biopsy (%):",    sum(GBM$Extent.of.Surgery == "Biopsy") / nrow(GBM), "\n")

summary(filter(GBM, Extent.of.Surgery == "Debulking")$Age)
summary(filter(GBM, Extent.of.Surgery == "Biopsy")$Age)

shapiro.test(filter(GBM, Extent.of.Surgery == "Debulking")$Age)
shapiro.test(filter(GBM, Extent.of.Surgery == "Biopsy")$Age)

wilcox.test(
  x = filter(GBM, Extent.of.Surgery == "Debulking")$Age,
  y = filter(GBM, Extent.of.Surgery == "Biopsy")$Age,
  alternative = "two.sided"
)

# Age by procedure
ggDensity.Age.Surgery <- ggplot(
  data = GBM, aes(x = Age, colour = Extent.of.Surgery, fill = Extent.of.Surgery)
) +
  geom_histogram(
    aes(y = after_stat(density)), binwidth = 2,
    position = "identity", alpha = 0.5
  ) +
  scale_colour_manual(values = c(COL_BIOPSY, COL_DEBULKING)) +
  scale_fill_manual(values = c(COL_BIOPSY, COL_DEBULKING)) +
  labs(
    x = "Age (Years)", y = "Density",
    colour = "Extent of Surgery", fill = "Extent of Surgery"
  ) +
  annotate("text", label = "p < 0.001", x = 33, y = 0.075) +
  geom_vline(
    xintercept = median(filter(GBM, Extent.of.Surgery == "Debulking")$Age),
    linetype = "dashed", colour = COL_DEBULKING, linewidth = 1
  ) +
  geom_vline(
    xintercept = median(filter(GBM, Extent.of.Surgery == "Biopsy")$Age),
    linetype = "dashed", colour = COL_BIOPSY, linewidth = 1
  ) +
  theme_pubr(legend = "bottom")

# --- 5.6 Intraoperative 5-ALA use (debulking cases only) ---
cat(
  "5-ALA use in debulking cases (%):",
  100 * sum(GBM$Intraoperative.5.ALA.Status == "Yes") /
    sum(GBM$Extent.of.Surgery == "Debulking"),
  "\n"
)

# --- 5.7 First-line treatment: chemoradiation vs other ---

crt_age_df <- data.frame(Age = GBM$Age, Tx = GBM$First.line.treatment)
crt_age_df$Tx[crt_age_df$Tx != "Radiotherapy with concomitant TMZ"] <- "Non Concomitant Chemoradiotherapy"

cat(
  "Concomitant chemoradiation (%):",
  100 * sum(GBM$First.line.treatment == "Radiotherapy with concomitant TMZ") / nrow(GBM),
  "\n"
)

summary(as.numeric(
  filter(GBM, Timing.to.chemoradiotherpay.since.surgery..DAYS. != "#NUM!")$Timing.to.chemoradiotherpay.since.surgery..DAYS.
))

summary(filter(crt_age_df, Tx == "Radiotherapy with concomitant TMZ")$Age)
summary(filter(crt_age_df, Tx == "Non Concomitant Chemoradiotherapy")$Age)

wilcox.test(
  x = filter(crt_age_df, Tx == "Radiotherapy with concomitant TMZ")$Age,
  y = filter(crt_age_df, Tx == "Non Concomitant Chemoradiotherapy")$Age,
  alternative = "two.sided"
)

# Age by treatment group
ggDensity.Age.CRT <- ggplot(
  data = crt_age_df, aes(x = Age, colour = Tx, fill = Tx)
) +
  geom_histogram(
    aes(y = after_stat(density)), binwidth = 2,
    position = "identity", alpha = 0.7
  ) +
  scale_colour_manual(
    values = c(COL_NCRT, COL_CRT),
    labels = c("Non Concomitant Chemoradiation", "Concomitant Chemoradiation")
  ) +
  scale_fill_manual(
    values = c(COL_NCRT, COL_CRT),
    labels = c("Non Concomitant Chemoradiation", "Concomitant Chemoradiation")
  ) +
  labs(x = "Age (Years)", y = "Density", colour = "First-line Treatment", fill = "First-line Treatment") +
  annotate("text", label = "p < 0.001", x = 33, y = 0.075) +
  geom_vline(
    xintercept = median(filter(crt_age_df, Tx == "Non Concomitant Chemoradiotherapy")$Age),
    linetype = "dashed", colour = COL_NCRT, linewidth = 1
  ) +
  geom_vline(
    xintercept = median(filter(crt_age_df, Tx == "Radiotherapy with concomitant TMZ")$Age),
    linetype = "dashed", colour = COL_CRT, linewidth = 1
  ) +
  theme_pubr(legend = "bottom")

# 2x2 demographics grid
ggarrange(
  ggDensity.Age.All, ggDensity.Age.Gender,
  ggDensity.Age.Surgery, ggDensity.Age.CRT,
  ncol = 2, nrow = 2
)


# =============================================================================
# 6. RADIOLOGICAL MRI FOLLOW-UP
# =============================================================================

follow_up_cols <- paste0("Duration.Between.Follow.Up.", c(1:16, "16.1"))

GBM_MRI_follow_up <- GBM[, c("Unit.Hospital.Number", "Age", "Extent.of.Surgery",
                               "Survival.Time..DAYS.", follow_up_cols)]
colnames(GBM_MRI_follow_up)[1:4] <- c("ID", "Age", "Surgery", "Survival")
GBM_MRI_follow_up$Tx <- crt_age_df$Tx

GBM_MRI_follow_up[GBM_MRI_follow_up == "0" | GBM_MRI_follow_up == "#NUM!"] <- "Not Applicable"

fu_col_idx <- which(colnames(GBM_MRI_follow_up) %in% follow_up_cols)
GBM_MRI_follow_up$No.Follow.Up <- sapply(
  seq_len(nrow(GBM_MRI_follow_up)),
  function(n) sum(GBM_MRI_follow_up[n, fu_col_idx] != "Not Applicable")
)

cat("Patients with ≥1 MRI follow-up (n):", sum(GBM_MRI_follow_up$No.Follow.Up != 0), "\n")
cat("Patients with ≥1 MRI follow-up (%):",
    100 * sum(GBM_MRI_follow_up$No.Follow.Up != 0) / nrow(GBM_MRI_follow_up), "\n")

# Mean inter-follow-up interval
GBM_MRI_follow_up$Mean.Int.Follow.Up <- sapply(
  seq_len(nrow(GBM_MRI_follow_up)),
  function(n) {
    vals <- as.numeric(GBM_MRI_follow_up[n, fu_col_idx][GBM_MRI_follow_up[n, fu_col_idx] != "Not Applicable"])
    if (length(vals) == 0) return(NA)
    mean(vals)
  }
)

# Summary by surgical procedure
for (proc in c("Debulking", "Biopsy")) {
  sub <- filter(GBM_MRI_follow_up, Surgery == proc)
  cat("\n--- MRI Follow-up:", proc, "---\n")
  cat("With follow-up (n):", sum(sub$No.Follow.Up != 0), "\n")
  cat("With follow-up (%):", 100 * sum(sub$No.Follow.Up != 0) / nrow(sub), "\n")
  print(summary(sub$Mean.Int.Follow.Up))
  print(confint(lm(sub$Mean.Int.Follow.Up ~ 1), level = 0.95))
}

# Survival-adjusted follow-up counts by procedure
for (proc in c("Debulking", "Biopsy")) {
  sub <- filter(GBM_MRI_follow_up, Surgery == proc)
  fit <- lm(sub$No.Follow.Up ~ sub$Survival)
  cat("\nSurvival-adjusted follow-up count [", proc, "] — fitted values:\n")
  print(summary(fit$fitted.values))
  print(confint(lm(fit$fitted.values ~ 1), level = 0.95))
}

# Summary by treatment group
for (tx in c("Radiotherapy with concomitant TMZ", "Non Concomitant Chemoradiotherapy")) {
  sub <- filter(GBM_MRI_follow_up, Tx == tx)
  cat("\n--- MRI Follow-up:", tx, "---\n")
  cat("With follow-up (n):", sum(sub$No.Follow.Up != 0), "\n")
  cat("With follow-up (%):", 100 * sum(sub$No.Follow.Up != 0) / nrow(sub), "\n")
  print(summary(sub$Mean.Int.Follow.Up))
  print(confint(lm(sub$Mean.Int.Follow.Up ~ 1), level = 0.95))
}


# =============================================================================
# 7. SURVIVAL ANALYSIS DATASET
# =============================================================================

GBM_survival <- data.frame(
  Age      = GBM$Age,
  Gender   = GBM$Gender,
  Surgery  = GBM$Extent.of.Surgery,
  Treatment = crt_age_df$Tx,
  OSDays   = GBM$Survival.Time..DAYS.,
  Status   = GBM$Alive.
)

GBM_survival$Status[GBM_survival$Status == "Y"] <- 0
GBM_survival$Status[GBM_survival$Status == "N"] <- 1
GBM_survival$Status <- as.numeric(GBM_survival$Status)

GBM_survival$OSMonths <- GBM_survival$OSDays / 30.4167

GBM_survival$Treatment[GBM_survival$Treatment == "Non Concomitant Chemoradiotherapy"] <- "No Chemoradiation"
GBM_survival$Treatment[GBM_survival$Treatment == "Radiotherapy with concomitant TMZ"] <- "Chemoradiation"

cat("Total deaths (events):", sum(GBM_survival$Status == 1), "\n")


# =============================================================================
# 8. KAPLAN-MEIER SURVIVAL ANALYSIS (UNADJUSTED)
# =============================================================================

# --- 8.1 Whole cohort ---
KM.OS.Whole <- survfit(Surv(OSMonths, Status) ~ 1, data = GBM_survival)
for (t in c(12, 24, 36, 48, 60)) summary(KM.OS.Whole, times = t)

ggKM.Whole <- KM.OS.Whole %>%
  ggsurvfit() +
  labs(x = "Months", y = "Probability of Overall Survival") +
  scale_x_continuous(breaks = seq(0, 60, 4)) +
  add_censor_mark() + add_confidence_interval() + add_risktable() +
  theme_classic2()

# --- 8.2 By surgical procedure ---
KM.OS.Surgery <- survfit(Surv(OSMonths, Status) ~ Surgery, data = GBM_survival)
survdiff(Surv(OSMonths, Status) ~ Surgery, data = GBM_survival)
summary(KM.OS.Surgery, times = 12)

ggKM.Surgery <- KM.OS.Surgery %>%
  ggsurvfit() +
  labs(x = "Months", y = "Probability of Overall Survival") +
  scale_x_continuous(breaks = seq(0, 60, 4)) +
  add_censor_mark() + add_confidence_interval() + add_risktable() +
  scale_colour_manual(labels = c("Biopsy", "Debulking"), values = c(COL_BIOPSY, COL_DEBULKING)) +
  scale_fill_manual(labels = c("Biopsy", "Debulking"),   values = c(COL_BIOPSY, COL_DEBULKING)) +
  theme_classic2() + theme_pubr(legend = "bottom")

# --- 8.3 By first-line treatment ---
KM.OS.Tx <- survfit(Surv(OSMonths, Status) ~ Treatment, data = GBM_survival)
survdiff(Surv(OSMonths, Status) ~ Treatment, data = GBM_survival)
for (t in c(12, 24, 36, 48, 60)) summary(KM.OS.Tx, times = t)

ggKM.Tx <- KM.OS.Tx %>%
  ggsurvfit() +
  labs(x = "Months", y = "Probability of Overall Survival") +
  scale_x_continuous(breaks = seq(0, 60, 4)) +
  add_censor_mark() + add_confidence_interval() + add_risktable() +
  scale_colour_manual(labels = c("Chemoradiation", "No Chemoradiation"), values = c(COL_CRT, COL_NCRT)) +
  scale_fill_manual(labels = c("Chemoradiation", "No Chemoradiation"),   values = c(COL_CRT, COL_NCRT)) +
  theme_classic2() + theme_pubr(legend = "bottom")

# --- 8.4 By age category ---
GBM_survival$Age.Category <- cut(
  GBM_survival$Age,
  breaks = c(-Inf, 50, 60, 70, Inf),
  labels = c("<50", "50~59", "60~69", "\u226570"),
  right  = FALSE
)

KM.OS.Age <- survfit(Surv(OSMonths, Status) ~ Age.Category, data = GBM_survival)
survdiff(Surv(OSMonths, Status) ~ Age.Category, data = GBM_survival)
summary(KM.OS.Age, times = 12)

ggKM.Age <- KM.OS.Age %>%
  ggsurvfit() +
  labs(x = "Months", y = "Probability of Overall Survival") +
  scale_x_continuous(breaks = seq(0, 60, 4)) +
  add_risktable() + add_confidence_interval() + add_censor_mark() +
  scale_color_npg(labels = c("Age <50", "Age \u226570", "50\u2264 Age \u226459", "60\u2264 Age \u226469")) +
  scale_fill_npg(labels  = c("Age <50", "Age \u226570", "50\u2264 Age \u226459", "60\u2264 Age \u226469")) +
  theme_classic2() + theme_pubr(legend = "bottom")

# 2x2 KM grid
ggarrange(ggKM.Whole, ggKM.Surgery, ggKM.Tx, ggKM.Age, ncol = 2, nrow = 2)


# =============================================================================
# 9. COX PROPORTIONAL HAZARDS MODEL (AGE AND SEX-ADJUSTED)
# =============================================================================

GBM_survival$Surgery <- as.factor(GBM_survival$Surgery)

cox_OS <- coxph(
  Surv(OSMonths, Status) ~ Surgery + Age + Gender,
  data = GBM_survival,
  x    = TRUE
)

summary(cox_OS)
cox_OS %>% tbl_regression(exponentiate = TRUE)


# =============================================================================
# 10. COVARIATE-ADJUSTED SURVIVAL CURVES
# =============================================================================

OS_adjSurgery <- adjustedsurv(
  data          = GBM_survival,
  variable      = "Surgery",
  ev_time       = "OSMonths",
  event         = "Status",
  method        = "direct",
  outcome_model = cox_OS,
  conf_int      = TRUE,
  bootstrap     = TRUE,
  n_boot        = 1000
)

plot(
  OS_adjSurgery,
  conf_int      = TRUE,
  censoring_ind = "lines",
  xlab          = "Months",
  ylab          = "Adjusted Probability of Overall Survival",
  legend.title  = "Surgery",
  conf_int_alpha = 0.3,
  line_size     = 0.5
)

adjusted_surv_quantile(OS_adjSurgery, conf_int = TRUE, use_boot = TRUE)

adjusted_curve_diff(
  OS_adjSurgery,
  times   = 12,
  group_1 = "Debulking",
  group_2 = "Biopsy",
  conf_int = TRUE
)

adjusted_curve_test(
  OS_adjSurgery,
  conf_level = 0.95,
  to         = max(GBM_survival$OSMonths, na.rm = TRUE)
)


# =============================================================================
# 11. MOLECULAR BIOMARKER PANEL
# =============================================================================

# Extract molecular variables (columns as per dataset structure)
GBM_molecular <- GBM[, c(1, 13, 15:32)]
GBM_molecular$MGMT.Promoter.Methylation[
  GBM_molecular$MGMT.Promoter.Methylation == "Not reported"
] <- "Not Reported"

# Reporting completeness by parameter
mol_cols <- c(4:17, 19)
Percentage.R <- data.frame(
  Molecular.Parameter = colnames(GBM_molecular)[mol_cols],
  Percentage.Reported = 100 * sapply(mol_cols, function(i) {
    sum(GBM_molecular[, i] != "Not Reported") / nrow(GBM_molecular)
  })
) %>% arrange(Percentage.Reported)

Order <- Percentage.R$Molecular.Parameter

# Stacked bar chart: reporting completeness
GBM_molecular_long <- pivot_longer(
  GBM_molecular,
  cols      = all_of(mol_cols),
  names_to  = "Molecular.Parameter",
  values_to = "Report"
)

GBM_molecular_long %>%
  count(Molecular.Parameter, Report) %>%
  ggplot() +
  geom_col(aes(x = Molecular.Parameter, y = n, fill = as.character(Report)),
           position = "fill") +
  coord_flip() +
  scale_x_discrete(
    limits = Order,
    name   = "Molecular Parameter",
    labels = c(
      "IDH.1.132H.Mutation"             = "IDH1 Mutation",
      "IDH.2.172H.Mutation"             = "IDH2 Mutation",
      "ATRX.Mutation"                   = "ATRX Mutation",
      "X1p.19q.Codeletion"              = "1p19q Co-deletion",
      "BRAF.Fusion.and.Gene.Mutation."  = "BRAF Mutation",
      "Histone.H3.3.K27M.Mutation"      = "Histone H3.3 Mutation",
      "H3C2.Mutation"                   = "H3C2 Mutation",
      "TP53.Mutation"                   = "TP53 Mutation",
      "EGFR.Amplification"              = "EGFR Amplification",
      "PDGFRA.Amplification"            = "PDGFRA Amplification",
      "MYC.Amplification"               = "MYC Amplification",
      "PTEN.Deletion"                   = "PTEN Deletion",
      "CDKN2A.Deletion"                 = "CDKN2A Deletion",
      "MGMT.Promoter.Methylation"       = "MGMT Promoter Methylation",
      "TERT.Promoter.Mutation"          = "TERT Promoter Mutation"
    )
  ) +
  ylab("Proportion of Cases Reported") +
  scale_fill_npg(name = "Status") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Full glioma panel completion (8 core markers)
core_markers <- c(4:9, 17, 19)
cat(
  "Cases with complete 8-marker glioma panel (n):",
  sum(sapply(seq_len(nrow(GBM_molecular)),
             function(n) sum(GBM_molecular[n, core_markers] != "Not Reported")) == 8),
  "\n"
)


# =============================================================================
# 12. MGMT PROMOTER METHYLATION ANALYSIS
# =============================================================================

# Build MGMT-specific survival dataset (exclude cases not reported)
MGMT_reported <- GBM_molecular$MGMT.Promoter.Methylation != "Not Reported"

GBM_MGMT_Survival <- data.frame(
  Age       = GBM$Age[MGMT_reported],
  Gender    = GBM$Gender[MGMT_reported],
  Surgery   = GBM$Extent.of.Surgery[MGMT_reported],
  Treatment = crt_age_df$Tx[MGMT_reported],
  MGMT      = GBM_molecular$MGMT.Promoter.Methylation[MGMT_reported],
  OSMonths  = GBM_survival$OSMonths[MGMT_reported],
  Status    = GBM$Alive.[MGMT_reported]
)

GBM_MGMT_Survival$Status[GBM_MGMT_Survival$Status == "Y"] <- 0
GBM_MGMT_Survival$Status[GBM_MGMT_Survival$Status == "N"] <- 1
GBM_MGMT_Survival$Status <- as.numeric(GBM_MGMT_Survival$Status)

GBM_MGMT_Survival$Age.Category <- cut(
  GBM_MGMT_Survival$Age,
  breaks = c(-Inf, 50, 60, 70, Inf),
  labels = c("<50", "50~59", "60~69", "\u226570"),
  right  = FALSE
)

# MGMT % methylation distribution
summary(as.numeric(
  GBM_molecular$Average...Methylation.Across.Tested.CpG.Sites[
    GBM_molecular$Average...Methylation.Across.Tested.CpG.Sites != "Not reported"
  ]
))

ggDensity.MGMT.Percentage <- ggplot(
  GBM_molecular,
  aes(x = as.numeric(Average...Methylation.Across.Tested.CpG.Sites))
) +
  geom_histogram(
    aes(y = after_stat(density)), binwidth = 2,
    position = "identity", alpha = 0.7, colour = "black", fill = "grey"
  ) +
  labs(x = "Percentage Methylation (%)", y = "Density") +
  theme_pubr()

# Helper function: KM curve for MGMT stratification
km_mgmt <- function(data, title_label) {
  fit <- survfit(Surv(OSMonths, Status) ~ MGMT, data = data)
  survdiff(Surv(OSMonths, Status) ~ MGMT, data = data)
  fit %>%
    ggsurvfit() +
    labs(x = "Months", y = "Probability of Overall Survival") +
    scale_x_continuous(breaks = seq(0, 60, 4)) +
    add_censor_mark() + add_risktable() +
    scale_colour_npg(labels = c("Equivocal", "Methylated", "Unmethylated")) +
    scale_fill_npg(labels   = c("Equivocal", "Methylated", "Unmethylated")) +
    theme_classic2() + theme_pubr(legend = "bottom") +
    annotate("text", label = title_label, x = 40, y = 0.85)
}

ggKM.MGMT.All      <- km_mgmt(GBM_MGMT_Survival, "Whole cohort; p = 0.7")
ggKM.MGMT.CRT      <- km_mgmt(filter(GBM_MGMT_Survival, Treatment == "Radiotherapy with concomitant TMZ"),
                               "With Chemoradiation; p = 0.4")
ggKM.MGMT.nCRT     <- km_mgmt(filter(GBM_MGMT_Survival, Treatment != "Radiotherapy with concomitant TMZ"),
                               "Without Chemoradiation; p = 0.6")
ggKM.MGMT.Over65   <- km_mgmt(filter(GBM_MGMT_Survival, Age >= 65), "Age \u226565; p = 0.6")
ggKM.MGMT.Under65  <- km_mgmt(filter(GBM_MGMT_Survival, Age <  65), "Age <65; p = 0.4")
ggKM.MGMT.Debulking <- km_mgmt(filter(GBM_MGMT_Survival, Surgery == "Debulking"), "Debulking; p = 0.3")
ggKM.MGMT.Biopsy    <- km_mgmt(filter(GBM_MGMT_Survival, Surgery == "Biopsy"),    "Biopsy; p = 0.3")

# MGMT panel (4x2)
ggarrange(
  ggKM.MGMT.All,      ggDensity.MGMT.Percentage,
  ggKM.MGMT.Debulking, ggKM.MGMT.Biopsy,
  ggKM.MGMT.CRT,       ggKM.MGMT.nCRT,
  ggKM.MGMT.Over65,    ggKM.MGMT.Under65,
  ncol = 2, nrow = 4
)

# Multivariable Cox model incorporating MGMT
cox.OS.MGMT <- coxph(
  Surv(OSMonths, Status) ~ Age.Category + Gender + Surgery + Treatment + MGMT,
  data = GBM_MGMT_Survival,
  x    = TRUE
)

summary(cox.OS.MGMT)
tbl_regression(cox.OS.MGMT, exponentiate = TRUE)


# =============================================================================
# 13. TERT PROMOTER MUTATION ANALYSIS
# =============================================================================

TERT_reported <- GBM_molecular$TERT.Promoter.Mutation.Loci..if.applicable. != "Not Reported"

GBM_TERT_Survival <- data.frame(
  Age      = GBM$Age[TERT_reported],
  Gender   = GBM$Gender[TERT_reported],
  Surgery  = GBM$Extent.of.Surgery[TERT_reported],
  Treatment = crt_age_df$Tx[TERT_reported],
  TERT     = GBM_molecular$TERT.Promoter.Mutation.Loci..if.applicable.[TERT_reported],
  OSMonths = GBM_survival$OSMonths[TERT_reported],
  Status   = GBM$Alive.[TERT_reported]
)

GBM_TERT_Survival$Status[GBM_TERT_Survival$Status == "Y"] <- 0
GBM_TERT_Survival$Status[GBM_TERT_Survival$Status == "N"] <- 1
GBM_TERT_Survival$Status <- as.numeric(GBM_TERT_Survival$Status)

# Exclude rare compound mutations (c.-124C>T and c.-146C>T)
GBM_TERT_Survival <- filter(GBM_TERT_Survival, TERT != "c.-124C>T and c.-146C>T")

# TERT locus distribution
cat("TERT Not Detected (n):", sum(GBM_molecular$TERT.Promoter.Mutation.Loci..if.applicable. == "Not Detected"), "\n")
cat("TERT c.-124C>T (n):",   sum(GBM_molecular$TERT.Promoter.Mutation.Loci..if.applicable. == "c.-124C>T"), "\n")
cat("TERT c.-146C>T (n):",   sum(GBM_molecular$TERT.Promoter.Mutation.Loci..if.applicable. == "c.-146C>T"), "\n")

KM.OS.TERT <- survfit(Surv(OSMonths, Status) ~ TERT, data = GBM_TERT_Survival)
survdiff(Surv(OSMonths, Status) ~ TERT, data = GBM_TERT_Survival)
summary(KM.OS.TERT, times = 12)

KM.OS.TERT %>%
  ggsurvfit() +
  labs(x = "Months", y = "Probability of Overall Survival") +
  scale_x_continuous(breaks = seq(0, 60, 4)) +
  add_censor_mark() +
  scale_colour_npg(labels = c("c.-124C>T", "c.-146C>T", "No Mutation")) +
  scale_fill_npg() +
  theme_classic2() + theme_pubr(legend = "bottom")

# =============================================================================
# END OF SCRIPT
# =============================================================================
