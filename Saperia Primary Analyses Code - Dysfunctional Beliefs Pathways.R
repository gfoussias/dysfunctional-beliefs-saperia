library(lavaan)
library(dplyr)
library(psych)
library(GPArotation)

dataset = read.csv("Baseline Variables 1935.csv")
subset = subset(subset,subset$Invalid=='0')


#DIMENSION REDUCTION#
# Subset the belief variables
belief_vars <- dataset[, c("zdpb", "zabs", "zrses", "zplex", "zismi", "zsaraq")]

# Replace missing data with variable means
belief_vars_imputed <- as.data.frame(lapply(belief_vars, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))

# 1. KMO and Bartlett's test
kmo_result <- KMO(cor(belief_vars_imputed))
print(kmo_result)

bartlett_result <- cortest.bartlett(cor(belief_vars_imputed), n = nrow(belief_vars_imputed))
print(bartlett_result)

# 2. PCA - first extract all components to check eigenvalues
pca_initial <- principal(belief_vars_imputed, nfactors = ncol(belief_vars_imputed), rotate = "none")
print(pca_initial$values)  # eigenvalues

# Determine number of factors with eigenvalue >1
n_factors <- sum(pca_initial$values > 1)
cat("Number of factors retained (eigenvalue > 1):", n_factors, "\n")

# Run final PCA with Varimax rotation
pca_final <- principal(belief_vars_imputed, nfactors = n_factors, rotate = "varimax", scores = TRUE)
print(pca_final$loadings)  # rotated factor loadings

# Extract regression-based factor scores
factor_scores <- pca_final$scores

factor_scores <- pca_final$scores

# Standardize the factor scores (z-scores)
zFactor <- as.numeric(scale(factor_scores))

# Add zFactor to the original dataset
dataset$zFactor <- zFactor

# 3. Exploratory Factor Analysis (EFA) with principal axis factoring and Varimax rotation
efa_final <- fa(belief_vars_imputed, nfactors = n_factors, fm = "pa", rotate = "varimax")
print(efa_final$loadings)  # rotated factor loadings
print(efa_final$values)    # eigenvalues

#PATH ANALYSIS W/ PCA FACTOR#
#OVERALL NEGATIVE SYMPTOMS
#model f1 - factor - sans - psp
modelf1 = '
  zsans ~ a*zFactor
  zpsp ~ b*zsans
  
  factor.through.sans:=a*b
   
'
# Estimate path model
path.analysisf1 = modelf1 %>%
  sem(data = dataset)
# Print results, including the fit statistics
path.analysisf1 %>%
  summary(fit.measures = TRUE, 
          standardized = TRUE,
          ci = T)

#DIMINISHED MOTIVATION
#model f1m - factor - dmot - psp
modelf1m = '
  zsans_mot ~ a*zFactor
  zpsp ~ b*zsans_mot
  
  factor.through.sans:=a*b
   
'
# Estimate path model
path.analysisf1m = modelf1m %>%
  sem(data = dataset)
# Print results, including the fit statistics
path.analysisf1m %>%
  summary(fit.measures = TRUE, 
          standardized = TRUE,
          ci = T)

#FULL SAMPLE - OVERALL NEGATIVE SYMPTOMS - INDIVIDUAL PATHS

### ---------
### 1. Define all models
### ---------
model1 <- '
  zsans ~ a*zdpb
  zpsp ~ b*zsans
  das.through.sans := a*b
'

model2 <- '
  zsans ~ a*zabs
  zpsp ~ b*zsans
  abs.through.sans := a*b
'

model3 <- '
  zsans ~ a*zrses
  zpsp ~ b*zsans
  rses.through.sans := a*b
'

model4 <- '
  zsans ~ a*zplex
  zpsp ~ b*zsans
  plex.through.sans := a*b
'

model5 <- '
  zsans ~ a*zismi
  zpsp ~ b*zsans
  ismi.through.sans := a*b
'

model6 <- '
  zsans ~ a*zsaraq
  zpsp ~ b*zsans
  saraq.through.sans := a*b
'

### ---------
### 2. Run direct-effect models (default SEs)
### ---------
direct_models <- list(
  DAS   = sem(model1, data = dataset),
  ABS   = sem(model2, data = dataset),
  RSES  = sem(model3, data = dataset),
  PLEX  = sem(model4, data = dataset),
  ISMI  = sem(model5, data = dataset),
  SARAQ = sem(model6, data = dataset)
)

### ---------
### 3. Run indirect-effect models (bootstrap SEs)
### ---------
boot_models <- list(
  DAS   = sem(model1, data = dataset, se = "bootstrap", bootstrap = 1000),
  ABS   = sem(model2, data = dataset, se = "bootstrap", bootstrap = 1000),
  RSES  = sem(model3, data = dataset, se = "bootstrap", bootstrap = 1000),
  PLEX  = sem(model4, data = dataset, se = "bootstrap", bootstrap = 1000),
  ISMI  = sem(model5, data = dataset, se = "bootstrap", bootstrap = 1000),
  SARAQ = sem(model6, data = dataset, se = "bootstrap", bootstrap = 1000)
)

### ---------
### 4. Functions to extract effects
### ---------
get_direct <- function(fit, name) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  direct <- subset(pe, op == "~")
  tibble(
    Model = name,
    EffectType = "Direct",
    Path = paste(direct$lhs, "~", direct$rhs),
    Estimate = direct$est.std,
    SE = direct$se,
    p = direct$pvalue
  )
}

get_indirect <- function(fit, name) {
  pe <- parameterEstimates(fit, standardized = TRUE, boot.ci.type = "perc")
  indirect <- subset(pe, op == ":=")
  tibble(
    Model = name,
    EffectType = "Indirect",
    Path = indirect$lhs,
    Estimate = indirect$est.std,
    SE = indirect$se,
    p = indirect$pvalue
  )
}

### ---------
### 5. Collect results
### ---------
direct_results   <- bind_rows(Map(get_direct, direct_models, names(direct_models)))
indirect_results <- bind_rows(Map(get_indirect, boot_models, names(boot_models)))

all_results <- bind_rows(direct_results, indirect_results)

### ---------
### 6. Apply selective FDR correction
### ---------
# Identify rows to FDR-correct
fdr_idx <- which(
  (all_results$EffectType == "Direct" & !grepl("zpsp ~ zsans", all_results$Path)) |
    (all_results$EffectType == "Indirect")
)

# Initialize FDR_p column with raw p-values
all_results$FDR_p <- all_results$p

# Replace only the selected rows with FDR-adjusted p-values
all_results$FDR_p[fdr_idx] <- p.adjust(all_results$p[fdr_idx], method = "fdr")

### ---------
### 7. View final table
### ---------
print(all_results)



#DIMINISHED MOTIVATION#
### ---------
### 1. Define all models (diminished motivation)
### ---------
model1m <- '
  zsans_mot ~ a*zdpb
  zpsp ~ b*zsans_mot
  das.through.sans_mot := a*b
'

model2m <- '
  zsans_mot ~ a*zabs
  zpsp ~ b*zsans_mot
  abs.through.sans_mot := a*b
'

model3m <- '
  zsans_mot ~ a*zrses
  zpsp ~ b*zsans_mot
  rses.through.sans_mot := a*b
'

model4m <- '
  zsans_mot ~ a*zplex
  zpsp ~ b*zsans_mot
  plex.through.sans_mot := a*b
'

model5m <- '
  zsans_mot ~ a*zismi
  zpsp ~ b*zsans_mot
  ismi.through.sans_mot := a*b
'

model6m <- '
  zsans_mot ~ a*zsaraq
  zpsp ~ b*zsans_mot
  saraq.through.sans_mot := a*b
'

### ---------
### 2. Run direct-effect models (default SEs)
### ---------
direct_models <- list(
  DAS   = sem(model1m, data = dataset),
  ABS   = sem(model2m, data = dataset),
  RSES  = sem(model3m, data = dataset),
  PLEX  = sem(model4m, data = dataset),
  ISMI  = sem(model5m, data = dataset),
  SARAQ = sem(model6m, data = dataset)
)

### ---------
### 3. Run indirect-effect models (bootstrap SEs)
### ---------
boot_models <- list(
  DAS   = sem(model1m, data = dataset, se = "bootstrap", bootstrap = 1000),
  ABS   = sem(model2m, data = dataset, se = "bootstrap", bootstrap = 1000),
  RSES  = sem(model3m, data = dataset, se = "bootstrap", bootstrap = 1000),
  PLEX  = sem(model4m, data = dataset, se = "bootstrap", bootstrap = 1000),
  ISMI  = sem(model5m, data = dataset, se = "bootstrap", bootstrap = 1000),
  SARAQ = sem(model6m, data = dataset, se = "bootstrap", bootstrap = 1000)
)

### ---------
### 4. Functions to extract effects
### ---------
get_direct_m <- function(fit, name) {
  pe <- parameterEstimates(fit, standardized = TRUE)
  direct <- subset(pe, op == "~")
  tibble(
    Model = name,
    EffectType = "Direct",
    Path = paste(direct$lhs, "~", direct$rhs),
    Estimate = direct$est.std,
    SE = direct$se,
    p = direct$pvalue
  )
}

get_indirect_m <- function(fit, name) {
  pe <- parameterEstimates(fit, standardized = TRUE, boot.ci.type = "perc")
  indirect <- subset(pe, op == ":=")
  tibble(
    Model = name,
    EffectType = "Indirect",
    Path = indirect$lhs,
    Estimate = indirect$est.std,
    SE = indirect$se,
    p = indirect$pvalue
  )
}

### ---------
### 5. Collect results
### ---------
direct_results_m   <- bind_rows(Map(get_direct_m, direct_models, names(direct_models)))
indirect_results_m <- bind_rows(Map(get_indirect_m, boot_models, names(boot_models)))

all_results_m <- bind_rows(direct_results_m, indirect_results_m)

### ---------
### 6. Apply selective FDR correction
### ---------
# Identify rows to FDR-correct (exclude zpsp ~ zsans_mot)
fdr_idx <- which(
  (all_results_m$EffectType == "Direct" & !grepl("zpsp ~ zsans_mot", all_results_m$Path)) |
    (all_results_m$EffectType == "Indirect")
)

# Initialize FDR_p column with raw p-values
all_results_m$FDR_p <- all_results_m$p

# Replace only the selected rows with FDR-adjusted p-values
all_results_m$FDR_p[fdr_idx] <- p.adjust(all_results_m$p[fdr_idx], method = "fdr")

### ---------
### 7. View final table
### ---------
print(all_results_m)

#EXPLORATORY OMNIBUS PATH MODELS#
#OVERALL NEGATIVE SYMPTOMS
#model 0 - beliefs - sans - psp
model0 = '
  zsans ~ a*zdpb + b*zabs + c*zrses + d*zplex + e*zismi + f*zsaraq
  zpsp ~ g*zsans
  
  zdpb ~~ zabs
  zdpb ~~ zrses
  zdpb ~~ zplex
  zdpb ~~ zismi
  zdpb ~~ zsaraq
  zabs ~~ zrses
  zabs ~~ zplex
  zabs ~~ zismi
  zabs ~~ zsaraq
  zrses ~~ zplex 
  zrses ~~ zismi
  zrses ~~ zsaraq
  zplex ~~ zismi
  zplex ~~ zsaraq
  zismi ~~ zsaraq

  das.through.sans:=a*g
  abs.through.sans:=b*g
  rses.through.sans:=c*g
  plex.through.sans:=d*g
  ismi.through.sans:=e*g
  saraq.through.sans:=f*g
  
'
# Estimate path model
path.analysis0 = model0 %>%
  sem(data = dataset)
# Print results, including the fit statistics
path.analysis0 %>%
  summary(fit.measures = TRUE, 
          standardized = TRUE,
          ci = TRUE)
# Print results for bootstrapped indirect effects
path.analysis0 <- sem(
  model0,
  data = dataset,
  se = "bootstrap",
  bootstrap = 1000 )

# Fit a regression model to calculate VIF for each belief variable
belief_model <- lm(zsans ~ zdpb + zabs + zrses + zplex + zismi + zsaraq, data = dataset)

# Calculate VIF for each predictor
vif(belief_model)

#DIMINISHED MOTIVATION
#model 0m - sans-mot - beliefs - psp
model0m = '
  zsans_mot ~ a*zdpb + b*zabs + c*zrses + d*zplex + e*zismi + f*zsaraq
  zpsp ~ g*zsans_mot
  
  zdpb ~~ zabs
  zdpb ~~ zrses
  zdpb ~~ zplex
  zdpb ~~ zismi
  zdpb ~~ zsaraq
  zabs ~~ zrses
  zabs ~~ zplex
  zabs ~~ zismi
  zabs ~~ zsaraq
  zrses ~~ zplex 
  zrses ~~ zismi
  zrses ~~ zsaraq
  zplex ~~ zismi
  zplex ~~ zsaraq
  zismi ~~ zsaraq

  das.through.sans:=a*g
  abs.through.sans:=b*g
  rses.through.sans:=c*g
  plex.through.sans:=d*g
  ismi.through.sans:=e*g
  saraq.through.sans:=f*g
  
'

# Estimate path model
path.analysis0m = model0m %>%
  sem(data = dataset)
# Print results, including the fit statistics
path.analysis0m %>%
  summary(fit.measures = TRUE, 
          standardized = TRUE,
          ci = TRUE)
# Print results for bootstrapped indirect effects
path.analysis0 <- sem(
  model0,
  data = dataset,
  se = "bootstrap",
  bootstrap = 1000 )
  
  
#CORRELATIONS B/W BELIEFS, NEG SX, FUNCTIONING, w/ FDR#
#FDR#
# Define belief variables (predictors)
belief_vars <- c("zdpb", "zabs", "zrses","zplex","zismi","zsaraq")

# Define sociodemographic, clinical, and other outcome variables
other_vars <- c("zsans", "zsans_mot", "zsans_exp", "zpsp")

# Extract only relevant variables
cor_data <- subset[, c(belief_vars, other_vars)]

# Initialize matrices to store results
cor_mat <- matrix(NA, nrow = length(belief_vars), ncol = length(other_vars),
                  dimnames = list(belief_vars, other_vars))
p_mat <- cor_mat

# Initialize matrix to store sample sizes
n_mat <- matrix(NA, nrow = length(belief_vars), ncol = length(other_vars),
                dimnames = list(belief_vars, other_vars))

# Loop through each pair and compute Pearson correlation and sample size
for (i in belief_vars) {
  for (j in other_vars) {
    x <- subset[[i]]
    y <- subset[[j]]
    complete_cases <- complete.cases(x, y)
    n_mat[i, j] <- sum(complete_cases)  # count non-missing pairs
    temp <- cor.test(x[complete_cases], y[complete_cases])
    cor_mat[i, j] <- temp$estimate
    p_mat[i, j] <- temp$p.value
  }
}

# FDR-corrected q-values
q_mat <- matrix(p.adjust(p_mat, method = "fdr"), nrow = nrow(p_mat), 
                dimnames = dimnames(p_mat))

# Significance matrix
sig_mat <- q_mat < 0.05

# Output results as a list
results <- list(
  correlations = cor_mat,
  unadjusted_p = p_mat,
  sample_sizes = n_mat,
  q_values = q_mat,
  significant = sig_mat
)