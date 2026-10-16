library(boot)
library(DescTools)
library(stringr)
set.seed(1)


fit_logistic <- function(phenotype, prs, data, covariates) {
  covar_str <- paste(covariates, collapse = " + ")
  null_formula <- as.formula(paste(phenotype, "~", covar_str))
  prs_formula <- as.formula(paste(phenotype, "~", prs, "+", covar_str))
  
  null_model <- glm(null_formula, family = binomial(), data = data, na.action = na.omit)
  prs_model <- glm(prs_formula, family = binomial(), data = data, na.action = na.omit)
  
  coef_summary <- summary(prs_model)$coefficients
  beta <- coef_summary[2, 1]
  se   <- coef_summary[2, 2]
  pval <- coef_summary[2, 4]
  
  z <- qnorm(0.975)
  beta_lower <- beta - z * se
  beta_upper <- beta + z * se
  
  OR <- exp(beta)
  OR_lower <- exp(beta_lower)
  OR_upper <- exp(beta_upper)
  
  null_r2 <- PseudoR2(null_model, which = r2_type)
  prs_r2  <- PseudoR2(prs_model, which = r2_type)
  inc_r2  <- prs_r2 - null_r2
  
  n_obs <- nobs(prs_model)
  
  result <- data.frame(
    Phenotype = phenotype,
    PRS = prs,
    N_obs = n_obs,
    Beta = beta,
    SE = se,
    P = pval,
    Beta_lower = beta_lower,
    Beta_upper = beta_upper,
    OR = OR,
    OR_lower = OR_lower,
    OR_upper = OR_upper,
    Null_R2 = null_r2,
    PRS_R2 = prs_r2,
    Inc_R2 = inc_r2
  )
  return(result)
}

boot_incR2 <- function(data, indices, formula_pair) {
  d <- data[indices, ]
  formulas <- unlist(strsplit(formula_pair, ";"))
  null_form <- as.formula(formulas[1])
  prs_form  <- as.formula(formulas[2])
  
  null_fit <- glm(null_form, family = binomial(), data = d, na.action = na.omit)
  prs_fit  <- glm(prs_form,  family = binomial(), data = d, na.action = na.omit)
  
  r2_null <- PseudoR2(null_fit, which = r2_type)
  r2_prs  <- PseudoR2(prs_fit,  which = r2_type)
  return(r2_prs - r2_null)
}


args <- commandArgs(T)
PRSfile <- args[1]
covfile <- args[2]
phenofile <- args[3]
trait <- args[4]


PRScore <- read.table(PRSfile, header=T)
covs <- read.table(covfile, header=T)
pheno <- read.table(phenofile, header=T)
df <- merge(PRScore, covs, by=c("FID", "IID"), all=F)
df <- merge(df, pheno, by=c("FID", "IID"), all=F)
df$PRS_scale <- scale(df$PRS)

prs_name <- "PRS_scale"
phenotype <- "Phenotype"
covariates <- c("Sex", "Age", 
                "PC1", "PC2", "PC3", "PC4", "PC5",
                "PC6", "PC7", "PC8", "PC9", "PC10")

r2_type <- "McFaddenAdj"
n_boot <- 1000
ci_type <- "perc"
result_file <- paste0(trait, "_increment_R2_result.csv")

df <- df[!(is.na(df[[phenotype]]) | is.na(df[[prs_name]])), ]
df[[phenotype]] <- as.factor(df[[phenotype]])
df[["Sex"]] <- as.factor(df[["Sex"]])

covar_str <- paste(covariates, collapse = " + ")
null_formula_str <- paste(phenotype, "~", covar_str)
prs_formula_str  <- paste(phenotype, "~", prs_name, "+", covar_str)
formula_pair <- paste(null_formula_str, prs_formula_str, sep = ";")

main_result <- fit_logistic(phenotype, prs_name, df, covariates)

boot_obj <- boot(data = df, statistic = boot_incR2, R = n_boot,
                 formula_pair = formula_pair)
saveRDS(boot_obj$t, paste0(trait,"_incrementalR2_boot_1000.rds"))
boot_ci <- boot.ci(boot_obj, type = ci_type, index = 1)

inc_r2_lower <- boot_ci$percent[4]
inc_r2_upper <- boot_ci$percent[5]

main_result$Inc_R2_lower <- inc_r2_lower
main_result$Inc_R2_upper <- inc_r2_upper
main_result$Phenotype <- trait

write.csv(main_result, file = result_file, row.names = FALSE, quote=F)
