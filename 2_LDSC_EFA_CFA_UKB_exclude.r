library(Matrix)
library(GenomicSEM)
library(psych)

##########################################################
######### estimate genetic correlation matrix
################################################################
traits = c("../gwas_data/cardiovascular_diseases/ldsc_munge_out/CAD_UKB_exclude.sumstats.gz",
		"../gwas_data/cardiovascular_diseases/ldsc_munge_out/Heart_Failure_UKB_exclude.sumstats.gz",
		"../gwas_data/cardiovascular_diseases/ldsc_munge_out/Stroke_UKB_exclude.sumstats.gz",
		"../gwas_data/cardiovascular_diseases/ldsc_munge_out/PAD_UKB_exclude.sumstats.gz",
		"../gwas_data/aging/ldsc_munge_out/FI.sumstats.gz",
		"../gwas_data/aging/ldsc_munge_out/Healthspan.sumstats.gz",
		"../gwas_data/aging/ldsc_munge_out/PhenoAge.sumstats.gz",
		"../gwas_data/aging/ldsc_munge_out/Longevity_rev.sumstats.gz")
sample.prev = c(0.31,0.08,0.13,0.05,NA,NA,NA,NA)
population.prev = c(0.08,0.02,0.02,0.02,NA,NA,NA,NA)
##the folder of LD scores
ld = "/home/haoruohan/software/ldsc/eur_w_ld_chr/"

##the folder of LD weights
wld = "/home/haoruohan/software/ldsc/eur_w_ld_chr/"
##name the traits
trait.names = c("CAD","HF","Stroke","PAD","FI","Healthspan","PhenoAge","Longevity")
##run LDSC
LDSCoutput = ldsc(traits=traits,
			sample.prev=sample.prev,
			population.prev=population.prev,
			ld=ld,
			wld=wld,
			trait.names=trait.names,
			ldsc.log="LDSCoutput_CVAging_8traits_UKB_exclude")
save(LDSCoutput,file="LDSCoutput_CVAging_8traits_UKB_exclude.RData")

S_LD = as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)

#####################################################
#### EFA
##########################################

EFA = factanal(covmat=S_LD, factors=1, rotation="promax")
print(EFA)
save(EFA, file=paste0("CVaging_EFA_promax_", i, "factors_UKB_exclude.RData"))

######################################## 
####CFA
########################################

model_single_factor = '
  CardioAge =~ NA*CAD + HF + Stroke + PAD + FI + Healthspan + PhenoAge + Longevity
  CardioAge ~~ 1*CardioAge
'
CFA = usermodel(
  covstruc = LDSCoutput,
  model = model_single_factor,
  estimation = "DWLS",
  CFIcalc = TRUE,
  std.lv = FALSE,
  imp_cov = TRUE
)
print(CFA)

save(CFA, file=paste0("CVaging_CFA_commonfactor_UKB_exclude.RData"))

