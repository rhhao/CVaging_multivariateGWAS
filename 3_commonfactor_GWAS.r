library(Matrix)
library(GenomicSEM)
library(psych)

load("LDSCoutput_CVaging_8traits_revCoded.RData")

reference.file = "1000G_Phase3_EUR_reference_maf0.005.txt"

save.dir = "./"

files = c("../gwas_data/aging/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_MetaAnalysis_formatted.gz",
        "../gwas_data/aging/healthspan_summary_formatted.gz",
        "../gwas_data/aging/PhenoAge_EUR_summary_statistics_formatted.gz",
        "../gwas_data/aging/Results_90th_percentile_reverse_coded_formatted.gz",
        "../gwas_data/cardiovascular_diseases/GCST90132314_buildGRCh37_formatted.gz",
        "../gwas_data/cardiovascular_diseases/FORMAT-METAL_Pheno1_EUR.formatted.gz",
        "../gwas_data/cardiovascular_diseases/GCST90104539_buildGRCh37_formatted.gz",
        "../gwas_data/cardiovascular_diseases/pad.vanzuydametal.2020.formatted.gz")

trait.names = c("FI","Healthspan","PhenoAge","Longevity","CAD","HF","Stroke","PAD")

se.logit = c(F, F, F, F, F, F, F, F)
OLS = c(T, T, T, T, T, T, T, T)
linprob = c(F, F, F, F, F, F, F, F)

N = c(175226, 300477, 34710, NA, NA, NA, NA, NA) #This is only needed for continuous outcomes or outcomes where linprob is TRUE
maf.filter = 0.01
keep.indel = FALSE
parallel = TRUE
cores = 8

CVaging_sumstat = sumstats(files=files,
                        ref=reference.file,
                        trait.names=trait.names,
                        se.logit=se.logit,
                        OLS=OLS,
                        linprob=linprob,
                        N=N,
                        maf.filter=maf.filter,
                        keep.indel=keep.indel,
                        parallel=parallel,
                        cores=cores)

save(CVaging_sumstat, file=paste(save.dir,"CVaging_commonfactor_sumstat.RData", sep=""))

#############
#load("CVaging_commonfactor_sumstat.RData")
Commonfactor = commonfactorGWAS(covstruc = LDSCoutput,
                                SNPs = CVaging_sumstat,
                                toler = 1e-40,
                                core = 20)

save(Commonfactor, file=paste(save.dir,"CVaging_commonfactor_multivariatesGWAS.RData", sep=""))

#calculate effective sample size
df_subset = subset(Commonfactor, Commonfactor$MAF >= .1 & Commonfactor$MAF <= .4)
N_eff = mean(1/((2*df_subset$MAF*(1-df_subset$MAF))*df_subset$SE^2))
cat("Effective sample size is:", N_eff, "\n")
