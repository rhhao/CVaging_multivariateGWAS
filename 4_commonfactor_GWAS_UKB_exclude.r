library(Matrix)
library(GenomicSEM)
library(psych)

load("LDSCoutput_CVAging_8traits_UKB_exclude.RData")

reference.file = "1000G_Phase3_EUR_reference_maf0.005.txt"

save.dir = "./"

files = c("../gwas_data/cardiovascular_diseases/cad.add.160614.website.formatted.txt.gz",
        "../gwas_data/cardiovascular_diseases/finngen_R12_I9_HEARTFAIL.formatted.gz",
        "../gwas_data/cardiovascular_diseases/finngen_R12_C_STROKE.formatted.gz",
        "../gwas_data/cardiovascular_diseases/MVP_R4.1000G_AGR.CircPVD.EUR.GIA.dbGaP.formatted.txt.gz",
        "../gwas_data/aging/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_MetaAnalysis_formatted.gz",
        "../gwas_data/aging/healthspan_summary_formatted.gz",
        "../gwas_data/aging/PhenoAge_EUR_summary_statistics_formatted.gz",
        "../gwas_data/aging/Results_90th_percentile_reverse_coded_formatted.gz")

trait.names = c("CAD","HF","Stroke","PAD","FI","Healthspan","PhenoAge","Longevity")

se.logit = c(F, T, T, T, F, F, F, F)
OLS = c(T, F, F, F, T, T, T, T)
linprob = c(F, F, F, F, F, F, F, F)

N = c(NA, NA, NA, NA, 175226, 300477, 34710, NA) #This is only needed for continuous outcomes or outcomes where linprob is TRUE
maf.filter = 0.01
keep.indel = FALSE
parallel = TRUE
cores = 8

CVaging_UKB_exclude_sumstat = sumstats(files=files,
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

save(CVaging_UKB_exclude_sumstat, file=paste(save.dir,"CVaging_commonfactor_UKB_exclude_sumstat.RData", sep=""))

Commonfactor_UKB_exclude = commonfactorGWAS(covstruc = LDSCoutput,
                                SNPs = CVaging_UKB_exclude_sumstat,
                                toler = 1e-50,
                                core = 20)
save(Commonfactor_UKB_exclude, file=paste(save.dir,"CVaging_commonfactor_UKB_exclude_multivariatesGWAS.RData", sep=""))

#calculate expected sample size
df_subset = subset(Commonfactor_UKB_exclude, Commonfactor_UKB_exclude$MAF >= .1 & Commonfactor_UKB_exclude$MAF <= .4)
N_eff_UKB_exclude = mean(1/((2*df_subset$MAF*(1-df_subset$MAF))*df_subset$SE^2))
cat("Effective sample size is:", N_eff_UKB_exclude, "\n")
