library(dplyr)
library(data.table)
library(mr.raps)
library(RadialMR)
library(MRPRESSO)
library(TwoSampleMR)

##################
add_snp = function(replace_snp_list, clump_out, exp_ma_data, outcome_ma_data, ref_path){
    additional_snps = c()
    for (snp in replace_snp_list){
        chr = clump_out[clump_out$SNP==snp,"CHR"] %>% as.character()
        plink_file = gsub("\\.bim","",list.files(pattern = paste0(".*\\.",chr,"\\D+.*bim$"), path = ref_path, ignore.case=T))
        system(paste0('plink --silent --bfile ', ref_path, '/', plink_file,' --r2 --ld-snp ', snp, ' --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out ', paste0("add_snp.",snp)))
        df = read.table(file = paste0("add_snp.",snp,".ld"),header = T)
        ld_snp = intersect(df$SNP_B, exp_ma_data$SNP)
        df2 = df[which(df$SNP_B %in% intersect(ld_snp, outcome_ma_data$SNP)),c("SNP_B","R2")]
        df2 = df2[order(df2$R2, decreasing = T), ]
        for(ld_snp in df2[, "SNP_B"]) {
            if (outcome_ma_data[which(outcome_ma_data$SNP == ld_snp),"BETA"] != 0){
                add_snp = ld_snp
                break
            }
        }
        if (exists("add_snp")){
            additional_snps = append(additional_snps, add_snp)
        }
        system(paste0("rm ", paste0("add_snp.",snp), "*"))
    }
    return(additional_snps)
}

radial_mr = function(exp_data, outcome_data){
    if (nrow(exp_data) <= 2){
        outliers = NULL
    }else{
        colnames(exp_data) = recode(colnames(exp_data), 'BETA'='beta.exposure', 'SE'='se.exposure', 'A1'='effect_allele.exposure', 'A2'='other_allele.exposure', 'MAF'='eaf.exposure', 'N'='samplesize.exposure')
        colnames(outcome_data) = recode(colnames(outcome_data), 'BETA'='beta.outcome', 'SE'='se.outcome', 'A1'='effect_allele.outcome', 'A2'='other_allele.outcome', 'MAF'='eaf.outcome', 'N'='samplesize.outcome')
        exp_data$'id.exposure' = "exposure"
        exp_data$exposure = "exposure"
        outcome_data$'id.outcome' = "outcome"
        outcome_data$outcome = "outcome"
        data = harmonise_data(exposure_dat = exp_data, outcome_dat = outcome_data)
        df = format_radial(data$beta.exposure,data$beta.outcome,data$se.exposure,data$se.outcome,data$SNP)
        fit1 = try(ivw_radial(df,0.05,1,0.0001), silent=T)
        if (!"try-error" %in% class(fit1)){
            res1 = ivw_radial(df,0.05,1,0.0001)
        }else{
            res1 = NA
        }
        fit2 = try(egger_radial(df,0.05,1), silent=T)
        if (!"try-error" %in% class(fit2)){
            res2 = egger_radial(df,0.05,1)
        }else{
            res2 = NA
        }
        try1 = try(res1$outliers$SNP, silent = T)
        try2 = try(res2$outliers$SNP, silent = T)
        if ("try-error" %in% class(try1) & "try-error" %in% class(try2)){
            outliers = NULL
        }else if("try-error" %in% class(try1) & !"try-error" %in% class(try2)){
            outliers = res2$outliers$SNP
        }else if(!"try-error" %in% class(try1) & "try-error" %in% class(try2)){
            outliers = res1$outliers$SNP
        }else{
            outliers = union(res1$outliers$SNP, res2$outliers$SNP)
        }
        return(outliers)
    }
}

calculate_F = function(exp_radial_data){
    R2 = sum((exp_radial_data$BETA)^2/((exp_radial_data$BETA)^2+(exp_radial_data$SE)^2*exp_radial_data$N))
    k = nrow(exp_radial_data)
    n = mean(exp_radial_data$N)
    F_statistic = (R2*(n-k-1))/((1-R2)*k)
    return(F_statistic)
}

MR = function(exp_radial_file, outcome_radial_file, exposure_id, outcome_id){
    exposure = read_exposure_data(
               filename = exp_radial_file,
               sep = "\t",
               snp_col = "SNP",
               beta_col = "BETA",
               se_col = "SE",
               effect_allele_col = "A1",
               other_allele_col = "A2",
               eaf_col = "MAF",
               pval_col = "P",
               samplesize_col = "N")

    outcome = read_outcome_data(
              filename = outcome_radial_file,
              sep = "\t",
              snp_col = "SNP",
              beta_col = "BETA",
              se_col = "SE",
              effect_allele_col = "A1",
              other_allele_col = "A2",
              eaf_col = "MAF",
              pval_col = "P",
              samplesize_col = "N")

    data = harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
    data$id.exposure = exposure_id
    data$id.outcome = outcome_id
    data$exposure = exposure_id
    data$outcome = outcome_id
    if (sum(data$beta.outcome == 0) != 0){
        data[which(data$beta.outcome == 0),"beta.outcome"] = 1e-10
    }
    if (sum(data$mr_keep) <= 2 & sum(data$mr_keep) >= 1){
        out = mr(data, method_list=c("mr_ivw", "mr_raps", "mr_wald_ratio"))
        out = generate_odds_ratios(out)
        out$id.exposure = exposure_id
        out$id.outcome = outcome_id
        mr_res = out[,c("id.exposure","id.outcome","method","nsnp","b","lo_ci","up_ci","or","or_lci95","or_uci95","pval")]
    }else if (sum(data$mr_keep) == 0){
        mr_res = NULL
    }else{
        out = mr(data, method_list=c("mr_weighted_median", "mr_ivw", "mr_raps", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_mode"))
        out = generate_odds_ratios(out)
        line1 = which(out$method == "MR Egger")
        out$lo_ci[line1] = out$b[line1] - (qt(0.05/2, out$nsnp[line1]-2, lower.tail = FALSE) * out$se[line1])
        out$up_ci[line1] = out$b[line1] + (qt(0.05/2, out$nsnp[line1]-2, lower.tail = FALSE) * out$se[line1])
        out$or_lci95[line1] = exp(out$lo_ci[line1])
        out$or_uci95[line1] = exp(out$up_ci[line1])
        line2 = which(out$method == "Weighted mode")
        out$lo_ci[line2] = out$b[line2] - (qt(0.05/2, out$nsnp[line2]-1, lower.tail = FALSE) * out$se[line2])
        out$up_ci[line2] = out$b[line2] + (qt(0.05/2, out$nsnp[line2]-1, lower.tail = FALSE) * out$se[line2])
        out$or_lci95[line2] = exp(out$lo_ci[line2])
        out$or_uci95[line2] = exp(out$up_ci[line2])
        out$id.exposure = exposure_id
        out$id.outcome = outcome_id
        mr_res = out[,c("id.exposure","id.outcome","method","nsnp","b","lo_ci","up_ci","or","or_lci95","or_uci95","pval")]
    }

    res = list(data = data, mr_res = mr_res)
    return(res)
}

sensitivity_test = function(harmonised_data){
    ###horizontal pleiotropy check
    bzx = harmonised_data$beta.exposure
    sebzx = harmonised_data$se.exposure
    bzy = harmonised_data$beta.outcome
    sebzy = harmonised_data$se.outcome
    df = data.frame(bzx,sebzx,bzy,sebzy)
    if(nrow(df) <= 3){
        MRPRESSO_RSSobs = "NA"
        MRPRESSO_P = "NA"
    }else{
        horizontal_pleiotropy_res = mr_presso(BetaOutcome = "bzy", BetaExposure = "bzx", SdOutcome= "sebzy", SdExposure = "sebzx", OUTLIERtest = TRUE, DISTORTIONtest= TRUE, data = df, NbDistribution = 1000 , SignifThreshold = 0.05)$`MR-PRESSO results`$`Global Test`
        MRPRESSO_RSSobs = horizontal_pleiotropy_res$RSSobs
        MRPRESSO_P = horizontal_pleiotropy_res$Pvalue
    }

    ###directional pleiotropy check
    directional_pleiotropy_res = mr_egger_regression(harmonised_data$beta.exposure, harmonised_data$beta.outcome, harmonised_data$se.exposure, harmonised_data$se.outcome)
    mr_egger_intercept = directional_pleiotropy_res$b_i
    mr_egger_intercept_se = directional_pleiotropy_res$se_i
    mr_egger_intercept_P = directional_pleiotropy_res$pval_i

    ###heterogeneity check
    hete = mr_heterogeneity(harmonised_data, method_list=c("mr_egger_regression","mr_ivw"))
    Cochrans_Q = hete[which(hete$method == "Inverse variance weighted"),"Q"]
    Cochrans_Q_P = hete[which(hete$method == "Inverse variance weighted"),"Q_pval"]

    sensitivity_res = data.frame("MRPRESSO_RSSobs" = MRPRESSO_RSSobs, "MRPRESSO_P" = MRPRESSO_P, "mr_egger_intercept" = mr_egger_intercept, "mr_egger_intercept_se" = mr_egger_intercept_se, "mr_egger_intercept_P" = mr_egger_intercept_P, "Cochrans_Q" = Cochrans_Q, "Cochrans_Q_P" = Cochrans_Q_P, check.names=F)

    return(sensitivity_res)
}

###########################

args = commandArgs(T)
exp_ma_file = args[1]
outcome_ma_file = args[2]
clump_file = args[3]
ref_path = args[4]

### 1. get input files

if (startsWith(outcome_ma_file, "CVaging")){
    outcome_id = "CVaging"
} else{
    outcome_id = strsplit(basename(outcome_ma_file), split="\\.")[[1]][1]
    outcome_id = sub("finngen_R12_", "", basename(outcome_id))
}
if (startsWith(exp_ma_file, "Cardio")){
    exposure_id = "CVaging"
} else {
    exposure_id = strsplit(basename(exp_ma_file), split="\\.")[[1]][1]
    exposure_id = sub("finngen_R12_", "", basename(exposure_id))
}
outcome_ma_data = fread(file = outcome_ma_file,header = T,data.table = F)

if (file.exists(paste0(getwd(),"/",exp_ma_file))){
    exp_ma_file = paste0(getwd(),"/",exp_ma_file)
}
if (file.exists(paste0(getwd(),"/",outcome_ma_file))){
    outcome_ma_file = paste0(getwd(),"/",outcome_ma_file)
}
if (file.exists(paste0(getwd(),"/",clump_file))){
    clump_file = paste0(getwd(),"/",clump_file)
}
if (file.exists(paste0(getwd(),"/",ref_path))){
    ref_path = paste0(getwd(),"/",ref_path)
}

outdir = paste0(dirname(clump_file), "/", exposure_id,"-",outcome_id)
if (!dir.exists(outdir)){
    dir.create(outdir)
}
setwd(outdir)

if (file.info(clump_file)$size == 0){
    file.create(paste0(exposure_id,".",outcome_id,".",cell_type,".MR.res"))
    file.create(paste0(exposure_id,".",outcome_id,".",cell_type,".MR.sensitivity"))
}else{
    exp_ma_data = fread(file = exp_ma_file, header = T, data.table = F)
    exp_ma_data = exp_ma_data[which(exp_ma_data$P < 5e-8), ]
    #colnames(exp_ma_data) = c("SNP","A1","A2","MAF","BETA","SE","P","N")
    clump = read.table(clump_file, header = T)
    snp_list = clump$SNP %>% as.character()

    ### 2. prepare exposure and outcome .ma files
    #Replace SNPs that do not exist in outcome with their proxy SNPs (r2>0.8) that exist.
    diff_snp = setdiff(snp_list, outcome_ma_data$SNP)
    if (length(diff_snp) != 0){
        additional_snps = add_snp(diff_snp, clump, exp_ma_data, outcome_ma_data, ref_path)
        snp_list = intersect(c(snp_list,additional_snps), outcome_ma_data[,"SNP"])
    }
    exp_data = exp_ma_data[match(snp_list, exp_ma_data$SNP) %>% na.omit,]
    outcome_data = outcome_ma_data[match(snp_list, outcome_ma_data$SNP) %>% na.omit,]
    rm(exp_ma_data, outcome_ma_data)

    ### 3. RadialMR filtering

    outliers = radial_mr(exp_data, outcome_data)
    if (!is.null(outliers)){
        snp_list = setdiff(snp_list, outliers)
    }
    if (length(snp_list) == 0){
        message(paste0("MR analysis stop!\nexposure ID: ", exposure_id, "\noutcome ID: ", outcome_id, "\nReason: All IVs are detected as outliers. Please check your data."))
        file.create(paste0(exposure_id,".",outcome_id,".MR.res"))
        file.create(paste0(exposure_id,".",outcome_id,".MR.sensitivity"))
    }else{
        exp_radial_data = exp_data[match(snp_list, exp_data$SNP) %>% na.omit,]
        outcome_radial_data = outcome_data[match(snp_list, outcome_data$SNP) %>% na.omit,]
        if (exposure_id == "CVaging") {
            exp_radial_data$N = 1783422
        }else{
            outcome_radial_data$N = 1783422
        }
        F_statistic = calculate_F(exp_radial_data)
        if (F_statistic < 10){
            message(paste0("MR analysis stop!\nexposure ID: ", exposure_id, "\noutcome ID: ", outcome_id, "\nReason: Abnormal bias detected from IVs (F_statistic < 10)!"))
            file.create(paste0(exposure_id,".",outcome_id,".",cell_type,".MR.res"))
            file.create(paste0(exposure_id,".",outcome_id,".",cell_type,".MR.sensitivity"))
        }else{
            write.table(file = paste0(exposure_id,".MR.radial.ma"), exp_radial_data, row.names = F, col.names = T,quote = F,sep = "\t")
            write.table(file = paste0(outcome_id,".MR.radial.ma"), outcome_radial_data, row.names = F, col.names = T,quote = F,sep = "\t")

            ### 4. MR analysis

            result = MR(exp_radial_file = paste0(exposure_id,".MR.radial.ma"), outcome_radial_file = paste0(outcome_id,".MR.radial.ma"), exposure_id = exposure_id, outcome_id = outcome_id)
            write.table(result$mr_res, file = paste0(exposure_id,".",outcome_id,".MR.res"), row.names = F, col.names = T,quote = F,sep = "\t")

            write.table(result$data$SNP[which(result$data$mr_keep == "TRUE")], file = paste0(exposure_id,".",outcome_id,".snps"), row.names = F, col.names = F,quote = F,sep = "\n") ### Output instrumental snps
            p1 = mr_scatter_plot(result$mr_res, result$data)

            ### 5. sensitivity test

            harmonised_data = result$data
            sensitivity = sensitivity_test(harmonised_data)
            pdf(paste0(exposure_id,".",outcome_id,".MR.plot.pdf"))
            print(p1)
            dev.off()
            sensitivity["F_statistic"] = F_statistic

            write.table(sensitivity, file = paste0(exposure_id,".",outcome_id,".MR.sensitivity"), row.names = F, col.names = T,quote = F,sep = "\t")
            system("rm *.MR.radial.ma")
        }
    }
}

