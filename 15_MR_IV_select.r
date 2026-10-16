library(dplyr)
library(data.table)
library(LDlinkR)
library(BiocParallel)

args = commandArgs(T)
options(warn = -1)

threads = 5
param = MulticoreParam(workers = threads, progressbar = TRUE)
################################################
query = function(snp_num, exp_ma_data, confounders, token){
   snp1 = exp_ma_data$SNP2[snp_num]
   snp2 = exp_ma_data$SNP[snp_num]
   scan_res = tryCatch(
   {
      LDtrait(snp1, pop="EUR", r2d = "r2", r2d_threshold = 0.8, win_size = 5e+05, token = token, file = FALSE, genome_build = "grch37", api_root = "https://ldlink.nih.gov/LDlinkRest")
   },
   error = function(e){
      return(NULL)}
   )
   if(!is.null(scan_res)){
      confounder_snps = scan_res[grep(paste0("\\d(",paste(confounders,collapse="|"),")\\d"),tolower(scan_res$GWAS_Trait)),]
      if(nrow(confounder_snps) != 0){
         return(snp2)
      }
   }
}

#---------------------------
outdir = args[1]
gwas_file = args[2]
token = args[3]

if (!dir.exists(outdir)){
    dir.create(outdir)
}
confounders = c("smoking","alcohol","education","college","university","smoke","drink","drinking")

exp_ma_data = fread(file = gwas_file, header = T, data.table = F, select = c("SNP", "P"))
exp_ma_data = exp_ma_data %>% unique()
exp_ma_data = exp_ma_data[exp_ma_data$P < 5e-8,]
if(nrow(exp_ma_data) == 0){
    confounder_snps=c()
    print("Notice: No availuable IVs!")
}else{
    if(sum(grepl("rs", exp_ma_data$SNP)) == 0){
        if(sum(grepl("chr", exp_ma_data$SNP)) == 0){
            exp_ma_data$SNP2 = lapply(exp_ma_data$SNP, function(x){paste0("chr",paste(strsplit(x, "[:_]")[[1]][1:2], collapse=":"))})
        }else{
            exp_ma_data$SNP2 = lapply(exp_ma_data$SNP, function(x){paste(strsplit(x, "[:_]")[[1]][1:2], collapse=":")})}
    }else{
        exp_ma_data$SNP2 = exp_ma_data$SNP
    }
    confounder_snps = bplapply(1:nrow(exp_ma_data), function(x) query(x, exp_ma_data=exp_ma_data, confounders=confounders, token=token), BPPARAM = param)
    confounder_snps = unlist(confounder_snps)
}
if(length(confounder_snps) != 0){
    exp_ma_data = exp_ma_data[-which(exp_ma_data$SNP %in% confounder_snps),]
}
snplist = exp_ma_data$SNP
write.table(unique(snplist), file = paste0(outdir, "/", gwas_file,".IV"), sep="\n", quote = F, col.names = F, row.names = F)
