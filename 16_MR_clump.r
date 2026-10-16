library(dplyr)
library(data.table)

Clump <- function(exp_ma_file,ref_path,duplication_path,clump_out){
    if(duplication_path != "None"){
        for (n in 1:22){
            dup_snp_file = list.files(pattern = paste0("chr",n,"(\\D+|$)"), path = duplication_path, ignore.case=T)
            plink_file = gsub("\\.bim","",list.files(pattern = paste0("chr",n,"\\D+.*bim$"), path = ref_path, ignore.case=T))
            system(paste0('plink --silent --bfile ', ref_path, '/', plink_file,' --clump ',exp_ma_file,' --clump-field P --clump-kb 10000 --clump-p1 5e-8 --clump-p2 0.01 --clump-r2 0.001 --out ',clump_out, '.chr', n,' --exclude ', duplication_path,'/',dup_snp_file))
        }
    }else{
        for (n in 1:22){
            plink_file = gsub("\\.bim","",list.files(pattern = paste0("\\D+",n,"\\..*bim$"), path = ref_path, ignore.case=T))
            system(paste0('plink --silent --bfile ', ref_path, '/', plink_file,' --clump ',exp_ma_file,' --clump-field P --clump-kb 10000 --clump-p1 5e-8 --clump-p2 0.01 --clump-r2 0.001 --out ',clump_out, '.chr', n))
        }
    }
}
#######################
args = commandArgs(T)

exp_ma_file = args[1]
snp_file = args[2]
outdir = args[3]
ref_path = args[4]
duplication_path = "None"

exp_ma_data = fread(file = exp_ma_file,header = T,data.table = F)
if (file.info(snp_file)$size == 0){
    file.create(paste0(outdir,'/',gsub("\\.IV", ".clumped", basename(snp_file))))
}else{
    snp_list = read.table(file = snp_file,header = F)[,1]
    match = match(snp_list, exp_ma_data$SNP) %>% na.omit
    exp_ma_data = exp_ma_data[match,]
    write.table(exp_ma_data,file = paste0(outdir,'/',gsub("\\.IV", ".ma", basename(snp_file))),row.names = F,col.names = T,sep = "\t",quote = F)
    Clump(exp_ma_file = paste0(outdir,'/',gsub("\\.IV", ".ma", basename(snp_file))), ref_path = ref_path, duplication_path = duplication_path, clump_out = paste0(outdir, '/', gsub("\\.IV", ".clump", basename(snp_file))))
    df = data.frame()
    files = list.files(outdir, pattern=paste0("^",gsub("\\.IV", "", basename(snp_file)),"\\.clump.chr.+\\.clumped$"))
    for (file in files){
        data = read.table(paste0(outdir,"/",file), header = T)
        df = rbind(df,data)
    }
    write.table(df, file = paste0(outdir,'/',gsub("\\.IV", ".clumped", basename(snp_file))), row.names = F,col.names = T,sep = "\t",quote = F)

    for (file in list.files(outdir, pattern=paste0("^",gsub("\\.IV", "", basename(snp_file)),"\\.clump\\.chr[1-9][0-9]*.*|","^",gsub("\\.IV", "", basename(snp_file)),"\\.ma$"))){
        file.remove(paste0(outdir, "/", file))
    }
}
