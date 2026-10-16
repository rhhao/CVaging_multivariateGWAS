library(dplyr)
args = commandArgs(T)
code = args[1]
output = args[2]


###
diagnose_date = read.table("ICD10_Date_2026_0307_long.txt",header=T)
birth_year = read.csv("34_Year_of_birth.csv", header=T)
birth_month = read.csv("52_Month_of_birth.csv",header=T)
df = merge(birth_year, birth_month, by="eid",all=F)
df$birth_date = paste0(df[,2],"-",df[,3],"-1")
Age = merge(diagnose_date,df,by="eid",all.x=T)
Age$diagnose_age=floor(as.numeric(as.Date(Age$Date)-as.Date(Age$birth_date))/365.25)
Age_filter = Age[which(Age[,4] %in% code), c("eid", "diagnose_age")]
Age_filter = unique(Age_filter)
Age_filter = Age_filter %>%
			group_by(eid) %>%
			filter(diagnose_age == min(diagnose_age)) %>%
			ungroup()
write.csv(Age_filter, output, row.names=F, quote=F)
