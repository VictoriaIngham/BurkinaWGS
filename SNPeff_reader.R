data = read.delim('/Volumes/Elements/Fst_ABBA/TEPs.txt')

library(dplyr)
library(Rfast)
data <- data %>%
  replace(data == "1|1", "1")

data <- data %>%
  replace(data == "0|1", "0.5")

data <- data %>%
  replace(data == "1|0", "0.5")

data <- data %>%
  replace(data == "0|0", "0")

pops = c('Tiassale','Bakaridjan','Gaoua','VK7','Tiefora','Tengrela','Banfora_Res','Banfora_Or','Banfora_Sus')

out=c()
for(i in 1:length(pops))
{
  subset = data[,grep(pops[i],colnames(data))]
  subset = as.matrix(subset)
  mode(subset) = 'numeric'
  total <- rowsums(subset)
  total = (total/ncol(subset)*100)
  out = cbind(out,total)
}

colnames(out) = pops

total_freq = rowsums(out)

out = cbind(data[,1:9],out)

out = out[which(total_freq >= 10),]

write.table(out,'/Users/vickyingham/Seafile/Papers/2023/WGS/Making tables/MODERATE_impact_frequencies.txt',sep='\t',row.names=F)
