ihs = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/Standardised iHS/VK7.txt')

h12 = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/H12/VK7.txt')

##Define a 10kb tolerance for position
tolerance = 10000/1e6

chr = c('AgamP4_2L','AgamP4_2R','AgamP4_3R','AgamP4_3L','AgamP4_X')

library(dplyr)
ihs <- ihs %>%
  readr::type_convert()
h12 <- h12%>%
  readr::type_convert()

library(fuzzyjoin)

#Combine H12 and iHHS by position with the tolerance defined above
combined = c()
for(i in 1:length(chr))
{
  is_true = which(ihs$pop == chr[i])
  chr_ihs = ihs[is_true,]
  
  is_true = which(h12$pop == chr[i])
  chr_h12 = h12[is_true,]
  
  out = fuzzy_join(chr_h12, chr_ihs, by = "pos",
             match_fun = ~ abs(.x - .y) < tolerance)

  out = subset(out,!duplicated(out$id.x))
  
  combined = rbind(combined,out)
}


##Split these into blocks using 100kb
combined$group = rep(0,nrow(combined))
for(i in 2:(nrow(combined))-1)
{
  if(as.numeric(combined[i+1,2])-as.numeric(combined[i,2]) > 0.1 )
  {
    combined[i,19] = 'SPLIT'
  }
}


write.table(combined,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/iHS and H12 overlap/VK7.txt',sep='\t',row.names=F)


##Investigate each block for peak##

groups = which(combined$group == 'SPLIT')

groups = c(0,groups,nrow(combined))

out_vector = c()
for(i in 1:(length(groups)-1))
{
  subset_group = combined[c((groups[i]+1):groups[i+1]),]
  
  peak_snp = subset_group[which(subset_group$ihh12 == max(subset_group$ihh12) ),]
  peak_snp = rbind(peak_snp,subset_group[which(subset_group$ihs == max(subset_group$ihs) ),])
  
  peak_snp$group = i
  
  
  out_vector = rbind(out_vector,peak_snp)
}

write.table(out_vector,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/iHS and H12 overlap/VK7_peak_SNPs.txt',sep='\t',row.names=F)

