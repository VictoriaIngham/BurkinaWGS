##File path and population of interest windows###
Gao =  c('ABBABABA_TiaBakGao1.csv','ABBABABA_VK7TieGao2.csv','ABBABABA_TenBanGao3.csv')
Tia =  c('ABBABABA_GaoBakTia4.csv','ABBABABA_VK7TieTia5.csv','ABBABABA_TenBanTia6.csv')
Bak = c('ABBABABA_TiaGaoBak7.csv','ABBABABA_VK7TieBak8.csv','ABBABABA_TenBanBak9.csv')
VK = c('ABBABABA_TiaBakVK710.csv','ABBABABA_GaoTieVK711.csv','ABBABABA_TenBanVK12.csv')
Tie = c('ABBABABA_TiaBakTie13.csv','ABBABABA_GaoVK7Tie14.csv','ABBABABA_TenBanTie15.csv')
Ten = c('ABBABABA_TiaBakTen16.csv','ABBABABA_GaoVK7Ten17.csv','ABBABABA_TieBanTen18.csv')
Ban = c('ABBABABA_TiaBakBan19.csv','ABBABABA_GaoVK7Ban20.csv','ABBABABA_TieTenBan21.csv')

###Create dataframes for each###
files = VK
names = gsub("(?:[^.]+\\.){0}([^.]+).*", "\\1", files)

names2 = substring(names,10)

all = c()
for(i in 1:length(files))
{
  filepath <- file.path("/Volumes/Elements/ABABABA/ABABABA",files[i])
  foi = read.csv(filepath,
                   colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"),
                   )
  foi$pop = as.factor(rep(names2[i],nrow(foi)))
  all = rbind(all,foi)
}
rm(foi)

all = subset(all,scaffold != 'AgamP4_Mt' )
all = subset(all,scaffold != 'AgamP4_UNKN' )
all = subset(all,scaffold != '3L0|0:1249vcf:3L' )

all$pop = as.factor(all$pop)

all$POS=floor(rowMeans(all[,c('start', 'end')], na.rm=TRUE))/1e6

colours = c("#56B4E9","#009E73","#CC79A7")
library(ggplot2)
ggplot(all) + 
  geom_line(aes(x=POS,y=fdM,color=pop),alpha=0.8) +
  ggtitle('VK7 D Statistic') +
  facet_wrap(~scaffold, scales = "free")  +
  ylim(c(-1,1)) +
  scale_color_manual(values=colours) +
  labs(x="Mb" , y='D', col="Population") +
  theme_bw() 


