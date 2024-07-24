##File path and population of interest windows###
files = list.files('/Volumes/Elements/Windowed_Fst/Burkina_pops',pattern='*.fst')

population = 'Gao'

###Create dataframes for each###
files = files[grep(population, files)]
names = gsub("(?:[^.]+\\.){0}([^.]+).*", "\\1", files)

all = c()
for(i in 1:length(files))
{
  filepath <- file.path("/Volumes/Elements/Windowed_Fst/Burkina_pops",files[i])
  foi = read.delim(filepath,
             colClasses=c("character","numeric","numeric","numeric","numeric","numeric"),
             sep = "\t")
  foi$pop = as.factor(rep(names[i],nrow(foi)))
  all = rbind(all,foi)
}
rm(foi)



###Hash out the population you're searching for###
all$pop2 = rep(NA,nrow(all))
all$pop2[grep('Bak',all$pop)] = 'Bakaridjan'
all$pop2[grep('Ban',all$pop)] = 'Banfora'
all$pop2[grep('Tia',all$pop)] = 'Tiassale'
all$pop2[grep('Tie',all$pop)] = 'Tiefora'
#all$pop2[grep('G',all$pop)] = 'Gaoua'
all$pop2[grep('VK',all$pop)] = 'VK7'
all$pop2[grep('Ten',all$pop)] = 'Tengrela'


all$pop2 = as.factor(all$pop2)

#fix distance to be in Mb
all$POS=floor(rowMeans(all[,c('BIN_START', 'BIN_END')], na.rm=TRUE))/1e6

all2 = subset(all,CHROM != 'AgamP4_Mt' )
all2 = subset(all2,CHROM != 'AgamP4_UNKN' )


group.colours = c(Banfora='#D41159',Bakaridjan='#0C7BDC',Tiassale='#E66100',Tiefora = '#5D3A9B',Gaoua = '#40B0A6',VK7 = '#000000',Tengrela='#1AFF1A')


library(ggplot2)
ggplot(all2) + 
  geom_smooth(aes(x=POS,y=WEIGHTED_FST,color=pop2),method='loess',span = 0.01) +
  ggtitle('Bakaridjan') +
  facet_wrap(~CHROM, scales = "free")  +
  labs(x="Mb" , y=expression(F[ST]), col="Population") +
  scale_color_manual(values=group.colours) +
  theme_bw() 
  