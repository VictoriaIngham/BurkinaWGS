##File path and population of interest windows###
files = list.files('/Users/vickyingham/Seafile/WGS/Metrics/TajimasD',pattern='*.txt')

###Create dataframes for each###
names = gsub("(?:[^.]+\\.){0}([^.]+).*", "\\1", files)

all = c()
for(i in 1:length(files))
{
  filepath <- file.path("/Users/vickyingham/Seafile/WGS/Metrics/TajimasD",files[i])
  foi = read.delim(filepath,
                   colClasses=c("character","numeric","numeric","numeric"),
                   sep = "\t",head=T)
  foi$pop = rep(names[i],nrow(foi))
  all = rbind(all,foi)
}
rm(foi)

all$pop[grep('bak',all$pop)] = 'Bakaridjan'
all$pop[grep('ban',all$pop)] = 'Banfora'
all$pop[grep('gou',all$pop)] = 'Gaoua'
all$pop[grep('ten',all$pop)] = 'Tengrela'
all$pop[grep('tie',all$pop)] = 'Tiefora'
all$pop[grep('tia',all$pop)] = 'Tiassale'
all$pop[grep('VK',all$pop)] = 'VK7'

group.colours = c(Banfora='#D41159',Bakaridjan='#0C7BDC',Tiassale='#E66100',Tiefora = '#5D3A9B',Gaoua = '#40B0A6',VK7 = '#000000',Tengrela='#1AFF1A')

all$position = (all$BIN_START)

all = droplevels(all[!all$CHROM == 'AgamP4_Mt',])
all = droplevels(all[!all$CHROM == 'AgamP4_UNKN',])
              
                 

library(ggplot2)
ggplot(all,aes(x=factor(CHROM),y=as.numeric(TajimaD),fill=pop)) + 
  geom_boxplot() +
  scale_fill_manual(values=group.colours)+ 
  labs(x="Chromosome" , y='Tajimas D', col="Population") +
  theme_linedraw() + guides(fill=guide_legend(title="Population"))

