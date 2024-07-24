library(RIdeogram)
library(chromoMap)

karyotype = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/Fst ABBA Overlap/agamp4_karyotype.txt',header=T)
#gene_density <- GFFex(input = "/Users/vickyingham/Seafile/WGS/VectorBase-55_AgambiaePEST.gff", karyotype = "/Users/vickyingham/Seafile/WGS/Analysis of vcf/Fst ABBA Overlap/agamp4_karyotype.txt", feature = "CDS", window = 500000)
features = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/Fst ABBA Overlap/Gaoua_features.txt',header=T)

#Gaoua
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Gaoua',anno_col = c("#40B0A6"),interactivity = F,
          left_margin = 90)

#Bakaridjan
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Bakaridjan',anno_col = c("#0C7BDC"),interactivity = F,
          left_margin = 90)      

#Tengrela
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Tengrela',anno_col = c("#1AFF1A"),interactivity = F,
          left_margin = 90)    


#Banfora
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Banfora',anno_col = c("#D41159"),interactivity = F,
          left_margin = 90) 

#Tiefora
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Tiefora',anno_col = c("#5D3A9B"),interactivity = F,
          left_margin = 90) 


#Tiassale
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='Tiassale',anno_col = c("#E66100"),interactivity = F,
          left_margin = 90) 

#VK7
chromoMap(list(karyotype),list(features),n_win.factor = 2,win.summary.display = T,title='VK7',anno_col = c("grey"),interactivity = F,
          left_margin = 90) 
