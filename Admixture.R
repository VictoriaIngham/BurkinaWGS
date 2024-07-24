

tbl=read.table("/Volumes/Elements/ADMIXTURE/admixture_out/admixture_file.6.Q")

tbl2 =rbind(tbl[1:4,],tbl[34:48,],tbl[5:9,],tbl[49:62,],tbl[10:14,],tbl[63:77,],tbl[15:19,],tbl[78:91,],tbl[20:24,],tbl[93:107,],tbl[108:117,],tbl[25:33,])

group.colours = c(Banfora='#D41159',Bakaridjan='#0C7BDC',Tiassale='#E66100',Tiefora = '#5D3A9B',Gaoua = '#40B0A6',VK7 = '#000000',TengrelaA='#FFC20A',TengrelaD='#FFC20A')

col = c('#E66100','#0C7BDC','#40B0A6','#000000','#5D3A9B','#FFC20A','#1AFF1A','#D41159')
barplot(t(as.matrix(tbl2)), col=col,
          xlab="Individual #", ylab="Ancestry", border=NA)