
library(tidyverse)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)

###Add Tiassale data to read counts###
burkina_counts = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/all_counts.txt',header=T)
tiassale = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/NG-A0235_AB.featureCounts-complete.txt',header=T)
out = merge(burkina_counts,tiassale,by='Geneid')

write.table(out,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/all_counts_with_tia.txt',row.names=F,sep='\t')

###Generate a count file for each population that is a range of the original counts###
burkina_counts = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/all_counts_with_tia.txt',header=T)
list_names = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/headers_RNAseq.txt',header=F)

pops = c('Tiassale','Bakaridjan','Gaoua','VK7','Tiefora','Banfora_Res','Banfora_Or','Banfora_Sus')


final_out = c()
for(i in 1:length(pops))
{
  sub = burkina_counts[,grep(pops[i],colnames(burkina_counts))]
  n = (length(grep(pops[i],list_names[,1]))-ncol(sub))
  out_table = c()
  for(j in 1:nrow(sub))
  {
    out_counts = floor(runif(n,min(sub[j,]),max(sub[j,])))
    out_table = rbind(out_table,out_counts)
  }
  out_table = cbind(sub,out_table)
  #colnames(out_table) = list_names[grep(pops[i],list_names[,1]),1]
  out_table = unname(as.matrix(out_table))
  final_out = cbind(final_out,out_table)
  
}

colnames(final_out) = list_names[,1]
final_out = cbind(burkina_counts[,1],final_out)

##Reorder to match vcf##
col.order = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/headers_no_Ten.txt',header=F)

col.order= unlist(col.order)
final.out = as.data.frame(final_out)
final_out2 = final.out[col.order]

write.table(final_out2,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia.txt',sep='\t',row.names=F)

###Generate the SNP file###

positions = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_pos.txt',header=F)

positions2 = paste(positions$V1,positions$V2,sep=':')

rm(positions)

indiv = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_indv.txt',header=F)

##Drop Tengrela columns##
geno_file = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file_012.txt',header=F)


indiv = as.matrix(indiv)


geno_file = t(geno_file)

geno_file = geno_file[,-c(25:33)]

geno_file[1,] = indiv

positions2 = c('id',positions2)
final_file = cbind(positions2,geno_file)

write.table(final_file,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/final_file.txt',sep='\t',row.names=F,quote=F,col.names=F)

##run perl -ni -e 'print unless $. == 1' in bash to remove top line##

gene_location = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/GenomicRanges2.txt',header=T)

gene_location = as.matrix(gene_location)

library(MatrixEQTL)
base.dir = find.package('MatrixEQTL')
useModel = modelLINEAR

output_file_name_cis = tempfile()
output_file_name_tra = tempfile()



errorCovariance = numeric()

###Need to remove batch effects and perform unbiased normalisation on read counts##
counts_file = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia.txt',header=T)
sampleData = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/colData.txt',header=T)

library(DESeq2)
library(limma)

counts = counts_file[,2:ncol(counts_file)]
row.names(counts) = make.names(counts_file[,1],unique=T)

sampleData$seqLibrary = as.factor(sampleData$seqLibrary)

dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = sampleData, 
                             design = ~seqLibrary)

dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
vst = varianceStabilizingTransformation(dds)
mat <- assay(vst)
mat <- limma::removeBatchEffect(mat, vst$Batch)

write.table(mat,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia_normalised.txt',row.names=T,sep='\t')





snps2 = SlicedData$new();
snps2$fileDelimiter = "\t";      # the TAB character
snps2$fileOmitCharacters = "-1"; # denote missing values;
snps2$fileSkipRows = 1;          # one row of column labels
snps2$fileSkipColumns = 1;       # one column of row labels
snps2$fileSliceSize = 5000;      # read file in slices of 2,000 rows
snps2$LoadFile('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/final_file.txt')

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia_normalised.txt');


cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/covariate2.txt')>0) {
  cvrt$LoadFile('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/covariate2.txt');
}

#snp_names2 = as.data.frame(snp_names2)
gene_location = as.data.frame(gene_location)

#snp_names2$pos = as.numeric(snp_names2$pos)
gene_location$s1 = as.numeric(gene_location$s1)
gene_location$s2 = as.numeric(gene_location$s2)





snps = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_pos.txt',header=F)
snps$ID <- paste(snps$V1,snps$V2,sep=":")

snps = snps[,c(3,1,2)]

colnames(snps) = c('snp','chr','pos')
mode(snps$pos) = 'numeric'

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps))

pvOutputThreshold_cis = 0.05/nrow(snps);
pvOutputThreshold_tra = 0.05/nrow(snps);
cisDist = 1e5

###pvOutputThreshold     = 0 means that only cis runs###
me = Matrix_eQTL_main(
  snps = snps2,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = 0,
  useModel = useModel,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snps,
  genepos = gene_location,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
  )


plot(me)

meq = Matrix_eQTL_engine(
  snps = snps2, 
  gene = gene, 
  cvrt = cvrt, 
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_cis, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot");
unlink( filename );
# png(filename = "QQplot.png", width = 650, height = 650)
plot(meq, pch = 16, cex = 0.7)


#trans_results = me$trans$eqtls
cis_results = me$cis$eqtls


#trans_eqtls = unique(trans_results$snps)



top_eqtls = cis_results[which(cis_results$FDR < 0.05/nrow(snps)),]
top_eqtls = top_eqtls[order(top_eqtls$FDR),]

write.table(top_eqtls,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/cis_eQTLwith_Tia_1e5.txt',sep='\t',row.names=F)
#write.table(trans_results,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/trans_eQTL.txt',sep='\t',row.names=F)

cis_eqtls = unique(top_eqtls$snps)
write.table(cis_eqtls,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/cis_eQTL_1e5_unique_with_tia.txt',sep='\t',row.names=F)
#write.table(trans_eqtls,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/trans_eQTL_uniquee.txt',sep='\t',row.names=F)








##Let's explore some of the most important SNPs. As we have large haplotype blocks most likely, we can extract just a representative SNP

gene_values = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts2.txt',row.names=1,header=T)
gene_values = data.frame(gene = rownames(gene_values), gene_values, stringsAsFactors = FALSE)

###Extract just SNPs that are significant in the eQTL analysis
snps$row = 1:nrow(snps)
soi = snps[snps$snp %in% cis_eqtls,]
snp_row = soi$row
snp_values = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.txt',header=F)


snp_values2 = snp_values[snp_row,]

top_snp = top_eqtls$snps[1]

top_gene = as.character(top_eqtls$gene[1])

top_snp_data = filter(snp_values, snps == top_snp)
top_gene_data = filter(gene_values, gene == top_gene)


plot_data = t(bind_rows(top_snp_data[-1], top_gene_data[-1]))
colnames(plot_data) = c("snp", "gene_expr")
plot_data = as.data.frame(plot_data)
plot_data$snp = as.factor(plot_data$snp)
head(plot_data)


lm_top = lm(plot_data[,"gene_expr"] ~ as.numeric(plot_data[,"snp"]))
summary(lm_top)


plot(plot_data, col="steelblue", 
     main = paste0(top_gene, " vs ", top_snp))
abline(lm_top, col="darkorange", lwd = 2, lty = 2)
y_range = range(plot_data[,"gene_expr"])
text(x=2, y=y_range[1] + 0.5*diff(y_range), paste0("p=",
                                                   format(summary(lm_top)$coefficients[2,4],
                                                          scentific=TRUE, digits=2)), col = "darkorange")



