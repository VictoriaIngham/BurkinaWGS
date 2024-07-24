###Permutation test for p-value cut-off###

##Scramble the samples to break down the relationship between phenotype and genotype###
genotype_file = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia_normalised.txt',header=T)

names = genotype_file[,1]
genotype_file = genotype_file[,-1]

for(i in 1:100)
{
  df2 <-cbind(names, genotype_file[sample(nrow(genotype_file)),])
  
  write.table(df2,paste('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/Permutation_counts/',i,'.txt',sep=""),row.names=F,sep='\t')
}








###Now to run the eQTL through###

positions = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_pos.txt',header=F)

positions2 = paste(positions$V1,positions$V2,sep=':')

rm(positions)

indiv = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_indv.txt',header=F)

gene_location = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/GenomicRanges2.txt',header=T)

gene_location = as.matrix(gene_location)
#snp_names2 = as.data.frame(snp_names2)
gene_location = as.data.frame(gene_location)

#snp_names2$pos = as.numeric(snp_names2$pos)
gene_location$s1 = as.numeric(gene_location$s1)
gene_location$s2 = as.numeric(gene_location$s2)

library(MatrixEQTL)
base.dir = find.package('MatrixEQTL')
useModel = modelLINEAR

output_file_name_cis = tempfile()
output_file_name_tra = tempfile()


errorCovariance = numeric()

###Need to remove batch effects and perform unbiased normalisation on read counts##
counts_file = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/generated_counts_with_Tia.txt',header=T)

##File path and population of interest windows###
files = list.files('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/Permutation_counts/',pattern='*.txt')


snps = fread('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/geno_file.012_pos.txt',header=F)

snps$ID <- paste(snps$V1,snps$V2,sep=":")

snps = snps[,c(3,1,2)]

colnames(snps) = c('snp','chr','pos')
mode(snps$pos) = 'numeric'

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.05);
#snps$RowReorder(maf>0.05);
#cat('SNPs before filtering:',nrow(snps))

pvOutputThreshold_cis = 0.05/nrow(snps);
pvOutputThreshold_tra = 0.05/nrow(snps);
cisDist = 1e5


snps2 = SlicedData$new();
snps2$fileDelimiter = "\t";      # the TAB character
snps2$fileOmitCharacters = "-1"; # denote missing values;
snps2$fileSkipRows = 1;          # one row of column labels
snps2$fileSkipColumns = 1;       # one column of row labels
snps2$fileSliceSize = 5000;      # read file in slices of 2,000 rows
snps2$LoadFile('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/final_file.txt')


cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/covariate2.txt')>0) {
  cvrt$LoadFile('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/covariate2.txt');
}

real = read.delim('/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/cis_eqtls_final.txt',header=F)

for(i in 2:length(files))
{
  filepath <- file.path("/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/Permutation_counts",files[i])
  

  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(filepath);
  
  
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
  
  #trans_results = me$trans$eqtls
  cis_results = me$cis$eqtls
  
  
  #trans_eqtls = unique(trans_results$snps)
  
  
  
  top_eqtls = cis_results[which(cis_results$FDR < 0.05/nrow(snps)),]
  top_eqtls = top_eqtls[,c(1,2,5)]
  #top_eqtls = top_eqtls[order(top_eqtls$FDR),]
  
  if(i==1)
  {
    df_out = merge(real,top_eqtls, by.x=c("V1", "V2"), by.y=c("snps", "gene"),all.x=T)
    colnames(df_out)[ncol(df_out)] = paste('perm',i,sep='')
  }
  else
  {
    df_out = merge(df_out,top_eqtls, by.x=c("V1", "V2"), by.y=c("snps", "gene"),all.x=T)
    colnames(df_out)[ncol(df_out)] = paste('perm',i,sep='')
  }
  
}



df_out[is.na(df_out)] = 1
df4 = df_out[,c(5,7:ncol(df_out))]


df5 = t(apply(df4, 1, rank))

p_vals = df5[,1]/100

write.table(p_vals,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/Permutation_counts/p_vals.txt',sep='\t',row.names=F)

library(qvalue)

q_obj = qvalue(p_vals,lambda=0)

real_out = cbind(real,q_obj$qvalues,q_obj$pvalues,q_obj$lfdr)
colnames(real_out) = c('snp','gene','statistic','pvalue','FDR','beta','perm_qvalue','perm_pvalue','perm_localFDR')

write.table(real_out,'/Users/vickyingham/Seafile/WGS/Analysis of vcf/eQTL/Permutation Test/Permutation_counts/permuted_final.txt',sep='\t',row.names=F)

