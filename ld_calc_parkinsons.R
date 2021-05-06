library(LDlinkR)
library(stringr)
library(rlist)
library(xlsx)
library(tidyverse)

pd_info = read.delim("pd_rs_gene.txt")

snps = as.vector(pd_info[,1])
genes = as.vector(pd_info[,2])

print(snps)
print(genes)

rs_snps = c()
gene_nums = c()

for(gene in genes){
  temp = str_extract(gene, "\\d*")
  gene_nums = append(gene_nums, temp, after = length(gene_nums))
}

for(snp in snps){
  temp = str_extract(snp, "rs[0-9]+")
  rs_snps = append(rs_snps, temp, after = length(rs_snps))
}

rs_snps = rs_snps[!is.na(rs_snps)]

all_ld <- list()
for (filename in list.files("ld_files")){
  all_ld[[length(all_ld)+1]] = read_table(paste0("ld_files/",filename))
}

final_table = do.call(rbind, all_ld)

ldData = filter(final_table, SNP_B %in% rs_snps)
newSNPs = ldData$SNP_A

print(newSNPs)