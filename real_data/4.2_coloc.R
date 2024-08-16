require(data.table)
require(openxlsx)
require(coloc)
require(snpStats)
require(locuscomparer)
require(foreach)
require(doParallel)

sites = fread("./data/queue/site.txt", header=F)$V1
snps = read.xlsx("./data/Novel_BMD_loci_qced.xlsx")$SNP


plot_and_coloc <- function(snp, site, gene) {
    fh1 = paste0("./result/bmd/2MB/", snp, "_", site, ".txt.gz")
    fh2 = paste0("./Sumstats/eQTL/data/cleaned/osteoclast/", gene, ".txt")

    try({
        pdf(paste0("./result/bmd/coloc_14/", snp, "_", site, "_", gene, ".pdf"), width = 8, height = 4)
        print(locuscompare(in_fn1=fh1, in_fn2=fh2, title1="GWAS", title2="eQTL", marker_col1="SNP", pval_col1="P", marker_col2="SNP", pval_col2="P", snp=snp))
        dev.off()
    }, silent = T)

    gwas = fread(fh1)
    eqtl = fread(fh2)
    if (sum(snp %in% eqtl$SNP)>1) {
        chr = ref$CHR[ref$SNP==snp]
        bp = ref$BP[ref$SNP==snp]
        eqtl = eqtl[-which(eqtl$SNP==snp & (eqtl$BP!=bp | eqtl$CHR!=chr)),]
    }
    input <- merge(eqtl, gwas, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
    input = input[!duplicated(input$SNP), ]

    if (nrow(input) > 0) {
        result <- coloc.abf(dataset1=list(snp=input$SNP, pvalues=input$P_gwas, type="cc", s=0.33, N=nrow(gwas)), dataset2=list(snp=input$SNP, pvalues=input$P_eqtl, type="quant", N=nrow(eqtl)), MAF=ifelse(input$EAF<0.5, input$EAF, 1-input$EAF))$summary
        data.frame(SNP=snp, SITE=site, GENE=gene, t(result))
    } else{
        c()
    }
}

registerDoParallel(cores = 32)
df <- rbindlist(foreach (snp = snps) %dopar% {
    genes = fread(paste0("./results/snp_GENE_list/", snp, ".txt"), header = T)$gene
    data.frame(SNP = rep(snp, length(genes)), GENE = genes)
})
nrow(df)
df = data.frame(SITE = rep(sites, each=nrow(df)), do.call(rbind, replicate(length(sites),df,simplify=F)))
nrow(df)
out <- rbindlist(foreach (d = split(df, 1:nrow(df))) %dopar% {
    plot_and_coloc(d$SNP, d$SITE, d$GENE)
})
nrow(out)
stopImplicitCluster()

write.xlsx(out, "./result/bmd/novel_loci_qced_coloc_14.xlsx", colNames=T, rowNames=F, quote=F, overwrite=T)