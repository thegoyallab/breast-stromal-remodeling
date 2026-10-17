suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- '/mnt/stk35data/BreastCancer_Census_v1'
SRC <- file.path(ROOT,'StageM13_PostIndependentReview_v0_1','results')
OUT <- file.path(ROOT,'StageM20B_MethodConcordance_v0_1','results')
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

jobs <- data.frame(
  contrast=c('FB1plusFB3','PV3','VEAT','VEATplusVEC'),
  counts=c('M13D_scanvi_FB1plusFB3_floorall_counts_v0_1.tsv.gz',
           'M13D_scanvi_PV3_floorall_counts_v0_1.tsv.gz',
           'M13D_scanvi_VEAT_floorall_counts_v0_1.tsv.gz',
           'M13D_scanvi_VEATplusVEC_floorall_counts_v0_1.tsv.gz'),
  meta=c('M13D_scanvi_FB1plusFB3_floorall_meta_v0_1.tsv',
         'M13D_scanvi_PV3_floorall_meta_v0_1.tsv',
         'M13D_scanvi_VEAT_floorall_meta_v0_1.tsv',
         'M13D_scanvi_VEATplusVEC_floorall_meta_v0_1.tsv'),
  stringsAsFactors=FALSE
)

sumlist <- list()
for(i in seq_len(nrow(jobs))){
  lab <- jobs$contrast[i]
  m <- read.delim(file.path(SRC,jobs$meta[i]),stringsAsFactors=FALSE,check.names=FALSE)
  m <- m[m$n_cells>=20 & m$role %in% c('reference_normal','core_primary'),,drop=FALSE]
  tab <- read.delim(gzfile(file.path(SRC,jobs$counts[i])),check.names=FALSE)
  sym <- tab$gene_symbol
  counts <- as.matrix(tab[,m$sample_id,drop=FALSE]); storage.mode(counts) <- 'integer'
  rn <- make.unique(ifelse(is.na(sym)|sym=='','UNNAMED',sym)); rownames(counts)<-rn
  grp <- factor(ifelse(m$role=='core_primary','Primary','Normal'),levels=c('Normal','Primary'))
  design <- model.matrix(~grp)

  y <- DGEList(counts=counts,genes=data.frame(gene_symbol=sym,stringsAsFactors=FALSE))
  keep <- filterByExpr(y,group=grp,min.count=10)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method='TMM')

  y2 <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y2,design,robust=TRUE)
  qlf <- glmQLFTest(fit,coef=2)
  et <- topTags(qlf,n=Inf,sort.by='none')$table
  et$gene_symbol <- y2$genes$gene_symbol
  et$edgeR_FDR <- p.adjust(et$PValue,'BH')
  et$gene_key <- rownames(et)

  v <- voom(y,design,plot=FALSE)
  vf <- lmFit(v,design)
  vf <- eBayes(vf,robust=TRUE,trend=FALSE)
  vt <- topTable(vf,coef=2,number=Inf,sort.by='none',adjust.method='BH')
  vt$gene_key <- rownames(vt)

  mm <- merge(
    et[,c('gene_key','gene_symbol','logFC','PValue','edgeR_FDR')],
    vt[,c('gene_key','logFC','P.Value','adj.P.Val')],
    by='gene_key',suffixes=c('_edgeR','_voom'),all=FALSE
  )
  names(mm)[names(mm)=='logFC_edgeR'] <- 'edgeR_log2FC'
  names(mm)[names(mm)=='PValue'] <- 'edgeR_P'
  names(mm)[names(mm)=='logFC_voom'] <- 'voom_log2FC'
  names(mm)[names(mm)=='P.Value'] <- 'voom_P'
  names(mm)[names(mm)=='adj.P.Val'] <- 'voom_FDR'
  mm$direction_concordant <- sign(mm$edgeR_log2FC)==sign(mm$voom_log2FC)
  mm$both_FDR05 <- mm$edgeR_FDR<0.05 & mm$voom_FDR<0.05
  mm$edgeR_only_FDR05 <- mm$edgeR_FDR<0.05 & mm$voom_FDR>=0.05
  mm$voom_only_FDR05 <- mm$voom_FDR<0.05 & mm$edgeR_FDR>=0.05
  write.table(mm,file.path(OUT,paste0('M20B_',lab,'_edgeR_voom_gene_concordance.tsv')),
              sep='\t',row.names=FALSE,quote=FALSE)

  rho <- suppressWarnings(cor(mm$edgeR_log2FC,mm$voom_log2FC,method='spearman',use='complete.obs'))
  pear <- suppressWarnings(cor(mm$edgeR_log2FC,mm$voom_log2FC,method='pearson',use='complete.obs'))
  ed_sig <- mm$edgeR_FDR<0.05
  vo_sig <- mm$voom_FDR<0.05
  union_sig <- ed_sig | vo_sig
  sumlist[[lab]] <- data.frame(
    contrast=lab,n_normal=sum(grp=='Normal'),n_primary=sum(grp=='Primary'),
    genes_tested=nrow(mm),spearman_logFC=rho,pearson_logFC=pear,
    direction_concordance_all=mean(mm$direction_concordant,na.rm=TRUE),
    direction_concordance_unionFDR05=if(any(union_sig)) mean(mm$direction_concordant[union_sig],na.rm=TRUE) else NA,
    edgeR_FDR05=sum(ed_sig),voom_FDR05=sum(vo_sig),both_FDR05=sum(ed_sig & vo_sig),
    edgeR_only_FDR05=sum(ed_sig & !vo_sig),voom_only_FDR05=sum(vo_sig & !ed_sig),
    jaccard_FDR05=if(sum(union_sig)>0) sum(ed_sig & vo_sig)/sum(union_sig) else NA
  )
}
S <- do.call(rbind,sumlist)
write.table(S,file.path(OUT,'M20B_edgeR_voom_concordance_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(S)
