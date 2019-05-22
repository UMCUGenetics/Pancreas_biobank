#####
#RNAseq
#####

outdir = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/RNAseq/"
metadata <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/fingerprint_pancreas_correct_match.xlsx", col_names = T,sheet = "Overview"))
names(metadata)
#Extracting some genes and positions from Biomart
#Using biomart
library(biomaRt)
library(edgeR)
library(gplots) 
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(cowplot)
library(org.Hs.eg.db)
library(RColorBrewer)
library(DESeq2)
library(DESeq)
library(ggpubr)
#BiocManager::install("DESeq", version = "3.8")

listMarts()
ensembl=useMart("ensembl")
#To use hg19/GrCh37
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
# grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
# listMarts(grch37)
ensembl<- useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
listFilters(ensembl)
listDatasets(ensembl) 
# #genes.with.id <- getBM(attributes=c("ensembl_gene_id", "external_gene_id"),values=gene_names, mart= ensembl)
# human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# cc <- getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position"),
#            filters = "hgnc_symbol", values = *,mart = human)
all.genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol','chromosome_name','start_position','end_position'),   mart = ensembl)
colnames(all.genes)
head(all.genes)


raw_counts_file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/RNAseq/Pancreas/raw_data/SU2C_raw_counts.txt"
RPKM_file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/RNAseq/Pancreas/raw_data/SU2C_RPKM.txt"
normalized_file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/RNAseq/Pancreas/raw_data/SU2C_normalized.txt"

countsTable = read.table(raw_counts_file, header = T)
RPKM = read.table(RPKM_file, header = T)
normalized = read.table(normalized_file, header = T)

data.frame(RNAseq_ID=colnames(countsTable)) %>% left_join(., dplyr::select(metadata,RNAseq_ID,SampleID_short), by = "RNAseq_ID") -> new_names
colnames(countsTable)
new_names$RNAseq_ID
colnames(countsTable) <- new_names$SampleID_short

data.frame(RNAseq_ID=colnames(normalized)) %>% left_join(., dplyr::select(metadata,RNAseq_ID,SampleID_short), by = "RNAseq_ID") -> new_names
colnames(normalized)
new_names$RNAseq_ID
colnames(normalized) <- new_names$SampleID_short

data.frame(RNAseq_ID=colnames(RPKM)) %>% left_join(., dplyr::select(metadata,RNAseq_ID,SampleID_short), by = "RNAseq_ID") -> new_names
colnames(RPKM)
new_names$RNAseq_ID
colnames(RPKM) <- new_names$SampleID_short

# ------- SAMPLE CORRELATION ---------
# discard genes with < 10 counts over all samples
expressed = which(rowSums(normalized) > 10)

# spearman correlation between samples based on expressed genes
cor_samples = cor(normalized[expressed,], method = 'spearman')
plot_cor = heatmap.2(cor_samples, scale = 'none', col = rev(heat.colors(256)))
pdf(paste(outdir, "correlation_samples.pdf", sep=""), 11, 11)
heatmap(cor_samples, scale = 'none', col = rev(heat.colors(256)),  margins=c(15,15))
dev.off()

head(df)
df <- normalized
df <- log2(df+1)
df$gene <- row.names(df)
row.names(df) <- NULL
df = melt(df, id.vars = "gene")
colnames(df) = c("gene", 'sample', "log2normcounts")
plot_counts = ggplot(df, aes(x = log2normcounts)) +
  geom_histogram(binwidth=0.2, colour="black", fill="lightblue") +
  scale_y_sqrt()

pdf(paste(outdir, "hist_nonzero_norm_counts.pdf", sep=""), 9, 6)
plot_counts
dev.off()



# plot RPKM histogram
df = melt(log2(RPKM))
colnames(df) = c("sample","log2RPKM")
plot_RPKM = ggplot(df, aes(x = log2RPKM)) +
  geom_histogram(binwidth=0.2, colour="black", fill="lightblue")

pdf(paste(outdir, "hist_nonzero_RPKM2.pdf", sep=""), 9, 6)
plot_RPKM
dev.off()


candidates <- c("AHCY","AMD1","AS3MT","ASMT","BHMT","CARM1","CBS","CHDH","COMT","CTH","DNMT1","DNMT3A","DNMT3B","DOT1L","EZH2","GAMT","GNMT","HNMT","ICMT","IL4I1","SETD8","MARS","MAT1A","MAT2A","MAT2B","FTSJ2","MSRB2","MSRB3","MTAP","MTFMT","MTHFR","MTR","MTRR","NNMT","PCMT1","PEMT","PNMT","PRMT1","PRMT2","PRMT3","PRMT5","PRMT6","PRMT7","RNMT","SETD7","SETDB1","SHMT1","SMS","SMYD2","SRM","SUV39H1","SUV39H2","TRDMT1")
candidates <- all.genes[which(all.genes$hgnc_symbol %in% candidates),]

df_heatmap = as.data.frame(t(normalized[which(row.names(normalized) %in% candidates$ensembl_gene_id),]))
df_heatmap  %>% tibble::rownames_to_column("SampleID_short") %>% left_join(., dplyr::select(metadata,SampleID_short,mtap_DNA_status,PMRT5_sens,PMRT5_sens_2), by = "SampleID_short") %>% as.data.frame()-> df_heatmap
row.names(df_heatmap) <- df_heatmap$SampleID_short
df_heatmap$SampleID_short <- NULL
#without MTAP
#df_heatmap$ENSG00000099810 <- NULL
df_heatmap<-df_heatmap[which(df_heatmap$mtap_DNA_status==1),]

ha01 = rowAnnotation(df = df_heatmap["mtap_DNA_status"],
                     width = unit(1, "cm"),
                     na_col = "white",
                     col = list(mtap_DNA_status=c("1"="dodgerblue3","0"="orangered")))

ha02 = rowAnnotation(df = df_heatmap["PMRT5_sens_2"],
                     width = unit(1, "cm"),
                     na_col = "white",
                     #col = list(PMRT5_sens=c("1"="brown","0"="goldenrod3","NA"="white")))
                     col = list(PMRT5_sens_2=c("YES"="olivedrab4","NO"="orangered3","NA"="grey")))

group <- df_heatmap[,grep("ENSG", names(df_heatmap),value = T)]
mapke = Heatmap(log2(group+1),
                cluster_columns = as.dendrogram(agnes(t(log2(group+1)))),
                clustering_distance_rows = "euclidean",#cluster_rows = hc.sample
                row_dend_width = unit(40, "mm"),
                row_names_side = "left",
                show_column_dend = FALSE,
                row_dend_side = "right",
                col = colorRamp2(c(0,5, 12), c("#FFFFC0","#54B3B0", "#1B2182")),
                column_names_gp =  gpar(fontsize = 12),
                row_names_gp = gpar(fontsize = 12))

plot_heatmap = paste(outdir,"heatmap_withouth_MTAP_without_MTAPdeleted.pdf", sep = "/")
#plot_heatmap = paste(outdir,"heatmap_without_MTAP.pdf", sep = "/")
pdf(plot_heatmap , useDingbats = F, width = 15, height = 14) 
mapke + ha01 + ha02
dev.off()


# ------ DIFFERENTIAL EXPRESSION ------ 

#select whitelist samples
whitelist_pos <- c("T3","T7","T21","T23")
whitelist_neg <- c("T4","T6","T20","T24")
whitelist <- c(whitelist_pos,whitelist_neg)

countsTable %>% dplyr::select(whitelist) -> countsTable_DE
colnames(countsTable_DE)

conds = factor(c(rep("Resistant",length(whitelist_pos)), rep("Sensitive",length(whitelist_neg))), levels = c("Resistant", "Sensitive"))

# ------ DIFFERENTIAL EXPRESSION ------ 

# Test for differential expression
# H0 : the measurements come from the same distribution and the gene is being expressed 
# at the same level across conditions

# general function for differential expression calculation
# annotate with gene symbols & description
diff_expr = function(counts=countsTable_DE, conds=conds, map=org.Hs.eg.db)
{
  samples = colnames(counts)
  # counts normalization
  cds = newCountDataSet(countsTable_DE, conds)
  cds = estimateSizeFactors(cds)
  normalized_counts = counts(cds, normalized=TRUE)
  # Estimate variance
  cds = estimateDispersions(cds)
  DE = nbinomTest(cds, levels(conds)[1], levels(conds)[2])
  # add normalized counts
  DE = merge(normalized_counts, DE, by.x= "row.names", by.y = "id")
  # recalculate log2foldchange with pseudocount to avoid Inf
  DE$log2FoldChange_pseudo= log2((DE$baseMeanB + 1) / (DE$baseMeanA+1))
  # find gene symbol and description
  cols = c("SYMBOL","GENENAME")
  DE_genes = select(map, keys=DE$Row.names, columns=cols, keytype="ENSEMBL")
  DE = merge(DE, DE_genes, by.x = "Row.names", by.y = "ENSEMBL")
  # remove duplicate genes (which arises during gene conversion)
  DE = DE[which(!duplicated(DE$Row.names)),]
  # sort result on adjusted pvalue
  DE = DE[order(DE$pval),]
  return(DE)
}

# DE MUT & WT
DE_overall = diff_expr(countsTable, conds, org.Hs.eg.db)
head(DE_overall[,-1])
nrow(DE_overall)
# Find interesting DE genes
threshold = 2
DE_overall_int = subset(DE_overall, DE_overall$pval < 0.05 & abs(DE_overall$log2FoldChange_pseudo) > threshold)
hist(DE_overall_int$pval)
heatmap.2(as.matrix(DE_overall_int[whitelist]), scale = 'none',col = rev(heat.colors(256)),   margins=c(5,5))
DE_overall_int$SYMBOL
write.table(DE_overall_int, file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/RNAseq/DE_overall_on_8_samples.txt", sep = "\t", quote = F, row.names = F)



# ------ make plots ------ 

load("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/MTAP_CDKN2A.RData")

MTAP_CDKN2A %>% dplyr::mutate_all(funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(minorallele = ifelse(minoralleleploidy_mean >= 0.5 & gene_perc_deleted < 2 , 100,
                                     ifelse(minoralleleploidy_mean >= 0.5 & gene_perc_deleted > 2 , gene_perc_deleted, 0)),
                majorallele = ifelse(majoralleleploidy_mean >= 0.5 & gene_perc_deleted < 2 , 100,
                                     ifelse(majoralleleploidy_mean >= 0.5 & gene_perc_deleted > 2 , (100-gene_perc_deleted), 0)))  %>% 
  dplyr::mutate(minorallele = ifelse(minoralleleploidy_min=="Inf" & copyNumber_CN_max>0.9 & 1/observedbaf>1.2 ,100,minorallele),
                majorallele = ifelse(minoralleleploidy_min=="Inf" & copyNumber_CN_max>0.9 & 1/observedbaf>1.2 ,100,majorallele)) %>% 
  dplyr::mutate('Minor allele' = 100-minorallele,
                'Major allele' = 100-majorallele) %>%
  left_join(., dplyr::select(metadata,sampleId,SampleID_short), by = "sampleId") -> MTAP_CDKN2A_overview


MTAP_CDKN2A_overview %>%
  dplyr::select(SampleID_short,sampleID,geneID,'Minor.allele','Major.allele') %>%
  melt() -> MTAP_CDKN2A_plot_table
MTAP_CDKN2A_plot_table$value2 = 100-MTAP_CDKN2A_plot_table$value
MTAP_CDKN2A_plot_table$geneID = factor(MTAP_CDKN2A_plot_table$geneID, levels = c("CDKN2A","MTAP"))
MTAP_CDKN2A_plot_table$variable = factor(MTAP_CDKN2A_plot_table$variable, levels = c('Major.allele','Minor.allele'))
MTAP_CDKN2A_plot_table$percent_genes = MTAP_CDKN2A_plot_table$value2/100
MTAP_CDKN2A_plot_table<- MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$SampleID_short!="T27"),]
MTAP_CDKN2A_plot_table$SampleID_short = factor(MTAP_CDKN2A_plot_table$SampleID_short, levels = rev(sample_order))


###add expression levels
candidates <- c("MTAP","CDKN2A")
candidates <- all.genes[which(all.genes$hgnc_symbol %in% candidates),]
candidates %>% dplyr::arrange(match(hgnc_symbol,c("MTAP","CDKN2A"))) -> candidates
candidates$ensembl_gene_id <- factor(candidates$ensembl_gene_id, levels = c(candidates$ensembl_gene_id))
candidates$hgnc_symbol <- factor(candidates$hgnc_symbol, levels = c(candidates$hgnc_symbol))
df_heatmap = as.data.frame(t(normalized[which(row.names(normalized) %in% candidates$ensembl_gene_id),]))
df_heatmap <- df_heatmap[as.vector(candidates$ensembl_gene_id)]
colnames(df_heatmap) <- as.vector(candidates$hgnc_symbol)
df_heatmap  %>% tibble::rownames_to_column("SampleID_short") ->  df_heatmap
##remove tumor only samples
df_heatmap  <- df_heatmap[which(df_heatmap$SampleID_short!="T29"),]
df_heatmap  <- df_heatmap[which(df_heatmap$SampleID_short!="T30"),]
df_heatmap  <- df_heatmap[which(df_heatmap$SampleID_short!="T31"),]
df_heatmap  <- df_heatmap[which(df_heatmap$SampleID_short!="T1"),]
df_heatmap  <- df_heatmap[which(df_heatmap$SampleID_short!="T2"),]
##remove no cancer sample
df_heatmap<- df_heatmap[which(df_heatmap$SampleID_short!="T27"),]

sample_order <- as.character(rev(df_heatmap[order(df_heatmap[,"MTAP"]),]$SampleID_short))
df_heatmap <- df_heatmap[which(df_heatmap$SampleID_short %in% MTAP_CDKN2A_plot_table$SampleID_short),]
#df_heatmap$SampleID_short = factor(df_heatmap$SampleID_short, levels = rev(levels(MTAP_CDKN2A_plot_table$SampleID_short)))
df_heatmap %>% 
  dplyr::mutate(MTAP = ifelse(MTAP <=10, 1, MTAP),
                CDKN2A = ifelse(CDKN2A <=10, 1, CDKN2A)) -> df_heatmap
#dplyr::filter(sampleId=="FR11123204" |sampleId=="FR11123647" ) %>%
df_heatmap$SampleID_short = factor(df_heatmap$SampleID_short, levels = sample_order)
df_heatmap %>%   melt() -> df_heatmap_plot_table
df_heatmap_plot_table$variable = factor(df_heatmap_plot_table$variable, levels = c("CDKN2A","MTAP"))


P1 <- ggplot(MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="CDKN2A"),],
             aes(x=geneID,
                 y=percent_genes,
                 fill=variable,
                 group=variable))+
  geom_bar(colour="black",
           stat="identity",
           position=position_dodge2(reverse = TRUE)) +
  coord_flip()+
  #facet_grid( SampleID_short ~. )+
  facet_wrap(~SampleID_short,  ncol=1, strip.position = "left")+
  theme_bw()  +
  #scale_fill_manual(values= c("#1b9e77","#7570b3"))+
  scale_fill_manual(values = c(alpha("#5E4FA2", 1), alpha("#5E4FA2", 0.5)))+#, alpha("#33A02C", 1), alpha("#33A02C", 0.5)))+
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1))+
  theme(panel.spacing = unit(0, "lines"), 
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 180),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom", 
        legend.title=element_blank(),
        panel.border = element_rect(colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ylab("CDKN2A")

P2 <- ggplot(MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="MTAP"),],
             aes(x=geneID,
                 y=percent_genes,
                 fill=variable,
                 group=variable))+
  geom_bar(colour="black",
           stat="identity",
           position=position_dodge2(reverse = TRUE)) +
  coord_flip()+
  #scale_y_reverse()+
  #facet_grid( SampleID_short ~. )+
  facet_wrap(~SampleID_short,  ncol=1, strip.position = "left")+
  theme_bw()  +
  #scale_fill_manual(values= c("#1b9e77","#7570b3"))+
  scale_fill_manual(values = c(alpha("#3288BD", 1), alpha("#3288BD", 0.5)))+#, alpha("#33A02C", 1), alpha("#33A02C", 0.5)))+
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1),
                     trans = "reverse")+
  theme(panel.spacing = unit(0, "lines"), 
        axis.ticks = element_blank(),
        strip.text.y = element_blank(),  ##remove labels 
        axis.text.x = element_text(size = 10),
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom", 
        legend.title=element_blank(),
        panel.border = element_rect(colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ylab("MTAP")




P3 <- ggplot(df_heatmap_plot_table, aes(x = SampleID_short, y = variable, fill = log2(value+1))) +
  geom_tile(color = "black") +
  coord_flip()+
  scale_fill_gradientn(name = "", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")+
  theme(panel.spacing = unit(0, "lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10,colour = "gray26"),
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom", 
        legend.title=element_text(angle = 0,size = 5),
        panel.border = element_rect(colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ylab("")


metadata <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/fingerprint_pancreas_correct_match.xlsx", col_names = T,sheet = "Overview"))
metadata %>% dplyr::select(SampleID_short,IC50,IC50_2,AUC) -> sensitivity_MTAP
sensitivity_MTAP <- sensitivity_MTAP[which(sensitivity_MTAP$SampleID_short %in% MTAP_CDKN2A_plot_table$SampleID_short),]
sensitivity_MTAP<- sensitivity_MTAP[which(sensitivity_MTAP$SampleID_short!="T27"),]
sensitivity_MTAP$SampleID_short = factor(sensitivity_MTAP$SampleID_short, levels = sample_order)
sensitivity_MTAP$IC50 <- as.numeric(sensitivity_MTAP$IC50)
sensitivity_MTAP$IC50_2 <- as.numeric(sensitivity_MTAP$IC50_2)
sensitivity_MTAP$AUC <- as.numeric(sensitivity_MTAP$AUC)
sensitivity_MTAP %>% 
  melt() -> sensitivity_MTAP_plot_table

P4 <- ggplot(sensitivity_MTAP_plot_table[which(sensitivity_MTAP_plot_table$variable=="IC50_2"),], 
             aes(x = SampleID_short, y = variable, fill = value)) +
  geom_tile(color = "black") +
  coord_flip()+
  scale_fill_gradientn(name = "", colours = brewer.pal(4,"YlGnBu"), na.value = "white")+
  theme(panel.spacing = unit(0, "lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10,colour = "gray26"),
        axis.text.y=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="bottom", 
        legend.title=element_text(angle = 0,size = 5),
        panel.border = element_rect(colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  xlab("") +
  ylab("")


pFigure2 = plot_grid(P1, P2,P3,P4 ,ncol=4, align="hv",axis="b", rel_heights = c(10,10,10,10),rel_widths = c(5,5,8,1.5),scale=c(1,1,1,1))
save_plot("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/Figure X - MTAP_CKN2A_bc_2.png", pFigure2, base_width = 10, base_height = 5)

dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/"
plot_name="CDKN2A_MTAP"
# draw it
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 10) 
plot(pFigure2)
dev.off()





###add main figure MTAP####
MTAP_CN<-MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="MTAP"&MTAP_CDKN2A_plot_table$variable=="Major.allele"),]
CDKN2A_CN<-MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="CDKN2A"&MTAP_CDKN2A_plot_table$variable=="Major.allele"),]
MTAP_CN_m<-MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="MTAP"&MTAP_CDKN2A_plot_table$variable=="Minor.allele"),]
CDKN2A_CN_m<-MTAP_CDKN2A_plot_table[which(MTAP_CDKN2A_plot_table$geneID=="CDKN2A"&MTAP_CDKN2A_plot_table$variable=="Minor.allele"),]

plot_1_major<-left_join(MTAP_CN,CDKN2A_CN, by = "SampleID_short")
plot_1_minor<-left_join(MTAP_CN_m,CDKN2A_CN_m, by = "SampleID_short")
plot_1<-left_join(plot_1_major,plot_1_minor, by = "SampleID_short")

plot_1 %>% dplyr::mutate(MTAP_major=ifelse(percent_genes.x.x>0.8,"present","absent"),
                         MTAP_minor=ifelse(percent_genes.x.y>0.8,"present","absent"),
                         CDKN2A_major=ifelse(percent_genes.y.x>0.8,"present","absent"),
                         CDKN2A_minor=ifelse(percent_genes.y.y>0.8,"present","absent")) %>% 
  dplyr::select(SampleID_short,MTAP_major,MTAP_minor,CDKN2A_major,CDKN2A_minor) -> plot_1
plot_1<-plot_1[match(sample_order, plot_1$SampleID_short),]
row.names(plot_1)<-plot_2$SampleID_short
plot_1$SampleID_short<-NULL
CN_plot<-Heatmap(t(plot_1), 
                 name = "allele", 
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE,
                 col =  c(adjustcolor( "chartreuse4", alpha.f = 0.5),"chartreuse4"),
                 rect_gp = gpar(col = "white", lwd = 2),
                 column_names_side = "bottom",
                 show_row_names = T)


dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/MTAP/"
plot_name="CDKN2A_MTAP"
# draw it
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 10) 
draw(CN_plot,heatmap_legend_side = "bottom")
dev.off()



cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
plot_2<-df_heatmap
length(sample_order)
length(plot_2$SampleID_short)
plot_2<-plot_2[match(sample_order, plot_2$SampleID_short),]
row.names(plot_2)<-plot_2$SampleID_short
plot_2$SampleID_short<-NULL

expr_plot = Heatmap(t(log10(plot_2+1)),
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE,
                    col = rev(cols),
                    na_col = "white",
                    rect_gp = gpar(col = "white", lwd = 2),
                    name = "log2normcounts",
                    column_names_side = "top",
                    show_row_names = F)

dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/MTAP/"
plot_name="CDKN2A_MTAP_expr_with_names_2"
# draw it
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 10) 
draw(expr_plot,heatmap_legend_side = "bottom")
dev.off()




