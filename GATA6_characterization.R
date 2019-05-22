library(SPARQL, quietly = T)
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

########################
########################
###GATA6, FOXA2, and HNF4A WNT10A
candidates <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/sato_et_al/mmc4.xlsx", col_names = T,sheet = "Sheet2"))
candidates <- candidates$names
#candidates <- overlap_sato_DE
#candidates <- c("GATA6")
candidates <- all.genes[which(all.genes$hgnc_symbol %in% candidates),]
candidates %>% dplyr::arrange(match(hgnc_symbol,c("GATA6","FOXA2","HNF4A","WNT10A"))) -> candidates
candidates$ensembl_gene_id <- factor(candidates$ensembl_gene_id, levels = c(candidates$ensembl_gene_id))
candidates$hgnc_symbol <- factor(candidates$hgnc_symbol, levels = c(candidates$hgnc_symbol))

df_heatmap = as.data.frame(t(normalized[which(row.names(normalized) %in% candidates$ensembl_gene_id),]))
df_heatmap = as.data.frame(t(RPKM[which(row.names(RPKM) %in% candidates$ensembl_gene_id),]))

df_heatmap <- df_heatmap[as.vector(candidates$ensembl_gene_id)]
colnames(df_heatmap) <- as.vector(candidates$hgnc_symbol)
#df_heatmap$GATA6 <- df_heatmap[order(df_heatmap[,"GATA6"],decreasing = T),]
df_heatmap  %>% tibble::rownames_to_column("SampleID_short") %>% left_join(., dplyr::select(metadata,SampleID_short,germline_sample_ID,sampleId,medium), by = "SampleID_short") %>% as.data.frame()-> df_heatmap
#remove paired samples
df_heatmap<-df_heatmap[which(df_heatmap$germline_sample_ID =="NA"),]

plotGATA6 <- function(data = df_heatmap, pdf.path = NULL){
  data$medium <- factor(data$medium, levels = c("NIC","W-E"))
  data$GATA6_log <- data$GATA6+1
  data$GATA6_log <- log10(data$GATA6_log)
  
  #check normality
  print("run shapiro to test for normality")
  print(shapiro.test(data$GATA6))
  #not normal distribution -> wilcoxon test
  stats_data <- as.data.frame(compare_means(GATA6 ~ medium, data = data,method = "wilcox.test"));print(stats_data)
  stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
  #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
  stats_data_plot <- stats_data_plot[,c("V1")]
  stats_data_plot <- lapply(stats_data_plot, as.character)
  stats_data_plot <- as.list(stats_data_plot)
  
  nn = data %>% group_by(medium) %>% tally()
  
  set_max <- data %>% group_by(medium) %>%
    dplyr::summarize(count = n(),mean_sign = mean(GATA6),median_sign = median(GATA6), SD_sign = sd(GATA6),max(GATA6), sum(median_sign,SD_sign)) 
  set_max <- as.data.frame(set_max)
  number  <- nrow(set_max)
  set_max <- max(set_max$`max(GATA6)`)
  set_max<-4
  Y_pos<-2000
  Y_pos<-3.5
  
  plot <- ggboxplot(data, x = "medium", y = "GATA6_log",
                    color = "black",fill = "medium",palette = c("#00AFBB", "#E7B800"))+ 
    #stat_compare_means(aes(label = ..p.signif..))+
    #stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
    stat_compare_means(method = "t.test", label.y = round(set_max),size = 4)+
    xlab(c(""))+
    ylab(c("GATA6 expression level \n log10 normalized "))+
    #scale_y_continuous(trans='log10',
    #                   breaks=10**(1:round(max(Y_pos))),
    #                   limits = c(0.4,10^round(max(Y_pos))),
    #                   labels = scales::comma)+
    scale_y_continuous(limits = c(1.8, 4))+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1),
    #      legend.position = "none")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=1.9, colour="grey20", size=4)
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,3,4)
    plot(plot)
    dev.off()}
}
plotGATA6(data=df_heatmap)
plotGATA6(data=df_heatmap,pdf.path=sprintf("%sGATA6_expression_comparison_2.pdf",outdir))

query_GATA6copynumber_SU2C <- function(cohort = Pancreas_organoid_lines,gene="GATA6"){
  df_out=NULL
  df_out=data.frame()
  print(cohort)
  
  geneID=gene
  for (sampleID in unique(cohort)) {
    #sampleID="FR11123716"
    print(sampleID)
    query <- sprintf(
      "
      
      PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
      PREFIX sampleid_pf: <http://sparqling-genomics/Sample/>
      PREFIX gene_pf: <http://sparqling-genomics/Gene/>
      
      SELECT 
      STRAFTER(STR(?sampleid),\"http://sparqling-genomics/Sample/\") AS ?sampleid
      STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\") AS ?gene
      ?mincopynumber
      ?maxcopynumber
      ?meancopynumber
      
      
      FROM <http://arne/genecopynumber>
      WHERE 
      { ?row col:sampleid sampleid_pf:%s .
      ?row col:gene gene_pf:%s .
      ?row col:mincopynumber ?mincopynumber.
      ?row col:maxcopynumber ?maxcopynumber.
      ?row col:meancopynumber ?meancopynumber.
      
      
      BIND(sampleid_pf:%s AS ?sampleid)
      BIND(gene_pf:%s AS ?gene)}
      ",sampleID,geneID,sampleID,geneID)
    query_out <- SPARQL(url = endpoint, curl_args = auth_options, query = query)$results
    df_out=rbind(df_out,query_out)
    df_out_final=rbind(df_out)
  }
  return(df_out_final)
  rm(df_out_1)
  rm(df_out_2)
}
query_GATA6copynumber_SU2C_data <- query_GATA6copynumber_SU2C(cohort = Pancreas_organoid_lines,gene="GATA6")
query_GATA6copynumber_SU2C_data %>% dplyr::mutate(sampleId=sampleid)%>% as.data.frame() -> query_GATA6copynumber_SU2C_data
dplyr::left_join(df_heatmap, query_GATA6copynumber_SU2C_data, by = "sampleId") %>% as.data.frame()-> df_heatmap

row.names(df_heatmap) <- df_heatmap$SampleID_short
df_heatmap$SampleID_short <- NULL
df_heatmap$germline_sample_ID <- factor(df_heatmap$germline_sample_ID, levels = sort(unique(df_heatmap$germline_sample_ID)))


number_of_colors <- length(unique(df_heatmap$germline_sample_ID))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_germline=sample(col_vector, number_of_colors)
names(col_germline) = unique(df_heatmap$germline_sample_ID)
col_germline["NA"] = "white"


ha01 = rowAnnotation(df = df_heatmap["medium"],
                     width = unit(0.5, "cm"),
                     na_col = "white",
                     col = list(medium=c("W-E"="paleturquoise4","NIC"="palevioletred3","NA"="white")))

ha02 = rowAnnotation(df = df_heatmap["germline_sample_ID"],
                     width = unit(0.5, "cm"),
                     na_col = "white",
                     col = list(germline_sample_ID=col_germline),
                     show_legend = FALSE)

group <- df_heatmap[,as.vector(candidates$hgnc_symbol)]
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
?Heatmap
mapke_1 = Heatmap(as.matrix(log2(group+1)),
                  cluster_rows = T, 
                  cluster_columns = T,
                  col = colorRamp2(c(0,8, 12), c("#4474BC","#E9F3D2", "#DE3928")),
                  name = "log2normcounts",
                  column_names_side = "top",
                  row_names_side = "left")
group_2 <- df_heatmap["meancopynumber"]
colnames(group_2) <- "GATA6"
mapke_2 = Heatmap(group_2,
                  cluster_rows = FALSE, 
                  cluster_columns = FALSE,
                  col = colorRamp2(c(0,2,4,8), c("royalblue","lightsteelblue1", "orangered","orangered4")),
                  na_col = "white",
                  name = "copy number",
                  column_names_side = "top",
                  show_row_names = FALSE)

ht_list = mapke_1 +mapke_2 + ha01 + ha02
ht_list = mapke_1 + ha01 + ha02

#draw(ht_list, heatmap_legend_side = "right")
plot_heatmap = paste(outdir,"heatmap_GATA6_D.pdf", sep = "/")
pdf(plot_heatmap , useDingbats = F, width = 6, height = 6) 
draw(ht_list, heatmap_legend_side = "right")
dev.off()




# ------ DIFFERENTIAL EXPRESSION ------ 

# DE WT VS MUTANT OVERALL (both liver & small intestine)

# Test for differential expression
# H0 : the measurements come from the same distribution and the gene is being expressed 
# at the same level across conditions

# general function for differential expression calculation
# annotate with gene symbols & description


# DE MUT & WT
#whitelist_NIC <- metadata[which(metadata$medium=="NIC"),]$SampleID_short
#whitelist_WE <- metadata[which(metadata$medium=="W-E"),]$SampleID_short
#whitelist <- c(whitelist_NIC,whitelist_WE)

whitelist_NIC <- metadata[which(metadata$medium=="NIC"&metadata$germline_sample_ID=="NA"),]$SampleID_short
whitelist_WE <- metadata[which(metadata$medium=="W-E"&metadata$germline_sample_ID=="NA"),]$SampleID_short
whitelist <- c(whitelist_NIC,whitelist_WE)
countsTable %>% dplyr::select(whitelist) -> countsTable_DE
colnames(countsTable_DE)
conds = factor(c(rep("NIC",length(whitelist_NIC)), rep("WE",length(whitelist_WE))), levels = c("NIC", "WE"))
DE_overall = diff_expr(countsTable_DE, conds, org.Hs.eg.db)
head(DE_overall[,-1])
nrow(DE_overall)
# Find interesting DE genes
threshold = 2
DE_overall_int = subset(DE_overall, DE_overall$pval < 0.05 & abs(DE_overall$log2FoldChange_pseudo) > threshold)
View(DE_overall_int)
nrow(DE_overall_int)
hist(DE_overall_int$pval)
DE_overall[which(DE_overall$SYMBOL=="GATA6"),]
heatmap.2(as.matrix(DE_overall_int[whitelist]), scale = 'none',col = rev(heat.colors(256)),   margins=c(5,5))
DE_overall_int$SYMBOL
write.table(DE_overall_int, file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/RNAseq/DE_medium_comp_excluding_paired_data.txt", sep = "\t", quote = F, row.names = F)




#intersect
overlap_sato_DE<-intersect(DE_overall_int$SYMBOL, candidates$hgnc_symbol)
length(overlap_sato_DE)

