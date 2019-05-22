###################
##################
##make somatic overview driver genes
##################
#!/usr/bin/env Rscript
library(SPARQL, quietly = T)
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(BSgenome)
#available.genomes()[1:5]
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(dplyr)
library(grid)
library(gridExtra)
library(VariantAnnotation)
library("readxl")
library(stringr)
library(dndscv)
library(openxlsx)
library(data.table)
library(forcats)


#compare SNVs
vcf_list <- list.files("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Somatic_VCFs/SNV_INDEL/test/", pattern = ".vcf.gz$", full.names = TRUE)
vcf_files_names <- sub("/Users/avanhoeck/surfdrive/projects/SU2C/Colorectal-Pancreas/Somatic_VCFs/SNV_INDEL/test//","",vcf_list)
vcf_files_names <- sub("_pass_filtered.vcf.gz","",vcf_files_names)
vcf_files_names=str_split_fixed(vcf_files_names, "_", 2)[,2]

pancreas_GRList = list()  
for(i in 1:length(vcf_list)){
  # Read vcf files with VariantAnnotation
  vcf_object = readVcf(vcf_list[i], "hg19")
  seqlevels(vcf_object) <- paste("chr", seqlevels(vcf_object), sep="")
  # predict coding effects
  #coding_GR = predictCoding(vcf_object, txdb, seqSource = Hsapiens)
  # store result in list
  pancreas_GRList[[i]] = vcf_object
}

names(pancreas_GRList) = vcf_files_names

length(rowRanges(pancreas_GRList$`FR11123814`))
length(intersect(rowRanges(pancreas_GRList$FR11123814),rowRanges(pancreas_GRList$FR11123965)))
# write the muts to data.frame
coding_df = data.frame()
for(i in 1:length(pancreas_GRList)){
  temp = pancreas_GRList[[i]]
  temp2 = rowRanges(temp)
  
  temp2 <- temp2[sapply(temp2$ALT, length) == 1]
  temp2$ALT <- unlist(temp2$ALT)
  
  temp2 <- as.data.frame(temp2)
  
  temp2 %>% dplyr::select(seqnames,start,end,REF,ALT)
  rownames(temp2) =NULL
  temp2['sampleId'] = names(pancreas_GRList[i])
  temp2 <- temp2 %>% dplyr::select(sampleId,seqnames,start,end,REF,ALT)
  coding_df = rbind(coding_df, temp2)
}
unique(coding_df)

mutations <- coding_df[,c("sampleId","seqnames","start", "REF","ALT")]
colnames(mutations) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations$chr <- substr(mutations$chr, 4, nchar(as.vector(mutations$chr)))
dndsout = dndscv(mutations)
Pancreas_gene_muts <- dndsout$annotmuts

cancer_genes <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/driver_genes/Comprehensive_Characterization _Cancer_Driver_Genes_and_Mutations/1-s2.0-S009286741830237X-mmc1.xlsx", sheet = "Table S1",  col_names = TRUE))
colnames(cancer_genes)
cancer_genes <- cancer_genes %>% dplyr::filter(Cancer== "PANCAN" | Cancer== "PAAD") %>% dplyr::select(Gene, Cancer,KEY,suppressor_or_oncogene)
#cancer_genes <- cancer_genes %>% dplyr::filter(Cancer== "PAAD") %>% dplyr::select(Gene, Cancer,KEY,suppressor_or_oncogene)
cancer_genes_2 <- as.data.frame(read.table("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/driver_genes/cancer_genome_interpreter/catalog_of_cancer_genes_latest/gene_MoA.tsv", sep = "\t", header = T))
cancer_genes_2 <- cancer_genes_2  %>% dplyr::select(gene, gene_MoA)
colnames(cancer_genes_2) <- c("Gene","Gene_MoA")
unique(cancer_genes$Gene[! cancer_genes$Gene %in% cancer_genes_2$Gene])
unique(cancer_genes_2$Gene[! cancer_genes_2$Gene %in% as.factor(cancer_genes$Gene)])
cancer_genes <- cancer_genes %>% dplyr::left_join(.,cancer_genes_2,by="Gene")
unique(cancer_genes$Gene)

query_cancer_gemes_SU2C_data <- function(gene_list=cancer_genes$Gene[1:5],sample_list=Pancreas_organoid_lines){
  
  gene_list=gene_list
  sample_list=sample_list
  ##STEP1:create gene graph
  create_gene_graph <- sprintf("INSERT DATA{GRAPH <http://arne/metadata/1> {\n <http://sparqling-genomics/Gene/%s> <http://arne/inGroup> <http://arne/CancerGenes> . }}", unique(gene_list))
  for (gene_query in create_gene_graph){
    #print(gene_query)
    SPARQL(url = endpoint, curl_args = auth_options, query = gene_query)
  }
  
  ##STEP2:create sample graph
  create_sample_graph <- sprintf("INSERT DATA{GRAPH <http://arne/metadata/2> {\n <http://sparqling-genomics/Sample/%s> <http://arne/inGroup> <http://arne/sample_list> . }}", unique(sample_list))
  #SPARQL(url = endpoint, curl_args = auth_options, query = create_graph_2[1])
  for (sample_query in create_sample_graph){
    SPARQL(url = endpoint, curl_args = auth_options, query = sample_query)
  }
  
  query_genecopynumber <- sprintf(
    "
    PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
    SELECT 
    STRAFTER(STR(?sampleid),\"http://sparqling-genomics/Sample/\") AS ?sampleid
    STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\") AS ?gene 
    ?mincopynumber
    ?maxcopynumber
    WHERE {
      GRAPH <http://arne/genecopynumber>
      {
        ?row col:sampleid ?sampleid.
        ?row col:gene ?gene .
        ?row col:mincopynumber ?mincopynumber.
        ?row col:mincopynumber ?maxcopynumber.
      }
      GRAPH <http://arne/metadata/1>
        {
          ?gene <http://arne/inGroup> <http://arne/CancerGenes> .
      }
      GRAPH <http://arne/metadata/2>
      {
          ?sampleid <http://arne/inGroup> <http://arne/sample_list> .
      }
    }
    ")
  query_out_genecopynumber <- SPARQL(url = endpoint, curl_args = auth_options, query = query_genecopynumber)$results
  query_genesomaticvariant <- sprintf(
    "
    PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
    PREFIX filter_pf: <http://sparqling-genomics/vcf2rdf/FilterItem/>
    SELECT 
    STRAFTER(STR(?sampleid),\"http://sparqling-genomics/Sample/\") AS ?sampleid
    STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\") AS ?gene
    STRAFTER(STR(?chromosome), \"http://rdf.biosemantics.org/data/genomeassemblies/hg19#\") AS ?chromosome
    ?position
    STRAFTER(STR(?ref), \"http://sparqling-genomics/Sequence/\") AS ?ref
    STRAFTER(STR(?alt), \"http://sparqling-genomics/Sequence/\") AS ?alt
    STRAFTER(STR(?clonality), \"http://sparqling-genomics/Clonality/\") AS ?clonality
    STRAFTER(STR(?trinucleotidecontext), \"http://sparqling-genomics/Sequence/\") AS ?trinucleotidecontext
    STRAFTER(STR(?type), \"http://sparqling-genomics/MutationType/\") AS ?type
    ?minoralleleploidy
    ?adjustedvaf
    ?adjustedcopynumber
    ?highconfidence
    ?effect
    ?loh
    ?filter
    
    
    WHERE {
    GRAPH <http://arne/somaticvariant>
    {
    ?row col:sampleid ?sampleid.
    ?row col:gene ?gene .
    ?row col:chromosome ?chromosome .
    ?row col:position ?position .
    ?row col:ref ?ref .
    ?row col:alt ?alt .
    ?row col:type ?type .
    ?row col:clonality ?clonality .
    ?row col:adjustedvaf ?adjustedvaf .
    ?row col:trinucleotidecontext ?trinucleotidecontext .
    ?row col:minoralleleploidy ?minoralleleploidy .
    ?row col:adjustedcopynumber ?adjustedcopynumber .
    ?row col:highconfidence ?highconfidence .
    ?row col:effect ?effect .
    ?row col:loh ?loh .
    ?row col:filter ?filter .
    
    }
    GRAPH <http://arne/metadata/1>
    {
    ?gene <http://arne/inGroup> <http://arne/CancerGenes> .
    }
    GRAPH <http://arne/metadata/2>
    {
    ?sampleid <http://arne/inGroup> <http://arne/sample_list> .
    }
    FILTER (?filter = filter_pf:PASS)
    }")
  query_out_somaticvariants <- SPARQL(url = endpoint, curl_args = auth_options, query = query_genesomaticvariant)$results
  ##final STEP1:remove graph
  clear_graph <- sprintf("CLEAR GRAPH <http://arne/metadata/1>")
  SPARQL(url = endpoint, curl_args = auth_options, query = clear_graph)
  clear_graph <- sprintf("CLEAR GRAPH <http://arne/metadata/2>")
  SPARQL(url = endpoint, curl_args = auth_options, query = clear_graph)
  
  nrow(query_out_somaticvariants)
  query_out <- full_join(query_out_genecopynumber, query_out_somaticvariants, by = c("sampleid" = "sampleid", "gene" = "gene")) %>% 
    dplyr::select(sampleid,gene,mincopynumber,maxcopynumber,adjustedcopynumber,adjustedvaf,minoralleleploidy,loh,chromosome,position,ref,alt,clonality,type,effect)
  return(query_out)
}
query_out <- query_cancer_gemes_SU2C_data(gene_list=cancer_genes$Gene,sample_list=Pancreas_organoid_lines)

SNPandIndels<-query_out %>% dplyr::filter(!is.na(effect))
mutations_S_I <- SNPandIndels[,c("sampleid","chromosome","position", "ref","alt")]
colnames(mutations_S_I) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations_S_I$chr <- substr(mutations_S_I$chr, 4, nchar(as.vector(mutations_S_I$chr)))
dndsout = dndscv(mutations_S_I)
mutations_dndscv <- dndsout$annotmuts
signif_genes = dndsout$sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
sel_cv = dndsout$sel_cv
diff_genes<-dndsout$sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
names(mutations_dndscv)
mutations_dndscv %>% dplyr::select(chr,pos,ref,mut,gene,ntchange,impact,pid) ->mutations_dndscv
names(mutations_dndscv) <- c("chromosome","position","ref","alt","gene_2","ntchange","impact","pid")
query_out$chromosome <- substr(query_out$chromosome, 4, nchar(as.vector(query_out$chromosome)))
names(query_out)
query_out<-full_join(query_out, mutations_dndscv, by = c("chromosome","position","ref","alt"))
unique(query_out$impact2)
query_out <- query_out %>% dplyr::mutate(impact2 = as.character(impact)) %>% dplyr::mutate(impact2=ifelse(str_sub(ntchange,-7,-1)=="inframe","INDEL_inframe",
                                                                                                          ifelse(str_sub(ntchange,-7,-1)=="frshift","INDEL",
                                                                                                                 ifelse(str_sub(ntchange,-3,-1)=="mnv","MNV",impact2))))


cancer_genes_list<- unique(cancer_genes$Gene)
df_sample= setNames(data.frame(matrix(ncol = length(cancer_genes_list)+1, nrow = 0)),c("sampleID",cancer_genes_list))
for (i in 1:length(unique(query_out$sampleid))){
  sampleID<-unique(query_out$sampleid)[i]
  #sampleID<-"FR11124502"
  query_out_single_sample <- dplyr::filter(query_out,sampleid == sampleID)
  unique(query_out$impact2)
  unique(query_out_single_sample$gene)
  #dplyr::filter(query_out_single_sample,gene == "APC")
  #dplyr::filter(query_out_single_sample,gene == "KRAS")
  df_sample_gene=NULL
  df_sample_gene= setNames(data.frame(matrix(ncol = length(cancer_genes_list)+1, nrow = 1)),c("sampleID",cancer_genes_list))
  for (j in 1:length(unique(query_out_single_sample$gene))){
    geneID<-unique(query_out_single_sample$gene)[j]
    #geneID<-"CDKN2A"
    #geneID<-"APC"
    #geneID<-"KRAS"
    #geneID<-"TCF12"
    query_out_single_sample_gene <- dplyr::filter(query_out_single_sample,gene == geneID)
    #check whether MOA of geneID is tsg or act
    MOA <- unique(dplyr::filter(cancer_genes,Gene == geneID)$Gene_MoA_modified)
    if ( MOA=="LoF") {
      #print("LOF")
      allel_1=NA
      allel_2=NA
      #check for somatic mutations
      if(any(!is.na(unique(query_out_single_sample_gene$impact2)))){
        mutations <- dplyr::filter(query_out_single_sample_gene,!is.na(impact2)) %>% distinct(position, .keep_all = TRUE)
        target <-c("Nonsense","INDEL","Missense","Essential_Splice","MNV","INDEL_inframe","Synonymous")
        #mutations <- mutations[order(match(mutations$impact2,target)),]
        mutations <- dplyr::filter(mutations,impact2==target[min(match(mutations$impact2,target))]) %>% dplyr::arrange(minoralleleploidy)
        allel_1=mutations$impact2[1]
        if(any(mutations$minoralleleploidy)<0.2){
          allel_2="LOH"
        }
        #remove Synonymous variants
        list_remove <- c("Synonymous")
        if (allel_1 %in% list_remove){
          allel_1=NA
          allel_2=NA
        }
        #check for deep deletion
      }else{
        if(max(query_out_single_sample_gene$mincopynumber)<0.2){
          allel_1="Deep deletion"
          allel_2="Deep deletion"
        }
      }
      df_sample_gene[geneID]<-paste(allel_2,allel_1,sep = "/")
    }
    if (MOA=="Act") {
      #print("Act")
      allel=NA
      #check for somatic mutations
      
      if(max(query_out_single_sample_gene$mincopynumber)>5){
        allel="Amplification"}
      else{if(any(!is.na(unique(query_out_single_sample_gene$impact2)))){
        mutations <- dplyr::filter(query_out_single_sample_gene,!is.na(impact2)) %>% distinct(position, .keep_all = TRUE)
        target <-c("Missense","MNV","INDEL_inframe","INDEL","Nonsense","Essential_Splice","Synonymous")
        #mutations <- mutations[order(match(mutations$impact2,target)),]
        mutations <- dplyr::filter(mutations,impact2==target[min(match(mutations$impact2,target))]) %>% dplyr::arrange(minoralleleploidy)
        allel=mutations$impact2[1]
        if(min(mutations$mincopynumber)>6){
          allel="Amplification"
        }
        #remove Synonymous variants
        list_remove <- c("Synonymous")
        if (allel %in% list_remove){
          allel=NA
        }
      }}
      df_sample_gene[geneID]<-allel
    }
    df_sample_gene["sampleID"]<-sampleID
  }
  df_sample=rbind(df_sample,df_sample_gene)
}
df_sample[which(names(df_sample)%in%c("sampleID","RNF43","CTNNB1","APC"))]
#remove sample T27
df_sample<-df_sample[which(df_sample$sampleID!="FR11123214"),]
row.names(df_sample)<-df_sample$sampleID
df_sample$sampleID<-NULL

#remove genes with no variants
df_sample<-df_sample[ ,-which(names(df_sample) %in% names(which(colSums(df_sample == "NA/NA")==25)))]
df_sample<-df_sample[ ,-which(names(df_sample) %in% names(which(colSums(is.na(df_sample))==25)))]
#create suppressor and drives dfs
LoFs <- unique(cancer_genes[which(cancer_genes$Gene_MoA_modified=="LoF"),]$Gene)
Acts <- unique(cancer_genes[which(cancer_genes$Gene_MoA_modified=="Act"),]$Gene)
LoFs <- df_sample[which(names(df_sample) %in%LoFs )]
Acts <- df_sample[which(names(df_sample) %in%Acts )]
LoFs_names<-names(sort(colSums(LoFs != "NA/NA"),decreasing = T))
Acts_names<-names(sort(colSums(!is.na(Acts)),decreasing = T))
LoFs<-dplyr::select(LoFs,rev(LoFs_names))
Acts<-dplyr::select(Acts,rev(Acts_names))
LoFs <- tibble::rownames_to_column(LoFs)
names(LoFs)[1]<-("sampleId")
Acts <- tibble::rownames_to_column(Acts)
names(Acts)[1]<-("sampleId")
plotting_data<-left_join(LoFs, Acts, by = "sampleId")
plotting_data<-left_join(plotting_data, dplyr::select(metadata,sampleId,SampleID_short), by = "sampleId")
plotting_data<-dplyr::select(plotting_data,c("SampleID_short",rev(LoFs_names),rev(Acts_names)))
plotting_data[plotting_data=="NA/NA"]<-NA
row.names(plotting_data) <- plotting_data$SampleID_short
plotting_data$SampleID_short<-NULL
#driver_names_order <-names(sort(colSums(!is.na(plotting_data)),decreasing = T))[1:30]
driver_names_selection_1<-unique(cancer_genes[which(cancer_genes$Cancer=="PAAD"),]$Gene)
driver_names_selection_2<-diff_genes$gene_name
driver_names_selection_Tuveson_3<-c("KRAS","TP53","CDKN2A","SMAD4","ARID1A","RNF43","GNAS","ATM","MYC","PIK3CA","PBRM1","PALB2","MAP2K1","KDM6A","ERBB2","BRCA1","AKT2")
driver_names_selection_Sato_4<-c("KRAS","CDKN2A","TP53","SMAD4","TGFBR2","GNAS","APC","AXIN2","RNF43","CTNNB1","ARID1A","ARID1B","ARID2","SMARCA4","BRCA1")
driver_names_selection<-unique(c(driver_names_selection_1,driver_names_selection_2,driver_names_selection_Tuveson_3,driver_names_selection_Sato_4))
driver_names_selection<- driver_names_selection[!driver_names_selection %in% c("CDKN2A.p16INK4a","CDKN2A.p14arf","U2AF1","SMARCA4", "ARID1B", "AKT2","PALB2")]
plotting_data<-dplyr::select(plotting_data,rev(driver_names_selection))
driver_names_order <-names(sort(colSums(!is.na(plotting_data)),decreasing = T))
plotting_data<-dplyr::select(plotting_data,rev(driver_names_order))
names(plotting_data)
plotting_data <- tibble::rownames_to_column(plotting_data)
names(plotting_data)[1]<-("sampleId")
plotting_data %>% 
  dplyr::arrange(factor(KRAS, levels = c("Missense","Amplification","MNV","NA")),
                 factor(TP53, levels = c("LOH/Nonsense","LOH/INDEL","LOH/Missense","Deep deletion/Deep deletion","NA/NA","NA/INDEL")),
                 factor(CDKN2A,levels = c("LOH/INDEL","LOH/Missense","Deep deletion/Deep deletion","LOH/Nonsense","NA/NA","NA/INDEL")),
                 rev(EEF2),rev(SMAD4)) -> plotting_data


####
#add tumor only samples to dataframe
####
tumoronly <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/driver_genes/singles_Pancreas/filtered/overview_revised.xlsx", col_names = T,sheet = "Sheet2"))


tumoronly <- tumoronly %>% dplyr::mutate(impact2 = as.character(TYPE)) %>% dplyr::mutate(impact2=ifelse(str_sub(impact2,1,5)=="misse","Missense",
                                                                                                        ifelse(str_sub(impact2,1,5)=="frame","INDEL",
                                                                                                               ifelse(str_sub(impact2,1,5)=="stop_","Nonsense",
                                                                                                                      ifelse(str_sub(impact2,1,5)=="synon","synonymous",impact2)))))
Genes_tumoronly <- unique(tumoronly$GENE)
tumoronly <- tumoronly%>% dplyr::select(sampleId, impact2, GENE) %>% 
  dplyr::mutate(sampleId = as.factor(sampleId),
                GENE= as.factor(GENE)) %>% 
  spread(., GENE, impact2)
tumoronly <- left_join(tumoronly, dplyr::select(metadata,sampleId,SampleID_short), by = "sampleId")
row.names(tumoronly) <- tumoronly$SampleID_short
tumoronly$sampleId <-  NULL
tumoronly$SampleID_short<-  NULL

to_merge_df <- plotting_data[1,]
row.names(to_merge_df) <- to_merge_df$sampleId
to_merge_df$sampleId <- NULL
one <- tibble::rownames_to_column(as.data.frame(t(to_merge_df)),var = "sampleId")
two <- tibble::rownames_to_column(as.data.frame(t(tumoronly)),var = "sampleId")
one <- left_join(one, two, by = "sampleId")
row.names(one) <- one$sampleId
one$sampleId <-  NULL
to_merge_df <- as.data.frame(t(one))
to_merge_df <- tibble::rownames_to_column(to_merge_df,var = "sampleId")
to_merge_df[1,] <- "grey"
to_merge_df[1,1] <- "blank"
row.names(to_merge_df) <- to_merge_df$sampleId
to_merge_df$sampleId <- NULL
to_merge_df <- to_merge_df[c("blank","T2","T1","T30","T31","T29"),]
to_merge_df <- tibble::rownames_to_column(to_merge_df,var = "sampleId")
row.names(to_merge_df) <- NULL
######


row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

modify_gene_names <- function(gene_names = names(plotting_data),df=plotting_data){
  
  percent <- function(x, digits = 0, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  
  data <- df
  n_samples<-nrow(data)
  data[!is.na(data)] <- 1
  data[is.na(data)] <- 0
  data <- mutate_all(data, function(x) as.numeric(as.character(x)))
  percentage <- percent(colSums(as.matrix(data), na.rm = FALSE)/n_samples)
  
  df2 <- data.frame(Gene=names(df),
                    percentage_ID=percentage)
  
  
  df2 <- left_join(df2, dplyr::select(cancer_genes,Gene,Gene_MoA_modified), by = "Gene")
  df2 <- df2[!duplicated(df2), ] 
  df2 %>%
    dplyr::mutate(MOA = ifelse(Gene_MoA_modified == "LoF", " ", 
                               ifelse(Gene_MoA_modified == "Act", "*","NA"))) -> df2
  df2$percentage_ID_2 <- with(df2, paste("(",percentage_ID,")", sep=""))
  df2$Gene_2 <- with(df2, paste(MOA,Gene, sep=""))
  cols<-c("Gene_2","percentage_ID_2")
  df2$unite <- apply( df2[ , cols ] , 1 , paste , collapse = " " )
  return(df2$unite)
}
#gene_names<-names(plotting_data)
gene_names<-modify_gene_names(plotting_data)

plotting_data <- tibble::rownames_to_column(plotting_data,var = "sampleId")
plotting_data <- rbind(plotting_data,to_merge_df)
row.names(plotting_data) <- plotting_data$sampleId
plotting_data$sampleId <- NULL

plotting_data[plotting_data=="NA/NA"]<-NA
plotting_data[plotting_data=="Deep deletion/Deep deletion"]<-"Deep deletion"
names(plotting_data)<-seq(1, ncol(plotting_data), by=1)
m <- as.matrix(plotting_data)
dt <- data.table(melt(m))
dt <- dt[, strsplit(as.character(value), "/"), by=list(Var1, Var2)]  # this expands "X/Y/Z" into three rows
dt[, shift:=(1:(.N))/.N - 1/(2 * .N) - 1/2, by=list(Var1, Var2)]
dt[, height:=1/.N, by=list(Var1, Var2)]
dt$test=dt$Var2+dt$shift
dt$test2=as.numeric(format(dt$test,digits =0))
dt$height2=as.numeric(1)


plot<- ggplot() + 
  geom_tile(data = dt,aes(Var1,y=test, fill=V1, height=height),
            color="grey", size=0.5) + xlab('Organoid sample') + ylab('Gene')+
  scale_y_continuous(breaks = seq(1, ncol(plotting_data), by=1), labels = gene_names)+
  scale_fill_manual(values=colorpalette_drivers,na.value = "white") +
  guides(fill=guide_legend(title="Effect"))+
  theme(axis.title.y=element_text(size=20,vjust=3),
        axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=20),
        axis.text.x = element_text(colour="grey20",size=20,angle=45,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=0.5),
        strip.text.x=element_text(size=15),
        strip.text.y=element_text(size=15,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'grey20',size = 2),
        strip.background=element_blank())+
  geom_tile(data = dt,aes(Var1,y=test2, height=height2),
            color="black", size=0.5,alpha=0)


save_plot("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/Figure X - driver_and_suppressor_genes_2.png", plot, base_width = 15, base_height = 10)
dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/"
plot_name="drivers_and_supressors"
# draw it
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 15, height = 10) 
plot(plot)
dev.off()



