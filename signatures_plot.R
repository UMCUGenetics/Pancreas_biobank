###
#mutational signatures
library("readxl")
library(reshape2)
library(stringr)
library(tidyr)
library(reshape)
library(plyr)
library(dplyr)
library(ggplot2) 
library(MutationalPatterns)
library(cowplot)

colorpalette = c("Signature.1" =  '#8dd3c7',
                 "Signature.2" =  '#ffffb3',
                 "Signature.3" =  '#bebada',
                 "Signature.4" =  '#fb8072',
                 "Signature.5" =  '#80b1d3',
                 "Signature.6" =  '#fdb462',
                 "Signature.7" =  '#b3de69',
                 "Signature.8" =  '#fccde5',
                 "Signature.9" =  '#d9d9d9',
                 "Signature.10" = '#ff1417' ,
                 "Signature.11" = '#ff6611' ,
                 "Signature.12" = '#c4ff00' ,
                 "Signature.13" = '#ff8844' ,
                 "Signature.14" = '#ffee55' ,
                 "Signature.15" = '#ffff99' ,
                 "Signature.16" = '#78FA37' ,
                 "Signature.17" = '#aacc22' ,
                 "Signature.18" = '#bbdd77' ,
                 "Signature.19" = '#c8cf82' ,
                 "Signature.20" = '#92a77e' ,
                 "Signature.21" = '#5599ee' ,
                 "Signature.22" = '#0088cc' ,
                 "Signature.23" = '#226688' ,
                 "Signature.24" = '#175279' ,
                 "Signature.25" = '#557777' ,
                 "Signature.26" = '#ddbb33' ,
                 "Signature.27" = '#d3a76d' ,
                 "Signature.28" = '#a9834b' ,
                 "Signature.29" = '#aa6688',
                 "Signature.30" = '#767676',
                 "Signature.A" = '#458B00' ,
                 "Signature.B" = '#D2691E' ,
                 "Signature.C" = '#6495ED' ,
                 "Signature.D" = '#A2CD5A' ,
                 "Signature.E" = '#CD3333' ,
                 "Signature.F" = '#7AC5CD' ,
                 "Signature.G" = '#009ACD' ,
                 "Signature.H" = '#CD2626' ,
                 "Signature.I" = '#FFB90F' ,
                 "Signature.J" = '#76EEC6' ,
                 "Signature.K" = '#EEB422' ,
                 "Signature.L" = '#97FFFF' ,
                 "Signature.M" = '#E9967A' ,
                 "Signature.N" = '#5F9EA0')
colorpalette_drivers = c("Amplification" =  '#b3de69',
                         "Deep deletion" =  'coral1',
                         "INDEL" =  'skyblue3',
                         "Synonymous" =  '#fb8072',
                         "NA" =  'white',
                         "Missense" =  '#fdb462',
                         "LOH" =  '#ffffb3',
                         "MNV" =  '#8dd3c7',
                         "Essential_Splice" =  '#d9d9d9',
                         "INDEL_inframe" = 'skyblue',
                         "Nonsense" = 'plum3' ,
                         "grey" = 'gray80')

COLORS_indels = c("cadetblue", "darkgoldenrod", "darkolivegreen4","firebrick","salmon", "sienna3", "slateblue", "deepskyblue3", "mistyrose4")
COLORS_indels = c("#8dd3c7","#ffffb3","slateblue","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#ff1417" ,"#ff6611" ,"#c4ff00" ,"#ff8844" ,"#ffee55" ,"#ffff99" ,"#78FA37" ,"#aacc22" ,"#bbdd77" ,"#c8cf82" ,"#92a77e" ,"#5599ee" ,"#0088cc" ,"#226688" ,"#175279" ,"#557777" ,"#ddbb33" ,"#d3a76d" ,"#a9834b" ,"#aa6688","#767676","#458B00" ,"#D2691E" ,"#6495ED" ,"#A2CD5A" ,"#CD3333" ,"#7AC5CD" ,"#009ACD" ,"#CD2626" ,"#FFB90F")


#metadata_panc_organoids <- as.data.frame(read_excel("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/wetransfer/all data pool 8-1-2019.xlsx", "Sheet1"),row.names = F)
metadata <- as.data.frame(read_xlsx("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/fingerprint_pancreas_correct_match.xlsx", col_names = T,sheet = "Overview"))

clinical <- data.frame(read.csv("~/surfdrive/Shared/Sig17/HMF_data/Metadata_DR10-update/DR-047_metadata.tsv",sep = "\t" ,na="NULL"))
DR47 <- data.frame(read.csv("~/surfdrive/Shared/Sig17/HMF_data/Metadata_DR10-update/DR47_samples.txt",sep = "\t" ,na="NULL",header = 1))
clinical_DR47 <- clinical[which(substr(sort(clinical$sampleId), 1, 12) %in% substr(sort(DR47$DR47_samples), 1, 12)),]
nrow(clinical_DR47)
clinical_DR47$preTreatments_dub=clinical_DR47$preTreatments
clinical_DR47$treatment_dub=clinical_DR47$treatment
max_fields_preTreatments=max(sapply(strsplit(as.character(clinical_DR47$preTreatments_dub),'/'),length))
max_fields_treatment=max(sapply(strsplit(as.character(clinical_DR47$treatment),'/'),length))
clinical_DR47=separate(data=clinical_DR47, col=preTreatments_dub, into=paste0("preTreatment_",1:as.numeric(max_fields_preTreatments)), sep="/")
clinical_DR47=separate(data=clinical_DR47, col=treatment_dub, into=paste0("treatment_",1:as.numeric(max_fields_treatment)), sep="/")

sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
cancer_signatures = as.matrix(cancer_signatures[,4:33])



###HMF_total#####
filenames=list.files(path="~/surfdrive/projects/HMF/Metadata_DR10-update/snv_HMF_update/", pattern = "_unstranded96.txt",full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x, header = T, sep = "\t", row.names = NULL)})
my_matrix=datalist[[1]]$X
datalist2 = lapply(datalist, function(x) { x["X"] <- NULL; x })
datalist3_ALL = do.call(cbind, datalist2)
colnames(datalist3_ALL)=sapply(strsplit(colnames(datalist3_ALL),"_"), function(x) x[2])
colClean <- function(x){ colnames(x) <- paste0(colnames(x), "_ALL"); x }
datalist3_ALL = colClean(datalist3_ALL)
datalist4_ALL = cbind(as.data.frame(my_matrix),datalist3_ALL)
colnames(datalist4_ALL)[1] = "matrix"
row.names(datalist4_ALL) = datalist4_ALL$matrix
datalist4_ALL$matrix=NULL

#make somatic file of fluoropyrimidines
##HMF data -> select all first biopsies
somatic_DR47=as.data.frame(t(datalist4_ALL))
row.names(somatic_DR47)=gsub("_ALL","",row.names(somatic_DR47))
somatic_DR47_firstbiopsies_IDs <- grep("T$", row.names(somatic_DR47), value=T)
somatic_DR47_firstbiopsies <- t(somatic_DR47[which(row.names(somatic_DR47) %in% somatic_DR47_firstbiopsies_IDs),])


###Pancreas_organoids_lines#####
filenames=list.files(path="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Somatic_VCFs/SNV_matrices/", pattern = "_unstranded96.txt",full.names=TRUE)
#filenames=list.files(path="~/surfdrive/Shared/Bastiaan/Sig17/HMF_data/cohort_analyse/nature_comm_pancreas_organoids/matrices/", pattern = "SNV_unstranded96.txt",full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x, header = T, sep = "\t", row.names = NULL)})
my_matrix=datalist[[1]]$X
datalist2 = lapply(datalist, function(x) { x["X"] <- NULL; x })
datalist3_pancreas_organoids_data = do.call(cbind, datalist2)
colClean <- function(x){ colnames(x) <- gsub("_PASS_SNV", "", colnames(x)); x }
datalist3_pancreas_organoids_data = colClean(datalist3_pancreas_organoids_data)
datalist4_pancreas_organoids_data = cbind(as.data.frame(my_matrix),datalist3_pancreas_organoids_data)
colnames(datalist4_pancreas_organoids_data)[1] = "matrix"
row.names(datalist4_pancreas_organoids_data) = datalist4_pancreas_organoids_data$matrix
datalist4_pancreas_organoids_data$matrix=NULL
names(datalist4_pancreas_organoids_data)=str_split_fixed(names(datalist4_pancreas_organoids_data), "_", 2)[,2]
datalist4_pancreas_organoids_data <- datalist4_pancreas_organoids_data[,which(names(datalist4_pancreas_organoids_data) %in% metadata$sampleId)]

ncol(datalist4_pancreas_organoids_data)

mut_matrix_cohort = cbind(somatic_DR47_firstbiopsies,datalist4_pancreas_organoids_data)
somatic_df=as.data.frame(t(mut_matrix_cohort))
somatic_df$mut_load=rowSums(somatic_df)
somatic_df_rel = sweep(somatic_df[,1:96],1,(rowSums(somatic_df[,1:96])),`/`)
colClean <- function(x){ colnames(x) <- paste0(colnames(x), "_rel"); x }
somatic_df_rel = colClean(somatic_df_rel)
somatic_df=cbind(somatic_df,somatic_df_rel)
somatic_df <- tibble::rownames_to_column(somatic_df, "sampleId")
rm(somatic_df_rel)
refit <- fit_to_signatures(mut_matrix_cohort, cancer_signatures)
refit = as.data.frame(t(refit$contribution))
refit_rel = sweep(refit[,1:30],1,(rowSums(refit[,1:30])),`/`)
colClean <- function(x){ colnames(x) <- paste0(colnames(x), "_rel"); x }
refit_rel = colClean(refit_rel)
refit=cbind(refit,refit_rel)
refit <- tibble::rownames_to_column(refit, "sampleId")
#refit_rel = sweep(refit[,1:49],1,(rowSums(refit[,1:49])),`/`)
#refit_rel = sweep(refit[,1:65],1,(rowSums(refit[,1:65])),`/`)
rm(refit_rel)
all_final=merge(somatic_df,refit,by="sampleId")
rm(refit)
length(all_final$sampleId)

#get IDs of only first biopsies
firstbiopsies_IDs <- grep("T$", clinical_DR47$sampleId, value=T)
clinical_DR47_firstbiopsies_IDs <- clinical_DR47[which(clinical_DR47$sampleId %in% firstbiopsies_IDs),]
somatic_clinical=dplyr::left_join(all_final,clinical_DR47_firstbiopsies_IDs,by="sampleId")
somatic_clinical$Cohort="NA"
somatic_clinical[grepl("^CPCT", somatic_clinical$sampleId, perl = TRUE),]$Cohort <- "Metastatic"
somatic_clinical[grepl("^DRUP", somatic_clinical$sampleId, perl = TRUE),]$Cohort <- "Metastatic"
somatic_clinical[grepl("^FR", somatic_clinical$sampleId, perl = TRUE),]$primaryTumorLocation <- "Pancreas"
somatic_clinical[grepl("^FR", somatic_clinical$sampleId, perl = TRUE),]$Cohort <- "Organoid"
SNVs_Pancreas = somatic_clinical %>% 
  filter(primaryTumorLocation == "Pancreas") %>% 
  dplyr::select(sampleId,Cohort,matches("Signature"),-ends_with("rel"))
SNVs_Pancreas$sampleId
#rownames(SNVs_Pancreas)<-SNVs_Pancreas$sampleId
nrow(SNVs_Pancreas)
ncol(SNVs_Pancreas)
##INDELS
#locate data in list
file.indel_list_colon_pancr <- list.files("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Somatic_VCFs/INDEL_matrices/", pattern = "indel.txt", full.names = TRUE)
file.indel_list_hmf <- list.files("~/surfdrive/projects/HMF/Metadata_DR10-update/indel_HMF_update/", pattern = ".txt", full.names = TRUE)
file.indel_list  <-c(file.indel_list_colon_pancr,file.indel_list_hmf)

#read samples in list
file.indel_list <- lapply(file.indel_list,function(i){
  read.table(i, header = TRUE, sep = "\t", row.names = 1)
})

# number names of samples
length(file.indel_list)
#combine samples in dataframe
df_indels <- bind_cols(file.indel_list)
colnames(df_indels) <- sub("*_PASS_indel.vcf.gz","",colnames(df_indels))
colnames(df_indels) <- sub("*_PASS_indel","",colnames(df_indels))
colnames(df_indels)=sapply(strsplit(colnames(df_indels),"_"), function(x) x[2])
df_indels = as.data.frame(t(as.data.frame(sapply(df_indels, as.numeric))))
colnames(df_indels) <- rownames(file.indel_list[[1]])
df_indels_pancreas <- df_indels[which(row.names(df_indels) %in% SNVs_Pancreas$sampleId),]
df_indels_pancreas$sampleId <- rownames(df_indels_pancreas)
rownames(df_indels_pancreas) <- NULL

##SVs
#locate data in list
file.SV_list_colon_pancr <- list.files("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Somatic_VCFs/SV/matrices/", pattern = "_sv.txt", full.names = TRUE)
file.SV_list_hmf <- list.files("~/surfdrive/projects/HMF/Metadata_DR10-update/sv_HMF_update/", pattern = "_mutmatsv.txt", full.names = TRUE)
file.SV_list  <-c(file.SV_list_colon_pancr,file.SV_list_hmf)

#read samples in list
file.SV_list <- lapply(file.SV_list,function(i){
  read.table(i, header = TRUE, sep = "\t", row.names = 1)
})

# number names of samples
length(file.SV_list)
#combine samples in dataframe
df_SVs <- bind_cols(file.SV_list)
colnames(df_SVs)=sapply(strsplit(colnames(df_SVs),"_"), function(x) x[2])
df_SVs = as.data.frame(t(as.data.frame(sapply(df_SVs, as.numeric))))
colnames(df_SVs) <- rownames(file.SV_list[[1]])
df_SVs_pancreas <- df_SVs[which(row.names(df_SVs) %in% SNVs_Pancreas$sampleId),]
df_SVs_pancreas$sampleId <- rownames(df_SVs_pancreas)
rownames(df_SVs_pancreas) <- NULL

Pancreas <- dplyr::inner_join(SNVs_Pancreas, df_indels_pancreas, by = "sampleId")
Pancreas <- dplyr::inner_join(Pancreas, df_SVs_pancreas, by = "sampleId")
# cluster samples based on eucledian distance between relative contribution
names(Pancreas)
Pancreas %>% dplyr::mutate(
  "repeat deletion 1bp"=repeat_deletion_len.1,
  "repeat deletion 2bp"=repeat_deletion_len.2,
  "repeat deletion >2bp"=repeat_deletion_len.3+repeat_deletion_len.4+repeat_deletion_len.5,
  "repeat insertion 1bp"=repeat_insertion_len.1,
  "repeat insertion 2bp"=repeat_insertion_len.2,
  "repeat insertion >2bp"=repeat_insertion_len.3+repeat_insertion_len.4+repeat_insertion_len.5,
  "micro hom deletion"=micro_hom_deletion_BiMH.1+micro_hom_deletion_BiMH.2+micro_hom_deletion_BiMH.3+micro_hom_deletion_BiMH.4+micro_hom_deletion_BiMH.5,
  "micro hom insertion"=micro_hom_insertion_BiMH.1+micro_hom_insertion_BiMH.2+micro_hom_insertion_BiMH.3+micro_hom_insertion_BiMH.4+micro_hom_insertion_BiMH.5,
  "no context 1bp"=no_context_len.1,
  "no context 2bp"=no_context_len.2,
  "no context >2bp"=no_context_len.3+no_context_len.4+no_context_len.5,
  "small deletions"=DEL_0e00_1e03_bp,
  "medium deletions"=DEL_1e03_1e04_bp+DEL_1e04_1e05_bp+DEL_1e05_1e06_bp,
  "large deletions"=DEL_1e06_1e07_bp+DEL_1e07_Inf_bp,
  "small duplications"=DUP_0e00_1e03_bp,
  "medium duplications"=DUP_1e03_1e04_bp+DUP_1e04_1e05_bp+DUP_1e05_1e06_bp,
  "large duplications"=DUP_1e06_1e07_bp+DUP_1e07_Inf_bp,
  "small inversions"=INV_0e00_1e03_bp,
  "medium inversions"=INV_1e03_1e04_bp+INV_1e04_1e05_bp+INV_1e05_1e06_bp,
  "large inversions"=INV_1e06_1e07_bp+INV_1e07_Inf_bp,
  "translocations"=TRA
) %>% dplyr::select(sampleId,Cohort,Signature.1:Signature.30,'repeat deletion 1bp':'translocations') -> Pancreas_transformed

##remove sample T27 due to low tumor percentage
Pancreas_transformed <- Pancreas_transformed[which(!Pancreas_transformed$sampleId=="FR11123214"),]

data_plot <- Pancreas_transformed#Pancreas_transformed_relative//Pancreas_transformed
data_plot %>% #Pancreas_transformed_relative
  dplyr::select(Signature.1:'translocations') %>% 
  #dplyr::select(Signature.1:Signature.30) %>% 
  #scale %>% # Scale the data
  dist(method="euclidean") %>% # calculate a distance matrix, 
  hclust(method = "ward.D2") -> dista

# order samples according to clustering
sample_order = Pancreas$sampleId[dista$order]
dista$height=log10(dista$height)-2
dhc <- as.dendrogram(dista) 
plot(dhc)
P0 <- ggplot(dhc, horiz = F,labels = FALSE) 

Pancreas_transformed %>% dplyr::select(sampleId,Signature.1:Signature.30) -> Pancreas_SNVplot
Pancreas_SNVplot$sampleId = factor(Pancreas_SNVplot$sampleId, levels = sample_order)
Pancreas_SNVplot_rel <- sweep(Pancreas_SNVplot[,2:31],1,(rowSums(Pancreas_SNVplot[,2:31])),`/`)
Pancreas_SNVplot_rel$sampleId <- Pancreas_SNVplot$sampleId
Pancreas_SNVplot_rel_m <- melt(Pancreas_SNVplot_rel,id.vars = "sampleId")

Pancreas_transformed %>% dplyr::select(sampleId,'repeat deletion 1bp':'no context >2bp') -> Pancreas_INDELplot
Pancreas_INDELplot$sampleId = factor(Pancreas_INDELplot$sampleId, levels = sample_order)
Pancreas_INDELplot_rel <- sweep(Pancreas_INDELplot[,2:12],1,(rowSums(Pancreas_INDELplot[,2:12])),`/`)
Pancreas_INDELplot_rel$sampleId <- Pancreas_INDELplot$sampleId
Pancreas_INDELplot_rel_m <- melt(Pancreas_INDELplot_rel,id.vars = "sampleId")

Pancreas_transformed %>% dplyr::select(sampleId,'small deletions':'translocations') -> Pancreas_SVplot
Pancreas_SVplot$sampleId = factor(Pancreas_SVplot$sampleId, levels = sample_order)
Pancreas_SVplot_rel <- sweep(Pancreas_SVplot[,2:11],1,(rowSums(Pancreas_SVplot[,2:11])),`/`)
Pancreas_SVplot_rel$sampleId <- Pancreas_SVplot$sampleId
Pancreas_SVplot_rel_m <- melt(Pancreas_SVplot_rel,id.vars = "sampleId")

Pancreas_transformed_relative <- dplyr::inner_join(Pancreas_SNVplot_rel, Pancreas_INDELplot_rel, by = "sampleId")
Pancreas_transformed_relative <- dplyr::inner_join(Pancreas_transformed_relative, Pancreas_SVplot_rel, by = "sampleId")



data_plot %>% dplyr::select(sampleId,Signature.1:Signature.30) -> SNVplot
SNVplot <- sweep(SNVplot[,2:31],1,(rowSums(SNVplot[,2:31])),`/`)
SNVplot$sampleId <- data_plot$sampleId
SNVplot$sampleId = factor(SNVplot$sampleId, levels = sample_order)
SNVplot <- melt(SNVplot,id.vars = "sampleId")
P1 <- ggplot(data=SNVplot, aes(x = sampleId, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=colorpalette) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(), 
    strip.text.x = element_text(size = 9.2), 
    legend.position="bottom", legend.title = element_blank()) + 
  #facet_grid(~cancerType, scales = "free_x", labeller = labeller(cancerType = display_cancer_types)) + 
  guides(fill = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))


data_plot %>% dplyr::select(sampleId,'repeat deletion 1bp':'no context >2bp') -> INDELplot
INDELplot <- sweep(INDELplot[,2:12],1,(rowSums(INDELplot[,2:12])),`/`)
INDELplot$sampleId <- INDELplot$sampleId
INDELplot$sampleId = factor(data_plot$sampleId, levels = sample_order)
INDELplot <- melt(INDELplot,id.vars = "sampleId")
P2 <- ggplot(data=INDELplot, aes(x = sampleId, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=COLORS_indels) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(), 
    strip.text.x = element_text(size = 9.2), 
    legend.position="bottom", legend.title = element_blank()) + 
  #facet_grid(~cancerType, scales = "free_x", labeller = labeller(cancerType = display_cancer_types)) + 
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

data_plot %>% dplyr::select(sampleId,'small deletions':'translocations') -> SVplot
SVplot <- sweep(SVplot[,2:11],1,(rowSums(SVplot[,2:11])),`/`)
SVplot$sampleId <- SVplot$sampleId
SVplot$sampleId = factor(data_plot$sampleId, levels = sample_order)
SVplot <- melt(SVplot,id.vars = "sampleId")
P3 <- ggplot(data=SVplot, aes(x = sampleId, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=COLORS_indels) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(), 
    strip.text.x = element_text(size = 9.2), 
    legend.position="bottom", legend.title = element_blank()) + 
  #facet_grid(~cancerType, scales = "free_x", labeller = labeller(cancerType = display_cancer_types)) + 
  guides(fill = guide_legend(nrow = 4)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))


Pancreas %>% dplyr::select(sampleId,Cohort) %>%
  mutate(score = ifelse(Cohort == "Organoid", 1, 1)) -> Pancreas_cohort
Pancreas_cohort$sampleId = factor(Pancreas_cohort$sampleId, levels = sample_order)

P4 <- ggplot(data=Pancreas_cohort, aes(x = sampleId, y = score)) +
  geom_bar(aes(fill = Cohort), stat = "identity", width=1) + ylab("") +
  scale_fill_manual(values=COLORS_indels) +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.background = element_blank(), 
    strip.text.x = element_text(size = 9.2), 
    legend.position="bottom", legend.title = element_blank()) + 
  #facet_grid(~cancerType, scales = "free_x", labeller = labeller(cancerType = display_cancer_types)) + 
  guides(fill = guide_legend(nrow = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1))

pFigure1 = plot_grid(P0,P1, P2, P3, P4, ncol=1, align="v", rel_heights = c(2,4, 4, 4, 1),scale=c(1.095,1,1,1,1), labels = c("A", "B", "C", "D","C"))
#pFigure1
save_plot("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/Figure X - somaticMutoverview_absolute.png", pFigure1, base_width = 14, base_height = 20)

dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/"
plot_name="somaticMutoverview_absolute"
# draw it
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 14, height = 20) 
plot(pFigure1)
dev.off()





