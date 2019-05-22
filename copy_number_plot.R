library("foreach")
library("reshape2")
library("ggplot2")
library("gridExtra")
library(grid)


copy_numbers <- read.csv("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/Results_CopyNumber.tsv", sep="\t")
colnames(copy_numbers) <- do.call(rbind, strsplit(colnames(copy_numbers), ".purple"))[,1]
Pancreas_organoid_lines <- c("FR11124502","FR11124096","FR11124060","FR11123764","FR11123208","FR11123647","FR11123965","FR11123814","FR11123204","FR11123568","FR11123640","FR11123200","FR11123782","FR11124713","FR11123939","FR11123104","FR11124064","FR11123256","FR11123775","FR11124045","FR11123963","FR11123524","FR11123171","FR11124161","FR11123214","FR11124015")
length(Pancreas_organoid_lines)
copy_numbers <- cbind(copy_numbers$chromArms,copy_numbers[which(names(copy_numbers) %in%Pancreas_organoid_lines) ])
names(copy_numbers)[1] <- c("chromArms")

write.table(copy_numbers, file = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/copy_numbers_Pancreas.tsv", sep="\t", col.names = TRUE, row.names = FALSE)
cn <- foreach(i = 2:ncol(copy_numbers), .combine = "cbind") %do% {
  x <- copy_numbers[,i] - copy_numbers[nrow(copy_numbers),i]
  return(data.frame(x))
}
colnames(cn) <- colnames(copy_numbers)[2:ncol(copy_numbers)]
rownames(cn) <- copy_numbers$chromArms

aneuploidy <- foreach(i=1:ncol(cn), .combine = "rbind") %do% {
  return(data.frame("Sample"=colnames(cn)[i], "AS"=length(which(cn[,i]!=0))))
}


#### ANEUPLOIDY PLOT HMF 

p_aneuploidy <- ggplot(aneuploidy, aes(x=AS)) + 
  geom_bar() +
  scale_x_continuous(breaks=seq(0,40,5)) +
  labs(x="Aneuploidy Score", y="Number of samples") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        text = element_text(size=20))



copy_numbers_b <- copy_numbers[,Pancreas_organoid_lines]
copy_numbers_b <- foreach(i = 1:ncol(copy_numbers_b), .combine = "cbind") %do% {
  x <- copy_numbers_b[,i] - copy_numbers_b[nrow(copy_numbers_b),i]
  return(data.frame(x))
}
colnames(copy_numbers_b) <- colnames(copy_numbers[,Pancreas_organoid_lines])
rownames(copy_numbers_b) <- copy_numbers$chromArms
  
copy_numbers_b_melted <- melt(cbind(sample=rownames(copy_numbers_b), copy_numbers_b))
copy_numbers_b_melted$sample <- factor(copy_numbers_b_melted$sample, levels = unique(copy_numbers_b_melted$sample))

P1 <- ggplot(copy_numbers_b_melted[copy_numbers_b_melted$sample!="genome_level",], aes(x=0, y=value)) + 
  geom_point(position = position_jitter(w = 0.05, h = 0),aes(color=variable)) + 
  facet_wrap(~sample, nrow=1, strip.position = "bottom") + 
  labs(x="Chromosome arm", y="Copy number change") +
  scale_y_continuous(breaks=c(seq(-5,5,1),10,15)) +
  theme(panel.spacing.x=unit(0.01, "lines"), axis.text.x = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA, color = "grey", linetype = "dotted"), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=10))

dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/"
plot_name="CopyNumber_gains_losses"
# draw it
pdf(sprintf("%s%s_B.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 3) 
grid.draw(P1)
dev.off()

P2 <- ggplot(copy_numbers_b_melted[copy_numbers_b_melted$sample!="genome_level",], aes(x=0, y=value)) + 
  geom_count() +
  scale_size_area()+
  facet_wrap(~sample, nrow=1, strip.position = "bottom") + 
  labs(x="Chromosome arm", y="Copy number change") +
  scale_y_continuous(breaks=c(seq(-5,5,1),10,15)) +
  theme(panel.spacing.x=unit(0.01, "lines"), axis.text.x = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA, color = "grey", linetype = "dotted"), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=10))


dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/"
plot_name="CopyNumber_gains_losses"
# draw it
pdf(sprintf("%s%s_A.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 3) 
grid.draw(P2)
dev.off()

P3 <- ggplot(copy_numbers_b_melted[copy_numbers_b_melted$sample!="genome_level",], aes(x=0, value)) + 
  geom_violin(fill="darkolivegreen3",colour = "darkolivegreen4", draw_quantiles = c(0.25, 0.5, 0.75)) + #, draw_quantiles = c(0.25, 0.5, 0.75)
  facet_wrap(~sample, nrow=1, strip.position = "bottom") + #, scale = "area", size = 0.1
  #scale_fill_manual(values=cancerTypeColours, guide=FALSE) +
  #scale_colour_manual(values=cancerTypeColours, guide=FALSE) +
  labs(x="Chromosome arm", y="Copy number change") +
  scale_y_continuous(breaks=c(seq(-5,5,1),10,15)) +
  theme(panel.spacing.x=unit(0.01, "lines"),axis.text.x = element_blank(), axis.ticks = element_blank(),  #panel.spacing.x=unit(0.01, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = NA, color = "grey", linetype = "blank"), 
        axis.line = element_line(colour = "black"),
        text = element_text(size=10))

dirpath="~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/"
plot_name="CopyNumber_gains_losses"
# draw it
pdf(sprintf("%s%s_C_1.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 3) 
grid.draw(P3)
dev.off()
