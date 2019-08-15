# Re-making PCA and cluster for methylation


setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Methylation/Diff_Meth/graphs")
library(methylKit)
library(grid)
library(readr)
library(ggplot2)
library(ggdendro)
library(ggfortify)
library(ggrepel)

# Get data from methylkit

file.listA <- list("m08_subsetted_final.txt","m19_subsetted_final.txt",
                   "m23_subsetted_final.txt","m37_subsetted_final.txt",
                   "q08_subsetted_final.txt","q19_subsetted_final.txt",
                   "q23_subsetted_final.txt","q37_subsetted_final.txt",
                   "w08_subsetted_final.txt","w19_subsetted_final.txt",
                   "w23_subsetted_final.txt","w37_subsetted_final.txt")

raw_data <- methRead(file.listA,
                     sample.id = list("m08","m19","m23","m37","q08","q19", "q23", "q37",
                                      "w08","w19","w23","w37"),
                     treatment = c(2,2,2,2,1,1,1,1,0,0,0,0),
                     assembly="bter_1.0", 
                     context="CpG")

meth_all_data <- unite(raw_data)


# Make a new dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)
labs$group <- c(rep("Queen", 4), rep("Worker",4),rep("Male",4))

labs$label<- c("Q23","Q37","Q19","Q08",
               "W08","W37","W19","W23",
               "M19","M23","M08","M37")


ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, size=0.8),
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.01, colour=group, size=3, fontface="bold"),
            show.legend = F)+
  scale_color_manual(values=c("#FF6600","#9900CC","#00CC66"))+
 # scale_color_manual(values=c("red","blue"))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.title = element_blank(),
        axis.text = element_blank())


# Make a new PCA plot

PCA_data <- PCASamples(meth_all_data, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)
PCA_data1$status <- c(rep("Male",4), rep("Queen",4), rep("Worker",4))

percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste( colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=status))+
  geom_point(size=5)+
  geom_text_repel(aes(label=sample), size=6,show.legend=FALSE, point.padding = 0.5, box.padding = 0.25)+
  scale_color_manual(values=c("#FF6600","#9900CC","#00CC66"))+
  theme_bw()+
  xlab(paste0("PC1: ",percentage[1],"variance")) +
  ylab(paste0("PC2: ",percentage[2],"variance")) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank())
