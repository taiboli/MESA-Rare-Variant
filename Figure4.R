library(tidyverse)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(data.table)
library(elementalist)
library(pheatmap)
library(gridExtra)


data_dir <- "Data/"
output_dir <- "Results/"

eval.data.watershed <- readRDS(paste0(data_dir, "MultiomicCombinedUpdated_Watershed_exact_N2pair_0.02_pvalthresh_0.05_evaluation_object.rds"))
# eval.data.rna.reduced <- readRDS("MultiomicWatershed/Results/Splicing-ReducedKeepN2pair_RIVER_N2pair_0.03_pvalthresh_0.05_evaluation_object.rds")
# eval.data.watershed <- readRDS("MultiomicWatershed/Results/CombinedPilot/Updated/MultiomicCombinedUpdated_Watershed_approximate_N2pair_0.02_pvalthresh_0.05_evaluation_object.rds")


plot.data <- data.frame(precision = eval.data.watershed$auc[[1]]$evaROC$GAM_precision,
                        recall = as.numeric(eval.data.watershed$auc[[1]]$evaROC$GAM_recall),
                        type = "RNA GAM")
tmp <- data.frame(precision = eval.data.watershed$auc[[1]]$evaROC$watershed_precision,
                  recall = as.numeric(eval.data.watershed$auc[[1]]$evaROC$watershed_recall),
                  type = "RNA Watershed")
plot.data <- rbind(plot.data, tmp)


tmp <- data.frame(precision = eval.data.watershed$auc[[2]]$evaROC$GAM_precision,
                  recall = as.numeric(eval.data.watershed$auc[[2]]$evaROC$GAM_recall),
                  type = "Protein GAM")
plot.data <- rbind(plot.data, tmp)

tmp <- data.frame(precision = eval.data.watershed$auc[[2]]$evaROC$watershed_precision,
                  recall = as.numeric(eval.data.watershed$auc[[2]]$evaROC$watershed_recall),
                  type = "Protein Watershed")
plot.data <- rbind(plot.data, tmp)


tmp <- data.frame(precision = eval.data.watershed$auc[[3]]$evaROC$GAM_precision,
                  recall = as.numeric(eval.data.watershed$auc[[3]]$evaROC$GAM_recall),
                  type = "Methylation GAM")
plot.data <- rbind(plot.data, tmp)

tmp <- data.frame(precision = eval.data.watershed$auc[[3]]$evaROC$watershed_precision,
                  recall = as.numeric(eval.data.watershed$auc[[3]]$evaROC$watershed_recall),
                  type = "Methylation Watershed")
plot.data <- rbind(plot.data, tmp)


tmp <- data.frame(precision = eval.data.watershed$auc[[4]]$evaROC$GAM_precision,
                  recall = as.numeric(eval.data.watershed$auc[[4]]$evaROC$GAM_recall),
                  type = "Splicing GAM")
plot.data <- rbind(plot.data, tmp)

tmp <- data.frame(precision = eval.data.watershed$auc[[4]]$evaROC$watershed_precision,
                  recall = as.numeric(eval.data.watershed$auc[[4]]$evaROC$watershed_recall),
                  type = "Splicing Watershed")
plot.data <- rbind(plot.data, tmp)





plot.data$algorithm <- ifelse(grepl("GAM", plot.data$type), "GAM", "Watershed")
plot.data$algorithm <- factor(plot.data$algorithm, levels = c("Watershed", "GAM"))

plot.data$color <- unlist(lapply(plot.data$type, function(x) strsplit(x, " ")[[1]][[1]]))

plot.data$color <- factor(plot.data$color, levels = c("RNA", "Methylation", "Splicing", "Protein"))




p4a <- ggplot(plot.data, aes(x = recall, y = precision, color = color, linetype = algorithm)) +
  geom_line(size = 1.1) +
  
  scale_color_manual(values = c('#02884C', '#89023E', '#818802', '#090288')) +
  scale_linetype_manual(values=c("solid", "dotted"))+
  xlab("Recall") + ylab("Precision") +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=16),
        axis.text  = element_text(size=12),
        panel.grid = element_blank(),
        # axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16), 
        legend.position = c(0.8, 0.75)) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.8, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect_round(radius = unit(8, "pt"))) 


# ggsave(p4a, filename = "MultiomicWatershed/FinalResults/Fig4a-MultiomicCombinedAUPRC_Updated.tiff", width = 9, height = 6)
# ggsave(p4a, filename = "MultiomicWatershed/FinalResults/Fig4a-MultiomicCombinedAUPRC.pdf", width = 8, height = 6)






# 
# 
# 
# 
auprc <- matrix(data = NA, ncol = 4, nrow = 2)
signal <- c("watershed", "GAM")
rownames(auprc) <- signal
colnames(auprc) <- c("RNA", "Protein", "Methylation", "Splicing")
for(i in 1:4){
  for(j in 1:2){
    auprc[j,i] <- eval.data.watershed$auc[[i]]$evaROC[[paste0(signal[[j]], "_pr_auc")]]
  }
}



auprc <- auprc[,c("RNA", "Methylation", "Splicing", "Protein")]
# 
# # RNA    Protein Methylation   Splicing
# # watershed 0.07896181 0.06848336   0.1121677 0.07489005
# # GAM       0.02649098 0.02680941   0.0606422 0.02740583
# 
# # pauprc <- tableGrob(auprc)
# 
# 
# # pauprc <- as.table(t(auprc))
# 
pauprc <- as.data.table(t(auprc))
setnames(pauprc, "watershed", "Watershed")
rownames(pauprc) <- colnames(auprc)
chosencols <- colnames(pauprc)
pauprc[, (chosencols) := lapply(.SD, function(x){sprintf(x, fmt = '%#.3f')}), .SDcols = chosencols]


# pauprc <- as.table(pauprc)


cols <- matrix(rep(c('#02884C', '#89023E', '#818802', '#090288'), 2), ncol=2, byrow=F)


tt <- ttheme_minimal(base_size = 14,
                     core=list(fg_params = list(col=cols),
                               bg_params = list(col=NA)),
                     # rowhead=list(bg_params = list(col=NA)),
                     # colhead=list(bg_params = list(col=NA)),
                     colhead = list(fg_params=list(cex = 1.0)),
                     rowhead = list(fg_params=list(cex = 1.0, fontface=3)))

mytable = tableGrob(pauprc, theme = tt)


p4a <- ggplot(plot.data, aes(x = recall, y = precision, color = color, linetype = algorithm)) +
  geom_line(size = 1.1) +
  annotation_custom(mytable, xmin = 0.25, xmax = 0.5, ymin = 0.5, ymax = 1.1) +
  scale_color_manual(values = c('#02884C', '#89023E', '#818802', '#090288')) +
  scale_linetype_manual(values=c("solid", "dotted"))+
  xlab("Recall") + ylab("Precision") +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=16),
        axis.text  = element_text(size=12),
        panel.grid = element_blank(),
        # axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.8, 0.75)) +
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.8, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_rect_round(radius = unit(8, "pt")))
ggsave(p4a, filename = paste0(output_dir, "Fig4a-MultiomicCombinedAUPRC_Updated.tiff"), width = 9, height = 6)




################################
## theta pair
################################


m <- matrix(NA, ncol = 4, nrow = 4)
rownames(m) <- c("RNA", "Protein", "Methylation", "Splicing")
colnames(m) <- c("RNA", "Protein", "Methylation", "Splicing")
upperm = lower.tri(m, diag = F)
m[upperm] = as.vector(eval.data.watershed$model_params$theta_pair)

m[upper.tri(m)] = t(m)[upper.tri(m)]

m <- m[c("RNA", "Methylation", "Splicing", "Protein"), c("RNA", "Methylation", "Splicing", "Protein")]

# m[is.na(m)] <- ""

m2 <- m
m2 <- matrix(sprintf(m2, fmt = "%#.3f"), ncol = 4, nrow = 4)
m2[which(m2 == "NA")] <- ""


p <- pheatmap(m,
              # annotation_col = patient.annotation,
              cluster_rows = F,
              cluster_cols = F,
              show_colnames = T,
              show_rownames = T,
              border_color = F,
              angle_col = 45,
              fontsize = 13,
              # display_numbers = T,
              # number_format = "%.3f",
              display_numbers = m2,
              fontsize_col = 16,
              fontsize_row = 16,
              fontsize_number = 14,
              legend = F)

ggsave(p, filename = paste0(output_dir, "/Fig4b-Theta_pair_updated.tiff"), width = 4.5, height = 4.5)
# ggsave(p, filename = "MultiomicWatershed/FinalResults//Fig4b-Theta_pair.pdf", width = 4.5, height = 3)







################################
## theta for features
################################

# all.annotations <- fread("MultiomicWatershed/CombinedData/WatershedCombinedUpdated-MethylDownsampled-Training.txt.gz",
#                          nrows = 2)
# cat(names(all.annotations)[3:75], file = "MultiomicWatershed/FinalScripts/Github/MESA-Rare-Variant/Data/Watershed_annotations.txt", sep = "\n")
 

annotation.names <- scan(paste0(data_dir, "/Watershed_annotations.txt"), what = character())


theta.data <- data.frame(eval.data.watershed$model_params$theta)
names(theta.data)[1:4] <- c("RNA", "Protein", "Methylation", "Splicing")
theta.data <- abs(theta.data)

theta.data$feature <- annotation.names



tmp <- apply(theta.data[,1:4], 1, max)

theta.data <- theta.data[which(abs(tmp) >= 0.05),]

theta.data$group <- "Group"


p1 <- ggplot(theta.data, 
             aes(x = RNA, y = reorder(feature, RNA))) +
  ylab("") + xlab("RNA Weights") +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_manual(values = c('#02884C')) +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=13),
        axis.text.x  = element_text(size=12),
        legend.position = "none",
        axis.text.y=  element_text(face=c(rep("plain", nrow(theta.data) - 5), rep("bold", 5)),
                                   size = 12))

p2 <- ggplot(theta.data, 
             aes(x = Protein, y = reorder(feature, Protein))) +
  ylab("") + xlab("Protein Weights") +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_manual(values = c('#090288')) +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=13),
        axis.text.x  = element_text(size=12),
        legend.position = "none",
        axis.text.y=  element_text(face=c(rep("plain", nrow(theta.data) - 5), rep("bold", 5)),
                                   size = 12))

p3 <- ggplot(theta.data, 
             aes(x = Methylation, y = reorder(feature, Methylation))) +
  ylab("") + xlab("Methylation Weights") +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_manual(values = c('#89023E')) +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=13),
        axis.text.x  = element_text(size=12),
        legend.position = "none",
        axis.text.y=  element_text(face=c(rep("plain", nrow(theta.data) - 5), rep("bold", 5)),
                                   size = 12))

p4 <- ggplot(theta.data, 
             aes(x = Splicing, y = reorder(feature, Splicing))) +
  ylab("") + xlab("Splicing Weights") +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_manual(values = c('#818802')) +
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=13),
        axis.text.x  = element_text(size=12),
        legend.position = "none",
        axis.text.y=  element_text(face=c(rep("plain", nrow(theta.data) - 5), rep("bold", 5)),
                                   size = 12))



p <- ggarrange(p1, p3, p4, p2,  ncol = 4, nrow = 1)



ggsave(p, filename = paste0(output_dir, "Fig4c-FeatureWeights_updated.tiff"), width = 19, height = 7)
# ggsave(p, filename = "MultiomicWatershed/FinalResults/Fig4c-FeatureWeights.pdf", width = 19, height = 8)






























