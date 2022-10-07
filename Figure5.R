library(tidyverse)
library(dplyr)
library(ggpubr)
library(elementalist)
library(data.table)
library(cowplot)
library(grid)
library(gridExtra)
library(rtracklayer)
library(ggvenn)



type_cols = c('#02884C','#89023E', '#818802',  '#090288')
names(type_cols) = c('RNA', 'Methylation', 'Splicing', 'Protein')


data_dir <- "Data/"
out_dir <- "Results/"


inputGenes <- readRDS(paste0(data_dir, "GenesDirectlyMeasured.RDS"))

absmax <- function(x) { x[which.max( abs(x) )]}
perc.rank <- function(x) trunc(rank(x))/length(x)





#################################################
## Accounting - Fig 5a
## 30.5M gene-RV-ind triplets 
#################################################

# rv.scores.summarized <- readRDS("MultiomicWatershed/Results/CombinedPilot/MultiomicWatershedUpdated-Prediction-AllData_IndividualGenotypeSummarized-filteredbygnomAD.RDS")
# 
# 
# posterior.count <- rv.scores.summarized %>%
#   group_by(Ind) %>%
#   summarise(RNAcount = sum(abs(RNA) >= 0.5),
#             Proteincount = sum(abs(Protein) >= 0.5),
#             RNAcountHigh = sum(abs(RNA) >= 0.9),
#             ProteincountHigh = sum(abs(Protein) >= 0.9),
#             MethylationCount = sum(abs(Methylation) >= 0.5),
#             MethylationCountHigh = sum(abs(Methylation) >= 0.9),
#             SplicingCount = sum(abs(Splicing) >= 0.5),
#             SplicingCountHigh = sum(abs(Splicing) >= 0.9))
# 
# plot.data <- posterior.count %>%
#   pivot_longer(!Ind, names_to = "Type", values_to = "Count")
# plot.data$Cutoff <- ifelse(grepl("High", plot.data$Type), 0.9, 0.5)
# plot.data$Signal <- ifelse(grepl("RNA", plot.data$Type), "RNA", 
#                            ifelse(grepl("Protein", plot.data$Type), "Protein",
#                                   ifelse(grepl("Methylation", plot.data$Type), "Methylation", "Splicing")))
# 
# plot.data$Signal <- factor(plot.data$Signal, levels = c("RNA", "Methylation", "Splicing", "Protein"))
# 
# plot.data$logCount <- log10(plot.data$Count + 1)
# 
# 
# ### remove gloal outliers 
# 
# rna.globaloutliers <- scan("MultiomicWatershed/GlobalOutliers/RNA.txt", what = character())
# protein.globaloutliers <- scan("MultiomicWatershed/GlobalOutliers/Protein.txt", what = character())
# splicing.globaloutliers <- scan("MultiomicWatershed/GlobalOutliers/Splicing.txt", what = character())
# methyl.globaloutliers <- scan("MultiomicWatershed/GlobalOutliers/Methylation.txt", what = character())
# 
# 
# 
# to.remove <- which(plot.data$Ind %in% rna.globaloutliers & plot.data$Signal == "RNA")
# plot.data <- plot.data[-to.remove,]
# 
# to.remove <- which(plot.data$Ind %in% protein.globaloutliers & plot.data$Signal == "Protein")
# plot.data <- plot.data[-to.remove,]
# 
# to.remove <- which(plot.data$Ind %in% methyl.globaloutliers & plot.data$Signal == "Methylation")
# plot.data <- plot.data[-to.remove,]
# 
# to.remove <- which(plot.data$Ind %in% splicing.globaloutliers & plot.data$Signal == "Splicing")
# plot.data <- plot.data[-to.remove,]
# 
# plot.data$countplus1 <- plot.data$Count + 1
# 
# # plot.data$Signal <- factor(plot.data, levels = c("RNA", "Methylation", "Splicing", "Protein"))


plot.data <- readRDS(paste0(data_dir, "FigureData-Figure5a.RDS"))


cutoff.labels <- c("Watershed posterior \u2265 0.5", "Watershed posterior \u2265 0.9")
names(cutoff.labels) <- c("0.5", "0.9")


p5a <- ggplot(plot.data, aes(Signal, countplus1)) + 
  geom_boxplot(aes(fill = Signal), alpha = 0.8) + 
  scale_fill_manual(values = type_cols) + 
  ylab('Number of rare variants per individual') + xlab('') + scale_y_log10() +
  
  facet_grid(.~Cutoff, scales = 'free_x', space = "free_x", labeller = labeller(Cutoff = cutoff.labels)) + 
  
  theme_pubclean() +
  theme(axis.title = element_text(face="bold", colour="black", size=16),
        axis.text  = element_text(size=14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_blank()) +
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"),
    strip.background = element_rect_round(radius = unit(8, "pt")),
    strip.text = element_text(size=16, face = "bold"),
    strip.placement = "inside")




# 
# p <- ggboxplot(plot.data, x = "Cutoff", y = "countplus1",
#                color = "Signal",
#                # palette = "jco",
#                xlab = "Watershed Cutoff", ylab = "RVs per individual") +
#   scale_color_manual(values = c('#02884C', '#89023E', '#818802', '#090288')) +
#   yscale("log10", .format = F) +
#   theme_pubclean() +
#   theme(axis.title = element_text(face="bold", colour="black", size=13),
#         axis.text  = element_text(size=12),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank())

# 
# 
# ggsave(p, filename = "MultiomicWatershed/FinalResults/Fig5a-CountWatershedByCutoff-ByRV-FilterGlobalOutliers.pdf",
#        width = 7, height = 4)


# plot.data %>% 
#   group_by(Type) %>%
#   summarise(median = median(Count),
#             high = quantile(Count, probs = 0.75),
#             low = quantile(Count, probs = 0.25))
# 
# 
# # Type                 median  high   low
# # <chr>                 <dbl> <dbl> <dbl>
# #   1 MethylationCount       65.5   114    47
# # 2 MethylationCountHigh   17      26    10
# # 3 Proteincount           52      91    37
# # 4 ProteincountHigh        1       2     0
# # 5 RNAcount               11      16     6
# # 6 RNAcountHigh            0       1     0
# # 7 SplicingCount           7      14     3
# # 8 SplicingCountHigh       0       1     0
# 
# ##### Only with omics measurements 
# 
# mesa.inds <- fread("Zscores/mesa.1319_samples.covariates.txt")
# mesa.inds <- mesa.inds$ID
# 
# tmp <- subset(plot.data, Ind %in% mesa.inds)
# 
# tmp %>% 
#   group_by(Type) %>%
#   summarise(median = median(Count),
#             high = quantile(Count, probs = 0.75),
#             low = quantile(Count, probs = 0.25))
# 
# # Type                 median  high   low
# # <chr>                 <dbl> <dbl> <dbl>
# #   1 MethylationCount         84   131    58
# # 2 MethylationCountHigh     19    29    13
# # 3 Proteincount             65   107    46
# # 4 ProteincountHigh          2     3     1
# # 5 RNAcount                 14    19     9
# # 6 RNAcountHigh              1     2     0
# # 7 SplicingCount             9    16     5
# # 8 SplicingCountHigh         0     1     0
# p <- ggboxplot(tmp, x = "Cutoff", y = "countplus1",
#                color = "Signal",
#                # palette = "jco",
#                xlab = "Watershed Cutoff", ylab = "RVs per individual") +
#   scale_color_manual(values = c('#02884C', '#89023E', '#818802', '#090288')) +
#   yscale("log10", .format = F) +
#   theme_pubclean() +
#   theme(axis.title = element_text(face="bold", colour="black", size=13),
#         axis.text  = element_text(size=12),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank())
# 
# 
# ggsave(p, filename = "MultiomicWatershed/FinalResults/Supp_Fig5a-CountWatershedByCutoff-DirectMeasuredIndsOnly-ByRV-FilterGlobalOutliers.pdf",
#        width = 7, height = 4)
# 


################################################
### Distribution of effect sizes in height
################################################

# watershed.summarized <- readRDS(paste0(data_dir, "MultiomicWatershedUpdated-Prediction-AllData_summarizedPosterior-FilteredBygnomAD.RDS"))
# 
# ## Gencode V35
# gencode <- readRDS(paste0(data_dir, "gencodeProcessed.RDS"))
# gencode$Gene <- substr(gencode$gene_id, 1, 15)
# abnormal.height.genes <- scan(paste0(data_dir, "abnormalHeightGenes.txt"),
#                               what = character())
# abnormal.height.genes.ensembl <- gencode$Gene[match(abnormal.height.genes,
#                                                     gencode$gene_name)]
# 
# ch = import.chain(paste0(data_dir, "hg19ToHg38.over.chain"))
# 
# convertNealeGWASToHg38Neale <- function(thisFile){
#   
#   this.gwas <- fread(thisFile)
#   this.gwas$Chr <- unlist(lapply(this.gwas$variant, function(x) strsplit(x, ":")[[1]][[1]]))
#   this.gwas$Pos <- unlist(lapply(this.gwas$variant, function(x) strsplit(x, ":")[[1]][[2]]))
#   this.granges <- makeGRangesFromDataFrame(df = this.gwas, seqnames.field = "Chr",
#                                            start.field = "Pos", end.field = "Pos",
#                                            ignore.strand = T, keep.extra.columns = T)
#   
#   seqlevelsStyle(this.granges) = "UCSC"  # necessary
#   this.gwas.hg38 = liftOver(this.granges, ch)
#   this.gwas.hg38 <- unlist(this.gwas.hg38)
#   
#   
#   this.gwas.hg38 <- as.data.frame(this.gwas.hg38)
#   setDT(this.gwas.hg38)
#   setnames(this.gwas.hg38, "seqnames", "Chr")
#   # this.gwas.hg38$Chr <- as.character(this.gwas.hg38$seqnames)
#   this.gwas.hg38$RV <- paste0(this.gwas.hg38$Chr, ":", this.gwas.hg38$start)
#   return(this.gwas.hg38)
#   
# }
# 
# # thisFile <- "C:/Users/litb9/Downloads/harmonized_imputed_gwas/NealeLab/50_raw.gwas.imputed_v3.both_sexes.tsv.gz"
# 
# thisFile <- paste0(data_dir, "50_irnt.gwas.imputed_v3.both_sexes.tsv.gz")
# this.gwas.hg38.height <- convertNealeGWASToHg38Neale(thisFile)
# 
# 
# transformNealeGWAS <- function(this.gwas){
#   this.gwas <- merge.data.table(this.gwas, watershed.summarized[,c("RV", "Gene", "geneRV", "RNA", "Protein", "Methylation", "Splicing")],
#                                 by = c("RV"), all = F)
#   this.gwas <- within(this.gwas, absbetaR <- perc.rank(abs(beta)))
#   this.gwas$zscore <- this.gwas$beta / this.gwas$se
#   this.gwas <- within(this.gwas, absZR <- perc.rank(abs(zscore)))
# }
# 
# this.gwas.hg38.height <- transformNealeGWAS(this.gwas.hg38.height)
# 
# 
# 
# thresholds <- c(0.5, 0.9)
# signals <- c("RNA", "Methylation", "Splicing", "Protein")
# 
# 
# tmp <- subset(this.gwas.hg38.height, Gene %in% abnormal.height.genes.ensembl)
# 
# plot.data.distribution <- data.frame(beta = tmp$absbetaR,
#                                      Group = paste0("Background\n", nrow(tmp)))
# 
# 
# setDT(tmp)
# 
# 
# for(this.s in signals){
#   for(this.tr in thresholds){
#     tmp2 <- subset(tmp, abs(get(this.s)) >= this.tr )
#     # tmp2 <- subset(tmp2, Gene %in% inputGenes[[this.s]])
#     if(nrow(tmp2) == 0) next
#     plot.data.distribution <- rbind(plot.data.distribution,
#                                     data.frame(beta = tmp2$absbetaR,
#                                                Group = paste0("|", this.s, "| >= ", this.tr, "\n", nrow(tmp2))))
#   }
# }
# 
# 
# 
# plot.data.distribution$color <- ifelse(grepl("RNA", plot.data.distribution$Group), "RNA",
#                                        ifelse(grepl("Protein", plot.data.distribution$Group), "Protein", 
#                                               ifelse(grepl("Methylation", plot.data.distribution$Group), "Methylation",
#                                                      ifelse(grepl("Splicing", plot.data.distribution$Group), "Splicing", "Background"))))
# 
# plot.data.distribution$color <- factor(plot.data.distribution$color, levels = c("Background", "RNA", "Methylation", "Splicing", "Protein"))


plot.data.distribution <- readRDS(paste0(data_dir, "FigureData-Figure5b.RDS"))
p5b <- ggboxplot(plot.data.distribution, x = "Group", y = "beta",
                 # label = "Count",
                 # color = "color", 
                 color = NA,
                 fill = NA,
                 # fill = "color",
                 add = "median",
                 ylab = "Height Effect Size Percentile") +
  geom_boxplot(aes(fill = color), alpha = 0.8) +
  ylim(c(0, 1)) +
  # stat_summary(fun.data=flabel, geom="text", vjust=-0.5, col="blue") +
  # scale_color_manual(values = c('#999999', 
  # '#02884C', '#89023E', '#818802', '#090288')) +
  scale_fill_manual(values = c('#999999',
                               '#02884C', '#89023E', '#818802', '#090288')) +
  theme_pubclean() +
  theme(axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=15),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  rotate_x_text(45)




#########################
### Fig 5c Venn Diagram
#########################
# 
# tmp <- subset(this.gwas.hg38.height, Gene %in% abnormal.height.genes.ensembl)
# 
# # set1 <- tmp$RV[which(abs(tmp$RNA) >= 0.5 & tmp$Gene %in% inputGenes[["RNA"]])]
# # set2 <- tmp$RV[which(abs(tmp$Methylation) >= 0.5 & tmp$Gene %in% inputGenes[["Methylation"]])]
# # set3 <- tmp$RV[which(abs(tmp$Splicing) >= 0.5 & tmp$Gene %in% inputGenes[["Splicing"]])]
# # set4 <- tmp$RV[which(abs(tmp$Protein) >= 0.5 & tmp$Gene %in% inputGenes[["Protein"]])]
# 
# set1 <- tmp$RV[which(abs(tmp$RNA) >= 0.5)]
# set2 <- tmp$RV[which(abs(tmp$Methylation) >= 0.5 )]
# set3 <- tmp$RV[which(abs(tmp$Splicing) >= 0.5 )]
# set4 <- tmp$RV[which(abs(tmp$Protein) >= 0.5 )]
# 
# set1 <- set1[!is.na(set1)]
# set2 <- set2[!is.na(set2)]
# set3 <- set3[!is.na(set3)]
# set4 <- set4[!is.na(set4)]
# 
# x <- list(RNA = set1, Methylation = set2, Splicing = set3, Protein = set4)

# names(x) <- c("Stage 1","Stage 2","Stage 3", "Stage4")

x <- readRDS(paste0(data_dir, "FigureData-Figure5c.RDS"))
p5c <- ggvenn(
  x, 
  fill_color = c('#02884C', '#89023E', '#818802', '#090288'),
  fill_alpha = 0.5,
  stroke_size = 0.5, stroke_alpha = 0.5,
  stroke_color = NA,
  set_name_size = 8,
  text_size = 6
)


#####################################
#### Fig 5d Multi-modal Effects
#####################################


# tmp <- subset(this.gwas.hg38.height, Gene %in% abnormal.height.genes.ensembl)
# # tmp <- within(tmp, absbetaR <- perc.rank(abs(beta)))
# setDT(tmp)
# 
# 
# this.threshold <- 0.5
# signals <- c("RNA", "Methylation", "Splicing", "Protein")
# 
# 
# ### Splicing multimodal
# this.signal <- "Splicing"
# plot.data <- data.frame(beta = tmp$absbetaR,
#                         Group = paste0("Background\n", nrow(tmp)),
#                         color = "Background")
# tmp2 <- subset(tmp, abs(get(this.signal)) >= this.threshold )
# # tmp2 <- subset(tmp2, Gene %in% inputGenes[[this.signal]])
# plot.data <- rbind(plot.data,
#                    data.frame(beta = tmp2$absbetaR,
#                               Group = paste0(this.signal, "\n", nrow(tmp2)),
#                               color = this.signal))
# 
# for(this.s in setdiff(signals, this.signal)){
#   tmp2 <- subset(tmp, abs(get(this.s)) >= this.threshold & abs(get(this.signal)) >= this.threshold )
#   # tmp2 <- subset(tmp2, Gene %in% inputGenes[[this.s]])
#   plot.data <- rbind(plot.data,
#                      data.frame(beta = tmp2$absbetaR,
#                                 Group = paste0(this.signal, " & ", this.s, "\n", nrow(tmp2)),
#                                 color = this.s))
# }
# 
# 
# plot.data$color <- factor(plot.data$color, levels = c("Background", "RNA", "Methylation", "Splicing", "Protein"))
# 
# plot.data$Signal <- "Splicing"
# 
# plot.data.splicing <- plot.data
# 
# 
# ### Protein multimodal
# this.signal <- "Protein"
# plot.data <- data.frame(beta = tmp$absbetaR,
#                         Group = paste0("Background\n", nrow(tmp)),
#                         color = "Background")
# tmp2 <- subset(tmp, abs(get(this.signal)) >= this.threshold )
# # tmp2 <- subset(tmp2, Gene %in% inputGenes[[this.signal]])
# plot.data <- rbind(plot.data,
#                    data.frame(beta = tmp2$absbetaR,
#                               Group = paste0(this.signal, "\n", nrow(tmp2)),
#                               color = this.signal))
# 
# for(this.s in setdiff(signals, this.signal)){
#   tmp2 <- subset(tmp, abs(get(this.s)) >= this.threshold & abs(get(this.signal)) >= this.threshold )
#   # tmp2 <- subset(tmp2, Gene %in% inputGenes[[this.s]])
#   plot.data <- rbind(plot.data,
#                      data.frame(beta = tmp2$absbetaR,
#                                 Group = paste0(this.signal, " & ", this.s, "\n", nrow(tmp2)),
#                                 color = this.s))
# }
# 
# 
# plot.data$color <- factor(plot.data$color, levels = c("Background", "RNA", "Methylation", "Splicing", "Protein"))
# 
# plot.data$Signal <- "Protein"
# 
# 
# plot.data <- rbind(plot.data.splicing, plot.data)
# 
# plot.data$Signal <- factor(plot.data$Signal, levels = c("Splicing", "Protein"))
# 
# 
# signal.labels <- c("Splicing + other signals", "Protein + other signals")
# names(signal.labels) <- c("Splicing", "Protein")
# 

plot.data <- readRDS(paste0(data_dir, "FigureData-Figure5d.RDS"))
p5d <- ggboxplot(plot.data, x = "Group", y = "beta",
               # label = "Count",
               # color = "color", 
               fill = NA,
               color = NA,
               add = "median",
               ylab = "Height Effect Size Percentile") +
  geom_boxplot(aes(fill = color), alpha = 0.8) +
  facet_grid(.~Signal, scales = 'free_x', space = "free_x", labeller = labeller(Signal = signal.labels)) + 
  ylim(c(0, 1)) +
  # stat_summary(fun.data=flabel, geom="text", vjust=-0.5, col="blue") +
  # scale_color_manual(values = c('#999999', 
  #                               '#02884C', '#89023E', '#818802', '#090288')) +
  scale_fill_manual(values = c('#999999',
                               '#02884C', '#89023E', '#818802', '#090288')) +
  theme_pubclean() +
  theme(axis.title.y = element_text(face="bold", colour="black", size=15),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=15),
        axis.text.y = element_text(size = 12),
        legend.position = "none") +
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"),
    strip.background = element_rect_round(radius = unit(8, "pt")),
    strip.text = element_text(size=16, face = "bold"),
    strip.placement = "inside") + 
  rotate_x_text(45)






############################# 
## Fig 5e Height comparison
## Gene based distribution 
#############################




chosen.genes <- c("FBN1", "CHSY1", "SAMHD1", "PDE11A")
# chosen.genes.ensembl <- gencode$Gene[match(chosen.genes, gencode$gene_name)]

signals <- c("RNA", "Methylation", "Splicing", "Protein")
type_cols = c('#02884C','#89023E', '#818802',  '#090288', '#999999')
names(type_cols) = c('RNA', 'Methylation', 'Splicing',  'Protein', 'Control')


data.summary <- readRDS(paste0(data_dir, "HeightAssociationResults.RDS"))

pList <- list()

plot.data <- readRDS(paste0(data_dir, "FigureData-Figure5e.RDS"))


for(this.index in 1:4){
  top.gene <- chosen.genes.ensembl[this.index]
  top.gene.name <- chosen.genes[this.index]

  top.gene.data <- subset(plot.data, title == chosen.genes[this.index])
  top.gene.data$absCombined <- abs(top.gene.data[[signals[this.index]]])
  top.gene.data$color <- ifelse(top.gene.data$category == "Control", "Control", signals[this.index])
  top.gene.data$title <- top.gene.name
  
  
  labels <- symnum(c(data.summary[[paste0(signals[this.index], "p")]][which(data.summary$Gene==top.gene)]),
                   corr = FALSE, cutpoints = c(0, .001, .01, .05, 1), symbols = c("***","**","*","NS"))
  p <- ggplot(top.gene.data, aes(x=category, y=absCombined, fill=color), alpha = 0.8) +
    # scale_fill_manual(values = get_palette(palette = "jco", k = 2)) +
    scale_fill_manual(values = type_cols) + 
    geom_violin(scale = 'width') +
    geom_point(position = position_jitterdodge(jitter.width = 0.85)) +
    facet_grid(. ~ title) + 
    ylim(c(min(top.gene.data$absCombined),max(top.gene.data$absCombined)*1.25)) +
    scale_x_discrete(labels=c("Control" = "Normal", "Outlier" = "Outlier")) +
    geom_signif(xmin = c("Control"),
                xmax = c("Outlier"),
                y_position = c(max(top.gene.data$absCombined)*1.1),
                annotations = labels,
                textsize = 6
    ) +
    # ggtitle(top.gene.name) +
    xlab("Height Category") +
    ylab(paste0(signals[this.index], " Posterior")) +
    theme_pubclean() +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_blank()) +
    theme(#strip.background = element_rect(colour = NULL, fill = "grey"),
      strip.background = element_rect_round(radius = unit(8, "pt")),
      strip.text = element_text(size=16, face = "bold"),
      strip.placement = "inside")
  
  pList[[length(pList) + 1]] <- p
  
  
}


p5e <- plot_grid(pList[[1]], pList[[2]], pList[[3]], pList[[4]], nrow = 1, align = "h")

### Trick to add common x axis title
x.grob <- textGrob("Height Category", 
                   gp=gpar(fontface="bold", fontsize=16))

#add to plot
p5e <- grid.arrange(arrangeGrob(p5e, bottom = x.grob))



### combine plots , label_x = c(0, -.05), scale = 0.95
first_row = plot_grid(p5a, p5b, labels=c('A', 'B'), nrow=1, align='h', label_size = 24, rel_widths = c(1, 1.5),
                      label_x = c(0, .05),
                      scale = 0.95)
second_row = plot_grid(p5c, p5d, labels=c('C', 'D'), label_size = 24, align = "h", rel_widths = c(1, 1.6),
                       label_x = c(0, -.05), scale = 0.95)
third_row = plot_grid(p5e, labels = c("E"), label_size = 24, scale = 0.9)
combined_plot = plot_grid(first_row, second_row, third_row, nrow=3, align='hv', scale=0.95, rel_heights = c(1, 1, 1.2))


ggsave(combined_plot, file=paste0(out_dir, 'Figure5.pdf'), width=18, height=22, units="in", device=cairo_pdf, dpi = 600)
ggsave(combined_plot, file=paste0(out_dir, 'Figure5.eps'), width=18, height=22, units="in")
ggsave(combined_plot, file=paste0(out_dir, 'Figure5.tiff'), width=18, height=22, units="in")


