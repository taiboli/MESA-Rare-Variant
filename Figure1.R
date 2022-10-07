library(tidyverse)
library(cowplot)
library(ggpubr)
library(elementalist)

data_dir = 'Data/'
out_dir = 'Results/'
type_cols = c('#89023E', '#818802', '#02884C', '#090288')
names(type_cols) = c('Methylation', 'Splicing', 'Expression', 'Protein')

# ### A - correlation of values over time - moved to supplement
# #cor_data = read.table(paste0(data_dir, 'corrected_value_correlations_across_exams_withSplicing.txt'), header=T)
# cor_data = read.table(paste0(data_dir, 'corrected_value_correlations_across_exams_geneLevel_withCpGsAndClusters.txt'), header=T)
# cor_data$Level = sapply(cor_data$Type, function(x) ifelse(x == 'CpG sites', 'CpG sites',
#                                                           ifelse(x == 'Splicing clusters', 'Splicing clusters', 'Gene level')))
# cor_data$Type = sapply(cor_data$Type, function(x) ifelse(x == 'CpG sites', 'Methylation',
#                                                          ifelse(x == 'Splicing clusters', 'Splicing', as.character(x))))
# cor_data$Type = factor(cor_data$Type, levels=c('Expression', 'Methylation', 'Splicing', 'Protein'))
# cor_data$Level = factor(cor_data$Level, levels=c('Gene level', 'CpG sites', 'Splicing clusters'))
# 
# cplot = ggplot(cor_data, aes(Gene_cor)) + geom_density(aes(fill=Type),alpha=0.5) +
#   theme_bw() + xlab('Correlation of measurements at the gene level across exams') +
#   scale_fill_manual(values=type_cols) +
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         legend.title=element_blank(),
#         legend.text=element_text(size=14),
#         legend.position='bottom') + facet_grid(.~Level, scales = 'free_x', space = "free_x") +
#   theme(strip.background = element_blank(), strip.text = element_text(size=12))
# 
# ggsave(cplot, file=paste0(out_dir, 'figureS1_correlations.pdf'), width=12, height=6, units="in")

### A - gene level replication of outliers across exams
rep_data = read_tsv(paste0(data_dir, 'outlier_replication_Z2_geneLevel_proportions_across_exams_withCpGsAndClusters.txt'))
rep_data = filter(rep_data, Type != 'CpG sites', Type != 'Splicing clusters')
rep_data$Type = factor(rep_data$Type, levels=c('Expression', 'Methylation', 'Splicing', 'Protein'))
rplot = ggplot(rep_data, aes(x=Threshold, y=Prop,Group=Type)) + geom_point(aes(color=Type),size=4) +
  geom_line(aes(color=Type)) + theme_bw() + xlab('Outlier threshold') +
  ylab('Proportion exam 1 gene-level\noutliers in exam 5') +
  scale_color_manual(values=type_cols) +
  facet_grid(.~Direction) + theme(strip.background = element_blank(),
                                  strip.text = element_text(size=14)) + 
  theme_pubclean() + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.title=element_text(size=14, face = "bold")) + 
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"), 
        strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold"),
        strip.placement = "inside", 
        legend.position = "bottom", 
        legend.title=element_blank(),
        legend.text=element_text(size=18, face = "bold"))

### B - proportion outlier events identified of all those tested 
prop_data = read_tsv(paste0(data_dir, 'proportion_outlier_events_across_thresholds_geneLevel_withCpGsAndClusters.txt'))
prop_data$Level = sapply(prop_data$Type, function(x) ifelse(x == 'CpG sites', 'CpG sites',
                                                            ifelse(x == 'Splicing clusters', 'Splicing\nclusters', 'Gene level')))
prop_data$Type = sapply(prop_data$Type, function(x) ifelse(x == 'CpG sites', 'Methylation',
                                                           ifelse(x == 'Splicing clusters', 'Splicing', x)))
prop_data$Type = factor(prop_data$Type, levels=c('Expression', 'Methylation', 'Splicing', 'Protein'))
prop_data$Level = factor(prop_data$Level, levels=c('Gene level', 'CpG sites', 'Splicing\nclusters'))
prop_data = prop_data %>% filter(Threshold == 3, Type == 'Splicing' | Direction != 'Both')
pplot = ggplot(prop_data, aes(x=Type, y=Prop, Group=Type)) + 
  geom_bar(stat='identity', aes(fill=Type, alpha=Direction)) +
  # theme_bw() + 
  xlab('') + guides(fill=F) +
  ylab('Proportion identified as gene-level outliers') +
  facet_grid(.~Level, scales = 'free_x', space = "free_x") +
  scale_alpha_manual(values=c(1,0.5,1)) +
  scale_fill_manual(values=type_cols) +
  theme_pubclean() +
  theme(axis.text.y=element_text(size=12),
        axis.text.x = element_text(size = 14),
        axis.title=element_text(size=14, face = "bold"),
        legend.title=element_blank(),
        legend.position=c(0.1,0.85),
        legend.text=element_text(size=14)) + 
  theme(#strip.background = element_rect(colour = "grey", fill = "white"), 
        strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold"),
        strip.placement = "inside")

### C - number of outliers identified per person
outlier_counts = read_tsv(paste0(data_dir, 'joint_outliers_per_individual_geneLevel_filteredMethyl_withCpGsAndClusters.txt'))
outlier_counts$Level = sapply(outlier_counts$Type, function(x) ifelse(x == 'CpG sites (603052)', 'CpG sites',
                                                                      ifelse(x == 'Splicing clusters (15407)', 'Splicing\nclusters', 'Gene level')))
outlier_counts$Level = factor(outlier_counts$Level, levels=c('Gene level', 'CpG sites', 'Splicing\nclusters'))
outlier_counts$Type = factor(outlier_counts$Type, levels=c("Expression (14290)", "Methylation (19919)", "Splicing (8211)", "Protein (1317)", "CpG sites (603052)", "Splicing clusters (15407)"))

outlier_counts$color = sapply(outlier_counts$Type, function(x) ifelse(grepl('CpG sites', x), 'Methylation',
                                                                      ifelse(grepl("Methylation", x), 'Methylation',
                                                                             ifelse(grepl("Splicing", x), "Splicing",
                                                                                    ifelse(grepl("Protein", x), "Protein", "Expression")))
                                                                      ))
outlier_counts$color = factor(outlier_counts$color, levels=c('Expression', 'Methylation', 'Splicing', 'Protein'))

outlier_counts$label = sapply(outlier_counts$Type, function(x) ifelse(grepl('CpG sites', x), 'CpG sites\n(603052)',
                                                                      ifelse(grepl("Methylation", x), 'Methylation\n(19919)',
                                                                             ifelse(grepl("Splicing clusters", x), "Splicing clusters\n(15407)",
                                                                                    ifelse(grepl("Splicing", x), "Splicing\n(8211)",
                                                                                           ifelse(grepl("Protein", x), "Protein\n(1317)", "Expression\n(14290)"))))
))

outlier_counts$label = factor(outlier_counts$label, levels=c("Expression\n(14290)", "Methylation\n(19919)", "Splicing\n(8211)", "Protein\n(1317)", "CpG sites\n(603052)", "Splicing clusters\n(15407)"))



options(scipen=10000)
oplot = ggplot(outlier_counts, aes(x=label, y=N+0.1)) + geom_boxplot(aes(fill = color)) + theme_bw() +
  scale_fill_manual(values=type_cols) +
  ylab('Number of joint outliers per individual') + xlab('') + scale_y_log10() +
  facet_grid(.~Level, scales = 'free_x', space = "free_x") +
  theme_pubclean() + 
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=12),
        axis.title=element_text(size=14, face = "bold")) + 
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"),
        strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold"),
        strip.placement = "inside", legend.position = "none")



first_row = plot_grid(rplot, labels=c('A'), nrow=1, label_size = 24)
second_row = plot_grid(oplot, pplot, labels=c('B', 'C'), ncol=2, align='hv', label_size = 24)
fig_time = plot_grid(first_row, second_row, nrow=2, align='hv')

ggsave(fig_time, file=paste0(out_dir, 'figure1.pdf'), width=18, height=11, units="in")

ggsave(fig_time, file=paste0(out_dir, 'figure1.eps'), width=18, height=11, units="in")
ggsave(fig_time, file=paste0(out_dir, 'Figure1.tiff'), width=18, height=11, units="in")
