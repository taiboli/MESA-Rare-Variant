library(tidyverse)
library(cowplot)
library(ggpubr)
library(elementalist)

data_dir = 'Data/'
out_dir = 'Results/'

## colors
type_cols = c('#89023E', '#818802', '#02884C', '#090288')
names(type_cols) = c('mOutliers', 'sOutliers', 'eOutliers', 'pOutliers')

### overall enrichments
exp_enrich_0.01 = read_tsv(paste0(data_dir, 'relativeRisk_eOutliers_peer30_freeze8_covs_age_ciseQTLs_globalRemoved_MAF0.01_updated.txt'))
prot_enrich_0.01 = read_tsv(paste0(data_dir, 'relativeRisk_pOutliers_peer30_freeze8_covs_age_cispQTLs0.25_globalRemoved_MAF0.01_ensemblIdFiltered.txt'))
splice_enrich_0.01 = read_tsv(paste0(data_dir, 'relativeRisk_sOutliers_globalRemoved_MAF0.01.txt'))
methyl_enrich_0.01 = read_tsv(paste0(data_dir, 'relativeRisk_collapsedGenes_PCandlincRNA_allSites_mOutliers_peer30_freeze8_covs_age_ciseQTLs_globalRemoved_MAF0.01_updated.txt'))
methyl_enrich_0.01$Category = sapply(methyl_enrich_0.01$Direction, function(x) ifelse(x == 'Hyper-methylated', 'Over',
                                                                                      ifelse(x == 'Hypo-methylated', 'Under', 'Combined')))
enrich_df = rbind(exp_enrich_0.01 %>% filter(ZT != 5) %>% mutate(Type = 'eOutliers'),
                  prot_enrich_0.01 %>% filter(ZT != 5) %>% mutate(Type = 'pOutliers'),
                  methyl_enrich_0.01 %>% filter(Threshold < 5) %>% mutate(ZT = Threshold, Type = 'mOutliers') %>% dplyr::select(Riskratio, Lower, Upper, Pval, ZT, Category, Type),
                  splice_enrich_0.01 %>% filter(ZT != 1) %>% filter(Controls == 'All') %>% dplyr::select(-Controls) %>% mutate(Type = 'sOutliers'))

enrich_df$ZT = factor(enrich_df$ZT, levels=c(2,3,4))
enrich_df$Type = factor(enrich_df$Type, levels=c('eOutliers', 'mOutliers', 'sOutliers', 'pOutliers'))

eplot = ggplot(enrich_df, aes(x=ZT, y=Riskratio, group=Type, color=Type)) + 
  geom_point(size=5) + geom_line() +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0) +
  facet_grid(.~Category, scales='free') +
  theme_bw() + geom_hline(yintercept=1, linetype='dashed') + xlab('|Z| Threshold') + ylab('Relative risk') +
  scale_color_manual(values=type_cols) + guides(color=F) +
  # theme_pubclean() + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, face = "bold")) + 
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"), 
    strip.background = element_rect_round(radius = unit(8, "pt"), colour = NA),
    strip.text = element_text(size=18, face = "bold"),
    strip.placement = "inside", 
    legend.position = "bottom", 
    legend.title=element_blank(),
    legend.text=element_text(size=18, face = "bold"))
  # theme(axis.text=element_text(size=15), axis.title=element_text(size=15),
  #       legend.title=element_blank(), legend.text=element_text(size=15),
  #       strip.background = element_blank(), legend.position='right',
  #       strip.text=element_text(size=15)) 






### annotations per category
all_single_rrs = read_tsv(paste0(data_dir, 'all_geneLevel_allSites_outlier_annotation_enrichments_MAF0-1_proteinExpUpdated.txt'))

protein_domain_single = read_tsv(paste0(data_dir, 'all_geneLevel_allSites_outlier_annotation_enrichments_MAF0-1_proteindomain.txt'))
all_single_rrs = rbind(all_single_rrs, protein_domain_single)



both_annot_rrs = filter(all_single_rrs, Direction == 'Both')
both_annot_rrs$Annotation = factor(both_annot_rrs$Annotation, levels=c('no_variant', 'non_coding', 'coding', 'conserved_noncoding', 'protein_domain', 'TSS', 'splice', 'stop', 'frameshift'))
annot_plotB = ggplot(both_annot_rrs %>% filter(!is.na(Upper), !is.na(Lower)), aes(x=Annotation, y=Riskratio, color=Type)) + 
  geom_point(size=4, position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, position=position_dodge(width = 0.5)) +
  xlab('') + ylab('Relative risk') +
  scale_color_manual(values=type_cols) +
  geom_hline(yintercept=1, linetype='dashed') + scale_y_log10() +
  theme_bw() +
  guides(color=F) +
  # theme_pubclean() + 
  theme(axis.text.x=element_text(size=15, hjust=1, angle=35),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=16, face = "bold"),
        title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position='bottom')

dir_annot_rrs = filter(all_single_rrs, Direction != 'Both')
dir_annot_rrs$Annotation = factor(dir_annot_rrs$Annotation, levels=c('no_variant', 'non_coding', 'coding', 'conserved_noncoding', 'protein_domain', 'TSS', 'splice', 'stop', 'frameshift'))

annot_plot = ggplot(dir_annot_rrs %>% filter(Type != 'mOutliers'), aes(x=Annotation, y=Riskratio, color=Type)) + 
  geom_point(size=4, position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, position=position_dodge(width = 0.5)) +
  xlab('') + ylab('Relative risk') +
  scale_color_manual(values=type_cols[c("eOutliers", "pOutliers")]) +
  geom_hline(yintercept=1, linetype='dashed') + scale_y_log10() +
  theme_bw() +
  theme(axis.text.x=element_text(size=14, hjust=1, angle=35),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=14),
        title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position='bottom') +
  facet_grid(.~Direction) + theme(strip.background = element_blank(), strip.text = element_text(size=16, face = "bold"))

ggsave(annot_plot, file='../Manuscript/Figures/figureS6_annotation_by_direction_updated.pdf', width=10, height=5, units="in")
ggsave(annot_plot, file='../Manuscript/Figures/figureS6_annotation_by_direction_updated.png', width=10, height=6.5, units="in")




### annotation enrichments per combination
all_combined_rrs = read_tsv(paste0(data_dir, 'all_geneLevel_allSites_overlap_outlier_annotation_enrichments_MAF0-1_proteinExpUpdated.txt'))
protein_domain_combined = read_tsv(paste0(data_dir, 'all_geneLevel_allSites_overlap_outlier_annotation_enrichments_MAF0-1_proteindomain.txt'))

all_combined_rrs = rbind(all_combined_rrs, protein_domain_combined)


ccols = c(type_cols, '#ff7f00')
names(ccols)[5] = 'Overlap'
all_combined_rrs$Annotation = factor(all_combined_rrs$Annotation, levels=c('no_variant', 'non_coding', 'coding', 'conserved_noncoding', 'protein_domain', 'TSS', 'splice', 'stop', 'frameshift'))
all_combined_rrs$Type = factor(all_combined_rrs$Type, levels=c('eOutliers', 'mOutliers', 'sOutliers', 'pOutliers', 'Overlap'))
all_combined_rrs$Group = factor(all_combined_rrs$Group, levels=c('eOutliers + mOutliers', 'eOutliers + sOutliers', 'eOutliers + pOutliers'))


overlap_plot = ggplot(all_combined_rrs %>% filter(!is.na(Upper), !is.na(Lower)), aes(x=Annotation, y=Riskratio, color = Type)) + 
  geom_point(size=4, position=position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0, position=position_dodge(width = 0.8)) + 
  xlab('') + ylab('Relative risk') + 
  scale_color_manual(values=ccols) + 
  geom_hline(yintercept=1, linetype='dashed') + scale_y_log10() +
  facet_grid(.~Group, scales='free') +
  # theme_pubclean() +
  theme_bw() +
  theme(#strip.background = element_rect(colour = NULL, fill = "grey"), 
    strip.background = element_rect_round(radius = unit(8, "pt"), colour = NA),
    strip.text = element_text(size=18, face = "bold")) +
  theme(axis.text=element_text(size=15),
        axis.text.x = element_text(hjust=1, angle=45),
        axis.title=element_text(size=16, face= "bold"),
        title=element_text(size=14),
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position='bottom')


### combine plots
first_row = plot_grid(eplot, annot_plotB, labels=c('A', 'B'), nrow=1, align='h', label_size = 24, rel_widths = c(1.2, 1), 
                      label_x = c(0, -.05), scale = 0.95)
second_row = plot_grid(overlap_plot, labels=c('C'), label_size = 24)
combined_plot = plot_grid(first_row, second_row, nrow=2, align='hv', scale=0.95, rel_heights = c(1,1.5))

# library(patchwork)
# combined_plot = (eplot + annot_plotB) / overlap_plot + plot_annotation(tag_levels = 'A')

ggsave(combined_plot, file=paste0(out_dir, '/figure3.eps'), width=12, height=12, units="in")
ggsave(combined_plot, file=paste0(out_dir, '/Figure3.tiff'), width=12, height=12, units="in")
