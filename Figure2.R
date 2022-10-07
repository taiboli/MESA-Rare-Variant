library(data.table)
library(tidyverse)
library(ggridges)
library(cowplot)
library(ggpubr)


data_dir = 'Data/'
type_cols = c('#89023E', '#818802', '#02884C', '#090288', 'grey')
names(type_cols) = c('mOutlier', 'sOutlier', 'eOutlier', 'pOutlier', 'non-outlier')

out_dir = "Results/"

outlier_zt = 3
matched_controls = FALSE

##### A - mOutlier overlap
### Includes methylation to expression, methylation to splicing, and methylation to protein
if (isTRUE(matched_controls)) {
  moutlier_exp = fread(paste0(data_dir, 'expression_zscores_matchedControl_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt'))
  moutlier_splice = fread(paste0(data_dir, 'splicing_zscores_matchedControls_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt')) 
  moutlier_prot = fread(paste0(data_dir, 'protein_zscores_matchedControls_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt')) 
} else {
  moutlier_exp = fread(paste0(data_dir, 'expression_zscores_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt')) 
  moutlier_splice = fread(paste0(data_dir, 'splicing_zscores_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt')) 
  moutlier_prot = fread(paste0(data_dir, 'protein_zscores_moutlier_collapsedMedianGene_allSites_Z', outlier_zt, '.txt')) 
}

moutlier_combined = rbind(moutlier_exp %>% select(Category, MeanZ) %>% mutate(Type = 'Expression'),
                          moutlier_splice %>% select(Category, MeanZ) %>% mutate(Type = 'Splicing'),
                          moutlier_prot %>% select(Category, MeanZ) %>% mutate(Type = 'Protein'))
moutlier_combined$Category = factor(moutlier_combined$Category, levels=c('non-outlier', 'mOutlier'))
moutlier_combined$Type = factor(moutlier_combined$Type, levels=c('Protein', 'Splicing', 'Expression'))
moutlier_combined$title = "mOutliers"


mo_plot = ggplot(moutlier_combined, aes(x = abs(MeanZ), y = Type)) + 
  geom_density_ridges(aes(fill=Category), alpha=0.6, scale=0.9) +
  facet_grid(.~title) +
  # scale_fill_manual(values=c('grey', type_cols[1])) +
  scale_fill_manual(values = type_cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .6))) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
  xlab('|Mean Z-score across exams|') + ylab('') +
  xlim(c(0,4)) + #ggtitle('mOutliers') +
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title=element_text(size=16, face = "bold"),
        plot.title=element_text(size=20, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position='none') + 
  theme(strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold", lineheight = 1.5),
        strip.placement = "inside")

wilcox.test(abs(filter(moutlier_exp, Category == 'mOutlier')$MeanZ), abs(filter(moutlier_exp, Category != 'mOutlier')$MeanZ), alternative='g')
wilcox.test(abs(filter(moutlier_prot, Category == 'mOutlier')$MeanZ), abs(filter(moutlier_prot, Category != 'mOutlier')$MeanZ), alternative='g') 
wilcox.test(abs(filter(moutlier_splice, Category == 'mOutlier')$MeanZ), abs(filter(moutlier_splice, Category != 'mOutlier')$MeanZ), alternative='g') 

##### B - sOutlier overlap 
### Includes splicing to expression, splicing to protein, splicing to methylation
if (isTRUE(matched_controls)) {
  soutlier_exp = fread(paste0(data_dir, 'expression_zscores_matchedControls_soutliers_Z', outlier_zt, '.txt'))
  soutlier_prot = fread(paste0(data_dir, 'protein_zscores_matchedControls_soutliers_Z', outlier_zt, '.txt')) 
  soutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_matchedControls_soutliers_Z', outlier_zt, '.txt')) 
} else {
  soutlier_exp = fread(paste0(data_dir, 'expression_zscores_soutliers_Z', outlier_zt, '.txt')) 
  soutlier_prot = fread(paste0(data_dir, 'protein_zscores_soutliers_Z', outlier_zt, '.txt')) 
  soutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_soutliers_Z', outlier_zt, '.txt'))
}

soutlier_combined = rbind(soutlier_exp %>% select(Category, MeanZ) %>% mutate(Type = 'Expression'),
                          soutlier_methyl %>% select(Category, MeanZ) %>% mutate(Type = 'Methylation'),
                          soutlier_prot %>% select(Category, MeanZ) %>% mutate(Type = 'Protein'))
soutlier_combined$Category = factor(soutlier_combined$Category, levels=c('non-outlier', 'sOutlier'))
soutlier_combined$Type = factor(soutlier_combined$Type, levels=c('Protein', 'Methylation', 'Expression'))
soutlier_combined$title = "sOutliers"

so_plot = ggplot(soutlier_combined, aes(x = abs(MeanZ), y = Type)) + 
  geom_density_ridges(aes(fill=Category), alpha=0.6, scale=0.9) +
  facet_grid(.~title) +
  scale_fill_manual(values=type_cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0, .5))) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
  xlab('|Mean Z-score across exams|') + ylab('') +
  xlim(c(0,4)) + #ggtitle('sOutliers') +
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title=element_text(size=16, face = "bold"),
        plot.title=element_text(size=20, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position='none') + 
  theme(strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold", lineheight = 1.5),
        strip.placement = "inside")

wilcox.test(abs(filter(soutlier_exp, Category == 'sOutlier')$MeanZ), abs(filter(soutlier_exp, Category != 'sOutlier')$MeanZ), alternative='g') 
wilcox.test(abs(filter(soutlier_prot, Category == 'sOutlier')$MeanZ), abs(filter(soutlier_prot, Category != 'sOutlier')$MeanZ), alternative='g') 
wilcox.test(abs(filter(soutlier_methyl, Category == 'sOutlier')$MeanZ), abs(filter(soutlier_methyl, Category != 'sOutlier')$MeanZ), alternative='g')

##### C - eOutlier overlap
### Includes expression to splicing, expression to protein
if (isTRUE(matched_controls)) {
  eoutlier_splice = fread(paste0(data_dir, 'splicing_zscores_matchedControls_eoutliers_Z', outlier_zt, '.txt'))
  eoutlier_prot = fread(paste0(data_dir, 'protein_zscores_matchedControls_eoutliers_Z', outlier_zt, '.txt')) 
  eoutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_matchedControls_eoutliers_Z', outlier_zt, '.txt')) 
} else {
  eoutlier_splice = fread(paste0(data_dir, 'splicing_zscores_eoutliers_Z', outlier_zt, '.txt'))
  eoutlier_prot = fread(paste0(data_dir, 'protein_zscores_eoutliers_Z', outlier_zt, '.txt')) 
  eoutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_eoutliers_Z', outlier_zt, '.txt'))
}

eoutlier_combined = rbind(eoutlier_splice %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Splicing'),
                          eoutlier_methyl %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Methylation'),
                          eoutlier_prot %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Protein'))
eoutlier_combined$Category = factor(eoutlier_combined$Category, levels=c('non-outlier', 'eOutlier'))
eoutlier_combined$Type = factor(eoutlier_combined$Type, levels=c('Protein', 'Splicing', 'Methylation'))
eoutlier_combined$Title = "eOutliers"

eoutlier_combined$label <- ifelse(eoutlier_combined$Direction == "Over", "Over-eOutliers", "Under-eOutliers")



eo_plot = ggplot(eoutlier_combined, aes(x = abs(MeanZ), y = Type)) + 
  geom_density_ridges(aes(fill=Category), alpha=0.6, scale=1) +
  facet_grid(.~label) +
  scale_fill_manual(values=type_cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .6))) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + xlab('|Mean Z-score across exams|') + ylab('') + 
  xlim(c(0,4)) + #ggtitle('eOutliers') +
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title=element_text(size=16, face = "bold"),
        plot.title=element_text(size=20, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position='none') + 
  theme(strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold", lineheight = 1.5),
        strip.placement = "inside")

dcurrent = 'Over'
wilcox.test(abs(filter(eoutlier_splice, Category == 'eOutlier', Direction == dcurrent)$MeanZ), abs(filter(eoutlier_splice, Category != 'eOutlier', Direction == dcurrent)$MeanZ), alternative='g') 
wilcox.test(abs(filter(eoutlier_prot, Category == 'eOutlier', Direction == dcurrent)$MeanZ), abs(filter(eoutlier_prot, Category != 'eOutlier', Direction == dcurrent)$MeanZ), alternative='g')
wilcox.test(abs(filter(eoutlier_methyl, Category == 'eOutlier', Direction == dcurrent)$MeanZ), abs(filter(eoutlier_methyl, Category != 'eOutlier', Direction == dcurrent)$MeanZ), alternative='g') 

##### D - pOutlier overlap
### Includes protein to expression, protein to splicing 
if (isTRUE(matched_controls)) {
  poutlier_splice = fread(paste0(data_dir, 'splicing_zscores_matchedControls_poutliers_Z', outlier_zt, '.txt'))
  poutlier_exp = fread(paste0(data_dir, 'expression_zscores_matchedControls_poutliers_Z', outlier_zt, '.txt')) 
  poutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_matchedControls_poutliers_Z', outlier_zt, '.txt')) 
} else {
  poutlier_splice = fread(paste0(data_dir, 'splicing_zscores_poutliers_Z', outlier_zt, '.txt')) 
  poutlier_exp = fread(paste0(data_dir, 'expression_zscores_poutliers_Z', outlier_zt, '.txt')) 
  poutlier_methyl = fread(paste0(data_dir, 'moutlier_collapsedMedianGene_allSites_zscores_poutliers_Z', outlier_zt, '.txt')) 
}

poutlier_combined = rbind(poutlier_splice %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Splicing'),
                          poutlier_methyl %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Methylation'),
                          poutlier_exp %>% select(Category, MeanZ, Direction) %>% mutate(Type = 'Expression'))
poutlier_combined$Category = factor(poutlier_combined$Category, levels=c('non-outlier', 'pOutlier'))
poutlier_combined$Type = factor(poutlier_combined$Type, levels=c('Expression', 'Methylation', 'Splicing'))

pLabeler <- c("Over-pOutliers", "Under-pOutliers")
names(pLabeler) <- c("Over", "Under")

poutlier_combined$label <- ifelse(poutlier_combined$Direction == "Over", "Over-pOutliers", "Under-pOutliers")

po_plot = ggplot(poutlier_combined, aes(x = abs(MeanZ), y = Type)) + 
  geom_density_ridges(aes(fill=Category), alpha=0.6, scale=0.9) +
  facet_grid(.~label) +
  scale_fill_manual(values=type_cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = expand_scale(mult = c(0, 0.5))) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + xlab('|Mean Z-score across exams|') + ylab('') + 
  xlim(c(0,4)) + #ggtitle('eOutliers') +
  theme(axis.text.y=element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title=element_text(size=16, face = "bold"),
        plot.title=element_text(size=20, hjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position='none') + 
  theme(strip.background = element_rect_round(radius = unit(8, "pt")),
        strip.text = element_text(size=16, face = "bold", lineheight = 1.5),
        strip.placement = "inside")



dcurrent = 'Over'
wilcox.test(abs(filter(poutlier_splice, Category == 'pOutlier', Direction == dcurrent)$MeanZ), abs(filter(poutlier_splice, Category != 'pOutlier', Direction == dcurrent)$MeanZ), alternative='g') 
wilcox.test(abs(filter(poutlier_exp, Category == 'pOutlier', Direction == dcurrent)$MeanZ), abs(filter(poutlier_exp, Category != 'pOutlier', Direction == dcurrent)$MeanZ), alternative='g') 
wilcox.test(abs(filter(poutlier_methyl, Category == 'pOutlier', Direction == dcurrent)$MeanZ), abs(filter(poutlier_methyl, Category != 'pOutlier', Direction == dcurrent)$MeanZ), alternative='g')

first_row = plot_grid(eo_plot, mo_plot,  labels=c('A', 'B'), nrow=1, align = "h", axis = "bt", rel_widths = c(2, 1), label_size = 24)
second_row = plot_grid(so_plot, po_plot, labels = c("C", "D"), nrow = 1, align = "h", axis = "bt", rel_widths = c(1,2), label_size = 24)

combined_plot = plot_grid(first_row, second_row, nrow = 2, align = "hv")

# grobs <- ggplotGrob(so_plot)$grobs
legend <- get_legend(
  so_plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.justification = "center")
)

combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))


ggsave(combined_plot, file=paste0(out_dir, 'Figure2.pdf'), width=16, height=10, units="in", device=cairo_pdf, dpi = 600)
ggsave(combined_plot, file=paste0(out_dir, 'Figure2.eps'), width=16, height=10, units="in")
ggsave(combined_plot, file=paste0(out_dir, 'Figure2.tiff'), width=16, height=10, units="in")

