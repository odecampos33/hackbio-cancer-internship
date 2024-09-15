
library(gplots)
library(tidyverse)
library(RColorBrewer)


# load dataset
url <- "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv"

gene_data <- read_csv(url)

head(gene_data)


# set genes as rownames
gene_data <- gene_data %>% column_to_rownames(var = '...1')


## column names of variables in gene data
names(gene_data)


## Visualise expression levels using heatmap and showing clusters

# Create heatmap
heatmap.2(as.matrix(gene_data), 
          #          col = diverging_palette, 
          trace = "none", 
          margins = c(8, 8), 
          cexRow = 0.8, 
          cexCol = 0.8,
          main = "Heatmap of Scaled Data")


# Scale the data
scaled_data <- scale(gene_data)

# Define color palette for heatmap
diverging_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
#heat map across both genes and samples 

# Create heatmap both samples and genes
heatmap.2(as.matrix(scaled_data), 
          col = diverging_palette, 
          trace = "none", 
          dendrogram = "both",  # Cluster both rows and columns
          scale = "row",        # Scale data by row
          margins = c(10, 10),  # Margins around the plot
          cexRow = 0.5,         # Size of row labels
          cexCol = 0.5,         # Size of column labels
          main = "Heatmap of Glioblastoma Data")


dev.off()

# Create heatmap across samples only 
heatmap.2(as.matrix(gene_data), 
          col = diverging_palette, 
          Colv=TRUE,
          Rowv=FALSE,
          trace = "none", 
          dendrogram = "col",  # Cluster columns
          scale = "row",        # Scale data by row
          margins = c(10, 10),  # Margins around the plot
          cexRow = 0.5,         # Size of row labels
          cexCol = 0.5,         # Size of column labels
          main = "Heatmap of Glioblastoma Data")


# Create heatmap across genes only 
heatmap.2(as.matrix(gene_data), 
          col = diverging_palette, 
          Colv=FALSE,
          Rowv=TRUE,
          trace = "none", 
          dendrogram = "row",  # Cluster rows
          scale = "row",        # Scale data by row
          margins = c(10, 10),  # Margins around the plot
          cexRow = 0.5,         # Size of row labels
          cexCol = 0.5,         # Size of column labels
          main = "Heatmap of Glioblastoma Data")

# cluster columns and rows (samples and genes) - sequential
heatmap.2(as.matrix(gene_data), Colv = T, Rowv = T, dendrogram = 'both', 
          col = RColorBrewer::brewer.pal(9, 'Blues'), 
          trace='none', scale='row',
          main="Clustering both Genes and Samples")


# creating metadata
# all samples with 2005 grouped in same group

metadata <- data.frame(
  row.names = colnames(gene_data), 
  groups=ifelse(grepl('2005', names(gene_data)), 'group1', 'group2' )
)

metadata$groups <- relevel(factor(metadata$groups), ref='group1')

metadata



## get group datasets
group1 <- gene_data[, which(metadata$groups == 'group1')]
group2 <- gene_data[, which(metadata$groups == 'group2')]


## function to get pvalues from t-test statistical test
run_differential_genes <- function(row){
  t.test(row[names(group1)], row[names(group2)])$p.value
}


# calculate pvalues, log fold change and mean difference
pvalues <- apply(cbind(group1, group2), 1, run_differential_genes)
logfoldchange <- log2(rowMeans(group1)) - log2(rowMeans(group2))
mean_diff <- rowMeans(group1) - rowMeans(group2)


## join mean_diff, fold change and pvalues of genes as dataframe
results <- as.data.frame(cbind(mean_diff, logfoldchange, pvalues))


## get differentially expressed genes
# selecting upregulated and downregulated genes
upregulated <- results %>% 
  filter(pvalues < 0.05, logfoldchange >= 1.5)

downregulated <- results %>% 
  filter(pvalues < 0.05, logfoldchange <= 1.5)


## show them
View(upregulated)
View(downregulated)


## write to file for gene enrichment analysis
downregulated %>% write.csv('downregulated_genes_new.csv', row.names = T)
upregulated %>% write.csv('upregulated_genes_new.csv', row.names = T)


## Visualisation of results

# volcano plots
results %>%
  # get negative log10 values
  mutate(neg_log_pval = -log10(pvalues)) %>%
  # group into up, down regulated and insignificant genes
  mutate(group = case_when(
    rownames(gene_data) %in% rownames(upregulated) ~ 'up',
    rownames(gene_data) %in% rownames(downregulated) ~ 'down',
    .default = 'insig')) %>%
  # visualise with ggplot2
  ggplot(aes(logfoldchange, neg_log_pval, color=group)) +
  geom_point(alpha=0.7) +
  theme_bw() +
  theme(legend.key = element_blank(), 
        legend.position = 'inside',
        legend.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        panel.grid = element_blank(),
        legend.position.inside = c(0.89, 0.83)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual(values=c('seagreen', 'black', 'coral'), 
                     labels=c('Down', 'Insig', 'Up')) +
  labs(title='Volcano Plot', x=expression('Log'[2]*'FC'), 
       y=expression('-Log'[10]*' pvalue'), color='Regulation Type') +
  scale_x_continuous(breaks=seq(-12.5,12.5,2.5), expand = c(0.1,0.1,0.05,0.1))



## Lollipop chart to show top biological pathways in upregulated and downregulated genes

# load the files
BP_upregulated <- read_csv(
  'Upregulated enrichment biological processes.csv', 
  show_col_types = FALSE)

BP_downregulated <- read_csv(
  'Downregulated enrichment biological processes.csv', 
  show_col_types = FALSE)


## visualise

## for downregulated genes
BP_downregulated %>%
  # filter by statistically significant result
  filter(`Enrichment FDR` < 0.05) %>%
  slice_min(`Enrichment FDR`, n=5) %>%
  mutate(`Enrichment FDR` = -log10(`Enrichment FDR`)) %>%
  mutate(Pathway = str_to_title(str_trim(str_remove(Pathway, 'GO:\\d+\\s')))) %>%
  ggplot(aes(x=`Enrichment FDR`, y=fct_reorder(Pathway, `Enrichment FDR`))) +
  geom_segment(aes(x=0, xend=`Enrichment FDR`, y=Pathway, yend=Pathway, 
                   color=`Fold Enrichment`), linewidth=1) +
  geom_point(aes(size=nGenes), color='darkgray') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        plot.title=element_text(face='bold'),
        axis.text.y = element_text(size=9),
        legend.title = element_text(size=9, face='bold'),
        panel.background = element_rect(fill='white')) +
  scale_color_gradient(low='blue', high='red') +
  scale_size(range = c(4,6)) +
  guides(color=guide_legend(title='Fold Enrichment'),scale='none') +
  labs(title='Biological Processes Pathway', 
       subtitle = 'Downregulated genes', y='', x=expression("-Log"[10]*" (pvalue)"))


# for upregulated genes
BP_upregulated %>%
  # filter by statistically significant result
  filter(`Enrichment FDR` < 0.05) %>%
  slice_min(`Enrichment FDR`, n=5) %>% 
  mutate(`Enrichment FDR` = -log10(`Enrichment FDR`)) %>%
  mutate(Pathway = str_to_title(str_trim(str_remove(Pathway, 'GO:\\d+\\s')))) %>%
  ggplot(aes(x=`Enrichment FDR`, y=fct_reorder(Pathway, `Enrichment FDR`))) +
  geom_segment(aes(x=0, xend=`Enrichment FDR`, y=Pathway, yend=Pathway, 
                   color=`Fold Enrichment`), linewidth =1) +
  geom_point(aes(size=nGenes), color='darkgray') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        plot.title=element_text(face='bold'),
        axis.text.y = element_text(size=9),
        legend.title = element_text(size=9, face='bold'),
        panel.background = element_rect(fill='white')) +
  scale_color_gradient(low='blue', high='red') +
  scale_size(range = c(4,5)) +
  guides(color=guide_legend(title='Fold Enrichment'),scale='none') +
  labs(title='Biological Processes Pathway', 
       subtitle = 'Upregulated genes', y='', x=expression("-Log"[10]*" (pvalue)"))

Collapse



