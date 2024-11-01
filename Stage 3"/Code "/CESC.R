## Load necessary libraries

# Install the packages if you don't have them
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
# BiocManager::install(c('EDASeq'))


# Load the libraries
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(tidyverse))


## Get TCGA projects
# GDCprojects <- getGDCprojects()
# getProjectSummary("TCGA-CESC")

# query database
query_TCGA <- GDCquery(project='TCGA-CESC',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       data.type = 'Gene Expression Quantification')

Output_query <- getResults(query_TCGA)


## download queried data
GDCdownload(query_TCGA)


## Prepare downloaded data
data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)


## Extract necessary samples
# Extracting tumor samples (Primary solid tumor)
tumor_samples <- TCGAquery_SampleTypes(data$barcode, typesample = 'TP')

# Extract normal samples (Solid Tissue Normal)
normal_samples <- TCGAquery_SampleTypes(data$barcode, typesample = 'NT')



# create clinical data
clinical_data <- data.frame(row.names = data$barcode, 
                            patientID = data$patient, 
                            age = data$age_at_index, 
                            race=data$race)


tumor_ids <- rownames(clinical_data) %in% tumor_samples # get tumour samples

## Create age groups to sample tumor samples by age group

# create age groups
age_groups <- cut(clinical_data$age[tumor_ids], 
                  breaks = c(0, 30, 40, 50, 60, 70, 100),
                  labels = c('20-30', '31-40', '41-50', '51-60', '61-70', '>70')
                  )
# add weights to sample observations from each group
probs <- as.vector((table(age_groups)))

# add age groups and weights to tumour info data frame

set.seed(302) # for reproducibility
tumor_info <- clinical_data |>
  filter(tumor_ids) |>
  dplyr::mutate(age_groups = age_groups, 
                prob = dplyr::case_when(age_groups == '20-30' ~ probs[1],
                                        age_groups == '31-40' ~ probs[3],
                                        age_groups == '41-50' ~ probs[3],
                                        age_groups == '51-60' ~ probs[4],
                                        age_groups == '61-70' ~ probs[5],
                                        age_groups == '>70' ~ probs[6])) |>
  #dplyr::group_by(age_groups) 
  dplyr::slice_sample(n=37, replace = F, weight_by = prob)


## Number of primary samples 37; Number of normal samples

# update tumor samples with randomly selected samples
tumor_samples <- rownames(tumor_info)


## Data Preprocessing
# preprocess data
data_preprocessed <- TCGAanalyze_Preprocessing(data, cor.cut = 0.4, 
                                               datatype = "unstranded",
                                               filename = "CERC_preprocess.png")


## normalise data by gene length and GC content
# normalising by gene length and GC content
data_preprocessed <- TCGAanalyze_Normalization(data_preprocessed, 
                                               geneInfo = geneInfoHT, 
                                               method='geneLength')


data_preprocessed <- TCGAanalyze_Normalization(data_preprocessed, 
                                               geneInfo = geneInfoHT, 
                                               method='gcContent')


## Filtering
# filter zero counts using quantile method (< 25% percentile)
data_preprocessed <- TCGAanalyze_Filtering(data_preprocessed, 
                                           method='quantile',
                                           qnt.cut = 0.25)


## Differential Gene Expression Analysis

# Performing differential expression analysis
dea_results <- TCGAanalyze_DEA(
  mat1 = data_preprocessed[, normal_samples],   
  mat2 = data_preprocessed[, tumor_samples],  
  Cond1type = "Normal",
  Cond2type = "Tumor",
  method = "glmLRT"                          
)


# View top differentially expressed genes
head(arrange(dea_results, FDR))

# Saving the results
write.csv(dea_results, "DEA_results_Tumor_vs_Normal.csv")


# get upregulated and downregulated genes
upregulated <- dea_results |>
  filter(FDR < 0.01, logFC > 5)

downregulated <- dea_results |>
  filter(FDR < 0.01, logFC < -5)


downregulated %>% write.csv('downregulated_genes.csv', row.names = T)
upregulated %>% write.csv('upregulated_genes.csv', row.names = T)


## Data Visualisation
# volcano plots
dea_results %>%
  mutate(neg_log_pval = -log10(FDR)) %>%
  mutate(group = factor(case_when(
    rownames(dea_results) %in% rownames(upregulated) ~ 'up',
    rownames(dea_results) %in% rownames(downregulated) ~ 'down',
    .default = 'insig'))) %>%
  mutate(group = relevel(group, ref='insig')) %>%
  ggplot(aes(logFC, neg_log_pval, color=group)) +
  geom_point(alpha=0.6) +
  theme_bw() +
  theme(legend.key = element_blank(), 
        legend.position = 'inside',
        legend.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        panel.grid = element_blank(),
        legend.position.inside = c(0.89, 0.83)) +
  geom_vline(xintercept = c(-5,5), linetype='dashed') +
  geom_hline(yintercept = -log10(0.01), linetype='dashed') +
  scale_color_manual(values=c('black', 'seagreen', 'firebrick'), 
                     labels=c('Not Significant', 
                              'Down Regulated' , 
                              'Up Regulated')) +
  labs(title='Volcano Plot', x=expression('LogFC'), 
       y=expression('-Log'[10]*' (FDR corrected P-values)'), color='') +
  scale_x_continuous(breaks=seq(-10,10,5))


TCGAVisualize_volcano(dea_results$logFC, dea_results$FDR, 
                      show="both", y.cut = -log10(0.01),
                      filename = 'volcano.pdf', x.cut=c(-5,5))


## Prepare metadata for ML
metadata <- data.frame(sampleIDs=c(tumor_samples, normal_samples), 
                       group=c(rep('Tumor', length(tumor_samples)),
                               rep('Normal', length(normal_samples))
                               )
                       )

write.csv(metadata, 'metadata.csv', row.names = F)


## Get preprocessed data
preprocessed_counts <- data.frame(
  data_preprocessed[, c(tumor_samples, normal_samples)]
)


raw_counts <- assay(data[, c(tumor_samples, normal_samples)])
raw_counts <- as.data.frame(raw_counts)



head(preprocessed_counts)
head(raw_counts)


# write to csv
write.csv(preprocessed_counts, 'preprocessed_counts.csv', row.names = T)
write.csv(raw_counts, 'raw_counts.csv', row.names = T)


## Enrichment Analysis

all_files <- list.files('enrichment/', full.names = T)
all_files



# Function to create lollipop chart
plot_lollipop_chart <- function(file_list, enrichment_type, gene_type, 
                                title=NULL, n=10){
  filenames <- file_list[str_detect(file_list, enrichment_type)]
  
  if (length(filenames) == 0) stop(paste(
    'enrichment type:', enrichment_type, 'file not found'), call.=F)
  
  gene_type_file <- filenames[str_detect(filenames, gene_type)]
  
  
  read_csv(gene_type_file, show_col_types = F) |>
    # filter by statistically significant result
    filter(`Enrichment FDR` < 0.05) |>
    slice_min(`Enrichment FDR`, n=n) |>
    select(Pathway, everything()) |>
    mutate(`Enrichment FDR` = -log10(`Enrichment FDR`)) |>
    mutate(Pathway = str_to_title(
      str_trim(str_remove_all(Pathway, 'GO:\\d+\\s|Path:\\Hsa\\d+')))) |>
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
    labs(title=title, 
         subtitle = paste(gene_type, 'regulated genes'), y='', 
         x=expression("-Log"[10]*" (pvalue)"))
  
}


# biological process
plot_lollipop_chart(all_files, 'Bio', 'Up', 'Biological Process Pathway')
plot_lollipop_chart(all_files, 'Bio', 'Down', 'Biological Process Pathway')


# Molecular Function
plot_lollipop_chart(all_files, 'Molecular', 'Up', 'Molecular Function Pathway')
plot_lollipop_chart(all_files, 'Molecular', 'Down', 'Molecular Function Pathway')


# Cellular localisation
plot_lollipop_chart(all_files, 'Cell', 'Up', 'Cellular Localisation')
plot_lollipop_chart(all_files, 'Cell', 'Down', 'Cellular Localisation')


# KEGG pathway
plot_lollipop_chart(all_files, 'KEGG', 'Up', 'KEGG Pathway')
plot_lollipop_chart(all_files, 'KEGG', 'Down', 'KEGG Pathway')
