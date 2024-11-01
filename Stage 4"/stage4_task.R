################ Differential Gene Expression #############
rm(list=ls())

# load libraries
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SummarizedExperiment))


# build a query to retrieve DNA expression data
query <- GDCquery(project = 'TCGA-LGG',
                  data.category = 'Transcriptome Profiling',
                  experimental.strategy = 'RNA-Seq',
                  workflow.type = 'STAR - Counts',
                  access = 'open',
                  data.type = 'Gene Expression Quantification')


# download methylated data
output_query <- getResults(query)

GDCdownload(query, files.per.chunk = 267)


# prepare data
## run to save prepared data as TCGA_LGG.rda
#expr_data <- GDCprepare(query, 
#                       summarizedExperiment = TRUE, 
#                      save=TRUE, 
#                     save.filename = 'TCGA_LGG.rda')


# load saved prepared data and assign to expr_data
load('TCGA_LGG.rda')

expr_data <- data

rm(data)


# creating clinical data
clinical_data <- data.frame(
  row.names = expr_data$barcode,
  patient_id = expr_data$patient,
  age = expr_data$age_at_index,
  gender = expr_data$gender,
  tumor_descriptor = expr_data$tumor_descriptor,
  n_mutations = expr_data$paper_Mutation.Count,
  tumor_grade = expr_data$paper_Grade,
  tumor_type = expr_data$paper_Histology,
  IDH_status = expr_data$paper_IDH.status,
  IDH_cluster_met = expr_data$paper_Pan.Glioma.DNA.Methylation.Cluster,
  IDH_cluster = expr_data$paper_Pan.Glioma.RNA.Expression.Cluster,
  patient_status = expr_data$vital_status
)

head(clinical_data)


## checking for missing records
anyNA.data.frame(clinical_data)

apply(is.na.data.frame(clinical_data),2, sum)


# drop samples without IDH status 
# replace missing values in tumor grade and type as unknown
clinical_data <- clinical_data |>
  drop_na(IDH_status) |>
  mutate(across(c(tumor_grade, tumor_type), \(x) as.character(x))) |>
  replace_na(list(tumor_grade = 'unk', tumor_type = 'unknown')) |>
  mutate(tumor_grade = factor(tumor_grade),
         tumor_type = factor(tumor_type))


# plot Age distribution
ggplot(clinical_data %>% select(age) %>% drop_na(), 
       aes(y=after_stat(density), x=age)) +
  geom_histogram(fill='steelblue', color='white', bins=40) +
  geom_density() +
  geom_vline(xintercept = mean(clinical_data$age, na.rm=T), 
             linetype='dashed', color='red') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size=12, face='bold')) +
  labs(title='Age Distribution', x='\nAge', y='Density\n') +
  scale_x_continuous(breaks=seq(0,100,10), expand = c(0.1,0.15,0.1,0.05))


## Data Preprocessing

# removing sample outliers with Pearson correlation less than 0.6
data_preprocessed <- TCGAanalyze_Preprocessing(
  expr_data, 
  cor.cut = 0.6, 
  datatype = 'unstranded',
  filename = 'LGG_preprocessing_result.png')


# normalising by gene length and GC content
data_preprocessed <- TCGAanalyze_Normalization(data_preprocessed, 
                                               geneInfo = geneInfoHT, 
                                               method='geneLength')


data_preprocessed <- TCGAanalyze_Normalization(data_preprocessed, 
                                               geneInfo = geneInfoHT, 
                                               method='gcContent')


# selected samples
data_preprocessed <- data_preprocessed[, rownames(clinical_data)]


# filter zero counts using quantile method (< 25% percentile)
data_preprocessed <- TCGAanalyze_Filtering(data_preprocessed, 
                                           method='quantile',
                                           qnt.cut = 0.25)


# Performing Differential gene expression

# by tumor grade
grade2 <- subset(clinical_data, tumor_grade == "G2")$tumor_grade
grade3 <- subset(clinical_data, tumor_grade == "G3")$tumor_grade

# by IDH status
IDH_wt <- subset(clinical_data, IDH_status == "WT")$IDH_status
IDH_mutant <- subset(clinical_data, IDH_status == "Mutant")$IDH_status


# Performing differential expression analysis
dea_results_IDH <- TCGAanalyze_DEA(
  mat1 = data_preprocessed[, IDH_wt],
  mat2 = data_preprocessed[, IDH_mutant],
  Cond1type = "WildType",
  Cond2type = "Mutant",
  method = "glmLRT",
  pipeline = 'edgeR'
)


# Performing differential expression analysis
dea_results_grade <- TCGAanalyze_DEA(
  mat1 = data_preprocessed[, grade2],
  mat2 = data_preprocessed[, grade3],
  Cond1type = "G2",
  Cond2type = "G3",
  method = "glmLRT",
  pipeline = 'edgeR'
)

# filter upregulated and downregulated genes

# IDH
upregulated_IDH <- dea_results_IDH |>
  filter(FDR != 0) |>
  filter(FDR < 0.05, logFC > 2)

downregulated_IDH <- dea_results_IDH |>
  filter(FDR != 0) |>
  filter(FDR < 0.05, logFC < -2)


# Tumor grade
upregulated_grade <- dea_results_grade |>
  filter(FDR != 0) |>
  filter(FDR < 0.05, logFC > 2)

downregulated_grade <- dea_results_grade |>
  filter(FDR != 0) |>
  filter(FDR < 0.05, logFC < -2)

# saving raw and preprocessed counts
preprocessed_counts <- data.frame(data_preprocessed)
raw_counts <- data.frame(assay(expr_data, 'unstranded'))

# write to csv
write.csv(preprocessed_counts, 'preprocessed_counts.csv', row.names = T)
write.csv(raw_counts, 'raw_counts.csv', row.names = T)

# save metadata
write.csv(clinical_data, 'metadata.csv', row.names=T)

# saving DGE results
write.csv(dea_results_IDH, 'DE_result_IDH_status.csv', row.names=T)
write.csv(dea_results_grade, 'DE_result_tumor_grade.csv', row.names=T)


# save upregulated and downregulated genes
upregulated_IDH |> write.csv('upregulated_IDH.csv', row.names = T)
downregulated_IDH |> write.csv('downregulated_IDH.csv', row.names = T)

upregulated_grade |> write.csv('upregulated_grade.csv', row.names = T)
downregulated_grade |> write.csv('downregulated_grade.csv', row.names = T)


# volcano plots
plot_volcano <- function(DE_result, upregulated, 
                         downregulated, title='Volcano Plot'){
  DE_result %>%
    mutate(neg_log_pval = -log10(FDR)) %>%
    mutate(group = factor(case_when(
      rownames(DE_result) %in% rownames(upregulated) ~ 'up',
      rownames(DE_result) %in% rownames(downregulated) ~ 'down',
      .default = 'insig'))) %>%
    mutate(group = relevel(group, ref='insig')) %>%
    filter(FDR != 0) |> 
    ggplot(aes(logFC, neg_log_pval, color=group)) +
    geom_point(alpha=0.6) +
    theme_classic() +
    theme(legend.key = element_blank(), 
          legend.position = 'inside',
          legend.title = element_text(face='bold', size=10),
          plot.title = element_text(face='bold'),
          panel.grid = element_blank(),
          legend.position.inside = c(0.89, 0.83)) +
    geom_vline(xintercept = c(-2,2), linetype='dashed') +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') +
    scale_color_manual(values=c('black', 'seagreen', 'firebrick'), 
                       labels=c('Not Significant', 
                                'Down Regulated' , 
                                'Up Regulated')) +
    labs(title=title, x=expression('LogFC'), 
         y=expression('-Log'[10]*' (FDR corrected P-values)'), color='') +
    scale_x_continuous(breaks=seq(-10,10,5))
}


plot_volcano(dea_results_IDH, upregulated_IDH, 
             downregulated_IDH, 
             'Volcano plot for IDH Status')

plot_volcano(dea_results_grade, upregulated_grade, 
             downregulated_grade, 
             'Volcano plot for Tumor Grade')


TCGAVisualize_volcano(dea_results_IDH$logFC, dea_results_IDH$FDR, 
                      show="both", y.cut = -log10(0.05),
                      filename = 'volcano_IDH_status.png', 
                      x.cut=c(-2,2))

TCGAVisualize_volcano(dea_results_grade$logFC, dea_results_grade$FDR, 
                      show="both", y.cut = -log10(0.05),
                      filename = 'volcano_tumor_grade.png', 
                      x.cut=c(-2,2))

rm(list=ls())


#################### Machine Learning ##############################

########## For IDH status #############
suppressPackageStartupMessages(library(tidyverse))
suppressMessages(library(kknn))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))


# load preprocessed data
counts_data <- read_csv('preprocessed_counts.csv', show_col_types = FALSE)
counts_data <- counts_data |> column_to_rownames(var='...1')

metadata <- read_csv('metadata.csv', show_col_types = F)
metadata <- metadata |> column_to_rownames(var='...1')

dea_results_IDH <- read.csv('DE_result_IDH_status.csv', row.names=1)

# clean variable names
names(counts_data) <- gsub(x=names(counts_data), 
                           pattern = '\\.', 
                           replacement = '-')


# merging counts data to metadata
prep_data <- t(counts_data) |>
  as.data.frame() |>
  rownames_to_column(var='sampleID') |>
  inner_join(metadata |> select(IDH_status) |> rownames_to_column(var='sampleID'), 
             by='sampleID') |>
  column_to_rownames('sampleID')


set.seed(12)
data.splits <- caret::createDataPartition(prep_data$IDH_status, 
                                          p=0.75, list=FALSE)

# feature selection (1st level)
remove_genes <- dea_results_IDH |>
  filter(FDR >= 0.05) |>
  rownames()

# remove genes that are not statistically significant
prep_data <- select(prep_data, -all_of(remove_genes)) |>
  mutate(IDH_status = factor(IDH_status, levels=c('WT', 'Mutant')))

# create train and test sets
train.data <- prep_data[data.splits, ]
test.data <- prep_data[-data.splits,]


## Standardisation scaler preprocessor
scale_preprocessor <- function(df_train, df_test=NULL){
  
  df_train_log <- log2(1+df_train)
  
  mean_vals <- colMeans(as.matrix(df_train_log))
  sd_vals <- matrixStats::colSds(as.matrix(df_train_log))
  
  df_train_log <- (t(df_train_log) - mean_vals) / sd_vals
  df_train_log <- data.frame(t(df_train_log))
  
  if (!is.null(df_test)){
    df_test_log <- log2(1+df_test)
    
    df_test_log <- (t(df_test_log) - mean_vals) / sd_vals
    df_test_log <- data.frame(t(df_test_log))
    
    return(list(df_train_log, df_test_log))
  } else{
    return (df_train_log)
  }
}


xtrain <- train.data %>% select(-IDH_status)
xtest <- test.data %>% select(-IDH_status)
ytrain <- train.data$IDH_status
ytest <- test.data$IDH_status

# scale train and test
scaled_data <- scale_preprocessor(xtrain, xtest)
scaled_train <- scaled_data[[1]]
scaled_test <- scaled_data[[2]]

## Feature selection using Lasso regression

# fit a lasso logistic regression model
# alpha = 1 (Lasso)
lasso <- glmnet(as.matrix(scaled_train), ytrain, lambda =2e-7,
                family=binomial, alpha=1, standardize = F)


# principal component analysis
pca <- prcomp(log2(1+select(prep_data, -IDH_status)), scale.=TRUE)

exp_var <- pca$sdev**2/sum(pca$sdev**2)

ggplot(pca$x, aes(PC1, PC2, color=metadata$IDH_status)) +
  geom_point() +
  theme_bw() +
  theme(legend.key = element_blank(), 
        legend.position = 'top',
        plot.title=element_text(face='bold'),
        panel.grid = element_blank()) +
  labs(title='Principal Component Analysis by IDH status', 
       color='IDH status', x=paste0('PC1 (', round(100*exp_var[1],1), '%)'),
       y=paste0('PC2 (', round(100*exp_var[2],1), '%)'))


# select coefficients with non-zero values
selected_features_df <- as.data.frame(as.matrix(coefficients(lasso))) |>
  arrange(s0) |>
  filter(s0 != 0) |>
  rownames_to_column('features') |>
  mutate(odds_ratio = exp(s0)) |>
  filter(features != '(Intercept)') |>
  left_join(dea_results_IDH |> select(gene_name) |> rownames_to_column('features'), 
            by='features') |>
  replace_na(list(gene_name = 'unknown gene'))

selected_features_df %>%
  mutate(s0 = abs(s0)) |>
  slice_max(s0, n=20) |>
  ggplot(aes(y=reorder(gene_name, s0), x=s0)) +
  geom_col(fill='firebrick', alpha=0.6) +
  geom_text(aes(label=round(odds_ratio,2)), hjust = 1.2, size=3.2, 
            fontface='bold', color='white', vjust=.2) +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title = 'Feature Coefficients', subtitle='(Top 20)', 
       x='Absolute Coefficients', y='')

ggsave('lasso_features_IDH.png', bg='white')

selected_features <- selected_features_df$features

## Modelling

# function for hyperparameter tuning
tune_parameters <- function(trainX, trainY, param_grid, 
                            cv=5, model_name='knn', 
                            seed=5, scale.=T, 
                            positive='Mutant'){
  set.seed(seed)
  folds <- createFolds(trainY, k=cv, list=T, returnTrain = TRUE)
  
  scores <- c()
  
  for (j in 1:nrow(param_grid)){
    params <- param_grid[j, ] 
    
    fold_scores = list(acc=c(), f1=c())
    
    # get data folds
    for (fold in folds){
      xtr <- trainX[fold, ]
      ytr <- trainY[fold]
      xte <- trainX[-fold, ]
      yte <- trainY[-fold]
      
      if (scale.){
        scale_xtr_xte <- scale_preprocessor(xtr, xte)
        xtr <- scale_xtr_xte[[1]]
        xte <- scale_xtr_xte[[2]]
      }
      
      # build model  
      if (model_name == 'knn'){
        model <- train.kknn(formula('target ~ .'), 
                            data=cbind(select(xtr, selected_features),
                                       target=ytr), 
                            scale=FALSE, ks=params)
        
        preds <- predict(model, xte)
        acc = mean(yte == preds)
        f1 <- MLmetrics::F1_Score(yte, preds, positive = positive)
        fold_scores$acc <- c(fold_scores$acc, acc)
        fold_scores$f1 <- c(fold_scores$f1, f1)
        
      } else if (model_name == 'rf'){
        model <- randomForest(formula('target ~ .'),
                              data=cbind(xtr[,selected_features], 
                                         target=ytr), importance = TRUE, 
                              replace=FALSE,
                              ntree=params$ntree, mtry=params$mtry, 
                              nodesize=params$nodesize) 
        
        preds <- predict(model, xte)
        acc = mean(yte == preds)
        f1 <- MLmetrics::F1_Score(yte, preds, positive = positive)
        fold_scores$acc <- c(fold_scores$acc, acc)
        fold_scores$f1 <- c(fold_scores$f1, f1)
      }
    }
    scores <- rbind(scores, c(acc.mean=mean(fold_scores$acc),
                              acc.sd=sd(fold_scores$acc),
                              f1.mean=mean(fold_scores$f1),
                              f1.sd=sd(fold_scores$f1)))
  }
  scores <- data.frame(scores)
  res <- arrange(cbind(param_grid, scores), desc(f1.mean))
  return (res)
}

confusion_matrix_plot <- function(predictions, actual, 
                                  title='Confusion Matrix'){
  as.data.frame(table(Predictions=predictions, Actual=actual)) |>
    ggplot(aes(y=Actual, x=Predictions, fill=Freq)) +
    geom_tile(alpha=0.9) +
    geom_text(aes(label=Freq), color='white', fontface='bold') +
    theme_light() +
    theme(panel.grid = element_blank(), 
          legend.position = 'none',
          axis.text.y = element_text(angle=90, hjust=0.5),
          plot.title = element_text(face='bold')) +
    labs(title = title, y='Actual', x='Predictions')
}


## Function to get model performance metrics

model_performance <- function(model, X, y, positive=NULL){
  model_predictions <- predict(model, X)
  
  acc <- mean(y == model_predictions)
  recall <- MLmetrics::Recall(y, model_predictions, positive = positive)
  f1 <- MLmetrics::F1_Score(y, model_predictions, positive = positive)
  specificity <- MLmetrics::Specificity(y, model_predictions, positive = positive)
  precision <- MLmetrics::Precision(y, model_predictions, positive = positive)
  
  res <- t(data.frame(Accuracy = acc, Recall = recall, F1=f1, 
                      Specificity=specificity, Precision=precision))
  data.frame(scores=100*res)
}

## kNN
set.seed(10)
tune_parameters(xtrain, ytrain, 
                expand.grid(ks=seq(3,21,2)), 
                cv=5, scale. = TRUE)

knn.model <- train.kknn(formula('IDH_status ~ .'), 
                        data=cbind(select(scaled_train, 
                                          selected_features),
                                   IDH_status=ytrain), 
                        scale=FALSE, ks=3)

knn_predictions <- predict(knn.model, scaled_test)

confusion_matrix_plot(knn_predictions, ytest, 'K-Nearest Neighbors')
ggsave('knn_confmat_IDH.png', dpi = 300)

model_performance(knn.model, scaled_test, ytest, 'Mutant')


## Random forest
set.seed(10)
rf_param <- expand.grid(ntree=c(100,150,200), 
                        mtry=c(3,4,5,7), 
                        nodesize=c(1))

tune_parameters(xtrain, ytrain, rf_param, 
                model_name = 'rf', scale.=FALSE,
                positive='Mutant')


set.seed(10)
rf <- randomForest(formula('IDH_status ~ .'),
                   data=cbind(scaled_train[,selected_features], 
                              IDH_status=ytrain), 
                   importance = TRUE, ntree=100, mtry=3)

# variable importance
rf_importance <- select(varImp(rf), 2)

rf_importance |>
  rownames_to_column(var='features') |>
  left_join(dea_results_IDH |> 
              select(gene_name) |> 
              rownames_to_column('features'), 
            by='features') |>
  replace_na(list(gene_name = 'unknown gene')) |>
  slice_max(Mutant, n=20) |>
  ggplot(aes(y=reorder(gene_name, Mutant), x=Mutant)) +
  geom_col(fill='firebrick', alpha=0.6) +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title = 'Feature Importance (Random Forest)', 
       subtitle='(Top 20)', x='Score', y='')

ggsave('rf_features_IDH.png', bg='white')

# random forest predictions
rf_predictions <- predict(rf, scaled_test)
confusion_matrix_plot(rf_predictions, ytest, 'Random Forest')
ggsave('rf_confmat_IDH.png', dpi = 300)

model_performance(rf, scaled_test, ytest, 'Mutant')

# merge confusion matrix for kNN and RF
ggpubr::ggarrange(confusion_matrix_plot(knn_predictions, ytest, 
                                        'Confusion Matrix (KNN)'),
                  confusion_matrix_plot(rf_predictions, ytest, 
                                        'Confusion Matrix (Random Forest)')
)
ggsave('confmat_IDH.png', dpi=300)


####### Tumor Grade ########

dea_results_grade <- read.csv('DE_result_tumor_grade.csv', row.names=1)

# merging counts data to metadata
prep_data <- t(counts_data) |>
  data.frame() |>
  rownames_to_column(var='sampleID') |>
  inner_join(metadata |> 
               select(tumor_grade, tumor_type) |> 
               rownames_to_column(var='sampleID'), 
             by='sampleID') |>
  column_to_rownames('sampleID')

# remove unknown tumor grades
prep_data <- prep_data |>
  filter(tumor_grade != 'unk')


metadata <- metadata |>
  filter(tumor_grade != 'unk')


remove_genes <- dea_results_grade |>
  filter(FDR >= 0.05) |>
  rownames()

# Feature selection

# remove genes that are not statistically significant
prep_data <- select(prep_data, -all_of(remove_genes)) |>
  mutate(tumor_grade = factor(tumor_grade))

# PCA by tumor grade and type
pca <- log2(1+select(prep_data, -tumor_grade, -tumor_type)) |>
  prcomp(scale. = TRUE)

# explained variance
exp_var <- pca$sdev**2/sum(pca$sdev**2)


ggplot(pca$x, aes(PC1, PC2, color=metadata$tumor_grade)) +
  geom_point() +
  theme_bw() +
  theme(legend.key = element_blank(), 
        legend.position = 'top',
        plot.title=element_text(face='bold'),
        panel.grid = element_blank()) +
  labs(title='Principal Component Analysis by Tumour Grade', 
       color = 'Tumour Grade')

ggplot(pca$x, aes(PC1, PC2, color=metadata$tumor_type)) +
  geom_point() +
  theme_bw() +
  theme(legend.key = element_blank(), 
        legend.position = 'top',
        plot.title=element_text(face='bold'),
        panel.grid = element_blank()) +
  labs(title='Principal Component Analysis by Tumour Type', 
       color = 'Tumour Type', 
       x=paste0('PC1 (', round(100*exp_var[1],1), '%)'),
       y=paste0('PC2 (', round(100*exp_var[2],1), '%)'))



# split data
set.seed(12)
data.splits <- caret::createDataPartition(prep_data$tumor_grade, 
                                          p=0.75, list=FALSE)

train.data <- prep_data[data.splits, ]
test.data <- prep_data[-data.splits,]

xtrain <- train.data |> 
  select(-tumor_grade, -tumor_type)

xtest <- test.data |> 
  select(-tumor_grade, -tumor_type)

ytrain <- train.data$tumor_grade
ytest <- test.data$tumor_grade


# scale data
scaled_data <- scale_preprocessor(xtrain, xtest)

scaled_train <- scaled_data[[1]]
scaled_test <- scaled_data[[2]]


## Feature selection using Lasso regression
# fit a lasso logistic regression model
# alpha = 1 (Lasso)
set.seed(0)
lasso <- glmnet(as.matrix(scaled_train), ytrain, 
                #lambda =2e-7,
                family=binomial, alpha=1, 
                standardize = F)

# fit lasso using min lasso penalty (lambda)
lasso <- glmnet(as.matrix(scaled_train), ytrain, 
                lambda = min(lasso$lambda),
                family=binomial, alpha=1, 
                standardize = F)

### Feature coefficients
# select coefficients with non-zero values
selected_features_df <- as.data.frame(as.matrix(coefficients(lasso))) |>
  arrange(s0) |>
  filter(s0 != 0) |>
  rownames_to_column('features') |>
  mutate(odds_ratio = exp(s0)) |>
  filter(features != '(Intercept)') |>
  left_join(dea_results_grade |> select(gene_name) |> rownames_to_column('features'), 
            by='features') |>
  replace_na(list(gene_name = 'unknown gene'))

selected_features_df %>%
  mutate(s0 = abs(s0)) |>
  slice_max(s0, n=20) |>
  ggplot(aes(y=reorder(gene_name, s0), x=s0)) +
  geom_col(fill='firebrick', alpha=0.6) +
  geom_text(aes(label=round(odds_ratio,2)), hjust = 1.2, size=3.2, 
            fontface='bold', color='white', vjust=.2) +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title = 'Feature Coefficients', subtitle='(Top 20)', 
       x='Absolute Coefficients', y='')

ggsave('lasso_features_grade.png', bg='white')



selected_features <- selected_features_df$features
# kNN

set.seed(5)
tune_parameters(xtrain, ytrain, 
                expand.grid(ks=seq(5,21,2)), 
                cv=5, scale.=TRUE,
                positive='G3')

knn.model <- train.kknn(formula('tumor_grade ~ .'), 
                        data=cbind(select(scaled_train, selected_features),
                                   tumor_grade=ytrain), 
                        scale=FALSE, ks=15)

knn_predictions <- predict(knn.model, select(scaled_test, selected_features))
confusion_matrix_plot(knn_predictions, ytest, 'Confusion Matrix (KNN)')
ggsave('knn_confmat_grade.png', dpi = 300)

# comparing model's predictions by tumor type

# tumor type vs grade
table(Tumour.Grade=ytest, 
      Tumor.Type=metadata[rownames(xtest), 'tumor_type'])

# predictions vs tumor grade
table(Predictions=knn_predictions, Tumour.Grade=ytest)

# predictions vs tumor type
table(Predictions=knn_predictions,
      Tumor_type=metadata[rownames(xtest), 'tumor_type'])

# Model performance
model_performance(knn.model, 
                  select(scaled_test, selected_features), 
                  y=ytest, positive='G3')


## Random Forest

# get optimal parameters
mtry = as.integer(seq(5,85,20)/100 * nrow(selected_features))


rf_param <- expand.grid(ntree=c(200,500,800,1000), 
                        mtry=mtry, 
                        nodesize=c(1,3,5,7))
set.seed(5)
tune_parameters(xtrain, ytrain, rf_param, 
                model_name = 'rf', scale. = FALSE)


set.seed(5)
rf <- randomForest(formula('tumor_grade ~ .'),
                   data=cbind(xtrain[,selected_features], 
                              tumor_grade=ytrain), 
                   replace=FALSE,
                   importance = TRUE, ntree=800, mtry=11, nodesize=1)


rf_importance <- select(varImp(rf), 2)

rf_importance |>
  rownames_to_column(var='features') |>
  left_join(dea_results_grade |> select(gene_name) |> rownames_to_column('features'), 
            by='features') |>
  replace_na(list(gene_name = 'unknown gene')) |>
  slice_max(G3, n=20) |>
  ggplot(aes(y=reorder(gene_name, G3), x=G3)) +
  geom_col(fill='firebrick', alpha=0.6) +
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title = 'Feature Importance (Random Forest)', 
       subtitle='(Top 20)', x='Score', y='')

ggsave('rf_features_grade.png', bg='white')


# random forest predictions
rf_predictions <- predict(rf, select(xtest, selected_features))

f <- subset(metadata, (rownames(metadata) %in% names(rf_predictions)))

table(predictions=rf_predictions, tumor_type=f$tumor_type)
table(predictions= rf_predictions, Actual=ytest)
table(tumor_grade=ytest, tumor_type=f$tumor_type)

# confusion matrix
confusion_matrix_plot(rf_predictions, ytest, 'Confusion Matrix (Random Forest)')
ggsave('rf_confmat_grade.png', dpi = 300)

# model performance
print(model_performance(rf, xtest, ytest))

ggpubr::ggarrange(confusion_matrix_plot(knn_predictions, 
                                        ytest, 
                                        'Confusion Matrix (KNN)'),
                  confusion_matrix_plot(rf_predictions, 
                                        ytest, 
                                        'Confusion Matrix (Random Forest)')
)

ggsave('confmat_grade.png', dpi=300)

# visualise model performance on tumor type and grade
data.frame(
  models = factor(c('knn', 'knn', 'rf', 'rf', 'actual', 'actual')),
  tumor_grades = factor(c('G2', 'G3', 'G2', 'G3', 'G2', 'G3')),
  astrocytoma = c(19, 27, 15,31,18,28),
  oligoastrocytoma = c(17, 9, 17, 9,15,11),
  oligodendroglioma = c(30, 11, 26, 15,20,21)
) |>
  pivot_longer(c(-models, -tumor_grades), 
               names_to = 'tumor_type', 
               values_to = 'predictions') |>
  ggplot(aes(models, predictions, fill=tumor_grades))+
  geom_col(alpha=0.7) +
  geom_text(aes(label=predictions), 
            size=3.3,
            position=position_stack(0.5)) +
  facet_grid(~tumor_type, labeller = as_labeller(function(x) str_to_title(x))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        plot.title = element_text(face='bold', size=12),
        legend.position = 'top',
        legend.key.size =  unit(0.15, 'in'),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank()) +
  scale_x_discrete(labels=c('Actual', 'KNN', 'RF')) +
  labs(title='Model Predictions by Tumor Grades and Types vs Actual values', 
       x='', y='Predictions', fill='Tumor Grade') +
  scale_fill_manual(values=c('cornflowerblue', 'lightcoral'))

ggsave('models_tumour_grades_types.png', dpi=300)

rm(list=ls())

