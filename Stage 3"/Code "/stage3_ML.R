
# install.packages(c('MLmetrics', 'matrixStats', 'kknn', 'glmnet'))

suppressPackageStartupMessages(library(tidyverse))
suppressMessages(library(kknn))
suppressPackageStartupMessages(library(glmnet))


# load preprocessed data
counts_data <- read_csv('preprocessed_counts.csv', show_col_types = FALSE)
counts_data <- counts_data |> column_to_rownames(var='...1')

metadata <- read_csv('metadata.csv', show_col_types = F)

dea_results <- read.csv('DEA_results_Tumor_vs_Normal.csv', row.names=1) 

head(counts_data)

dim(counts_data)

# clean variable names
names(counts_data) <- gsub(x=names(counts_data), 
                           pattern = '\\.', 
                           replacement = '-')


## Visualise differential genes by tissue sample types

# randomly selecting statistically significant genes with FDR < 1% & |logFC| > 2
set.seed(10)
selected_ids <- sample(1:nrow(filter(dea_results, FDR < 0.05, abs(logFC) > 2)), 25)
selected_genes <- rownames(counts_data)[selected_ids]



## ----fig.height=9, fig.width=12----------------------------------------------------------------------------------------------------------
t(counts_data[selected_genes, ]) |>
  as.data.frame() |>
  rownames_to_column(var='ID') |>
  dplyr::inner_join(metadata, by=c('ID' = 'sampleIDs')) |>
  mutate(group = factor(group)) |>
  pivot_longer(cols=c(-ID,-group), values_to = 'expr', names_to = 'genes') |>
  select(-ID) |>
  mutate(expr = log2(1+expr)) |>
  ggplot(aes(group, expr, fill=group)) +
  geom_boxplot() +
  facet_wrap(~genes, nrow = 5, ncol = 5, scales='free_y') +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid = element_blank(), 
        plot.title = element_text(face='bold'),
        strip.text = element_text(face='bold')) +
  labs(title= 'Gene expression levels of selected genes', x='Tissue Type', 
       y=expression('Expression Levels '*'(Log'[2]*')'))


# plot gene expression level distribution in each sample
## ----fig.height=6, fig.width=12----------------------------------------------------------------------------------------------------------
counts_data |>
  rownames_to_column(var='ID') |>
  pivot_longer(cols=-ID, names_to = 'samples', values_to = 'exprs') |>
  mutate(samples = str_extract(samples, 'TCGA-\\w{2}-\\w{4}')) |>
  ggplot(aes(samples, log2(1+exprs))) +
  geom_boxplot(fill='skyblue') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  scale_x_discrete(labels=paste0('s', seq(1,40), sep=' '))



# summary statistics of genes (log2 expression values)
counts_data |>
  rownames_to_column(var='ID') |>
  pivot_longer(cols=-ID, names_to = 'samples', values_to = 'exprs') |>
  mutate(exprs = log2(1+exprs)) |>
  summarize(mean_val = mean(exprs), 
            median_val = median(exprs),
            sd_val = sd(exprs), 
            .by = ID) |>
  arrange(desc(mean_val)) |>
  head(10)


## feature selection
# remove not statistically significant genes
remove_genes <- dea_results |>
  filter(FDR >= 0.05) |>
  rownames()

# droppping statistically insignificant genes
counts_data <- counts_data[-which(remove_genes %in% rownames(counts_data)),]


# Data Preparation
# transpose for ML (merging target class)
prep_data <- t(counts_data) |>
  as.data.frame() |>
  rownames_to_column('sampleIDs') |>
  dplyr::inner_join(metadata, by='sampleIDs') |>
  mutate(group = factor(group), 
         group = relevel(group, ref='Normal')) |>
  column_to_rownames('sampleIDs')


length(remove_genes)


## Target class distribution
ggplot(metadata, aes(group)) + 
  geom_bar(fill='firebrick', alpha=0.6) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title='Target class distribution', x='Tumor type', y='Frequency') +
  scale_y_continuous(breaks=seq(0,40,5))


## Data splitting into train and test data

# split into train and test data
set.seed(20)

split_data <- function(df, p=0.75){
  # get majority and minority class
  majority_class <- prep_data |> 
  filter(group == 'Tumor') 
  
  minority_class <- prep_data |> 
  filter(group == 'Normal') 
  
  # split majority and minority class into train and test data
  majority_train <- majority_class |>
  slice_sample(n=as.integer(p*nrow(majority_class)))
  
  majority_test <- majority_class |> 
    filter(!(rownames(majority_class)  %in% rownames(majority_train)))
  
  # sample 2 from minority class
  minority_train <- minority_class |>
  slice_sample(n=2)
  
  minority_test <- minority_class |> 
    filter(!(rownames(minority_class)  %in% rownames(minority_train)))
  
  
  train <- bind_rows(majority_train, minority_train)
  test <- bind_rows(majority_test, minority_test)
  
  return(list(train=train, test=test))
}

train_test_data <- split_data(prep_data)


## Get train and test data
train_data <- train_test_data$train
test_data <- train_test_data$test

xtrain <- train_data %>% select(-group)
xtest <- test_data %>% select(-group)
ytrain <- train_data$group
ytest <- test_data$group


## Standardisation scaler preprocessor
scale_preprocessor <- function(df_train, df_test=NULL){
  df_train_log <- log2(1+df_train)
  
  mean_vals <- rowMeans(as.matrix(t(df_train_log)))
  sd_vals <- matrixStats::rowSds(as.matrix(t(df_train_log)))
  
  df_train_log <- (t(df_train_log) - mean_vals) / sd_vals
  
  if (!is.null(df_test)){
    df_test_log <- log2(1+df_test)
    df_test_log <- (t(df_test_log) - mean_vals) / sd_vals
    
    return(list(as.data.frame(t(df_train_log)), 
                as.data.frame(t(df_test_log))
                )
           )
  } else{
    return (as.data.frame(t(df_train_log)))
  }
}


scaled_data <- scale_preprocessor(xtrain, xtest)

scaled_train <- scaled_data[[1]]
scaled_test <- scaled_data[[2]]


# check for genes with missing values as a result of scaling
# 5 variables have missing values in them
which(scaled_train |>
  apply(2, function(x) sum(is.na(x))) > 0)


remove_na_variables <- c("ENSG00000268799", "ENSG00000269138", 
                         "ENSG00000273177", "ENSG00000275508", 
                         "ENSG00000279273")

# remove genes with all missing values
scaled_train <- select(scaled_train, -all_of(remove_na_variables))
scaled_test <- select(scaled_test, -all_of(remove_na_variables))


## Feature selection using Lasso regression
# fit a lasso logistic regression model
# alpha = 1 (Lasso)
lasso <- glmnet(as.matrix(scaled_train), ytrain, lambda =2e-7,
                family=binomial, alpha=1, standardize = F)


# predictions on train
lasso_probabilities <- predict(lasso, as.matrix(scaled_train), type = 'response')
lasso_predictions <- factor(ifelse(lasso_probabilities > 0.5, 'Tumor', 'Normal'))
lasso_predictions <- relevel(lasso_predictions, ref='Normal')


### Feature coefficients
## ----fig.width=8-------------------------------------------------------------------------------------------------------------------------
# select coefficients with non-zero values
selected_features <- as.data.frame(as.matrix(coefficients(lasso))) |>
  arrange(s0) |>
  filter(s0 != 0) |>
  rownames_to_column('features') |>
  mutate(odds_ratio = exp(s0)) |>
  filter(features != '(Intercept)')

selected_features %>%
  ggplot(aes(y=reorder(features, odds_ratio), x=odds_ratio)) +
  geom_col(fill='red', alpha=0.5) +
  geom_text(aes(label=round(s0,2)), hjust = -.05, size=3.4, fontface='bold') +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(face='bold')) +
  labs(title = 'Feature Coefficients', x='Odds Ratio', y='')

ggsave('lasso_features.png', bg='white')


## Modelling
s
# kNN model
knn_model <- train.kknn(group ~ ., data = cbind(
  select(scaled_train, all_of(selected_features$features)), 
  group=ytrain), kmax = 5)


knn_predictions <- predict(knn_model, scaled_test)


as.data.frame(table(Predictions=knn_predictions, Actual=ytest)) |>
  ggplot(aes(y=Actual, x=Predictions, fill=Freq)) +
  geom_tile(alpha=0.9) +
  geom_text(aes(label=Freq), color='white', fontface='bold') +
  theme_light() +
  theme(panel.grid = element_blank(), 
        legend.position = 'none',
        axis.text.y = element_text(angle=90, hjust=0.5),
        plot.title = element_text(face='bold')) +
  labs(title = 'Confustion Matrix', y='Actual', x='Predictions')

ggsave('knn_conf_mat.png', bg='white')


## Function to get model performance metrics

model_performance <- function(model, X, y){
  model_predictions <- predict(model, X)
  
  acc <- mean(y == model_predictions)
  recall <- MLmetrics::Recall(y, model_predictions, positive = 'Tumor')
  f1 <- MLmetrics::F1_Score(y, model_predictions, positive = 'Tumor')
  specificity <- MLmetrics::Specificity(y, model_predictions, positive = 'Tumor')
  precision <- MLmetrics::Precision(y, model_predictions, positive = 'Tumor')
  
  res <- t(data.frame(Accuracy = acc, Recall = recall, F1=f1, 
                      Specificity=specificity, Precision=precision))
  data.frame(scores=100*res)
}


print(model_performance(knn_model, scaled_test, ytest))

print(model_performance(knn_model, scaled_train, ytrain))