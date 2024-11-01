
# **Identifying Key Genetic Biomarkers in Cervical Cancer Using Machine Learning and Differential Gene Expression Analysis**


## **Introduction to Cervical Cancer**

Cervical cancer is a significant global health challenge and one of the leading causes of cancer deaths in women, especially in the developing world (Bedell et al., 2020). Globally, it is the fourth most common cancer type and the fourth cause of cancer-related death among women (Bray et al., 2024). It is mainly caused by persistent infection with high-risk strains of human papillomavirus (HPV) (Banzola et al., 2018). If untreated, HPV can lead to cervical intraepithelial neoplasia (CIN), eventually progressing to cervical cancer (Balasubramaniam et al., 2019). Fortunately, HPV vaccination and regular screening programs, such as Pap smears and HPV tests, have played crucial roles in reducing cervical cancer incidence and mortality worldwide (Bedell et al., 2020).

## **Methodology**

### **Dataset Description and Preprocessing Steps**

Cervical cancer data for this study was sourced from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga). The dataset comprised 309 samples: 304 primary tumor samples, 3 normal samples, and the remainder were metastasis cases. We selected primary tumor and normal tissue samples for analysis, narrowing the data down to 40 samples by selecting all 3 normal samples and randomly choosing 37 tumor samples stratified by age group.

### **Data Preparation and Preprocessing**

Data preprocessing and normalization were conducted using the TCGAbiolinks library in R. The steps included

- **Filtering genes:** Genes with a correlation less than 0.4 across samples were removed.
- **Normalization:** Expression levels were normalized based on gene length and GC content.
- **Exclusion of low-expression genes:** Genes with expression levels below the 25th percentile were excluded.

### **Biomarker Discovery and Machine Learning**

#### **Differential Gene Expression (DGE) Analysis**

We performed DGE analysis using the likelihood ratio method, identifying differentially expressed genes based on the criteria: False discovery rate (FDR) < 0.01 and Absolute log fold change > 5. Genes with FDR < 0.01 and LogFC < -5 were classified as downregulated, while genes with FDR < 0.01 and LogFC > 5 were classified as upregulated. The preprocessed data, consisting of 30,753 genes, was used for subsequent machine learning predictions (Fig 1).

#### **Pathway Enrichment Analysis**
We conducted pathway enrichment analysis to determine the biological roles of the upregulated and downregulated genes. This included analyzing the molecular functions, cellular localization, and biochemical pathways in which these genes are involved.


#### **Machine Learning**

__Data Splitting__

Due to class imbalance in the dataset, we allocated 25% of the majority class and one sample from the minority class to the test set, while the remaining samples were used for training.

__Feature Selection__

- __Filtering:__ Genes with FDR > 0.05 were removed, reducing the number of features from 30,753 to 4,177.
- __Lasso Logistic Regression:__ A Lasso regression model was used to select features by applying a penalty which shrank unimportant features to exactly zero. A total of 12 genes selected for modeling based on non-zero coefficient values (Fig 6).

__Modeling__

A k-nearest neighbors (k-NN) algorithm (k=5) was used to train a model for predicting cervical tissue types. The model's performance was evaluated using accuracy, recall, F1 score, precision, and specificity.

## **Results and Interpretation**

Fig 1 shows a volcano plot showing upregulated (red) and downregulated genes (green) at the chosen cutoffs.

![volcano_plot](imgs/volcano.png)

_**Fig 1**: Volcano plot showing upregulated, down regulated and not statistically significant genes_


### **Biological Processes**

The results showed that most of the upregulated genes were involved in the cell cycle pathway, which plays a significant role in the pathogenesis of cervical cancer. Other enriched processes included the muscular system, nervous system signaling, and cellular differentiation pathways (Fig 2).

Conversely, downregulated genes were mainly linked to striated muscle cell differentiation, muscle structure development, and processes such as muscular contraction, heart rate regulation, and synaptic signaling (Fig 3).

![upregulated_biological_process](enrichment/upgene_bio_proc.png)

_**Fig 2**: Biological Processes (Upregulated)_


![downregulated_biological_process](enrichment/downgene_bio_proc.png)

_**Fig 3**: Biological Processes (downregulated)_

More info about enriched pathways for down and up regulated genes can be found [here](enrichment/Down%20regulated%20enrichment%20Biological%20Processes.csv) and [here](enrichment/Upregulated%20enrichment%20GO%20Biological%20Process.csv)


### **Molecular Function**

- Upregulated genes: Involved in processes such as cell division, mitotic regulation, and differentiation (Fig 4).
- Downregulated genes: Associated with muscle development, contraction, and signaling pathways (Fig 5).

Figure 4 and 5 show the molecular functions these gene products play in the human organism. Other molecular functions for [downregulated genes](enrichment/Down%20regulated%20enrichment%20GO%20Molecular%20function.csv) and [upregulated genes](enrichment/Upregulated%20enrichment%20GO%20Molecular%20Function.csv)

![upregulated_molecular_function](enrichment/upgene_molec_func.png)

_**Fig 4**: Molecular Function (Upregulated)_

![downregulated_molecular_function](enrichment/downgene_molec_func.png)

_**Fig 5**: Molecular Function (Downregulated)_


### **Biochemical Pathways**

The top 5 biochemical pathways enriched by downregulated and upregulated genes are shown below.

**Table 1: KEGG Pathways (Downregulated genes)**

Enrichment FDR |nGenes	|Pathway Genes|Fold Enrichment |Pathway
--------------:|-------:|------------:|---------------:|:----------------
0.002725386	   |10	    |362	      |5.138795311	   |Neuroactive ligand-receptor interaction
0.008437578	   |3	    |23	          |24.2640509	   |Renin-angiotensin system
0.008437578	   |7	    |221	      |5.892175257	   |cAMP signaling pathway
0.03256974	   |5	    |149	      |6.242429203	   |Adrenergic signaling in cardiomyocytes
0.041900311	   |5	    |166	      |5.603144284	   |cGMP-PKG signaling pathway


**Table 2: KEGG Pathways (Upregulated genes)**

Enrichment FDR |nGenes	|Pathway Genes|Fold Enrichment |Pathway
--------------:|-------:|------------:|---------------:|:----------------
2.03e-10	   |16	    |126	      |11.13227513	   |Cell cycle
3.37e-07	   |15	    |186	      |7.069892473	   |Alcoholism
2.82e-06	   |12	    |135	      |7.792592593	   |Systemic lupus erythematosus
1.37e-05	   |11	    |131	      |7.361323155	   |Oocyte meiosis
5.38e-05	   |8	    |73	          |9.607305936	   |p53 signaling pathway
5.38e-05	   |12	    |189	      |5.566137566	   |Neutrophil extracellular trap formation
0.0027576	   |10	    |202	      |4.339933993	   |Viral carcinogenesis


### **Model Performance**

Model performance on the test set showed a perfect score in all metrics evaluated (Tables 3 & 4). The confusion matrix (Fig 7) shows that all test samples in both tumour and normal classes were correctly classified by the kNN model.

__Table 3: Model Performance (Train data)__

Metric      | Score (%)
------------|---------
Accuracy    | 100			
Recall      | 100			
F1          | 100			
Specificity	| 100			
Precision	| 100


__Table 4: Model Performance (Test data)__

Metric      | Score (%)
------------|---------
Accuracy    | 100			
Recall      | 100			
F1          | 100			
Specificity	| 100			
Precision	| 100

![confusion_matrix](imgs/knn_conf_mat.png)

_**Fig 7**: Confusion Matrix (test set)_

## **Conclusion and Future Research**

This study highlights the significant role of the cell cycle in cervical cancer pathogenesis, as revealed by the upregulated genes. The successful identification of 12 key genes via Lasso regression demonstrates the potential for accurate machine learning-based tissue classification. Future research can focus on validating these biomarkers in larger cohorts and exploring their potential as therapeutic targets. Additionally, further investigation into the pathways associated with downregulated genes could uncover novel insights into cervical cancer progression and treatment.

## **REFERENCES**
1. Balasubramaniam, S. D., Balakrishnan, V., Oon, C. E., & Kaur, G. (2019). Key Molecular Events in Cervical Cancer Development. Medicina, 55(7). https://doi.org/10.3390/MEDICINA55070384

2. Banzola, I., Mengus, C., Wyler, S., Hudolin, T., Manzella, G., Chiarugi, A., Boldorini, R., Sais, G., Schmidli, T. S., Chiffi, G., Bachmann, A., Sulser, T., Spagnoli, G. C., & Provenzano, M. (2018). Expression of indoleamine 2,3-dioxygenase induced by IFN-γ and TNF-α as potential biomarker of prostate cancer progression. Frontiers in Immunology, 9(MAY). https://doi.org/10.3389/fimmu.2018.01051

3. Bedell, S. L., Goldstein, L. S., Goldstein, A. R., & Goldstein, A. T. (2020). Cervical Cancer Screening: Past, Present, and Future. Sexual Medicine Reviews, 8(1), 28–37. https://doi.org/10.1016/j.sxmr.2019.09.005

4. Bray, F., Laversanne, M., Sung, H., Ferlay, J., Siegel, R. L., Soerjomataram, I., & Jemal, A. (2024). Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. CA: A Cancer Journal for Clinicians, 74(3), 229–263. https://doi.org/10.3322/caac.21834

5. Rosati, D., Palmieri, M., Brunelli, G., Morrione, A., Iannelli, F., Frullanti, E., & Giordano, A. (2024). Differential gene expression analysis pipelines and bioinformatic tools for the identification of specific biomarkers: A review. Computational and Structural Biotechnology Journal, 23(October 2023), 1154–1168. https://doi.org/10.1016/j.csbj.2024.02.018

6. Shaon, M. S. H., Karim, T., Shakil, M. S., & Hasan, M. Z. (2024). A comparative study of machine learning models with LASSO and SHAP feature selection for breast cancer prediction. Healthcare Analytics, 6(June). https://doi.org/10.1016/j.health.2024.100353


### __Group Members__

>- Chigozie Nkwocha
>- Chaima Ben Mohamed
>- Charlotte Chinwendu Iwuji
>- Opeyemi De Campos
>- Reem Atawia

