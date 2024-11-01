## **Classifying IDH Status of Lower-Grade Gliomas Using Gene Expression Data and a Machine Learning Approach**

> ## **Group Members**
> Chigozie Nkwocha \
> Chaima Ben Mohamed \
> Charlotte Chinwendu Iwuji \
> Opeyemi De Campos \
> Reem Atawia

 ## 1. **Introduction to Gliomas, IDH Staus and their Significance**
Gliomas are the most common primary tumors found in the brain and spinal cord (Chen et al. 2017). They account for about 80% of all brain tumors (Li et al. 2022). Diffuse gliomas are a type of brain tumor that originates from the glial cells, which support and protect neurons (Neumaier, Zlatopolskiy, and Neumaier 2023; Yang et al. 2022). Adult diffuse gliomas are classified and graded based on histological features, including subtypes such as oligodendroglioma, oligoastrocytoma, astrocytoma, and glioblastoma, ranging from grade II to IV (Ceccarelli et al. 2016).
They are also usually classified based on the type of glial cells involved (e.g., astrocytes in astrocytomas) and their genetic mutations, particularly in the IDH gene(Louis et al. 2021). The IDH (Isocitrate Dehydrogenase) status is crucial for classifying diffuse gliomas, particularly astrocytomas and oligodendrogliomas. IDH mutations help distinguish glioma subtypes in the WHO classification, dividing tumors into IDH-mutant and IDH wild-type (Louis et al. 2021). Gliomas with IDH mutations generally have a better prognosis due to slower growth and improved survival. These mutations also produce an oncometabolite, altering cellular metabolism and driving tumorigenesis. Additionally, IDH mutations are potential therapeutic targets, with inhibitors being explored in clinical trials to offer new treatments (Louis et al. 2021).

## 2. **Methodology**
### 2.1 **Dataset Description and Preprocessing Steps**

Expression data for lower grade gliomas (LGG) cancer data for this study was sourced from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/ccg/research/genome-sequencing/tcga). The dataset comprised 534 samples: 516 primary tumor cases and 18 recurrent cases. Only primary cases were selected.

### 2.2 **Data Preparation and Preprocessing**

Data preprocessing and normalization were conducted using the TCGAbiolinks library in R. The steps included

- **Filtering genes:** Genes with a correlation less than 0.6 across samples were removed.
- **Normalization:** Expression levels were normalized based on gene length and GC content.
- **Exclusion of low-expression genes:** Genes with expression levels below the 25th percentile were excluded.

### 2.3 **Biomarker Discovery and Machine Learning**

#### 2.3.1 **Differential Gene Expression (DGE) Analysis**

We performed DGE analysis using the likelihood ratio method, identifying differentially expressed genes based on the criteria: False discovery rate (FDR) < 0.05 and Absolute log fold change > 2. Genes with FDR < 0.05 and LogFC < -2 were classified as downregulated, while genes with FDR < 0.05 and LogFC > 2 were classified as upregulated. The preprocessed data, consisting of about 26,000 genes, was used for subsequent machine learning predictions (Fig 1). DGE analysis was based on the tumor grade (Grade 2 vs Grade 3) and IDH status (Wild type vs Mutant).

#### 2.3.2 **Pathway Enrichment Analysis**
We conducted pathway enrichment analysis to determine the biological roles of the upregulated and downregulated genes. This included analyzing the molecular functions, cellular localization, and biochemical pathways in which these genes are involved.

#### 2.3.3 **K-Means Clustering**
We employed K-means clustering, which is an unsupervised machine learning algorithm that generates clusters using the cluster’s object mean value (Ikotun et al., 2023). The algorithm generated 4 clusters to group the data into distinct clusters based on gene expression profiles. We cross-referred the identified clusters with the metadata to assess the corresponding samples' IDH status (Wild Type or Mutant).

![Figure 1](https://github.com/Chygos/hackbio-cancer-internship/blob/main/stage4/stage4_report/imgs/K-means%20Clustering.jpg?raw=true)

__Figure 1:__ K-means Clusters

The clustering algorithm has grouped samples based on their gene expression profiles, with a distinction in IDH status within the clusters. However, the overlap between clusters indicates that the algorithm did not distinctly separate the samples.  Additionally, samples from both the Mutant and Wild Type groups were present across all four clusters which shows that the clustering failed to classify Wild Type and Mutant samples as anticipated. This result contrasts with the findings from a related study.

#### 2.3.34 **Machine Learning**

__Data Splitting__

The preprocessed dataset was split into train and test sets, with 75% for model development and the other 25% for model evaluation.

__Feature Selection__

- __Filtering:__ Statistically insignificant genes from DGE analysis (FDR > 0.05) were removed reducing the features to about 24,000 genes.
- __Lasso Logistic Regression:__ A Lasso regression model was used to select features by applying a penalty which shrank unimportant features to exactly zero. Genes with non-zero coefficient values were selected for modelling.

__Modeling__

The k-nearest neighbors (k-NN) and random forest models were used to train a model for predicting IDH status and tumor grades, respectively. The model's performance was evaluated using accuracy, recall, F1 score, precision, and specificity. To select optimal parameters, hyperparameter tuning with a 5-fold cross-validation was used and the set of parameters with the highest accuracy (IDH status) and F1 score (tumor grade) was selected. A confusion matrix representing each model's predictions on the test set and the actual class was visualized.

## 3. **Results**

### 3.1 **IDH Status**

___Table 1: Model Performance on test set___

Model        | Accuracy | Precision | Recall | F1  | Specificity
:------------|---------:|----------:|-------:|----:|-----------:
KNN          |99.2	     |100        |99.04   |99.52| 100
Random Forest|99.2	     |100        |99.04   |99.52| 100

![confmat_IDH](imgs/confmat_IDH.png)

__Figure 2:__ Confusion matrix on the test set (Tumour Grade)

From Figure 2, both models could almost perfectly predict IDH status based on a set of genes.

### 3.2 **Tumor Grade**
___Table 2: Model Performance on test set___

Model        | Accuracy | Precision | Recall  | F1   | Specificity
:------------|---------:|----------:|--------:|-----:|-----------:
KNN          |61.95	    |68.09      |53.33    |59.81 | 71.70
Random Forest|63.72	    |60.34      |66.72    |63.06 | 61.67


![confmat_grade](imgs/confmat_grade.png)

__Figure 3:__ Confusion matrix on the test set (Tumour Grade)

From Table 2 and Figure 3, both models have problems distinguishing both tumor grades. We further examined the number of samples in each tumor type that the models in the test data correctly detected. The result can be found in Figure 4. Both models had problems differentiating tumor grades in the oligodendroglioma tumor type than in any other tumor type.

![tumour_grade_type_vs_actual](imgs/models_tumour_grades_types.png)

__Figure 4:__ Models Predictions of Tumor types and Grades vs Actual values

### **Results for Gene Enrichment Analysis**

The gene enrichment analysis of IDH-mutant genes revealed significant upregulation of immune-related biological processes, such as the **regulation of immune system processes** and **response to bacterial molecules**, highlighting the immune system's role in IDH mutation biology. In terms of cellular components, IDH-related proteins were enriched in the **extracellular space** and **plasma membrane**, suggesting their involvement in membrane signaling and extracellular interactions. The analysis also identified **oligosaccharide binding** and **ethanol binding** as key molecular functions. For downregulated IDH analysis, there were "no significant genes" for KEGG, Biological processes, Cellular components, and Molecular function.

![IDH Upregulated Biological Processes](https://github.com/Chygos/hackbio-cancer-internship/blob/main/stage4/stage4_report/imgs/IDH_BPP_UP.png?raw=true)

__Figure 5:__ IDH Upregulated Biological processes

![IDH Upregulated Cellular Localisation](https://github.com/user-attachments/assets/40a82e49-2645-4e8a-9989-12c4f38739c1)


__Figure 6:__ IDH Upregulated Cellular localisation


![IDH Molecular Function](https://github.com/Chygos/hackbio-cancer-internship/blob/main/stage4/stage4_report/imgs/IDH_MFP_UP.png?raw=true)

__Figure 7__ IDH Upregulated Pathwway for Molecular Function


![Upregulated genes for all pathways](https://github.com/user-attachments/assets/66b2d2a6-83d2-4c0f-ae16-45e3faf52de4)

__Figure 8__: Upregulated genes for all Pathways

## **Conclusion**

While the models used in this study performed exceptionally well in distinguishing IDH-mutant from IDH wild-type samples, achieving high accuracy, precision, and specificity, they struggled with accurately classifying tumor grades, highlighting challenges in differentiating between grade II and grade III gliomas. These results underscore the importance of IDH status in glioma classification and the potential for improving tumor grade prediction with further refinement of feature selection methods and model training. Continued research into the integration of molecular biomarkers, alongside advanced machine learning techniques, could enhance glioma diagnosis and prognosis.

## **REFERENCES**

- Ceccarelli, Michele et al. 2016. “Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.” Cell 164(3): 550–63.

- Chen, Ricky, Matthew Smith-Cohn, Adam L. Cohen, and Howard Colman. 2017. “Glioma Subclassifications and Their Clinical Significance.” Neurotherapeutics 14(2): 284–97.

- Li, Yurong, Qin Qin, Yumeng Zhang, and Yuandong Cao. 2022. “Noninvasive Determination of the IDH Status of Gliomas Using MRI and MRI-Based Radiomics: Impact on Diagnosis and Prognosis.” Current Oncology 29(10): 6893–6907.

- Louis, David N. et al. 2021. “The 2021 WHO Classification of Tumors of the Central Nervous System: A Summary.” Neuro-oncology 23(8): 1231–51. https://pubmed.ncbi.nlm.nih.gov/34185076/ (October 6, 2024).

- Neumaier, Felix, Boris D. Zlatopolskiy, and Bernd Neumaier. 2023. “Mutated Isocitrate Dehydrogenase (MIDH) as Target for PET Imaging in Gliomas.” Molecules 28(7).

- Rosati, Diletta et al. 2024. “Differential Gene Expression Analysis Pipelines and Bioinformatic Tools for the Identification of Specific Biomarkers: A Review.” Computational and Structural Biotechnology Journal 23(October 2023): 1154–68. https://doi.org/10.1016/j.csbj.2024.02.018.

- Shaon, Md Shazzad Hossain, Tasmin Karim, Md Shahriar Shakil, and Md Zahid Hasan. 2024. “A Comparative Study of Machine Learning Models with LASSO and SHAP Feature Selection for Breast Cancer Prediction.” Healthcare Analytics 6(June).

- Yang, Keyang et al. 2022. “Glioma Targeted Therapy: Insight into Future of Molecular Approaches.” Molecular Cancer 21(1). https://doi.org/10.1186/s12943-022-01513-z.
- Ikotun, A.M. et al. (2023) “K-means clustering algorithms: A comprehensive review, variants analysis, and advances in the era of big data,” Information sciences, 622, pp. 178–210. Available at: https://doi.org/10.1016/j.ins.2022.11.139.
