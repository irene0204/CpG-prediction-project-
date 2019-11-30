# Alzheimer's Disease-related CpG sites prediction project pipeline 
This repository contains all components of the pipeline for predicting novel Alzheimer's Disease (AD)-associated CpG sites beyond 450K array across the whole human genome, including experimental set construction, features collection/processing, features selection and ensemble learning, for each of the AD-associated trait of interest. 

## Tools
* Python 3.6
* R 3.4
* Amazon Elastic Compute Cloud (AWS EC2)

## Prerequites
The following input files are needed:
* CSV files with summary level data from the ROSMAP study. For each trait, the file includes CpG ID, F statistics and p-values (null hypothesis: AD samples have the same methylation level as control samples) for *CpG sites whose methylation level was measured using Illumina 450K array (i.e. 450K sites)*. `ROSMAP.csv`
* a TXT file with the whole human genome spread across 200 base-pair intervals `wins.txt`
* BED files with window IDs of *all CpG sites across the whole human genome (i.e. WGBS sites)* and values of the 1806 features used in our previously published work on DIVAN.  `DIVAN_features.bed`
* TSV.GZ files with genomic locations of WGBS sites and CADD scores `CADD.tsv.gz`
* TSV.BGZ files with genomic locations of WGBS sites and DANN scores `DANN.tsv.bgz`
* TAB.BGZ files with genomic locations of WGBS sites and EIGEN scores `EIGEN.tab.bgz`
* BED.GZ files with genomic locations of WGBS sites and GWAVA scores `GWAVA.bed.gz`
* BED files with window IDs of WGBS sites and RNA-sequencing read counts data `RNASEQ.bed`
* BED files with window IDs of WGBS sites and ATAC-sequencing read counts data `ATACSEQ.bed`
* BED files with genomic locations of WGBS sites and WGBS read counts data `wgbs_readcounts.bed`
* a TXT file with genomic locations of transcription start sites (tss) `tss.txt`



## Running the pipeline 

**0) Preliminary step: hg38 to hg19 conversion**

The genomic coordinates in this pipeline are 1-indexed/in hg19. The original WGBS datasets are 0-indexed/in hg38 and therefore need to be converted. This conversion can be completed by running the WGBS_allsites_preprocess.py script. The number of studied WGBS sites is reduced from `approximatly 28 million` to `approximatly 26 million` after this conversion due to using LiftOver (inconsistent chromosome, multiple conversion results, etc). 

``` 
WGBS_allsites_preprocess.py ${wgbs_readcounts.bed} ${wins.txt}
```
The file `${wgbs_readcounts.bed}`contains the 0-indexed/hg38 genomic locations of WGBS sites.

This step generates `all_wgbs_sites_winid.csv`,which contains the genomic locations in both hg38 and hg19 and window IDs for WGBS sites. 




**1) Asssign feature values to WGBS sites**

To prepare for future prediction of AD-associated WGBS sites, we first assign all (2256) feature values to WGBS sites. 

By running the WGBS_all_sites_feature_preprocess.py script, we can processes features of WGBS sites in 14 batches of 2 million for the consideration of memory limit. 

``` 
WGBS_all_sites_feature_preprocess.py ${all_wgbs_sites_winid.csv} ${DIVAN_features.bed} ${CADD.tsv.gz} \
    ${DANN.tsv.bgz} ${EIGEN.tab.bgz} ${GWAVA.bed.gz} ${RNASEQ.bed} ${ATACSEQ.bed} \
    ${wgbs_readcounts.bed} ${tss.txt}
```

The features are processed as follows:

* 1806 features in the DIVAN study are constructed to cover the entire human genome in 200 base-pair resolution and assigned to each site by matching the window ID
* CADD, DANN, EIGEN, GenoCanyon and GWAVA scores are assigned to each site by matching genomic location
* RNA-sequencing, ATAC-sequencing and WGBS readcounts data are assigned to each site by matching window ID
* Distance to the nearest tss are found for each site
* Summary of all features 

  | Feature source         | Number   | 
  | -------------          |:--------:| 
  | REMC DNase             | 73       |
  | REMC Histone           | 735      | 
  | ENCODE DNase           | 80       |   
  | ENCODE FAIRE           | 31       |   
  | ENCODE TF(HAIB)        | 292      |   
  | ENCODE TF(SYDH)        | 279      |   
  | ENCODE Histone         | 267      |   
  | ENCODE RNA Polymerase  | 49       |   
  | ENCODE RNA-seq         | 243      |   
  | ENCODE ATAC-seq        | 66       |   
  | GenoCaynon             | 1        |   
  | Eigen                  | 4        |   
  | DANN                   | 2        |   
  | CADD                   | 2        |   
  | GWAVA                  | 4        | 
  | Distance to nearest TSS| 1        | 
  | WGBS                   | 127      | 
  | Total                  | 2256     | 
  
  


This step generates HDF5 files for all batches of WGBS sites and their feature values:
``` 
all_features_0_2000000.h5
....
```

**2) Assign feature values to all 450K sites**

By running the all450k_feature_preprocess.py script, we assign all (2256) feature values to 450K sites.

``` 
all450k_feature_preprocess.py ${all_450k_sites_winid.csv} ${DIVAN_features.bed} ${CADD.tsv.gz} \
    ${DANN.tsv.bgz} ${EIGEN.tab.bgz} ${GWAVA.bed.gz} ${RNASEQ.bed} ${ATACSEQ.bed} \
    ${wgbs_readcounts.bed} ${tss.txt}
```
This step generates a HDF5 file for 450K sites and their feature values:
``` 
all_450k_features.h5
```

**3) Experimental set construction for each trait**

For furture model training purpose, we constructed a experimental set for each trait. In each set, we include positive sites (signficantly associated with AD) and negative sites (not significantly associated with AD). The inclusion criteria are as follows:

a) Select positive sites whose p-values are below trait-specific threshold

b) For each selected positive site, select 10 negative sites that:

* have p-values greater than 0.4
* have the same methylation status (either hyper- and hypo-) as the positive site
* have the closest β-values as the positive site 

The above process is achieved by running the AD_sites_selection.py script. 

``` 
AD_sites_selection.py ${ROSMAP.csv} ${wins.txt} \
    --amyloid_positive_threshold 0.00005 \
    --cerad_positive_threshold 0.00001 \
    --ceradaf_positive_threshold 0.00005 \
    --tangles_positive_threshold 0.0000005 \
    --cogdec_positive_threshold 0.00003 \
    --gpath_positive_threshold 0.00001 \
    --braak_positive_threshold 0.00005       
```

This step outputs 7 CSV files for 7 traits: 
``` 
all_sites_winid.csv
```
which contains the selected experimental set with columns: CpG ID, chromosome, coordinate, p-value, β-value, label (0 for negative sites/1 for positive sites) and window ID. 

and 1 CSV file for all 450K sites:
``` 
all_450k_sites_winid.csv
```
which contains all 450k sites with columns: CpG ID, chromosome, coordinate, p-value, β-value and window ID. 


**4) Assign feature values to the experimental set for each trait**

By running the all_features_preprocess.py script, we assign feature values to the constructed experimental dataset for each trait.

``` 
all_features_preprocess.py ${all_sites_winid.csv} ${DIVAN_features.bed} ${CADD.tsv.gz} \
    ${DANN.tsv.bgz} ${EIGEN.tab.bgz} ${GWAVA.bed.gz} ${RNASEQ.bed} ${ATACSEQ.bed} \
    ${wgbs_readcounts.bed} ${tss.txt}
```


This step generates 7 HDF5 files for 7 traits:
``` 
all_features.h5
```
which contains all feature values of the experimental set for each trait.


**5) Features selection for each trait**

Considering the number of features is greater than the number of CpG sites in the experimental set, we need to perform features selection for each trait to avoid overfitting. 

The feature selection process is achived by running the feature_selection.py script. 

``` 
feature_selection.py ${all_features.h5}
```
The features are selected as follows for each trait: 

* Split training/testing data on 9:1 ratio 
* Select top 100 significant features by fitting the training data using random forest, xgboost, logistic regression and support vector classifier (SVC) with linear kernel, respectively. 
* For each selected feature, summarize the number of classifiers that select it, n (0<n≤4) 
* Keep features with n≥2
* For each of the kept feature, perform Wilcoxon rank-sum test and calculate the p-values under the null hypothesis: AD samples have the same feature values as control samples
* From top to bottom, sort selected features first by desceding n and then by ascending p-value 
* Select top ranked features for each trait 

This step outputs a CSV file and a HDF5 file for each trait:

1)`feature_stats.csv`, which contains information of the top ranked feaures, including feature name, p-value, and n

2)`selected_features.h5`, which contains the training and testing set with the values of top ranked features assigned and the labels for training and testing set.

**6) Model hyper-parameters tuning and model selection for each trait**

We use 4 base classifiers, random forest, xgboost, logistic regression with L2-regularization and SVC with linear kernel.The optimal  hyper-paramters for each classifier and best combination of base classifiers are selected by running the ModelSelectionTuning.py script.

``` 
ModelSelectionTuning.py ${selected_features.h5}
```

The optimal hyper-parameters for each base classifier and and the best combination of base classifiers are selected as follows:

* Use the training set (generated from step 5) with 3-fold cross-validation to select the optimal hyper-parameters for each base classifier
* Use the training/testing set (generated from step 5) with 10-fold cross-validation to evaluate all possible combinations of base classifiers
* Calculates average AUC, F1-score, etc across all 10 folds 
* Select the best combination of base classifiers 


**7) Predict 450K sites and WGBS sites for each trait**

We use the selected ensemble model from step 6 to make predictions on CpG sites within and beyond the 450K array, which is achieved by running the WGBS_prediction.py script.


``` 
WGBS_prediction.py ${selected_features.h5} ${all_features_0_2000000.h5} ${all_450k_features.h5}
```

In this step, we need to retrain the base classifiers in the ensemble model with the entire experimental set to obtain the optimal hyper-parameters and save the retrained model. However if the retrain process has already been done, we load the saved retrained model directly. Subsequently we use the retrained model to predict the probabilities of 450K sites/WGBS sites being positive. CpG sites are ranked in descendig orders of these probabilities and top 500 sites are selected as candidates for experimental validation. 

This step generates 2 CSV files, a HDF5 file and a pkl file for each trait:

1)`pred_positive_500.csv`, which contains the top 500 sites and their probabilities of being positive; 

2)`top500_nearest_450k.csv`, which contains the 450K sites within 5k up/downstream of top 500 predicted sites with their CpG ID and genomic location;

3)`pred_probs.h5`, which all WGBS sites and their probabilities of being positive; 

4)the saved retrained model in `prediction_model.pkl`;

and a HDF5 file `pred_probs_450k.h5`, which contains all 450K sites and their probabilities of being positive,


**8) Combine results for all AD traits**

Since we have 7 models trained from datasets of 7 traits, we take both the average and the weighted average of the 7 predicted probabilities as the predicted probability of a CpG site being AD-associated. The weights for each trait is determined based on the F1 score. 

The above process can be achieved by running WGBS_alltraits_prediction_AD.py.

``` 
WGBS_alltraits_prediction_AD.py ${10fold_test_results.pkl} ${pred_probs.h5} \
        ${pred_probs_450k.h5} ${tss.txt}
```
The file `${10fold_test_results.pkl}` is an intermediate output from step 6), which contains F1 score of the ensemble model for each trait. 

This step generates 2 CSV files for WGBS sites:

1)`common_top500_mean_nearest_450k.csv`, which contains the top 500 sites with highest unweighted average of 7 probabilities, and their distances to the nearest tss; 

2)`common_top500_weighted_nearest_450k.csv`,  which contains the top 500 sites with highest weighted average of 7 probabilities, and their distances to the nearest tss; 

and a csv file for 450K sites:
`450kwithpredictedprob.csv`, which contains the unweighted and weighted average of 7 probabilities of all 450K sites. 
