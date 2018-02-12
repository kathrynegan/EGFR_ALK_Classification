# EGFR_ALK_Classification
===========================

Updated by kathrynegan

A hybrid model to classify EGFR and ALK test use, results, and methods from free text pathology reports 

### This system classifies EGFR and ALK molecular tests from the free text of pathology reports.
_None of the SVM models needed for end to end classification are included in this repository_

Dependencies = python 2.7.13; scipy and sklearn 0.18.1

#### There are four principal algorithms
-------------------------------------
each needs its own folder to house svm models and feature sets

* __reported__: this is a combination of a keyword filter and an SVM that classifies all reports in the input as either "Reported" or "NotReported" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Reported" is considered a positive classification and "NotReported" a negative)

* __insufficient__: this is an SVM that classifies all reports classified as "NotReported" by the previous reported algorithm as either "Insufficient" or "Unknown" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Insufficient" is considered a positive classification and "Unknown" a negative)

* __positive__: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either "Positive" or "Negative" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, "Positive" is considered a positive classification and "Negative" a negative)

* __method__: this is an SVM that classifies all reports classified as "Reported" by the previous reported algorithm as either the standard testing methodology or "OTHER" for each test
    (for the purposes of sensitivity/specificity metrics and for consumption by downstream algorithms, the standard testing methodology ("MutationalAnalysis" for EGFR and "FISH" for ALK) is considered a positive classification and "OTHER" a negative)
    

- run.py is the main script to run the end to end classification pipeline
- utils/gentest_classifier.py reads in models and feature mappings, uses the svm models learned in training to classify one instance at a time
- utils/vectorizer.py creates a vector for a given pathology report

Vector creation and classification pipeline are run for both EGFR and ALK tests


#### The References folder contains
Internal validation performance as well as a project overview is reported in the attached abstract (microsoft word doc) "Validation of Natural Language Processing (NLP) for Automated Ascertainment of EGFR and ALK Tests in SEER Cases of Non-Small Cell Lung Cancer (NSCLC)"
As well as classifcation_architecture.png; a visual represetation of the classification architecture
