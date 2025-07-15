**This folder contains an example R script and data files for the study titled "Effects of Oversampling on Ecological Niche Modeling Using Imbalanced Data."**



---------------------------------------------------------------------------------------------------------------

This folder contains example R scripts developed and modified by the author for the present study, along with the corresponding sample data files. 

Descriptions of each R script (C1–C5; format: R) and sample data file (formats: CSV and RData) are provided below.

\# (format) file name.



**# (R) C1. Data preprocessing function**

-1.1. Outlier detection using Hampel filter

-1.2. Variable selection

-1.3. Environmentally filter 

-1.4. Data partition



**# (R) C2. Oversampling function**

\- Sampling With Replacement (SWR), SMOTE, BLSMOTE, ADASYN and DBSMOTE are available for use.



**# (R) C3. Statistical analysis function**

\- Data arrangement function for statistical test (Common function)

\- 3.1. MRPP tests with restricted permutation (multivariate test)

\- 3.2. Restricted permutation test (univariate test)



**# (R) C4. Ecological niche modeling function**  

\- MaxEnt, RFC, RFR, SVC and SVR are available for use.



**# (R) C5. Sensitivity analysis**

\- 5.1. Variable importance 

\- 5.2. Variable importance-based distance between models



**# (CSV) Sample data.C1\_D1.raw data** 

\- Example of collected biological-environmental data

\- Used in code C1



**# (CSV) Sample data.C1\_D2.sp1\_Preprocessed data**

\- Example of the file after preprocessing in code C1

\- Focused on Baetis fuscatus (i.e., sp1) among five species

\- Used in section 1.4. 'Data partition' in C1



**# (RData) Sample data.C2\_sp1\_p0.20-after partition**

\- Example of the preprocessed file from code C1

\- Presence proportion adjusted to 0.20 (20%) (i.e., p0.20)

\- Used in code C2



**# (RData) Sample data.C3-4\_sp1\_p0.20\_after oversampling**

\- Example of the dataset after oversampling in code C2

\- SMOTE applied

\- Two repetition conducted (i.e., k1 and k2)

\- Used in codes C3 and C4



**# (RData) Sample data.C5\_ENM result\_Max**

\- Output generated through code C4

\- Contains data used for ENM development, the developed ENM, evaluation results, and parameter settings

\- Maximum entropy model applied (i.e., Max)

\- Used in code C5



**# (CSV) Sample data.C5\_sp1\_0.20\_Max**

\- Results of variable importance analysis using ENM

\- Includes ENMs without oversampling, and with SMOTE and ADASYN applied

\- Contains variable importance values from five repeated ENM runs (i.e., k1–k5)

\- Used in section 5.2. 'Variable importance-based distance between models' in C5





