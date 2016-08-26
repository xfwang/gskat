# The "gskat" Package

This package implements a family based association test via GEE Kernel Machine (SKAT) score test. It has functions to test an association between SNP/SNV sets and binary/continuous phenotypes.

References: 
-------------------------
Wang X. et al. GEE-Based SNP Set Association Test for Continuous and Discrete Traits in Family-Based Association Studies. Genetic Epidemiology (2013) 37: 778-786

Wang X. et al. Rare variant association test in family based sequencing studies. Breifings in Bioinformatics (In press).


Installation Instructions
-------------------------
Start R and then install the 'gskat' package with the following commands.
```
#install.packages("devtools")
library(devtools)
install_github("gskat",username="xfwang")
```


-------------------------
zhang.zhenyu@stonybrook.edu

xuefeng.wang@stonybrook.edu



-------------------------
Q&A:
Gskat now have function for testing without any covariates.
A new function called gskat_seq_opt2 has just been added for testing without covaraites. 

Does gskat provides a function to do null model fitting to speed up the large scale test?
We are working on this in the new package.

-------------------------
Updates:



