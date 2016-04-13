# The "gskat" Package

This package implements a family based association test via GEE Kernel Machine (SKAT) score test. It has functions to test an association between SNP/SNV sets and binary/continuous phenotypes.


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
Q and A:
Does gskat have function for testing without any covariates?
A new function called gskat_seq_opt2 has just been added for testing without covaraites. 

Does gskat provides a function to do null model fitting to speed up the large scale test?
We are working on this in the new package.



