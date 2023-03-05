Conditional analysis default code
=======================================

The lead SNP in each iteration will be used as the covariate in the next iteration, 
until we do not find any significant SNP.

Execution
-------------

Execute the script 
Conditional_eQTL.sh

Data
--------

See the folder "Data" for sample dataset.
For genotype and SNPs, we have provided information for only chromosome 1.

User needs to check the DICE paper (Schmiedel et al. Cell 2018) regarding processing gene expression, genotype and covariates, as an input to eQTL studies.

Output
----------

Outputs for individual iterations are generated within the folders *Iteration_\** under the specified output folder.

