Conditional eQTL analysis
=======================================

Developed by: Sourya Bhattacharyya

La Jolla Institute for Immunology

Documents
-------------

For details of conditional eQTL analysis, see the following papers and documentation:

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6075455/

http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

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

