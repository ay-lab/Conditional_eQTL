#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=40GB
#PBS -l walltime=20:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

CodeExec='Conditional_eQTL.R'

covariatefile='Data/genotypes.clean.eigenvec'
geneexprfile='Data/GeneExpr.txt'
logGeneExpr=0	## gene expression values are not log transformed
genotypefile='Data/snps_chr1.txt'
snplocfile='Data/snpsloc1.txt'
genelocfile='Data/geneloc_hg37.txt'

Rscript ${CodeExec} --GenotypeFile ${genotypefile} --ExprFile ${geneexprfile} --SNPLocFile ${snplocfile} --GeneLocFile ${genelocfile} --OutDir ${MatrixeQTL_OutDir} --CovariateFile ${covariatefile} --logGeneExpr ${logGeneExpr}

