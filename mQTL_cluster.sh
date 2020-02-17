#!/bin/bash
#$ -cwd -V
#$ -pe smp 20
#$ -l h_vmem=40G
#$ -e ~/log
#$ -o ~/log

module add apps/R/3.4.0
Rscript --vanilla /home/b1017249/WORKING_DATA/matrixEQTL/meQTL_Run_All.R
