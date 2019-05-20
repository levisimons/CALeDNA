#Script to run a zeta diversity analysis on the cluster.
#!/bin/sh
#SBATCH -t100:00:00
#SBATCH -n1
#SBATCH -c1
#SBATCH --mem=120000m
#SBATCH -pSCCWRP
cd /home/rcf-40/alsimons/cmb/CALeDNA
Rscript Zeta4FamilyMapCluster.R
