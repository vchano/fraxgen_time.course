#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#              Script for running FastQC in rawdata              #
#                                                                #
##################################################################

#-------------------------------------------------------------------------------
#
#    Forest Genetics and Forest Tree Breeding
#    Büsgenweg 
#    Georg-August-Universtät Göttingen
#    https://github.com/vchano/
#
# Licence: GNU General Public Licence Version 3
#-------------------------------------------------------------------------------

#SBATCH -p gailing
#SBATC -A all
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=FRAX.TEST
#SBATCH --output=FASTQC.TEST_%J_out.txt
#SBATCH --error=FASTQC.TEST_%J_err.txt
#SBATC --exclusive
#SBATCH --ntasks-per-socket 24
#SBATC --qos long
#SBATC --begin=2021-12-03T23:30:00
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=vmchano.gaug@gmail.com


# Load modules

module load fastqc/0.11.4
module load anaconda3/2020.11

# Set PATHS to working directories
INPUT=$SCRATCH/FRAXGEN/0.RAWDATA/TEST
OUTPUT=$SCRATCH/FRAXGEN/1.QC/RAWDATA_TEST

# Run FastQC

fastqc $INPUT/*.fastq.gz --outdir $OUTPUT --threads 48

# Run MultiQC
source activate multiqc
multiqc $OUTPUT/
conda deactivate

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo 'FastQC test and MultiQC of FRAXGEN test samples'
echo 'That´s all Folks!'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/1.fraxgen.test.fastqc.rawdata.sh
#-------------------------------------------------------------------------------
