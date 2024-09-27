#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#                   Script for running HTSEQ                     #
#                                                                #
##################################################################

#-------------------------------------------------------------------------------
#
#    Forest Genetics and Forest Tree Breeding
#    Büsgenweg 2
#    Georg-August-Universtät Göttingen
#    https://github.com/vchano/
#
# Licence: GNU General Public Licence Version 3
#-------------------------------------------------------------------------------

#SBATCH -p gailing
#SBATC -A all
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --job-name=HTS1.FRAX.TEST
#SBATCH --output=HTSEQ_TEST_FRAXGEN_%J_out.txt
#SBATCH --error=HTSEQ_TEST_FRAXGEN_%J_err.txt
#SBATCH --ntasks-per-socket 24
#SBATCH --qos long
#SBATCH --time=120:00:00
#SBATC --begin=2022-11-25T01:00:00
#SBATC --exclusive
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=vmchano.gaug@gmail.com

# Load modules

module purge
export PATH="/opt/sw/rev/21.12/haswell/gcc-9.3.0/anaconda3-2021.05-e3srav/bin:$PATH"
export PATH="/opt/sw/rev/21.12/haswell/gcc-9.3.0/anaconda3-2021.05-e3srav/bin/conda:$PATH"
module load anaconda3/2021.05
source activate htseq

#module load htseq/2.0.2
#module load anaconda3/2020.11
##module load python/3.8.6
##export PATH=$HOME/.conda/envs/pysam/bin:$PATH
#source activate pysam

# Set PATHS to working directories
INPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/TEST
OUTPUT=$SCRATCH/FRAXGEN/4.COUNTING/TEST
REF=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3

# Running HTSEQ and processing count files 

list='p1194sResistantFar3Control-1 p1194sResistantFar3Control-2 p1194sResistantFar3Control-3 p1194sResistantFar3Inoculated-4 p1194sResistantFar3Inoculated-5 p1194sResistantFar3Inoculated-6 p1194sSusceptibleUW1Control-7 p1194sSusceptibleUW1Control-8 p1194sSusceptibleUW1Control-9 p1194sSusceptibleUW1Inoculated-10 p1194sSusceptibleUW1Inoculated-11 p1194sSusceptibleUW1Inoculated-12'

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT/${sample}_ht2.fix.psrt.bam $REF > $OUTPUT/${sample}_counts.v2.txt
done

paste $OUTPUT/*_counts.v2.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/test_count_table.v2.txt
cp -R $OUTPUT/test_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/

paste $OUTPUT/*_counts.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/test_count_table.txt
cp -R $OUTPUT/test_count_table.txt $HOME/FRAXGEN/4.FRAXEN.GEA/

#-------------------------------------------------------------------------------

# End
echo '****************************************************************'
echo 'Mapping with HISAT2 of FRAXGEN test reads to Reference'
echo 'That´s all Folks!'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/4.fraxgen.test.htseq.counting.sh
#-------------------------------------------------------------------------------