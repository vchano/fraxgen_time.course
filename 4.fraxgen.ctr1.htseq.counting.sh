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
#SBATCH --job-name=HTS1.FRAX.CTR1
#SBATCH --output=HTSEQ_CTR1_FRAXGEN_%J_out.txt
#SBATCH --error=HTSEQ_CTR1_FRAXGEN_%J_err.txt
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
#module purge 
#module load htseq/0.6.1
#module load anaconda3/2020.11
##module load python/3.8.6
##export PATH=$HOME/.conda/envs/pysam/bin:$PATH
#source activate pysam

# Set PATHS to working directories
INPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP1_CTR
OUTPUT=$SCRATCH/FRAXGEN/4.COUNTING/EXP1_CTR
REF=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3

# Running HTSEQ and processing count files 

list='p1258s49 p1258s50 p1258s51 p1258s52 p1258s57 p1258s58 p1258s59 p1258s60 p1258s61 p1258s68 p1258s69 p1258s71 p1258s72 p1258s79 p1258s80 p1258s81 p1258s82 p1258s87 p1258s88 p1258s89 p1258s90 p1258s96 p1258s97 p1258s98 p1258s99 p1258s106 p1258s108 p1258s109 p1258s110 p1258s118 p1258s119 p1258s120 p1258s121 p1258s127 p1258s128 p1258s129 p1258s130 p1258s135 p1258s136 p1258s137'

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT/${sample}_ht2.fix.psrt.bam $REF > $OUTPUT/${sample}_counts.v2.txt
done

paste $OUTPUT/*_counts.v2.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=$25=$27=$29=$31=$33=$35=$37=$39=$41=$43=$45=$47=$49=$51=$53=$55=$57=$59=$61=$63=$65=$67=$69=$71=$73=$75=$77=$79=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/ctr1_count_table.v2.txt
cp -R $OUTPUT/ctr1_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/



#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo '** Thats all Folks! **'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/4.fraxgen.ctr1.htseq.counting.sh
#-------------------------------------------------------------------------------