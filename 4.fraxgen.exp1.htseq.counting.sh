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
#SBATCH --job-name=HTS2.FRAX.EXP1
#SBATCH --output=HTSEQ_EXP1_FRAXGEN_%J_out.txt
#SBATCH --error=HTSEQ_EXP1_FRAXGEN_%J_err.txt
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
INPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP1_EXP
OUTPUT=$SCRATCH/FRAXGEN/4.COUNTING/EXP1_EXP
REF=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3


# Running HTSEQ and processing count files 

list='p1258s53 p1258s54 p1258s55 p1258s56 p1258s62 p1258s63 p1258s64 p1258s67 p1258s73 p1258s74 p1258s75 p1258s77 p1258s78 p1258s83 p1258s84 p1258s85 p1258s86 p1258s91 p1258s92 p1258s93 p1258s94 p1258s100 p1258s103 p1258s104 p1258s105 p1258s111 p1258s114 p1258s116 p1258s122 p1258s124 p1258s125 p1258s126 p1258s131 p1258s132 p1258s133 p1258s134 p1258s139 p1258s141 p1258s142 p1258s143'

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT/${sample}_ht2.fix.psrt.bam $REF > $OUTPUT/${sample}_counts.v2.txt
done

paste $OUTPUT/*_counts.v2.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=$25=$27=$29=$31=$33=$35=$37=$39=$41=$43=$45=$47=$49=$51=$53=$55=$57=$59=$61=$63=$65=$67=$69=$71=$73=$75=$77=$79=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/exp1_count_table.v2.txt
cp -R $OUTPUT/exp1_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/


#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo '** Thats all Folks! **'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/4.fraxgen.exp1.htseq.counting.sh
#-------------------------------------------------------------------------------