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

#SBATCH -p medium
#SBATCH -A all
#SBATCH -n 6
#SBATCH -N 1
#SBATCH --job-name=HTS4.FRAX.EXP2
#SBATCH --output=HTSEQ_EXP2_FRAXGEN_%J_out.txt
#SBATCH --error=HTSEQ_EXP2_FRAXGEN_%J_err.txt
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
INPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP2_EXP
OUTPUT=$SCRATCH/FRAXGEN/4.COUNTING/EXP2_EXP
REF=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3

# Running HTSEQ and processing count files 

list='p1258s150 p1258s151 p1258s154 p1258s155 p1258s163 p1258s164 p1258s166 p1258s168 p1258s177 p1258s178 p1258s179 p1258s180 p1258s188 p1258s189 p1258s191 p1258s192 p1258s199 p1258s200 p1258s202 p1258s203 p1258s211 p1258s214 p1258s217 p1258s225 p1258s226 p1258s227 p1258s229 p1258s236 p1258s237 p1258s238 p1258s249 p1258s251 p1258s252 p1258s253 p1258s260 p1258s261 p1258s262 p1258s264 p1258s274 p1258s275 p1258s276 p1258s277 p1258s285 p1258s286 p1258s287 p1258s290 p1258s297 p1258s298 p1258s299 p1258s300 p1258s309 p1258s310 p1258s312 p1258s315 p1258s323 p1258s324 p1258s325 p1258s326 p1258s334 p1258s335 p1258s336 p1258s338'

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT/${sample}_ht2.fix.psrt.bam $REF > $OUTPUT/${sample}_counts.v2.txt
done

paste $OUTPUT/*_counts.v2.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=$25=$27=$29=$31=$33=$35=$37=$39=$41=$43=$45=$47=$49=$51=$53=$55=$57=$59=$61=$63=$65=$67=$69=$71=$73=$75=$77=$79=$81=$83=$85=$87=$89=$91=$93=$95=$97=$99=$101=$103=$105=$107=$109=$111=$113=$115=$117=$119=$121=$123=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/exp2_count_table.v2.txt
cp -R $OUTPUT/exp2_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/


#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo '** Thats all Folks! **'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/4.fraxgen.exp2.htseq.counting.sh
#-------------------------------------------------------------------------------
