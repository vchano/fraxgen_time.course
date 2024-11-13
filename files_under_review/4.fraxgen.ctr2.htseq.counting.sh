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
#SBATCH --job-name=HTS3.FRAX.CTR2
#SBATCH --output=HTSEQ_CTR2_FRAXGEN_%J_out.txt
#SBATCH --error=HTSEQ_CTR2_FRAXGEN_%J_err.txt
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
INPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP2_CTR
OUTPUT=$SCRATCH/FRAXGEN/4.COUNTING/EXP2_CTR
REF=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3


# Running HTSEQ and processing count files 

list='p1258s145 p1258s147 p1258s148 p1258s149 p1258s156 p1258s157 p1258s160 p1258s161 p1258s169 p1258s170 p1258s171 p1258s173 p1258s181 p1258s183 p1258s184 p1258s185 p1258s193 p1258s194 p1258s195 p1258s198 p1258s205 p1258s206 p1258s208 p1258s210 p1258s221 p1258s222 p1258s223 p1258s230 p1258s231 p1258s232 p1258s235 p1258s243 p1258s244 p1258s245 p1258s247 p1258s254 p1258s255 p1258s257 p1258s259 p1258s267 p1258s270 p1258s271 p1258s272 p1258s280 p1258s281 p1258s282 p1258s284 p1258s292 p1258s293 p1258s295 p1258s296 p1258s303 p1258s304 p1258s306 p1258s308 p1258s316 p1258s317 p1258s320 p1258s321 p1258s329 p1258s330 p1258s331'

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT/${sample}_ht2.fix.psrt.bam $REF > $OUTPUT/${sample}_counts.v2.txt
done

paste $OUTPUT/*_counts.v2.txt | awk '{$3=$5=$7=$9=$11=$13=$15=$17=$19=$21=$23=$25=$27=$29=$31=$33=$35=$37=$39=$41=$43=$45=$47=$49=$51=$53=$55=$57=$59=$61=$63=$65=$67=$69=$71=$73=$75=$77=$79=$81=$83=$85=$87=$89=$91=$93=$95=$97=$99=$101=$103=$105=$107=$109=$111=$113=$115=$117=$119=$121=$123=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT/ctr2_count_table.v2.txt
cp -R $OUTPUT/ctr2_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/


#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo '** Thats all Folks! **'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/4.fraxgen.ctr2.htseq.counting.sh
#-------------------------------------------------------------------------------
