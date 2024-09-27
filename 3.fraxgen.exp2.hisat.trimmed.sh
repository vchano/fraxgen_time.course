#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#                   Script for running HISAT2                    #
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
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=FRAX4.EXP2.HS2
#SBATCH --output=HISAT2_EXP2_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_EXP2_FRAXGEN_%J_err.txt
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
module load hisat2/2.1.0
module load cufflinks/2.2.1
module load samtools/1.9

# Set PATHS to working directories
INPUT=$SCRATCH/FRAXGEN/2.TRIMMED
OUTPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP2_EXP
REF=$SCRATCH/FRAXGEN/0.REF.GENOME

# Indexing Fraxinus excelsior genome
#gunzip -c $REF/genome_GCA_019097785.1_FRAX_001_PL_genomic.fna.gz > $REF/genome_GCA_019097785.1_FRAX_001_PL_genomic.fna
# hisat2-build -p 96 $REF/genome_GCA_019097785.1_FRAX_001_PL_genomic.fna $REF/Fexcelsior.indexed.genome

# Conversion of gff3 file in gtf to extract exons and splice sites. This is done with gffread from cufflinks

# gffread $REF/annot_Bart_edited_new.gff3 -T -o $REF/Fexcelsior.annot.gtf

# Extraction of splice sites and exons
#hisat2_extract_splice_sites.py $REF/Fexcelsior.annot.gtf > $REF/Fexcelsior_splice_sites.txt
#hisat2_extract_exons.py $REF/Fexcelsior.annot.gtf > $REF/Fexcelsior_exons.txt

# Running HISAT2 and processing output sam files 

list='p1258s150 p1258s151 p1258s154 p1258s155 p1258s163 p1258s164 p1258s166 p1258s168 p1258s177 p1258s178 p1258s179 p1258s180 p1258s188 p1258s189 p1258s191 p1258s192 p1258s199 p1258s200 p1258s202 p1258s203 p1258s211 p1258s214 p1258s217 p1258s225 p1258s226 p1258s227 p1258s229 p1258s236 p1258s237 p1258s238 p1258s249 p1258s251 p1258s252 p1258s253 p1258s260 p1258s261 p1258s262 p1258s264 p1258s274 p1258s275 p1258s276 p1258s277 p1258s285 p1258s286 p1258s287 p1258s290 p1258s297 p1258s298 p1258s299 p1258s300 p1258s309 p1258s310 p1258s312 p1258s315 p1258s323 p1258s324 p1258s325 p1258s326 p1258s334 p1258s335 p1258s336 p1258s338'

for sample in ${list}
do

#hisat2 -p 48 --score-min L,0,-0.6 -x $REF/Fexcelsior.indexed.genome -1 $INPUT/${sample}_R1_p.fastq.gz -2 $INPUT/${sample}_R2_p.fastq.gz -S $OUTPUT/${sample}_ht2.sam --summary-file $OUTPUT/${sample}_summary.txt --un-gz $OUTPUT/OTHER_FILES/ --al-gz $OUTPUT/OTHER_FILES/ --un-conc-gz $OUTPUT/OTHER_FILES/ --al-conc-gz $OUTPUT/OTHER_FILES/
#
## Convert SAM to BAM
#samtools view -S -b $OUTPUT/${sample}_ht2.sam > $OUTPUT/${sample}_ht2.bam
# SORT -n sorts bam by name
samtools sort -l 1 -@48 -n -o $OUTPUT/${sample}_ht2.nsrt.bam -T $OUTPUT/${sample}_tmp $OUTPUT/${sample}_ht2.bam # Sort SAM by name and save in BAM
# FIXMATE adds ms and MC tags for later deduplication with markdup
samtools fixmate -@48 -O bam,level=1 -m $OUTPUT/${sample}_ht2.nsrt.bam $OUTPUT/${sample}_ht2.nsrt.fix.bam # Fill mate coordinates and save BAM
# SORT sorts by chromosome position
samtools sort -l 1 -@48 -o $OUTPUT/${sample}_ht2.fix.psrt.bam -T $OUTPUT/${sample}_tmp $OUTPUT/${sample}_ht2.nsrt.fix.bam # Sort SAM by chr position and save in BAM
# Mark duplicates in position sorted bam
samtools markdup -@48 -r -s -O bam,level=1 $OUTPUT/${sample}_ht2.fix.psrt.bam $OUTPUT/${sample}_ht2.fix.psrt.dedup.bam # Mark dups in coordinate sorted after fixmate
done
rm $OUTPUT/*_ht2.sam

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
echo '** Thats all Folks! **'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.exp2.hisat.trimmed.sh
#-------------------------------------------------------------------------------