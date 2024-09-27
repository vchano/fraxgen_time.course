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

#SBATCH -p gailing
#SBATC -A all
#SBATCH -n 48
#SBATCH -N 1
#SBATCH --job-name=FRAX1.CTR1.HS2
#SBATCH --output=HISAT2_CTR1_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_CTR1_FRAXGEN_%J_err.txt
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
OUTPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP1_CTR
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

list='p1258s49 p1258s50 p1258s51 p1258s52 p1258s57 p1258s58 p1258s59 p1258s60 p1258s61 p1258s68 p1258s69 p1258s71 p1258s72 p1258s79 p1258s80 p1258s81 p1258s82 p1258s87 p1258s88 p1258s89 p1258s90 p1258s96 p1258s97 p1258s98 p1258s99 p1258s106 p1258s108 p1258s109 p1258s110 p1258s118 p1258s119 p1258s120 p1258s121 p1258s127 p1258s128 p1258s129 p1258s130 p1258s135 p1258s136 p1258s137'

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
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.ctr1.hisat.trimmed.sh
#-------------------------------------------------------------------------------