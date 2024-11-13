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
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=FRAX1.TEST.HS2
#SBATCH --output=HISAT2_TEST_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_TEST_FRAXGEN_%J_err.txt
#SBATCH --ntasks-per-socket 24
#SBATC --qos long
#SBATCH --time=48:00:00
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
INPUT=$SCRATCH/FRAXGEN/2.TRIMMED/TEST
OUTPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/TEST
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

list='p1194sResistantFar3Control-1 p1194sResistantFar3Control-2 p1194sResistantFar3Control-3 p1194sResistantFar3Inoculated-4 p1194sResistantFar3Inoculated-5 p1194sResistantFar3Inoculated-6 p1194sSusceptibleUW1Control-7 p1194sSusceptibleUW1Control-8 p1194sSusceptibleUW1Control-9 p1194sSusceptibleUW1Inoculated-10 p1194sSusceptibleUW1Inoculated-11 p1194sSusceptibleUW1Inoculated-12'

for sample in ${list}
do

hisat2 -p 96 --score-min L,0,-0.2 -x $REF/Fexcelsior.indexed.genome -U $INPUT/${sample}_p.fastq.gz -S $OUTPUT/${sample}_ht2.sam --summary-file $OUTPUT/${sample}_summary.txt 

# Convert SAM to BAM
samtools view -S -b $OUTPUT/${sample}_ht2.sam > $OUTPUT/${sample}_ht2.bam
# SORT -n sorts bam by name
samtools sort -l 1 -@96 -n -o $OUTPUT/${sample}_ht2.nsrt.bam -T $OUTPUT/${sample}_tmp $OUTPUT/${sample}_ht2.bam # Sort SAM by name and save in BAM
# FIXMATE adds ms and MC tags for later deduplication with markdup
samtools fixmate -@96 -O bam,level=1 -m $OUTPUT/${sample}_ht2.nsrt.bam $OUTPUT/${sample}_ht2.nsrt.fix.bam # Fill mate coordinates and save BAM
# SORT sorts by chromosome position
samtools sort -l 1 -@96 -o $OUTPUT/${sample}_ht2.fix.psrt.bam -T $OUTPUT/${sample}_tmp $OUTPUT/${sample}_ht2.nsrt.fix.bam # Sort SAM by chr position and save in BAM
# Mark duplicates in position sorted bam
samtools markdup -@96 -r -s -O bam,level=1 $OUTPUT/${sample}_ht2.fix.psrt.bam $OUTPUT/${sample}_ht2.fix.psrt.dedup.bam # Mark dups in coordinate sorted after fixmate
done
rm $OUTPUT/*_ht2.sam


#-------------------------------------------------------------------------------

# End
echo '****************************************************************'
echo 'Mapping with HISAT2 of FRAXGEN test reads to Reference'
echo 'That´s all Folks!'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.test.hisat.trimmed.sh
#-------------------------------------------------------------------------------
