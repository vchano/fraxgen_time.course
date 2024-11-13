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
#SBATCH --job-name=FRAX2.EXP1.HS2
#SBATCH --output=HISAT2_EXP1_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_EXP1_FRAXGEN_%J_err.txt
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
OUTPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP1_EXP
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

list='p1258s53 p1258s54 p1258s55 p1258s56 p1258s62 p1258s63 p1258s64 p1258s67 p1258s73 p1258s74 p1258s75 p1258s77 p1258s78 p1258s83 p1258s84 p1258s85 p1258s86 p1258s91 p1258s92 p1258s93 p1258s94 p1258s100 p1258s103 p1258s104 p1258s105 p1258s111 p1258s114 p1258s116 p1258s122 p1258s124 p1258s125 p1258s126 p1258s131 p1258s132 p1258s133 p1258s134 p1258s139 p1258s141 p1258s142 p1258s143'

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
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.exp1.hisat.trimmed.sh
#-------------------------------------------------------------------------------
