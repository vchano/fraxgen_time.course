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
#SBATCH --job-name=FRAX3.CTR2.HS2
#SBATCH --output=HISAT2_CTR2_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_CTR2_FRAXGEN_%J_err.txt
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
OUTPUT=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/EXP2_CTR
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

list='p1258s145 p1258s147 p1258s148 p1258s149 p1258s156 p1258s157 p1258s160 p1258s161 p1258s169 p1258s170 p1258s171 p1258s173 p1258s181 p1258s183 p1258s184 p1258s185 p1258s193 p1258s194 p1258s195 p1258s198 p1258s205 p1258s206 p1258s208 p1258s210 p1258s221 p1258s222 p1258s223 p1258s230 p1258s231 p1258s232 p1258s235 p1258s243 p1258s244 p1258s245 p1258s247 p1258s254 p1258s255 p1258s257 p1258s259 p1258s267 p1258s270 p1258s271 p1258s272 p1258s280 p1258s281 p1258s282 p1258s284 p1258s292 p1258s293 p1258s295 p1258s296 p1258s303 p1258s304 p1258s306 p1258s308 p1258s316 p1258s317 p1258s320 p1258s321 p1258s329 p1258s330 p1258s331'
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
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.ctr2.hisat.trimmed.sh
#-------------------------------------------------------------------------------