#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#                   Script for running BOWTIE2                   #
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
#SBATCH --job-name=FRAX.HS2
#SBATCH --output=HISAT2_FRAXGEN_%J_out.txt
#SBATCH --error=HISAT2_FRAXGEN_%J_err.txt
#SBATCH --ntasks-per-socket 24
#SBATCH --qos long
#SBATCH --time=120:00:00
#SBATC --begin=2022-11-25T01:00:00
#SBATCh --exclusive
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

#list='p1258s49 p1258s50 p1258s51 p1258s52 p1258s53 p1258s54 p1258s55 p1258s56 p1258s57 p1258s58 p1258s59 p1258s60 p1258s61 p1258s62 p1258s63 p1258s64 p1258s67 p1258s68 p1258s69 p1258s71 p1258s72 p1258s73 p1258s74 p1258s75 p1258s77 p1258s78 p1258s79 p1258s80 p1258s81 p1258s82 p1258s83 p1258s84 p1258s85 p1258s86 p1258s87 p1258s88 p1258s89 p1258s90 p1258s91 p1258s92 p1258s93 p1258s94 p1258s96 p1258s97 p1258s98 p1258s99 p1258s100 p1258s103 p1258s104 p1258s105 p1258s106 p1258s108 p1258s109 p1258s110 p1258s111 p1258s114 p1258s116 p1258s118 p1258s119 p1258s120 p1258s121 p1258s122 p1258s124 p1258s125 p1258s126 p1258s127 p1258s128 p1258s129 p1258s130 p1258s131 p1258s132 p1258s133 p1258s134 p1258s135 p1258s136 p1258s137 p1258s139 p1258s141 p1258s142 p1258s143 p1258s145 p1258s147 p1258s148 p1258s149 p1258s150 p1258s151 p1258s154 p1258s155 p1258s156 p1258s157 p1258s160 p1258s161 p1258s163 p1258s164 p1258s166 p1258s168 p1258s169 p1258s170 p1258s171 p1258s173 p1258s177 p1258s178 p1258s179 p1258s180 p1258s181 p1258s183 p1258s184 p1258s185 p1258s188 p1258s189 p1258s191 p1258s192 p1258s193 p1258s194 p1258s195 p1258s198 p1258s199 p1258s200 p1258s202 p1258s203 p1258s205 p1258s206 p1258s208 p1258s210 p1258s211 p1258s214 p1258s217 p1258s221 p1258s222 p1258s223 p1258s225 p1258s226 p1258s227 p1258s229 p1258s230 p1258s231 p1258s232 p1258s235 p1258s236 p1258s237 p1258s238 p1258s243 p1258s244 p1258s245 p1258s247 p1258s249 p1258s251 p1258s252 p1258s253 p1258s254 p1258s255 p1258s257 p1258s259 p1258s260 p1258s261 p1258s262 p1258s264 p1258s267 p1258s270 p1258s271 p1258s272 p1258s274 p1258s275 p1258s276 p1258s277 p1258s280 p1258s281 p1258s282 p1258s284 p1258s285 p1258s286 p1258s287 p1258s290 p1258s292 p1258s293 p1258s295 p1258s296 p1258s297 p1258s298 p1258s299 p1258s300 p1258s303 p1258s304 p1258s306 p1258s308 p1258s309 p1258s310 p1258s312 p1258s315 p1258s316 p1258s317 p1258s320 p1258s321 p1258s323 p1258s324 p1258s325 p1258s326 p1258s329 p1258s330 p1258s331 p1258s334 p1258s335 p1258s336 p1258s338'
list='p1258s49 p1258s50 p1258s51 p1258s52 p1258s57 p1258s58 p1258s59 p1258s60 p1258s61 p1258s68 p1258s69 p1258s71 p1258s72 p1258s79 p1258s80 p1258s81 p1258s82 p1258s87 p1258s88 p1258s89 p1258s90 p1258s96 p1258s97 p1258s98 p1258s99 p1258s106 p1258s108 p1258s109 p1258s110 p1258s118 p1258s119 p1258s120 p1258s121 p1258s127 p1258s128 p1258s129 p1258s130 p1258s135 p1258s136 p1258s137 p1258s138'

for sample in ${list}
do

hisat2 -p 48 --score-min L,0,-0.6 -x $REF/Fexcelsior.indexed.genome -1 $INPUT/${sample}_R1_p.fastq.gz -2 $INPUT/${sample}_R2_p.fastq.gz -S $OUTPUT/${sample}_ht2.sam --summary-file $OUTPUT/${sample}_summary.txt --un-gz $OUTPUT/OTHER_FILES/ --al-gz $OUTPUT/OTHER_FILES/ --un-conc-gz $OUTPUT/OTHER_FILES/ --al-conc-gz $OUTPUT/OTHER_FILES/

# Convert SAM to BAM
samtools view -S -b $OUTPUT/${sample}_ht2.sam > $OUTPUT/${sample}_ht2.bam
# SORT -n sorts bam by name
samtools sort -l 1 -@48 -n -o $OUTPUT/${sample}_ht2.nsrt.bam -T $OUTPUT/tmp/${sample}_tmp $OUTPUT/${sample}_ht2.bam # Sort SAM by name and save in BAM
# FIXMATE adds ms and MC tags for later deduplication with markdup
samtools fixmate -@48 -O bam,level=1 -m $OUTPUT/${sample}_ht2.nsrt.bam $OUTPUT/${sample}_ht2.nsrt.fix.bam # Fill mate coordinates and save BAM
# SORT sorts by chromosome position
samtools sort -l 1 -@48 -o $OUTPUT/${sample}_ht2.fix.psrt.bam -T $OUTPUT/tmp/${sample}_tmp $OUTPUT/${sample}_ht2.nsrt.fix.bam # Sort SAM by chr position and save in BAM
# Mark duplicates in position sorted bam
samtools markdup -@48 -r -s -O bam,level=1 $OUTPUT/${sample}_ht2.fix.psrt.bam $OUTPUT/${sample}_ht2.fix.psrt.dedup.bam # Mark dups in coordinate sorted after fixmate
done
rm $OUTPUT/*_ht2.sam


#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/3.fraxgen.hisat.trimmed.sh
#-------------------------------------------------------------------------------