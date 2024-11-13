#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#   Script for running Trimmomatic and FastQC in trimmed data    #
#                                                                #
##################################################################

#-------------------------------------------------------------------------------
#
#    Forest Genetics and Forest Tree Breeding
#    Büsgenweg 
#    Georg-August-Universtät Göttingen
#    https://github.com/vchano/
#
# Licence: GNU General Public Licence Version 3
#-------------------------------------------------------------------------------

#SBATCH -p gailing
#SBATC -A all
#SBATCH -n 46
#SBATCH -N 1
#SBATCH --job-name=FRAXGEN.2.3.4
#SBATCH --output=2.3.4_%J_out.txt
#SBATCH --error=2.3.4_%J_err.txt
#SBATC --exclusive
#SBATC -q short
#SBATCH --ntasks-per-socket 24
#SBATCH --qos long
#SBATC --begin=2021-12-03T23:30:00
#SBATCH --time=96:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=vmchano.gaug@gmail.com



# Load modules

module purge 
module load fastqc/0.11.4
module load anaconda3/2021.05
module load trimmomatic/0.36
module load hisat2/2.1.0
module load cufflinks/2.2.1
module load samtools/1.9
module load anaconda3/2021.05


# Set PATHS to working directories
#INPUT=$SCRATCH/FRAXGEN/0.RAWDATA/RSEQ/RSEQ2
#TRIMMED=$SCRATCH/FRAXGEN/2.TRIMMED/RSEQ/RSEQ2
#FQC_RAW=$SCRATCH/FRAXGEN/1.QC/RAWDATA/RSEQ/RSEQ2
#FQC_TRIM=$SCRATCH/FRAXGEN/1.QC/TRIMMED/RSEQ/RSEQ2
#OUTPUT=$SCRATCH/FRAXGEN/1.QC/TRIMMED/RSEQ/RSEQ2

# Run FastQC on rawdata
#fastqc $INPUT/*.fastq.gz --outdir $FQC_RAW --threads 48
## Run MultiQC
#source activate multiqc
#multiqc $FQC_RAW/ -o $FQC_RAW
#conda deactivate

# Run Trimmomatic on rawdata

#list='p1258s49 p1258s50 p1258s51 p1258s52 p1258s53 p1258s54 p1258s55 p1258s56 p1258s57 p1258s58 p1258s59 p1258s60 p1258s61 p1258s62 p1258s63 p1258s64 p1258s67 p1258s68 p1258s69 p1258s71 p1258s72 p1258s73 p1258s74 p1258s75 p1258s77 p1258s78 p1258s79 p1258s80 p1258s81 p1258s82 p1258s83 p1258s84 p1258s85 p1258s86 p1258s87 p1258s88 p1258s89 p1258s90 p1258s91 p1258s92 p1258s93 p1258s94 p1258s96 p1258s97 p1258s98 p1258s99 p1258s100 p1258s103 p1258s104 p1258s105 p1258s106 p1258s108 p1258s109 p1258s110 p1258s111 p1258s114 p1258s116 p1258s118 p1258s119 p1258s120 p1258s121 p1258s122 p1258s124 p1258s125 p1258s126 p1258s127 p1258s128 p1258s129 p1258s130 p1258s131 p1258s132 p1258s133 p1258s134 p1258s135 p1258s136 p1258s137 p1258s139 p1258s141 p1258s142 p1258s143 p1258s145 p1258s147 p1258s148 p1258s149 p1258s150 p1258s151 p1258s154 p1258s155 p1258s156 p1258s157 p1258s160 p1258s161 p1258s163 p1258s164 p1258s166 p1258s168 p1258s169 p1258s170 p1258s171 p1258s173 p1258s177 p1258s178 p1258s179 p1258s180 p1258s181 p1258s183 p1258s184 p1258s185 p1258s188 p1258s189 p1258s191 p1258s192 p1258s193 p1258s194 p1258s195 p1258s198 p1258s199 p1258s200 p1258s202 p1258s203 p1258s205 p1258s206 p1258s208 p1258s210 p1258s211 p1258s214 p1258s217 p1258s221 p1258s222 p1258s223 p1258s225 p1258s226 p1258s227 p1258s229 p1258s230 p1258s231 p1258s232 p1258s235 p1258s236 p1258s237 p1258s238 p1258s243 p1258s244 p1258s245 p1258s247 p1258s249 p1258s251 p1258s252 p1258s253 p1258s254 p1258s255 p1258s257 p1258s259 p1258s260 p1258s261 p1258s262 p1258s264 p1258s267 p1258s270 p1258s271 p1258s272 p1258s274 p1258s275 p1258s276 p1258s277 p1258s280 p1258s281 p1258s282 p1258s284 p1258s285 p1258s286 p1258s287 p1258s290 p1258s292 p1258s293 p1258s295 p1258s296 p1258s297 p1258s298 p1258s299 p1258s300 p1258s303 p1258s304 p1258s306 p1258s308 p1258s309 p1258s310 p1258s312 p1258s315 p1258s316 p1258s317 p1258s320 p1258s321 p1258s323 p1258s324 p1258s325 p1258s326 p1258s329 p1258s330 p1258s331 p1258s334 p1258s335 p1258s336 p1258s338'
#list='p1258s243 p1258s249 p1258s259 p1258s275 p1258s280 p1258s281 p1258s282 p1258s284 p1258s304 p1258s321 p1258s325 p1258s334'
#list='p1258s112 p1258s138 p1258s218 p1258s332'
#list='p1361s88 p1361s93 p1361s94 p1361s98 p1361s104 p1361s105 p1361s109'
#for sample in ${list}
#do
#
#java -jar /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/TRIMMOMATIC/0.36/trimmomatic-0.36.jar PE -threads 48 -phred33 $INPUT/${sample}*_R1_001.fastq.gz $INPUT/${sample}*_R2_001.fastq.gz $TRIMMED/${sample}_R1_p.fastq.gz $TRIMMED/UNPAIRED/${sample}_R1_u.fastq.gz $TRIMMED/${sample}_R2_p.fastq.gz $TRIMMED/UNPAIRED/${sample}_R2_u.fastq.gz ILLUMINACLIP:/scratch/users/chano/TGC_WES/0.APPS/adapters.fa:2:30:10 CROP:49 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:20
#
#done

# Run FastQC on trimmed data
#fastqc $TRIMMED/*_p.fastq.gz --outdir $FQC_TRIM --threads 48
## Run MultiQC
#source activate multiqc
#multiqc $FQC_TRIM/ -o $FQC_TRIM 
#conda deactivate

# Set PATHS to working directories
INPUT2=$SCRATCH/FRAXGEN/2.TRIMMED/RSEQ/RSEQ2
OUTPUT2=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/RSEQ/RSEQ2
REF=$SCRATCH/FRAXGEN/0.REF.GENOME
list='p1361s88 p1361s93 p1361s94 p1361s98 p1361s104 p1361s105 p1361s109'
for sample in ${list}
do

hisat2 -p 46 --score-min L,0,-0.2 -x $REF/Fexcelsior.indexed.genome -1 $INPUT2/${sample}_R1_p.fastq.gz -2 $INPUT2/${sample}_R2_p.fastq.gz -S $OUTPUT2/${sample}_ht2.sam --summary-file $OUTPUT2/${sample}_summary.txt --un-gz $OUTPUT2/OTHER_FILES/ --al-gz $OUTPUT2/OTHER_FILES/ --un-conc-gz $OUTPUT2/OTHER_FILES/ --al-conc-gz $OUTPUT2/OTHER_FILES/

# Convert SAM to BAM
samtools view -S -b $OUTPUT2/${sample}_ht2.sam > $OUTPUT2/${sample}_ht2.bam
# SORT -n sorts bam by name
samtools sort -l 1 -@46 -n -o $OUTPUT2/${sample}_ht2.nsrt.bam -T $OUTPUT2/${sample}_tmp $OUTPUT2/${sample}_ht2.bam # Sort SAM by name and save in BAM
# FIXMATE adds ms and MC tags for later deduplication with markdup
samtools fixmate -@46 -O bam,level=1 -m $OUTPUT2/${sample}_ht2.nsrt.bam $OUTPUT2/${sample}_ht2.nsrt.fix.bam # Fill mate coordinates and save BAM
# SORT sorts by chromosome position
samtools sort -l 1 -@46 -o $OUTPUT2/${sample}_ht2.fix.psrt.bam -T $OUTPUT2/${sample}_tmp $OUTPUT2/${sample}_ht2.nsrt.fix.bam # Sort SAM by chr position and save in BAM
# Mark duplicates in position sorted bam
samtools markdup -@46 -r -s -O bam,level=1 $OUTPUT2/${sample}_ht2.fix.psrt.bam $OUTPUT2/${sample}_ht2.fix.psrt.dedup.bam # Mark dups in coordinate sorted after fixmate
done
#rm $OUTPUT2/*_ht2.sam


# Load HTSEQ
export PATH="/opt/sw/rev/21.12/haswell/gcc-9.3.0/anaconda3-2021.05-e3srav/bin:$PATH"
export PATH="/opt/sw/rev/21.12/haswell/gcc-9.3.0/anaconda3-2021.05-e3srav/bin/conda:$PATH"
source activate htseq

# Set PATHS to working directories
INPUT3=$SCRATCH/FRAXGEN/3.MAPPING.HISAT2/RSEQ/RSEQ2
OUTPUT3=$SCRATCH/FRAXGEN/4.COUNTING/RSEQ/RSEQ2
REF2=$SCRATCH/FRAXGEN/0.REF.GENOME/annot_Bart_edited_new.gff3

for sample in ${list}
do
htseq-count -f bam --idattr=ID --type=mRNA --stranded=no $INPUT3/${sample}_bt2.fix.psrt.bam $REF2 > $OUTPUT3/${sample}_counts.rseq2.txt
done

paste $OUTPUT3/*_counts.rseq2.txt | awk '{$3=$5=$7=$9=""; print $0}' | sed -e "s/\  \b/\ /g" > $OUTPUT3/rseq2_count_table.v2.txt
cp -R $OUTPUT3/rseq2_count_table.v2.txt $HOME/FRAXGEN/4.FRAXEN.GEA/


#-------------------------------------------------------------------------------

# End
echo '****************************************************************'
echo 'Trimming, qual checking, mapping with HISAT2, and counting with HTSEQ of FRAXGEN RSEQ2'
echo 'That´s all Folks!'
exit 0

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/2_3_4.fraxgen.trimm.fastqc.hisat2.sh
#-------------------------------------------------------------------------------
