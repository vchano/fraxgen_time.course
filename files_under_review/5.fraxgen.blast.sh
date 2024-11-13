#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#               Script for blast candidate sequences             #
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
#SBATC --ntasks-per-socket 24
#SBATC --exclusive
#SBATC --qos long
#SBATCH --time=48:00:00
#SBATC --begin=2023-05-05T18:00:00
#SBATCH --job-name=ULMI.BLAST
#SBATCH --output=BLAST.ULMI.UNIGENES_%J_out.txt
#SBATCH --error=BLAST.ULMI.UNIGENES_%J_err.txt
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=vmchano.gaug@gmail.com
#SBATC -a 0-79%79

# loading modules
module purge
module load blast-plus/2.11.0-py3.9 

# Setting PATHS to working directories

FASTA=$SCRATCH/FRAXGEN/0.REF.GENOME/Fraxinus_excelsior.gene_models_prot.faa
DB=$SCRATCH/NCBI.REF.SEQ
OUTPUT=$SCRATCH/FRAXGEN/5.ANNOTATION/BLASTP

#wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.protein.faa.gz
#gunzip *.protein.faa.gz -k (esto te permite guardar los archivos comprimidos, pero puedes no hacerlo sin la opción -k)
#cat *.protein.faa > refseq_plants.fasta

#makeblastdb -dbtype prot -in $DB/refseq_plants.fasta -out $DB/refseq_plants

blastp -query $FASTA -db $DB/refseq_plants -evalue 0.00001 -out $OUTPUT/fraxgen.bart.genome -outfmt 14 -max_target_seqs 50 -num_threads 48


echo "That's all folks"

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/5.fraxgen.blast.sh
#-------------------------------------------------------------------------------
