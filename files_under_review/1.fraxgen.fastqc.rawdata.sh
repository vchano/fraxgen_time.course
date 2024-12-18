#!/bin/bash
#-------------------------------------------------------------------------------

##################################################################
#                                                                #
#              Script for running FastQC in rawdata              #
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
#SBATCH -n 96
#SBATCH -N 1
#SBATCH --job-name=FRAX.FQC1_2
#SBATCH --output=FASTQC1.2_%J_out.txt
#SBATCH --error=FASTQC1.2_%J_err.txt
#SBATC --exclusive
#SBATCH --ntasks-per-socket 24
#SBATC --qos long
#SBATC --begin=2021-12-03T23:30:00
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=vmchano.gaug@gmail.com


# Load modules

module load fastqc/0.11.4
module load anaconda3/2020.11

# Set PATHS to working directories
INPUT=$SCRATCH/FRAXGEN/0.RAWDATA/RSEQ
OUTPUT=$SCRATCH/FRAXGEN/1.QC/RAWDATA/RSEQ

# Run FastQC

# this is not need as we are using regular expression
# Five samples failed during exome-capture sequencing (WES): P001_WA01, P001_WB01, P007_WA07, P007_WD07 and P007_WD08  
#list='P001_WA02 P001_WA03 P001_WA04 P001_WA05 P001_WA06 P001_WA07 P001_WA08 P001_WA09 P001_WA10 P001_WA11 P001_WA12 P001_WB02 P001_WB03 P001_WB04 P001_WB05 P001_WB06 P001_WB07 P001_WB08 P001_WB09 P001_WB10 P001_WB11 P001_WB12 P001_WC01 P001_WC02 P001_WC03 P001_WC04 P001_WC05 P001_WC06 P001_WC07 P001_WC08 P001_WC09 P001_WC10 P001_WC11 P001_WC12 P001_WD01 P001_WD02 P001_WD03 P001_WD04 P001_WD05 P001_WD06 P001_WD07 P001_WD08 P001_WD09 P001_WD10 P001_WD11 P001_WD12 P001_WE01 P001_WE02 P001_WE03 P001_WE04 P001_WE05 P001_WE06 P001_WE07 P001_WE08 P001_WE09 P001_WE10 P001_WE11 P001_WE12 P001_WF01 P001_WF02 P001_WF03 P001_WF04 P001_WF05 P001_WF06 P001_WF07 P001_WF08 P001_WF09 P001_WF10 P001_WF11 P001_WF12 P001_WG01 P001_WG02 P001_WG03 P001_WG04 P001_WG05 P001_WG06 P001_WG07 P001_WG08 P001_WG09 P001_WG10 P001_WG11 P001_WG12 P001_WH01 P001_WH02 P001_WH03 P001_WH04 P001_WH05 P001_WH06 P001_WH07 P001_WH08 P001_WH09 P001_WH10 P001_WH11 P001_WH12 P002_WA01 P002_WA02 P002_WA03 P002_WA04 P002_WA05 P002_WA06 P002_WA07 P002_WA08 P002_WA09 P002_WA10 P002_WA11 P002_WA12 P002_WB01 P002_WB02 P002_WB03 P002_WB04 P002_WB05 P002_WB06 P002_WB07 P002_WB08 P002_WB09 P002_WB10 P002_WB11 P002_WB12 P002_WC01 P002_WC02 P002_WC03 P002_WC04 P002_WC05 P002_WC06 P002_WC07 P002_WC08 P002_WC09 P002_WC10 P002_WC11 P002_WC12 P002_WD01 P002_WD02 P002_WD03 P002_WD04 P002_WD05 P002_WD06 P002_WD07 P002_WD08 P002_WD09 P002_WD10 P002_WD11 P002_WD12 P002_WE01 P002_WE02 P002_WE03 P002_WE04 P002_WE05 P002_WE06 P002_WE07 P002_WE08 P002_WE09 P002_WE10 P002_WE11 P002_WE12 P002_WF01 P002_WF02 P002_WF03 P002_WF04 P002_WF05 P002_WF06 P002_WF07 P002_WF08 P002_WF09 P002_WF10 P002_WF11 P002_WF12 P002_WG01 P002_WG02 P002_WG03 P002_WG04 P002_WG05 P002_WG06 P002_WG07 P002_WG08 P002_WG09 P002_WG10 P002_WG11 P002_WG12 P002_WH01 P002_WH02 P002_WH03 P002_WH04 P002_WH05 P002_WH06 P002_WH07 P002_WH08 P002_WH09 P002_WH10 P002_WH11 P002_WH12 P003_WA01 P003_WA02 P003_WA03 P003_WA04 P003_WA05 P003_WA06 P003_WA07 P003_WA08 P003_WA09 P003_WA10 P003_WA11 P003_WA12 P003_WB01 P003_WB02 P003_WB03 P003_WB04 P003_WB05 P003_WB06 P003_WB07 P003_WB08 P003_WB09 P003_WB10 P003_WB11 P003_WB12 P003_WC01 P003_WC02 P004_WA01 P004_WA02 P004_WA03 P004_WA04 P004_WA05 P004_WA06 P004_WA07 P004_WA08 P004_WA09 P004_WA10 P004_WA11 P004_WA12 P004_WB01 P004_WB02 P004_WB03 P004_WB04 P004_WB05 P004_WB06 P004_WB07 P004_WB08 P004_WB09 P004_WB10 P004_WB11 P004_WB12 P004_WC01 P004_WC02 P004_WC03 P004_WC04 P004_WC05 P004_WC06 P004_WC07 P004_WC08 P004_WC09 P004_WC10 P004_WC11 P004_WC12 P004_WD01 P004_WD02 P004_WD03 P004_WD04 P004_WD05 P004_WD06 P004_WD07 P004_WD08 P004_WD09 P004_WD10 P004_WD11 P004_WD12 P004_WE01 P004_WE02 P004_WE03 P004_WE04 P004_WE05 P004_WE06 P004_WE07 P004_WE08 P004_WE09 P004_WE10 P004_WE11 P004_WE12 P004_WF01 P004_WF02 P004_WF03 P004_WF04 P004_WF05 P004_WF06 P004_WF07 P004_WF08 P004_WF09 P004_WF10 P004_WF11 P004_WF12 P004_WG01 P004_WG02 P004_WG03 P004_WG04 P004_WG05 P004_WG06 P004_WG07 P004_WG08 P004_WG09 P004_WG10 P004_WG11 P004_WG12 P004_WH01 P004_WH02 P004_WH03 P004_WH04 P004_WH05 P004_WH06 P004_WH07 P004_WH08 P004_WH09 P004_WH10 P004_WH11 P004_WH12 P005_WA01 P005_WA02 P005_WA03 P005_WA04 P005_WA05 P005_WA06 P005_WA07 P005_WA08 P005_WA09 P005_WA10 P005_WA11 P005_WA12 P005_WB01 P005_WB02 P005_WB03 P005_WB04 P005_WB05 P005_WB06 P005_WB07 P005_WB08 P005_WB09 P005_WB10 P005_WB11 P005_WB12 P005_WC01 P005_WC02 P005_WC03 P005_WC04 P005_WC05 P005_WC06 P005_WC07 P005_WC08 P005_WC09 P005_WC10 P005_WC11 P005_WC12 P005_WD01 P005_WD02 P005_WD03 P005_WD04 P005_WD05 P005_WD06 P005_WD07 P005_WD08 P005_WD09 P005_WD10 P005_WD11 P005_WD12 P005_WE01 P005_WE02 P005_WE03 P005_WE04 P005_WE05 P005_WE06 P005_WE07 P005_WE08 P005_WE09 P005_WE10 P005_WE11 P005_WE12 P005_WF01 P005_WF02 P005_WF03 P005_WF04 P005_WF05 P005_WF06 P005_WF07 P005_WF08 P005_WF09 P005_WF10 P005_WF11 P005_WF12 P005_WG01 P005_WG02 P005_WG03 P005_WG04 P005_WG05 P005_WG06 P005_WG07 P005_WG08 P005_WG09 P005_WG10 P005_WG11 P005_WG12 P005_WH01 P005_WH02 P005_WH03 P005_WH04 P005_WH05 P005_WH06 P005_WH07 P005_WH08 P005_WH09 P005_WH10 P005_WH11 P005_WH12 P006_WA01 P006_WA02 P006_WA03 P006_WA04 P006_WA05 P006_WA06 P006_WA07 P006_WA08 P006_WA09 P006_WA10 P006_WA11 P006_WA12 P006_WB01 P006_WB02 P006_WB03 P006_WB04 P006_WB05 P006_WB06 P006_WB07 P006_WB08 P006_WB09 P006_WB10 P006_WB11 P006_WB12 P006_WC01 P006_WC02 P006_WC03 P006_WC04 P006_WC05 P006_WC06 P006_WC07 P006_WC08 P006_WC09 P006_WC10 P006_WC11 P006_WC12 P006_WD01 P006_WD02 P006_WD03 P006_WD04 P006_WD05 P006_WD06 P006_WD07 P006_WD08 P006_WD09 P006_WD10 P006_WD11 P006_WD12 P006_WE01 P006_WE02 P006_WE03 P006_WE04 P006_WE05 P006_WE06 P006_WE07 P006_WE08 P006_WE09 P006_WE10 P006_WE11 P006_WE12 P006_WF01 P006_WF02 P006_WF03 P006_WF04 P006_WF05 P006_WF06 P006_WF07 P006_WF08 P006_WF09 P006_WF10 P006_WF11 P006_WF12 P006_WG01 P006_WG02 P006_WG03 P006_WG04 P006_WG05 P006_WG06 P006_WG07 P006_WG08 P006_WG09 P006_WG10 P006_WG11 P006_WG12 P006_WH01 P006_WH02 P006_WH03 P006_WH04 P006_WH05 P006_WH06 P006_WH07 P006_WH08 P006_WH09 P006_WH10 P006_WH11 P006_WH12 P007_WA01 P007_WA02 P007_WA03 P007_WA04 P007_WA05 P007_WA06 P007_WA08 P007_WA09 P007_WA10 P007_WA11 P007_WA12 P007_WB01 P007_WB02 P007_WB03 P007_WB04 P007_WB05 P007_WB06 P007_WB07 P007_WB08 P007_WB09 P007_WB10 P007_WB11 P007_WB12 P007_WC01 P007_WC02 P007_WC03 P007_WC04 P007_WC05 P007_WC06 P007_WC07 P007_WC08 P007_WC09 P007_WC10 P007_WC11 P007_WC12 P007_WD01 P007_WD02 P007_WD03 P007_WD04 P007_WD05 P007_WD06 P007_WD09 P007_WD10 P007_WD11 P007_WD12 P007_WE01 P007_WE02 P007_WE03 P007_WE04 P007_WE05 P007_WE06 P007_WE07 P007_WE08 P007_WE09 P007_WE10 P007_WE11 P007_WE12 P007_WF01 P007_WF02 P007_WF03 P007_WF04 P007_WF05 P007_WF06 P007_WF07 P007_WF08 P007_WF09 P007_WF10 P007_WF11 P007_WF12 P007_WG01 P007_WG02 P007_WG03 P007_WG04 P007_WG05 P007_WG06 P007_WG07 P007_WG08 P007_WG09 P007_WG10 P007_WG11 P007_WG12 P007_WH01 P007_WH02 P007_WH03 P007_WH04 P007_WH05 P007_WH06 P007_WH07 P007_WH08 P007_WH09 P007_WH10 P007_WH11 P007_WH12 P008_WA01 P008_WA02 P008_WA03 P008_WA04 P008_WA05 P008_WA06 P008_WA07 P008_WA08 P008_WA09 P008_WA10 P008_WA11 P008_WA12 P008_WB01 P008_WB02 P008_WB03 P008_WB04 P008_WB05 P008_WB06'
#for sample in ${list}
#do
#
#done

fastqc $INPUT/*.fastq.gz --outdir $OUTPUT --threads 96

# Run MultiQC
source activate multiqc
multiqc $OUTPUT/
conda deactivate

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
# Excute this script as sbatch $SCRATCH/FRAXGEN/0.APPS/1.fraxgen.fastqc.rawdata.sh
#-------------------------------------------------------------------------------
