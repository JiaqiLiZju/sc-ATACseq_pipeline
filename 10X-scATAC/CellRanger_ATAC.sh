#BSUB -q normal
#BSUB -J scATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

# cellranger-atac mkfastq --id=tiny-bcl \
#                      --run=cellranger-atac-tiny-bcl-1.0.0 \
#                      --samplesheet=cellranger-atac-tiny-bcl-samplesheet-1.0.0.csv

# REF=/share/home/guoguoji/tools/CellRanger_ATAC_Reference_Earthworm/Earthworm
# cellranger-atac count --id=EW1 --sample=Sample --reference=$REF --fastqs=`pwd`

REF=/share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed/Smed
cellranger-atac count --id=Smed1 --sample=Sample --reference=$REF --fastqs=`pwd`
