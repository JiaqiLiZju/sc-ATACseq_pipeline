#BSUB -q normal
#BSUB -J scATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

REF=/share/home/guoguoji/tools/CellRanger_ATAC_Reference_Earthworm/Earthworm

# cellranger-atac mkfastq --id=tiny-bcl \
#                      --run=cellranger-atac-tiny-bcl-1.0.0 \
#                      --samplesheet=cellranger-atac-tiny-bcl-samplesheet-1.0.0.csv

# cellranger-atac mkfastq \
#     --id=EW1 \
#     --run=../BCL/211009_E00559_0224_AHGGCHCCX2/Data/Intensities/BaseCalls/L005/ \
#     --csv=bcl.csv

cellranger-atac count --id=Sample --sample Sample --reference=$REF --fastqs ./FASTQ
