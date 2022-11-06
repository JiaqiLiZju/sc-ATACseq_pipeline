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

cat Smed.cellranger.csv
# library_id,fragments,cells
# Smed1,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed1/fragments.tsv.gz,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed1/singlecell.csv
# Smed2,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed2/fragments.tsv.gz,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed2/singlecell.csv
# Smed3,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed3/fragments.tsv.gz,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed3/singlecell.csv
# Smed4,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed4/fragments.tsv.gz,/share/home/guoguoji/RAWDATA/Jiaqili/Smed/Smed/Smed4/singlecell.csv

cellranger-atac aggr --id=Smed --csv=Smed.cellranger.csv --normalize=depth --reference=$REF

# cellranger-atac reanalyze --id=Smed_reanalysis \
#                             --params=Smed_reanalysis.csv \
#                             --peaks=Smed/outs/peaks.bed \
#                             --fragments=Smed/outs/fragments.tsv.gz \
#                             --reference=$REF