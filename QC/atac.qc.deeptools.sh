# sort bam
samtools sort -o Aligned.sorted.bam Aligned.out.bam
samtools index Aligned.sorted.bam

# Insert QC
samtools view $outdir/Aligned.sorted.bam |awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' | sort | uniq | cut -f 2 >$outdir/fragment.length.txt

# peak-calling
bedtools bamtobed -i Aligned.sorted.bam >Aligned.sorted.bed
macs2 callpeak -g hs \ # effective genome size
    --nomodel --shift 100 --extsize 200 \
    --qval 5e-2 -B --SPMR \
    -f BED -t Aligned.sorted.bed \
    -n Aligned.sorted --outdir ./peaks/ &>log.macs2

# macs2 callpeak -g mm \
#     --nomodel --shift 37 --ext 73 \
#     --qval 1e-2 -B --SPMR --call-summits \
#     -f BED -t Aligned.sorted.bed \
#     -n Aligned.sorted --outdir ./peaks/ &>log.macs2
# -rw-rw-r-- 1 ggj ggj 1.3M 5月  16 10:02 Aligned.sorted_peaks.narrowPeak
# -rw-rw-r-- 1 ggj ggj 1.4M 5月  16 10:02 Aligned.sorted_peaks.xls
# -rw-rw-r-- 1 ggj ggj 910K 5月  16 10:02 Aligned.sorted_summits.bed

# peak number
# peaks=$(cat peaks/Aligned.sorted_summits.bed |wc -l)
peaks=$(tail -n 1 peaks/Aligned.sorted_summits.bed |cut -f 4)
echo '==> peaks:' $peaks >peaks.txt

# FRiP
Reads=$(bedtools intersect -a Aligned.sorted.bed -b ./peaks/Aligned.sorted_peaks.narrowPeak |wc -l|awk '{print $1}')
totalReads=$(wc -l Aligned.sorted.bed|awk '{print $1}')
echo $Reads $totalReads >>peaks.txt
echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%' >>peaks.txt

# TSS enrichment
bamCoverage -p 5 --normalizeUsing RPKM -b Aligned.sorted.bam -o Aligned.sorted.bw &>log.bamCoverage
# bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw & 

## both -R and -S can accept multiple files 
computeMatrix reference-point -p 15 \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ./Mix.Aligned.Refseq.bed \
    -S ./Aligned.sorted.bw \
    --skipZeros \
    -o matrix1_test_TSS.gz \
    --outFileSortedRegions regions1_test_genes.bed
## both plotHeatmap and plotProfile will use the output from computeMatrix
plotHeatmap -m matrix1_test_TSS.gz  -out test_Heatmap.png
plotHeatmap -m matrix1_test_TSS.gz  -out test_Heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m matrix1_test_TSS.gz  -out test_Profile.png
plotProfile -m matrix1_test_TSS.gz  -out test_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

## gene body
computeMatrix scale-regions -p 15 \
    -R ./Mix.Aligned.Refseq.bed \
    -S ./Aligned.sorted.bw \
    -b 1000 -a 1000  \
    --skipZeros -o matrix1_test_body.gz
# plotHeatmap -m matrix1_test_body.gz  -out ExampleHeatmap1.png
plotHeatmap -m matrix1_test_body.gz  -out test_body_Heatmap.png
plotProfile -m matrix1_test_body.gz  -out test_body_Profile.png
plotProfile -m matrix1_test_body.gz -out test_Body_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720

# show the prediction using USUC & IGV genome browser
# pyGenomeTracks
# pip install pyGenomeTracks
make_tracks_file --trackFiles \
    ./Aligned.sorted.bw \
    ./Mix.Human.Refseq.bed \
    -o tracks.ini

# track-plot
pyGenomeTracks --tracks tracks.ini --region HUMAN_8:0-25000000 --outFileName tracks_zoom.png --width 50 --fontSize 5 --dpi 300
pyGenomeTracks --tracks tracks.ini --region chr8:0-25000000 --outFileName tracks_zoom.png --width 50 --fontSize 5 --dpi 300

# TSS enrichment
python fraglist.py Aligned.snap barcode_fragment_Aligned.bed
cat barcode_fragment_Aligned.bed |cut -f 4 |uniq >barcode.idx
head Mix.Refseq.bed
ATACCellTSS -bed barcode_fragment_Aligned.bed -xgi barcode.idx -tss ./Mix.Refseq.bed -out barcode_tssEnrich.txt