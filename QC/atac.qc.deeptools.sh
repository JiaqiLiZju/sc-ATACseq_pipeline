# sort bam
samtools sort -f -m 10G -@ 20 human.bam human.sorted.bam
samtools index human.sorted.bam

# Insert QC
samtools view $outdir/human.sorted.bam |awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' | sort | uniq | cut -f 2 >$outdir/fragment.length.txt

# peak-calling
bedtools bamtobed -i human.sorted.bam >human.sorted.bed
macs2 callpeak -g hs --nomodel \
    --shift -100 --extsize 200 \
    -t human.sorted.bed \
    -n human.sorted --outdir ./peaks/ &>log.macs2
# -rw-rw-r-- 1 ggj ggj 1.3M 5月  16 10:02 human.sorted_peaks.narrowPeak
# -rw-rw-r-- 1 ggj ggj 1.4M 5月  16 10:02 human.sorted_peaks.xls
# -rw-rw-r-- 1 ggj ggj 910K 5月  16 10:02 human.sorted_summits.bed

# peak number
peaks=$(cat peaks/human.sorted_summits.bed |wc -l)
echo '==> peaks:' $peaks >peaks.txt

# FRiP
Reads=$(bedtools intersect -a human.sorted.bed -b ./peaks/human.sorted_peaks.narrowPeak |wc -l|awk '{print $1}')
totalReads=$(wc -l human.sorted.bed|awk '{print $1}')
echo $Reads $totalReads >>peaks.txt
echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%' >>peaks.txt

# TSS enrichment
bamCoverage -p 5 --normalizeUsing RPKM -b human.sorted.bam -o human.sorted.bw &>log.bamCoverage
# bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw & 

## both -R and -S can accept multiple files 
computeMatrix reference-point -p 15 \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R ./Mix.Human.Refseq.bed \
    -S ./human.sorted.bw \
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
    -R ./Mix.Human.Refseq.bed \
    -S ./human.sorted.bw \
    -b 1000 -a 1000  \
    --skipZeros -o matrix1_test_body.gz
# plotHeatmap -m matrix1_test_body.gz  -out ExampleHeatmap1.png
plotHeatmap -m matrix1_test_body.gz  -out test_body_Heatmap.png
plotProfile -m matrix1_test_body.gz  -out test_body_Profile.png
plotProfile -m matrix1_test_body.gz -out test_Body_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720
