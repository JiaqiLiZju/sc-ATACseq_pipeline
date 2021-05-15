#ls  *.bam  |xargs -i samtools index {} 
ls *last.bam |while read id;do
nohup bamCoverage -p 5 --normalizeUsingRPKM -b $id -o ${id%%.*}.last.bw &
done
# nohup bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw & 

samtools sort -f -m 10G -@ 20 human.bam human.sorted.bam
samtools index human.sorted.bam
bamCoverage -p 5 --normalizeUsing RPKM -b human.sorted.bam -o human.sorted.bw &>log.bamCoverage
# bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw & 

## both -R and -S can accept multiple files 
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
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
computeMatrix scale-regions  -p 15  \
    -R ~/project/atac/mm10_Refgene/Refseq.bed  \
    -S ~/project/atac/deeptools_result/*.bw  \
    -b 10000 -a 10000  \
    --skipZeros -o matrix1_test_body.gz
# plotHeatmap -m matrix1_test_body.gz  -out ExampleHeatmap1.png
plotHeatmap -m matrix1_test_body.gz  -out test_body_Heatmap.png
plotProfile -m matrix1_test_body.gz  -out test_body_Profile.png
plotProfile -m matrix1_test_body.gz -out test_Body_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720