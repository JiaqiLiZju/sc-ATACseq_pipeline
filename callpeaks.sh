#BSUB -q normal
#BSUB -J CHATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

bedtools bamtobed -i Aligned.out.bam >H.sorted.bed

samtools sort -o H.sorted.bam Aligned.out.bam
samtools index H.sorted.bam

bedtools bamtobed -i H.sorted.bam >H.sorted.bed

macs2 callpeak -g 32e9 \
    --nomodel --shift 100 --extsize 200 \
    --qval 5e-2 -B --SPMR \
    -f BED -t merged.sorted.bed \
    -n H.sorted --outdir ./peaks/ &>log.macs2
    
peaks=$(tail -n 1 peaks/H.sorted_summits.bed |cut -f 4)
echo '==> peaks:' $peaks >peaks.txt

# FRiP
Reads=$(bedtools intersect -a H.sorted.bed -b ./peaks/H.sorted_peaks.narrowPeak |wc -l|awk '{print $1}')
totalReads=$(wc -l H.sorted.bed|awk '{print $1}')
echo $Reads $totalReads >>peaks.txt
echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%' >>peaks.txt

# add PM
for i in {10..35}
do
snaptools snap-del --snap-file ./COL${i}/bowtie_out/H.snap --session-name PM && snaptools snap-add-pmat --snap-file ./COL${i}/bowtie_out/H.snap --peak-file ./H.sorted_peaks.bed &
done