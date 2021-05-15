# macs2 callpeak -t 2-cell-1.bed  -g mm --nomodel --shift -100 --extsize 200  -n 2-cell-1 --outdir ../peaks/
bedtools bamtobed -i human.sorted.bam >human.sorted.bed
macs2 callpeak -t human.sorted.bed -g mm --nomodel --shift  -100 --extsize 200  -n ${id%%.*} --outdir ./

# fraction of reads in called peak regions
cd ~/project/atac/peaks
ls *narrowPeak|while  read id;
do
echo $id
bed=$(basename $id "_peaks.narrowPeak").bed
#ls  -lh $bed 
Reads=$(bedtools intersect -a $bed -b $id |wc -l|awk '{print $1}')
totalReads=$(wc -l $bed|awk '{print $1}')
echo $Reads  $totalReads 
echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%'
done
