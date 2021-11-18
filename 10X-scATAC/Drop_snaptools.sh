#BSUB -q normal
#BSUB -J Smed-ATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

reference=/share/home/guoguoji/tools/BWA-Reference-EarthWorm/Earthworm.fasta

dropseq_root=/share/home/guoguoji/tools/Drop-seq_tools-1.12/
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
bwa_exec=/share/home/guoguoji/tools/bwa-0.7.15/bwa

sample_name=$(basename `pwd`)
outdir=bowtie_out
tmpdir=tmp

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

if [ ! -d $tmpdir ]; then
    mkdir $tmpdir
fi

# gunzip -c ${sample_name}_H_R1.fastq.gz >$tmpdir/H_R1.fastq
# gunzip -c ${sample_name}_H_R2.fastq.gz >$tmpdir/H_R2.fastq
gunzip -c Undetermined_S0_L008_R1_001.fastq.gz >$tmpdir/H_R1.fastq
gunzip -c Undetermined_S0_L008_R2_001.fastq.gz >$tmpdir/H_R2.fastq

# pre-alignment tag and trim
java -jar ${picard_jar} FastqToSam F1=$tmpdir/H_R1.fastq F2=$tmpdir/H_R2.fastq O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name 

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R1.txt \
    BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell_R1.bam COMPRESSION_LEVEL=0

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R2.txt \
    BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0

# barcode qname
samtools view $tmpdir/unaligned_tagged_Cell.bam -H > $tmpdir/unaligned_tagged_Cell.snap.sam
samtools view $tmpdir/unaligned_tagged_Cell.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); printf "%s:%s\n", td["CB"], $0 } } }' >> $tmpdir/unaligned_tagged_Cell.snap.sam

samtools view $tmpdir/unaligned_tagged_Cell.snap.bam | cut -f 1 | head
java -Xmx100g -jar ${picard_jar} SamToFastq INPUT=$tmpdir/unaligned_tagged_Cell.snap.sam FASTQ=$tmpdir/unaligned_R1.fastq SECOND_END_FASTQ=$tmpdir/unaligned_R2.fastq #READ1_TRIM=20

# Step 3. Alignment
snaptools align-paired-end \
	--input-reference=${reference} \
	--input-fastq1=$tmpdir/unaligned_R1.fastq  \
	--input-fastq2=$tmpdir/unaligned_R2.fastq  \
	--output-bam=$outdir/Aligned.out.bam \
	--aligner=$(basename $bwa_exec) \
	--path-to-aligner=$(dirname $bwa_exec) \
	--min-cov=0 \
	--num-threads=16 \
	--if-sort=True \
	--tmp-folder=$tmpdir \
	--overwrite=TRUE

# Insert QC
samtools view $outdir/Aligned.out.bam |awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' | sort | uniq | cut -f 2 >$outdir/fragment.length.txt

# Step 4. Pre-processing.
samtools view $outdir/Aligned.out.bam -H  |cut -f 2,3 |grep -v VN |grep -v bwa |grep -v samtools |awk '{printf"%s\t%s\n", substr($1,4,length($1)-2), substr($2,4,length($2)-2)}' >$outdir/mix.chrom.sizes

snaptools snap-pre \
	--input-file=$outdir/Aligned.out.bam \
	--output-snap=$outdir/H.snap \
	--genome-name=hg19 \
	--genome-size=$outdir/mix.chrom.sizes \
	--overwrite=True \
	--min-mapq=30 \
	--min-flen=0 \
	--max-flen=1000 \
	--keep-chrm=True \
	--keep-single=True \
	--keep-secondary=True \
	--keep-discordant=True \
	--max-num=50000 \
	--min-cov=10 \
	--verbose=True

# add bmat
snaptools snap-add-bmat \
    --snap-file=$outdir/H.snap \
    --bin-size-list 500 1000 5000 10000 \
    --verbose=True

snaptools snap-add-pmat \
    --snap-file=H.snap \
	--peak-file=peaks/WOC.sorted_summits.bed \
	--verbose=True

# rm tmpdir if final output exist
# if [ -f $outdir/H.snap ]; then
# 	rm -r $tmpdir;
# fi


samtools sort -o WOC1_barcode_sort.bam WOC1_barcode.bam
samtools index WOC1_barcode_sort.bam
bamCoverage --bam WOC1_barcode_sort.bam -o WOC1_barcode_sort.SeqDepthNorm.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 1300000000

