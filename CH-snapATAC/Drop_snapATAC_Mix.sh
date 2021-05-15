#!/bin/bash
cd /media/ggj/Files/scATAC/snapATAC/COL10

dropseq_root=~/jiaqiLi/general_global_soft/sc-software/Drop-seq_tools-1.12/
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
reference=/home/ggj/jiaqiLi/mount-1/ggj/jiaqiLi/general_global_soft/database_genome/Human-Mouse-Merged/hg19_mm10_transgenes.fasta
HY_barcode_file=/media/ggj/Files/scATAC/snapATAC/lig_768_bc.pickle2
RT_barcode_file=/media/ggj/Files/scATAC/snapATAC/RT_384_bc.pickle2
outdir=.
tmpdir=tmp
mkdir $tmpdir

# filter linker
~/jiaqiLi/general_global_soft/sc-software/bbmap/bbduk2.sh \
    in=${sample_name}_H_R1.fastq in2=${sample_name}_H_R2.fastq \
    outm=H_R1_linker1.fastq outm2=H_R2_linker1.fastq \
    fliteral=AAGCAGTGGTATCAACGCAGAGT k=23 skipr2=t hdist=3 -Xmx80g

~/jiaqiLi/general_global_soft/sc-software/bbmap/bbduk2.sh \
    in=H_R1_linker1.fastq in2=H_R2_linker1.fastq \
    outm=H_R1_linker2.fastq outm2=H_R2_linker2.fastq \
    fliteral=AGATGTGTATAAGAGACAG k=19 skipr2=t hdist=3 -Xmx80g  

# pre-alignment tag and trim
java -jar ${picard_jar} FastqToSam F1=H_R1_linker2.fastq F2=H_R2_linker2.fastq  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name 

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary.txt \
    BASE_RANGE=1-10:34-43 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary.txt \
    BASE_RANGE=1-10:34-43 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0

# barcode correct
samtools view -h $tmpdir/unaligned_tagged_Cell.bam > $tmpdir/unaligned_tagged_Cell.sam

python ../../CHseq_correctBC.py $tmpdir $tmpdir $HY_barcode_file $RT_barcode_file unaligned_tagged_Cell.correct.sam
samtools view -bS  $tmpdir/unaligned_tagged_Cell.correct.sam > $tmpdir/unaligned_tagged_Cell.correct.bam

# barcode qname
samtools view $tmpdir/unaligned_tagged_Cell.bam -H > $tmpdir/merged.header.sam
cat <( cat $tmpdir/merged.header.sam ) \
<( samtools view $tmpdir/unaligned_tagged_Cell.correct.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -bS - > $tmpdir/unaligned_tagged_Cell.snap.bam
# samtools view $tmpdir/unaligned_tagged_Cell.snap.bam | cut -f 1 | head
java -Xmx100g -jar ${picard_jar} SamToFastq INPUT=$tmpdir/unaligned_tagged_Cell.snap.bam READ1_TRIM=62 FASTQ=$tmpdir/unaligned_R1.fastq SECOND_END_FASTQ=$tmpdir/unaligned_R2.fastq

# Step 3. Alignment
snaptools align-paired-end \
	--input-reference=${reference} \
	--input-fastq1=$tmpdir/unaligned_R1.fastq \
	--input-fastq2=$tmpdir/unaligned_R2.fastq \
	--output-bam=$tmpdir/Aligned.out.bam \
	--aligner=bwa \
	--path-to-aligner=/usr/bin/ \
	--min-cov=0 \
	--num-threads=16 \
	--if-sort=True \
	--tmp-folder=./ \
	--overwrite=TRUE

# Fragment QC
samtools view $tmpdir/Aligned.out.bam |awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' | sort | uniq | cut -f 2 >fragment.length.txt

# FilterBAM
${dropseq_root}/FilterBAM I=$tmpdir/Aligned.out.bam O=$tmpdir/mouse.demo.bam REF_SOFT_MATCHED_RETAINED=MOUSE
${dropseq_root}/FilterBAM I=$tmpdir/Aligned.out.bam O=$tmpdir/human.demo.bam REF_SOFT_MATCHED_RETAINED=HUMAN

cat mouse.demo.header.sam |cut -f 2,3 |grep -v VN |grep -v bwa|awk '{printf"%s\t%s\n", substr($1,4,length($1)-2), substr($2,4,length($2)-2)}' >mix.chrom.sizes

# Step 4. Pre-processing.
snaptools snap-pre \
    --input-file=$tmpdir/mouse.demo.bam \
	--output-snap=$tmpdir/mouse.demo.snap \
	--genome-name=mm10 \
	--genome-size=./mix.chrom.sizes \
	--overwrite=True \
	--min-mapq=30 \
	--min-flen=0 \
	--max-flen=1000 \
	--keep-chrm=TRUE \
	--keep-single=TRUE \
	--keep-secondary=True \
    --keep-discordant=True \
	--max-num=1000000 \
	--min-cov=10 \
	--verbose=True

snaptools snap-pre \
    --input-file=$tmpdir/human.demo.bam \
	--output-snap=$tmpdir/human.demo.snap \
	--genome-name=hg19 \
	--genome-size=./mix.chrom.sizes \
	--min-mapq=30 \
	--min-flen=0 \
	--max-flen=1000 \
	--keep-chrm=TRUE \
	--keep-single=TRUE \
	--keep-secondary=True \
    --keep-discordant=True \
	--overwrite=True \
	--max-num=1000000 \
	--min-cov=10 \
	--verbose=True

snaptools snap-add-bmat \
    --snap-file=$tmpdir/mouse.demo.snap \
    --bin-size-list 5000 10000 \
    --verbose=True

snaptools snap-add-bmat \
    --snap-file=$tmpdir/human.demo.snap \
    --bin-size-list 5000 10000 \
    --verbose=True
