#BSUB -q normal
#BSUB -J CHATAC
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

reference=/share/home/guoguoji/tools/STAR_Reference_Axolotl/AmexG_v3.0.0.fa
HY_barcode_file=/share/home/guoguoji/tools/CHseq/barcode/lig_768_bc.pickle2
RT_barcode_file=/share/home/guoguoji/tools/CHseq/barcode/RT_384_bc.pickle2
barcode_file=/share/home/guoguoji/tools/CHseq/barcode/barcode_29K.txt

bbduk2=/share/home/guoguoji/tools/bbmap/bbduk2.sh
dropseq_root=/share/home/guoguoji/tools/Drop-seq_tools-1.12/
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar
bwa_exec=/share/home/guoguoji/tools/bwa-0.7.15/bwa
correctBC_script=/share/home/guoguoji/JiaqiLi/tools/CHATAC/CHseq_correctBC.py

sample_name=$(basename `pwd`)
outdir=bowtie_out
tmpdir=tmp

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

if [ ! -d $tmpdir ]; then
    mkdir $tmpdir
fi

for i in $(ls *fastq.gz);do gunzip $i;done

cp ${sample_name}_H_R1.fastq $tmpdir/H_R1_linker2.fastq
cp ${sample_name}_H_R2.fastq $tmpdir/H_R2_linker2.fastq

# filter linker
# ${bbduk2} in=${sample_name}_H_R1.fastq in2=${sample_name}_H_R2.fastq \
#     outm=$tmpdir/H_R1_linker1.fastq outm2=$tmpdir/H_R2_linker1.fastq \
#     fliteral=AAGCAGTGGTATCAACGCAGAGT k=23 skipr2=t hdist=3 -Xmx80g

# ${bbduk2} in=$tmpdir/H_R1_linker1.fastq in2=$tmpdir/H_R2_linker1.fastq \
#     outm=$tmpdir/H_R1_linker2.fastq outm2=$tmpdir/H_R2_linker2.fastq \
#     fliteral=AGATGTGTATAAGAGACAG k=19 skipr2=t hdist=3 -Xmx80g  

# pre-alignment tag and trim
java -jar ${picard_jar} FastqToSam F1=$tmpdir/H_R1_linker2.fastq F2=$tmpdir/H_R2_linker2.fastq O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name 

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R1.txt \
    BASE_RANGE=1-20 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=true DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=H.bam OUTPUT=$tmpdir/unaligned_tagged_Cell_R1.bam COMPRESSION_LEVEL=0

${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary_R2.txt \
    BASE_RANGE=1-20 BASE_QUALITY=10 BARCODED_READ=1 TAG_BARCODED_READ=false DISCARD_READ=false TAG_NAME=CB NUM_BASES_BELOW_QUALITY=1 \
	INPUT=$tmpdir/unaligned_tagged_Cell_R1.bam OUTPUT=$tmpdir/unaligned_tagged_Cell.bam COMPRESSION_LEVEL=0

# barcode correct
samtools view -h $tmpdir/unaligned_tagged_Cell.bam > $tmpdir/unaligned_tagged_Cell.sam
python $correctBC_script $tmpdir $tmpdir $HY_barcode_file $RT_barcode_file unaligned_tagged_Cell.sam
samtools view -bS $tmpdir/unaligned_tagged_Cell.corrected.sam > $tmpdir/unaligned_tagged_Cell.correct.bam

# barcode qname
samtools view $tmpdir/unaligned_tagged_Cell.correct.bam -H > $tmpdir/unaligned_tagged_Cell.snap.sam
samtools view $tmpdir/unaligned_tagged_Cell.correct.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); printf "%s:%s\n", td["CB"], $0 } } }' >> $tmpdir/unaligned_tagged_Cell.snap.sam

# samtools view $tmpdir/unaligned_tagged_Cell.snap.bam | cut -f 1 | head
java -Xmx100g -jar ${picard_jar} SamToFastq INPUT=$tmpdir/unaligned_tagged_Cell.snap.sam FASTQ=$tmpdir/unaligned_R1.fastq SECOND_END_FASTQ=$tmpdir/unaligned_R2.fastq #READ1_TRIM=20

echo -e "adaMis=2\ntrim=20,0,0,0" >$tmpdir/soapnuke.config
# java -jar trimmomatic.jar PE -phred33 $tmpdir/unaligned_R1.fastq $tmpdir/unaligned_R2.fastq -baseout unaligned.fq.gz ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50
/share/home/guoguoji/tools/SOAPnuke/SOAPnuke filter -1 $tmpdir/unaligned_R1.fastq -2 $tmpdir/unaligned_R2.fastq -o $tmpdir -C unaligned_1P.fq.gz -D unaligned_2P.fq.gz -f CTGTCTCTTATACACATCT -r CTGTCTCTTATACACATCT -J -c $tmpdir/soapnuke.config

# Step 3. Alignment
snaptools align-paired-end \
	--input-reference=${reference} \
	--input-fastq1=$tmpdir/unaligned_1P.fq.gz  \
	--input-fastq2=$tmpdir/unaligned_2P.fq.gz  \
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
	--max-num=5000 \ # num_core_barcode
	--min-cov=10 \
	--verbose=True

# add bmat
snaptools snap-add-bmat \
    --snap-file=$outdir/H.snap \
    --bin-size-list 5000 10000 20000 40000 50000 100000 \
    --verbose=True

# rm tmpdir if final output exist
if [ -f $outdir/H.snap ]; then
	gzip ${sample_name}_H_R1.fastq;
	gzip ${sample_name}_H_R2.fastq;
	rm -r $tmpdir;
fi