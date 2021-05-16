#BSUB -q normal
#BSUB -J axo118
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 36

bwa_exec=/share/home/guoguoji/tools/bwa-0.7.15/bwa
reference=/share/home/guoguoji/tools/BWA-Reference-Human-Mouse-Merged/hg19_mm10_transgenes.fasta

# dict
# samtools dict $reference

# chrom size
# fetchChromSizes mm10.fa > mm10.chrom.sizes

# bowtie index
snaptools index-genome \
	--input-fasta=${reference} \
	--output-prefix=${reference} \
    --aligner=$(basename $bwa_exec) \
	--path-to-aligner=$(dirname $bwa_exec) \
	--num-threads=10 \
	&>log.reference.bowtie
