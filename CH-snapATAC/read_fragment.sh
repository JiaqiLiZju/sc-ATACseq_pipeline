samtools view Cerebellum_62216.bam |awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print $1"\t"abs($9)}' | sort | uniq | cut -f 2 >fragment.length.txt
