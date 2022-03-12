#BSUB -q normal
#BSUB -J cellranger
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

### Earthworm.config
# {
#     organism: "Earthworm"
#     genome: ["Earthworm"]
#     input_fasta: ["/share/home/guoguoji/tools/STAR_Reference_Qiu/Earthworm.fasta"]
#     input_gtf: ["/share/home/guoguoji/tools/STAR_Reference_Qiu/Earthworm.gtf"]
# }
#     non_nuclear_contigs: ["chrM"]
#     input_motifs: "/path/to/jaspar/motifs.pfm"

# cd /share/home/guoguoji/tools/CellRanger_ATAC_Reference_Earthworm
# cellranger-atac mkref --config=./Earthworm.config

# cd /share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed
gffread gencode.v19.annotation.gff3 -T -o gencode.v19.gtf

{
    organism: "Smed"
    genome: ["Smed"]
    input_fasta: ["/share/home/guoguoji/tools/CellRanger_ATAC_Reference_dd_Smed_v4/dd_Smes_g4.fasta"]
    input_gtf: ["/share/home/guoguoji/tools/CellRanger_ATAC_Reference_dd_Smed_v4/smes_v2_hconf_SMESG.gtf"]
}

cellranger-atac mkref --config=Smed.config

# (base) [guoguoji@gpu05 CellRanger_ATAC_Reference_Smed]$ rm -r Smed && cellranger-atac mkref --config=Smed.contig

# >>> Creating reference for Smed <<<

# Creating new reference folder at /share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed/Smed
# ...done

# Writing genome FASTA file into reference folder...
# ...done

# Indexing genome FASTA file...
# ...done

# Writing genes GTF file into reference folder...


# mkref has FAILED
# Error building reference package
# Error while parsing GTF file /share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed/smed.gtf
# Invalid GTF annotation on line 8984:
# ['scaffold11194', 'ALN', 'transcript', '49735', '50682', '1433', '-', '.', 'gene_id "gene01135"; transcript_id "mRNA01755"']
# End position of feature = 50682 > contig scaffold11194 length = 28674

# Please fix your GTF and start again.

# (base) [guoguoji@gpu05 CellRanger_ATAC_Reference_Smed]$ rm -r Smed && cellranger-atac mkref --config=Smed.contig

# >>> Creating reference for Smed <<<

# Creating new reference folder at /share/home/guoguoji/tools/CellRanger_ATAC_Reference_dd_Smed_v4/Smed
# ...done

# Writing genome FASTA file into reference folder...
# ...done

# Indexing genome FASTA file...
# ...done

# Writing genes GTF file into reference folder...
# ...done

# Writing genome metadata JSON file into reference folder...
# Computing hash of genome FASTA file...
# ...done

# Computing hash of genes GTF file...
# ...done

# ...done

# Generating bwa index (may take over an hour for a 3Gb genome)...
# [bwa_index] Pack FASTA... 4.54 sec
# [bwa_index] Construct BWT for the packed sequence...
# [BWTIncCreate] textLength=1547878984, availableWord=120914104
# [BWTIncConstructFromPacked] 10 iterations done. 99999992 characters processed.
# [BWTIncConstructFromPacked] 20 iterations done. 199999992 characters processed.
# [BWTIncConstructFromPacked] 30 iterations done. 299999992 characters processed.
# [BWTIncConstructFromPacked] 40 iterations done. 399999992 characters processed.
# [BWTIncConstructFromPacked] 50 iterations done. 499999992 characters processed.
# [BWTIncConstructFromPacked] 60 iterations done. 599999992 characters processed.
# [BWTIncConstructFromPacked] 70 iterations done. 699999992 characters processed.
# [BWTIncConstructFromPacked] 80 iterations done. 799999992 characters processed.
# [BWTIncConstructFromPacked] 90 iterations done. 899511384 characters processed.
# [BWTIncConstructFromPacked] 100 iterations done. 990672200 characters processed.
# [BWTIncConstructFromPacked] 110 iterations done. 1071692152 characters processed.
# [BWTIncConstructFromPacked] 120 iterations done. 1143698872 characters processed.
# [BWTIncConstructFromPacked] 130 iterations done. 1207694632 characters processed.
# [BWTIncConstructFromPacked] 140 iterations done. 1264570200 characters processed.
# [BWTIncConstructFromPacked] 150 iterations done. 1315117368 characters processed.
# [BWTIncConstructFromPacked] 160 iterations done. 1360039768 characters processed.
# [BWTIncConstructFromPacked] 170 iterations done. 1399962904 characters processed.
# [BWTIncConstructFromPacked] 180 iterations done. 1435442664 characters processed.
# [BWTIncConstructFromPacked] 190 iterations done. 1466973128 characters processed.
# [BWTIncConstructFromPacked] 200 iterations done. 1494993448 characters processed.
# [BWTIncConstructFromPacked] 210 iterations done. 1519893976 characters processed.
# [BWTIncConstructFromPacked] 220 iterations done. 1542021560 characters processed.
# [bwt_gen] Finished constructing BWT in 223 iterations.
# [bwa_index] 467.09 seconds elapse.
# [bwa_index] Update BWT... 2.73 sec
# [bwa_index] Pack forward-only FASTA... 2.69 sec
# [bwa_index] Construct SA from BWT and Occ...

# 152.20 sec
# [main] Version: 0.7.17-r1188
# [main] CMD: bwa index /share/home/guoguoji/tools/CellRanger_ATAC_Reference_dd_Smed_v4/Smed/fasta/genome.fa
# [main] Real time: 630.353 sec; CPU: 629.254 sec
# done

# Writing TSS and transcripts bed file...


# mkref has FAILED
# Error building reference package
# Invalid gene annotation input: in GTF
# /share/home/guoguoji/tools/CellRanger_ATAC_Reference_dd_Smed_v4/Smed/genes/genes.gtf.gz
# records for gene_id SMESG000025367.1 are not contiguous in the file
