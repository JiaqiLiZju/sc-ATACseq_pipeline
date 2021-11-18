#BSUB -q normal
#BSUB -J scATAC
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
cellranger-atac mkref --config=./Earthworm.config

# cd /share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed

{
    organism: "Smed"
    genome: ["Smed"]
    input_fasta: ["/share/home/guoguoji/tools/BWA_Reference_Smed/SmedAsxl_genome_v1.1.nt"]
    input_gtf: ["/share/home/guoguoji/tools/CellRanger_ATAC_Reference_Smed/smed.gtf"]
}

cellranger-atac mkref --config=Smed.contig
