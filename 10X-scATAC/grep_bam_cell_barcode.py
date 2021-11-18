import pysam
from sys import argv

sam_fname, barcode_fname = argv[1:]

barcode_list = []
with open(barcode_fname, "r") as fh:
    for line in fh:
        barcode_list.append(line.strip())

samfile = pysam.AlignmentFile(sam_fname, "rb")
pairedreads = pysam.AlignmentFile("WOC1_barcode.bam", "wb", template=samfile)
for read in samfile:
    barcode = read.qname.split(":")[0]
    if barcode in barcode_list:
        pairedreads.write(read)

pairedreads.close()
samfile.close()
