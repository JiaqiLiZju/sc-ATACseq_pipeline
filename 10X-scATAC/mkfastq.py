import pysam
from copy import deepcopy

# filename = "H_R1_test.fastq"
# with pysam.FastxFile(filename) as fin, open("FASTQ/EW_S1_L001_R2_001.fastq", mode='w') as fbc, open("FASTQ/EW_S1_L001_R1_001.fastq", mode='w') as fr1:
#     for idx, entry in enumerate(fin):
#         entry_bc = deepcopy(entry)
#         entry_bc.name = "sra." + str(idx+1) + " " + entry_bc.name
#         entry_bc.comment = "length=16"
#         entry_bc.quality = entry_bc.quality[:16]
#         entry_bc.sequence = entry_bc.sequence[:16]
#         fbc.write(str(entry_bc) + '\n')

#         entry.name = "sra." + str(idx+1) + " " + entry.name
#         entry.quality = entry.quality[16:]
#         entry.sequence = entry.sequence[16:]
#         fr1.write(str(entry) + '\n')

# filename = "H_R2_test.fastq"
# with pysam.FastxFile(filename) as fin, open("FASTQ/EW_S1_L001_R3_001.fastq", mode='w') as fr2:
#     for idx, entry in enumerate(fin):
#         entry.name = "sra." + str(idx+1) + " " + entry.name
#         fr2.write(str(entry) + '\n')

#rev_comp
def _rev_comp(seq):
    trantab = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(trantab)[::-1]


filename = "Test_R2.fastq"
with pysam.FastxFile(filename) as fin, open("FASTQ/Sample_S1_L001_R2_001.fastq", mode='w') as fbc, open("FASTQ/Sample_S1_L001_R3_001.fastq", mode='w') as fr1:
    for idx, entry in enumerate(fin):
        entry_bc = deepcopy(entry)
        name = entry_bc.name.replace("/2", "")

        entry_bc.name = name
        # entry_bc.comment = "length=16"
        entry_bc.quality = entry_bc.quality[-16:][::-1]
        entry_bc.sequence = _rev_comp(entry_bc.sequence[-16:])
        fbc.write(str(entry_bc) + '\n')

        entry.name = name
        entry.quality = entry.quality[:-16]
        entry.sequence = entry.sequence[:-16]
        fr1.write(str(entry) + '\n')

filename = "Test_R1.fastq"
with pysam.FastxFile(filename) as fin, open("FASTQ/Sample_S1_L001_R1_001.fastq", mode='w') as fr2:
    for idx, entry in enumerate(fin):
        name = entry.name.replace("/1", "")
        entry.name = name
        fr2.write(str(entry) + '\n')
