import h5py
from sys import argv

snap_file, out_file = argv[1:]
f = h5py.File(snap_file, "r")

frag_list = []
barcode_id = 0
for barcode in f["BD"]["name"]:  
    _chroms = f["FM"]["fragChrom"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
    _chroms = [item.decode() for item in _chroms];
    _start = f["FM"]["fragStart"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
    _len = f["FM"]["fragLen"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
    _barcode = [barcode.decode()] * len(_chroms)

    frag_list_uniq = set(zip(_chroms, _start, _start + _len, _barcode)); # remove duplicated fragments
    frag_list.extend(list(frag_list_uniq))
    barcode_id += 1

# from snap import *
# from snap_pre import *

# frag_list = [];
# for read_group in group_reads_by_barcode('COL10_human.bam', 'bam'):
#   for (read1, read2, is_paired, is_secondary) in pairReadsByName(read_group):
#     if is_paired:
#       frag = readPairToFragment(read1, read2, is_secondary);
#     else: # supplementary alignments or unmated reads;
#       frag = readToFragment(read1, is_secondary);
#     barcode = frag.qname.split(":")[0].upper();
#     frag_list.append((read1, read2, frag.chrom, frag.pos, frag.pos+frag.flen, barcode))

  # frag_list_uniq = set(frag_list); # remove duplicated fragments
  # just in case the only fragments are chrM
  # if len(frag_list_uniq) == 0: continue;

def data_write_csv(file_name, datas):
    import codecs, csv
    file_csv = codecs.open(file_name,'w+','utf-8')
    writer = csv.writer(file_csv, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    for data in datas:
        writer.writerow(data)
    print("ok")

data_write_csv(out_file, frag_list)
