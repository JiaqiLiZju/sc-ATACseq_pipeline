from snap import *
from snap_pre import *

frag_list = [];
for read_group in group_reads_by_barcode('COL10_Aligned.out.bam', 'bam'):
    for (read1, read2, is_paired, is_secondary) in pairReadsByName(read_group):
        if is_paired:
            frag = readPairToFragment(read1, read2, is_secondary);
            logstr = ','.join((str(frag.qname), str(frag.chrom), str(frag.pos), str(frag.flen), str(frag.mapq), str(frag.is_single), str(frag.is_secondary), str(frag.is_proper_pair), str(is_paired)))
            print(logstr)        

        else: # supplementary alignments or unmated reads;
            frag = readToFragment(read1, is_secondary);
            logstr = ','.join((str(frag.qname), str(frag.chrom), str(frag.pos), str(frag.flen), str(frag.mapq), str(frag.is_single), str(frag.is_secondary), str(frag.is_proper_pair), str(is_paired)))
            print(logstr)        

        barcode = frag.qname.split(":")[0].upper();
        frag_list.append((frag.chrom, frag.pos, frag.pos+frag.flen, barcode))
        print('\t'.join((str(frag.chrom), str(frag.pos), str(frag.pos+frag.flen), barcode)), file=sys.stderr)

frag_list = set(frag_list); # remove duplicated fragments just in case the only fragments are chrM
# if len(frag_list_uniq) == 0: continue;

def data_write_csv(file_name, datas):
    import codecs, csv
    file_csv = codecs.open(file_name,'w+','utf-8')
    writer = csv.writer(file_csv, delimiter='\t', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    for data in datas:
        writer.writerow(data)
    print("ok")

data_write_csv("COL10_Aligned.bam.bed", frag_list)
