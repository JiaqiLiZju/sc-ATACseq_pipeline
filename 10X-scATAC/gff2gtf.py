import sys

chrom_size = {}
for line in open("SmedAsxl_genome_v1.1.nt.fai", 'r'):
    chrom, size = line.strip().split("\t")[:2]
    chrom_size[chrom] = int(size)

inFile = open(sys.argv[1],'r')
for line in inFile:
    #skip comment lines that start with the '#' character
    if line[0] != '#':
        #split line into columns by tab
        data = line.strip().split('\t')
        if len(data) < 3: continue
        end, chrom = int(data[4]), data[0]

        if data[2] == "gene" and end > chrom_size[chrom]: print(chrom, file=sys.stderr); FLAG=True; continue
        if data[2] == "gene" and end < chrom_size[chrom]: FLAG=False
        if FLAG: continue

        #if the feature is a gene 
        if data[2] == "gene":
            # get the id
            GENE_ID = data[-1].split('ID=')[-1].split(';')[0]

            #modify the last column
            data[-1] = 'gene_id "' + GENE_ID + '";'# transcript_id "' + GENE_ID + '"'
            #data.append(GENE_ID)
            #print out this new GTF line
            print( '\t'.join(data[:9]))
            
        #if the feature is anything else
        else:
            if data[2] == "mRNA": 
                data[2] = "transcript"
                Trans_ID = data[-1].split('ID=')[-1].split(';')[0]

            else:
                # get the parent as the ID
                # ID = data[-1].split('Parent=')[-1].split(';')[0].replace("model", "TU")
                Trans_ID = data[-1].split('Parent=')[-1].split(';')[0]#.replace("mRNA", "gene")

            #modify the last column
            data[-1] = 'gene_id "' + GENE_ID + '"; transcript_id "' + Trans_ID + '"'
            #data.append(GENE_ID)
            #print out this new GTF line
            print( '\t'.join(data[:9]))
