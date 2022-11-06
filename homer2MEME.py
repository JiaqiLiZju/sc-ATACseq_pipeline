import os
import glob
import pandas as pd

out_str = ''
d = []
for fname in glob.glob("./*motif*"):
    for line in open(fname, 'r'):
        if line.startswith(">"):
            if len(d) > 0:
                out_str = out_str + '\n' + motif + '\n' + \
                pd.DataFrame(d, columns=["A:", "C:", "G:", "T:"]) \
                    .T \
                    .to_csv(sep='\t', header=False)
            
            motif = line.rstrip().split("\t")[0].replace(">", "") + "-Homer"
            d = []

        else:
            lines = line.rstrip().split("\t")
            d.append([float(i) for i in lines])

out_str = out_str + '\n' + motif + '\n' + \
            pd.DataFrame(d, columns=["A:", "C:", "G:", "T:"]) \
                .T \
                .to_csv(sep='\t', header=False)

with open("tmp_transfac.txt", 'w') as out_fh:
    out_fh.write(out_str)

os.system("uniprobe2meme tmp_transfac.txt 1> out.meme 2>log.txt")