import pysam
import pandas as pd
import numpy as np
from collections import OrderedDict

barcode_human_cnt = OrderedDict()
barcode_mouse_cnt = OrderedDict()

fh = pysam.AlignmentFile("Aligned.out.sam", "rb")
for _read in fh:
    if _read.has_tag('CB'):
        barcode = _read.get_tag('CB')
        reference = _read.reference_name
        if barcode not in barcode_human_cnt.keys():
            barcode_human_cnt[barcode] = 0
        if barcode not in barcode_mouse_cnt.keys():
            barcode_mouse_cnt[barcode] = 0
        if "HUMAN" in reference:
            barcode_human_cnt[barcode] += 1
        if "MOUSE" in reference:
            barcode_mouse_cnt[barcode] += 1
    else:
        continue

fh.close()

barcode_human = pd.DataFrame(barcode_human_cnt, index=['cnt_human']).T
barcode_mouse = pd.DataFrame(barcode_mouse_cnt, index=['cnt_mouse']).T
read_cnt = pd.merge(barcode_human, barcode_mouse, left_index=True, right_index=True)
read_cnt.to_csv("./read_cnt.csv")