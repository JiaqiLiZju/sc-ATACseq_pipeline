"""
Created on 13 11 2019

@author: feilijiang

edit from UMI_barcode_attach_gzipped_with_dic.py of  junyue
"""
import subprocess
import sys
#from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial
import pickle
#import re



def correct_barcodes(input_folder, output_folder,barcodepath,samfile):
	# load the barcode_list
	print("Load barcode whitelist dictionary...")
	RT_bc1 = open(barcodepath+"/barcode1_96_bc.pickle2","rb")
	barcode_list1 = pickle.load(RT_bc1)
	RT_bc1.close()
	RT_bc2 = open(barcodepath+"/barcode2_96_bc.pickle2","rb")
	barcode_list2 = pickle.load(RT_bc2)
	RT_bc2.close()
	RT_bc3 = open(barcodepath+"/barcode3_96_bc.pickle2","rb")
	barcode_list3 = pickle.load(RT_bc3)
	RT_bc3.close()
	RT_bc4 = open(barcodepath+"/RT_384_bc.pickle2","rb")
	barcode_list4 = pickle.load(RT_bc4)
	RT_bc4.close()
	print("Start to correct barcodes...") 
	# open star.sam file and the output file
	samfile=input_folder + "/"+ samfile
	output_file=output_folder+"/"+"unaligned_tagged_filtered_corrected.sam"
	f1 = open(samfile)
	f2 = open(output_file,"w")	
	total_line = 0 
	filtered_line = 0	
	line=f1.readline()
	while(line):
		total_line+=1
		if(line[0]=="@"):
			f2.write(line)
			line=f1.readline()
		else:
			# find the 3 barcodes and hy barcode
			barcode=line[line.rfind("CB:Z:"):]
			ch1=line[:line.rfind("CB:Z:")]
			ch2=line[line.rfind("RG:Z:"):]
			barcode1=barcode[5:11]
			barcode2=barcode[11:17]
			barcode3=barcode[17:23]
			barcode4=barcode[23:33]
			# first check if the hy barcode match with the barcode
			if barcode1 in barcode_list1 and barcode2 in barcode_list2 and barcode3 in barcode_list3 and barcode4 in barcode_list4:
				filtered_line+=1
				rtbarcode_match1=barcode_list1[barcode1]
				rtbarcode_match2=barcode_list2[barcode2]
				rtbarcode_match3=barcode_list3[barcode3]
				rtbarcode_match4=barcode_list4[barcode4]
				line_write=ch1+"CB:Z:"+rtbarcode_match1+rtbarcode_match2+rtbarcode_match3+rtbarcode_match4+'\t'+ch2
				f2.writelines(line_write)
				line=f1.readline()
			else: 
				line=f1.readline()
	f1.close()
	f2.close()
	print("total line: %f, filtered line: %f, filter rate: %f"%(total_line, filtered_line, float(filtered_line) / float(total_line)))
	com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
	print(com_message)



if __name__ == "__main__":
    input_folder = sys.argv[1]
    #sampleID = sys.argv[2]
    output_folder = sys.argv[2]
    barcodepath = sys.argv[3]
    samfile = sys.argv[4]
    #mismatch_rate = sys.argv[6]
    correct_barcodes(input_folder,output_folder,barcodepath,samfile)
