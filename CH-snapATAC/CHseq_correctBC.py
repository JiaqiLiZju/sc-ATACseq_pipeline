# import subprocess
import sys
# import gzip
# from multiprocessing import Pool
# from functools import partial
import pickle
# import re


def correct_barcodes(input_folder, output_folder, HY_barcode_file, RT_barcode_file,samfile):
	# load the HY_barcode_list and HY_barcode_list
	print("Load Hybridization barcode dictionary...")
	# generate the ligation barcode list
	barcodes = open(HY_barcode_file, "rb")
	HY_barcode_list = pickle.load(barcodes)
	barcodes.close()
	print("Load RT barcode dictionary...")    
	# generate the RT barcode list:
	barcodes = open(RT_barcode_file, "rb")
	RT_barcode_list = pickle.load(barcodes)
	barcodes.close()
	print("Start to correct barcodes...") 
	# open star.sam file and the output file
	alignedsam=input_folder + "/"+ samfile
	output_file=output_folder+"/"+"unaligned_tagged_Cell.corrected.sam"
	#mismatch_rate = int(mismatch_rate)
	f1 = open(alignedsam)
	f2 = open(output_file,"w")	
	#f3 = open("/home/ggj/NEW/YF/1101_tissues/Mnu1020-Micelibnu/COL15/test/only_hy1.sam","wb")
	#f4 = open("/home/ggj/NEW/YF/1101_tissues/Mnu1020-Micelibnu/COL15/test/hyrt1.sam","wb")
	total_line = 0 
	filtered_line = 0	
	#while(line1[0]!="@")
	line=f1.readline()
	while(line):
		total_line+=1
		if(line[0]=="@"):
			f2.write(line)
			#f3.write(line)
			#f4.write(line)
			line=f1.readline()
		else:
			# find the rt and hy barcode
			barcode=line[line.rfind("CB:Z:"):]
			barcode1=line[:line.rfind("CB:Z:")]
			barcode2=line[line.rfind("CB:Z:"):]
			hybarcode_tmp=barcode[5:15]
			rtbarcode_tmp=barcode[15:25]
			# first check if the hy barcode match with the barcode
			if rtbarcode_tmp in RT_barcode_list and hybarcode_tmp in HY_barcode_list:
				filtered_line+=1
				rtbarcode_match=RT_barcode_list[rtbarcode_tmp]
				hybarcode_match=HY_barcode_list[hybarcode_tmp]
				line_write=barcode1+"CB:Z:"+hybarcode_match+rtbarcode_match+'\t'+barcode2[26:]
				f2.writelines(line_write)
				line=f1.readline()
				#line_write=re.sub(rtbarcode_tmp,rtbarcode_match,line)
				#f3.writelines(line_write)
				#line_new=re.sub(hybarcode_tmp,hybarcode_match,line_write)
				#f4.writelines(line_new)
				#strinfo = re.compile(rtbarcode_tmp)
				#line_write = strinfo.sub(rtbarcode_match,line)
				#strinfo = re.compile(hybarcode_tmp)
				#line_write = strinfo.sub(hybarcode_match,line_write)
				#f2.writelines(line_write)
				#rtbarcode_match=RT_barcode_list[rtbarcode_tmp]
				#line=line.replace(rtbarcode_tmp,rtbarcode_match)
				# first check if the rt barcode match with the barcode
				#if hybarcode_tmp in HY_barcode_list:
				#	filtered_line+=1
				#	hybarcode_match=HY_barcode_list[hybarcode_tmp]
				#	strinfo1 = re.compile(hybarcode_tmp)
				#	line_write = strinfo1.sub(hybarcode_match,line_write)
					#hybarcode_match=HY_barcode_list[hybarcode_tmp]
					#line=line.replace(hybarcode_tmp,hybarcode_match)
				#	f2.write(line_write)
				#	line=f1.readline()
				#else:
				#	line=f1.readline()
			else: 
				line=f1.readline()
	f1.close()
	f2.close()
	#f3.close()
	#f4.close()
	print("total line: %f, filtered line: %f, filter rate: %f"%(total_line, filtered_line, float(filtered_line) / float(total_line)))
	com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
	print(com_message)
	
#def read_barcodefiles(HY_barcode_file,RT_barcode_file):
#	print("Load Hybridization barcode dictionary...")
#	 # generate the ligation barcode list
#    barcodes = open(HY_barcode_file, "rb")
#    HY_barcode_list = pickle.load(barcodes)
#    barcodes.close()
#    print("Load RT barcode dictionary...")    
#    # generate the RT barcode list:
#    barcodes3 = open(RT_barcode_file, "rb")
#    RT_barcode_list = pickle.load(barcodes3)
#    barcodes3.close()   
#    # parallele for the functions
#    p = Pool(processes = int(core))
#    func = partial(correct_barcodes, input_folder = input_folder, output_folder=output_folder, HY_barcode_list = HY_barcode_list, RT_barcode_list=RT_barcode_list, mismatch_rate = 1)
#    result=p.map()
#    p.close()
#    p.join()
#    com_message = '''~~~~~~~~~~~~~~~correct barcodes done~~~~~~~~~~~~~~~~~~'''
#    print(com_message)	
	
if __name__ == "__main__":
    input_folder = sys.argv[1]
    #sampleID = sys.argv[2]
    output_folder = sys.argv[2]
    HY_barcode_file = sys.argv[3]
    RT_barcode_file = sys.argv[4]
    samfile = sys.argv[5]
    #mismatch_rate = sys.argv[6]
    correct_barcodes(input_folder,output_folder, HY_barcode_file, RT_barcode_file,samfile)
