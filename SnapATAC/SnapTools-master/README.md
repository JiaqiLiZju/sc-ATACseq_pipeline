[![Build Status](https://travis-ci.org/r3fang/SnapTools.svg?branch=master)](https://travis-ci.org/r3fang/SnapTools)

## SnapTools
A module for working with snap files in Python.

## Introduction
snap (Single Nucleus Accessibility Profile) file is a hierarchically structured hdf5 file that is specially designed for storing single nucleus ATAC-seq datasets. A snap file (version 4) contains the following sessions: header (HD), cell-by-bin accessibility matrix (AM), cell-by-peak matrix (PM), cell-by-gene matrix (GM), barcode (BD) and fragment (FM). 

* HD session contains snap-file version, created date, alignment and reference genome information. 
* BD session contains all unique barcodes and corresponding meta data. 
* AM session contains cell-by-bin matrices of different resolutions (or bin sizes). 
* PM session contains cell-by-peak count matrix. 
* GM session contains cell-by-gene count matrix. 
* FM session contains all usable fragments for each cell. Fragments are indexed for fast search. 
* Detailed information about snap file can be found [here](https://github.com/r3fang/SnapTools/blob/master/docs/snap_format.docx).

## Requirements 
* Python (both python2 and python3)
* pysam
* h5py
* numpy
* pybedtools

## Latest News
* add `snap-del` which deletes a session in snap file.
* version 1.4.+ now supports both python2 and python3.

## Quick Install 
Install snaptools from PyPI

```
$ pip install snaptools
```

Install snaptools from source code

```bash
$ git clone https://github.com/r3fang/snaptools.git
$ cd snaptools
$ pip install -e .
$ ./bin/snaptools
usage: snaptools [-h]  ...

Program: snaptools (A module for working with snap files in Python)
Version: 1.4.1
Contact: Rongxin Fang
E-mail:  r4fang@gmail.com

optional arguments:
  -h, --help        show this help message and exit

functions:

    dex-fastq       De-multicomplex fastq file.
    index-genome    Index reference genome.
    align-paired-end
                    Align paired-end reads.
    align-single-end
                    Align single-end reads.
    snap-pre        Create a snap file from bam or bed file.
    snap-add-bmat   Add cell x bin count matrix to snap file.
    snap-add-pmat   Add cell x peak count matrix to snap file.
    snap-add-gmat   Add cell x gene count matrix to snap file.
    snap-del        Delete a session.
```

## Example

**Step 1. Download test example**

```
$ wget http://renlab.sdsc.edu/r3fang/share/SnapTools/snaptools_test.tar.gz
$ tar -xf snaptools_test.tar.gz
$ cd snaptools_test/
$ gunzip mm10.fa.gz
```

**Step 2. Index Reference Genome (Optional)**. 
Index the reference genome before alingment if you do not have one. (skip this step if you already have indexed genome). Here we show how to index the genome using BWA. User can choose different `--aligner `. 

```bash
$ which bwa
/opt/biotools/bwa/bin/bwa
$ snaptools index-genome	\
	--input-fasta=mm10.fa	\
	--output-prefix=mm10	\
    --aligner=bwa	\
	--path-to-aligner=/opt/biotools/bwa/bin/	\
	--num-threads=5
```

**Step 3. Alignment**. 
We next align reads to the corresponding reference genome using snaptools with following command. After alignment, reads are sorted by the read names which allows for grouping reads according to the barcode (`--if-sort`). User can mutiple CPUs to speed up this step (`--num-threads`).

```bash
$ snaptools align-paired-end	\
	--input-reference=mm10.fa	\
	--input-fastq1=demo.R1.fastq.gz	\
	--input-fastq2=demo.R2.fastq.gz	\
	--output-bam=demo.bam	\
	--aligner=bwa	\
	--path-to-aligner=/opt/biotools/bwa/bin/	\
	--read-fastq-command=zcat	\
	--min-cov=0	\
	--num-threads=5	\
	--if-sort=True	\
	--tmp-folder=./	\
	--overwrite=TRUE                     
```

**Step 4. Pre-processing**.         
After alignment, we converted pair-end reads into fragments and for each fragment, we check the following attributes: 1) mapping quality score MAPQ; 2) whether two ends are appropriately paired according to the alignment flag information; 3) fragment length. We only keep the properly paired fragments whose MAPQ (`--min-mapq`) is greater than 30 with fragment length less than 1000bp (`--max-flen`). Because the reads have been sorted based on the names, fragments belonging to the same cell (or barcode) are naturally grouped together which allows for removing PCR duplicates. After alignment and filtration, we generated a snap-format (Single-Nucleus Accessibility Profiles) file that contains meta data, cell-by-bin count matrices of a variety of resolutions, cell-by-peak count matrix. Detailed information about snap file can be found in here. 

```bash
$ snaptools snap-pre  \
	--input-file=demo.bam  \
	--output-snap=demo.snap  \
	--genome-name=mm10  \
	--genome-size=mm10.chrom.size  \
	--min-mapq=30  \
	--min-flen=0  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=TRUE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=1000000  \
	--min-cov=100  \
	--verbose=True
```

This command creates two files `demo.snap` and `demo.snap.qc` which contains the library quality control metrics as shown below.

```bash
$ cat demo.snap.qc

Total number of unique barcodes:             3217
TN - Total number of fragments:              576676
UM - Total number of uniquely mapped:        540307
SE - Total number of single ends:            0
SA - Total number of secondary alignments:   1
PE - Total number of paired ends:            540306
PP - Total number of proper paired:          539772
PL - Total number of proper frag len:        539772
US - Total number of usable fragments:       539772
UQ - Total number of unique fragments:       537336
CM - Total number of chrM fragments:         0
```

**Step 5. Cell-by-Bin Matrix**.            
Using generated snap file, we next create the cell-by-bin matrix. Snap file allows for storing cell-by-bin matrices of different resolutions. In the below example, three cell-by-bin matrices are created with bin size of 5,000 and 10,000. The cell-by-bin matrices will be added to `demo.snap` without creating another file. Same with `snap-add-pmat` and `snap-add-gmat`.

```bash
$ snaptools snap-add-bmat  \
	--snap-file=demo.snap  \
	--bin-size-list 5000 10000  \
	--verbose=True
```
