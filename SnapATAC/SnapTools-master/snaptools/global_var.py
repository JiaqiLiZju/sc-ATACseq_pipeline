MIN_BIN_SIZE=1000
MIN_FRAGMENT_LEN=0
MAX_FRAGMENT_LEN=2000
MIN_MAPQ=0
MAGIC_STRING="SNAP"

# this is a predfined genome list
GENOMELIST=["hg38", "hg19", "hg18", "hg17", "hg16", "mm10", "mm9", "mm8", "mm7", 
            "vicPac2", "vicPac1", "xenLae2", "allMis1", "dasNov3", "gadMor1",
            "papAnu2", "papHam1", "bisBis1", "panPan2", "panPan1", "aptMan1",
            "melUnd1", "otoGar3", "felCat8", "felCat5", "felCat4", "felCat3",
            "felCat3", "galGal5", "galGal4", "galGal3", "galGal2", "panTro5",
            "panTro4", "panTro3", "panTro2", "panTro1", "criGri1", "criGriChoV1",
            "manPen1", "latCha1", "bosTau8", "bosTau7", "bosTau6", "bosTau4",
            "bosTau3", "bosTau2", "macFas5", "canFam3", "canFam2", "canFam1", 
            "turTru2", "loxAfr3", "loxAfr1", "calMil1", "musFur1", "fr3",
            "fr2", "fr1", "nomLeu3", "nomLeu2", "nomLeu1", "aquChr2", "rhiRox1",
            "gorGor5", "gorGor4", "gorGor3", "chlSab2", "cavPor3", "eriEur2",
            "eriEur1", "equCab2", "equCab1", "dipOrd1", "petMar3", "petMar2",
            "petMar1", "braFlo1", "anoCar2", "anoCar1", "galVar1", "triMan1",
            "calJac3", "calJac1", "oryLat2", "geoFor1", "pteVam1", "myoLuc2",
            "balAcu1", "micMur2", "micMur1", "hetGla2", "hetGla1", "oreNil2",
            "monDom5", "monDom4", "monDom1", "ponAbe2", "chrPic1", "ailMel1",
            "susScr11", "susScr3", "susScr2", "ochPri3", "ochPri2", "ornAna2",
            "ornAna1", "nasLar1", "oryCun2", "rn6", "rn5", "rn4", "rn3",
            "rheMac8", "rheMac3", "rheMac2", "proCap1", "oviAri3", "oviAri1",
            "sorAra2", "sorAra1", "choHof1", "speTri2", "saiBol1", "gasAcu1",
            "tarSyr2", "tarSyr1", "sarHar1", "echTel2", "echTel1", "tetNig2",
            "tetNig1", "nanPar1", "tupBel1", "melGal5", "melGal1", "macEug2",
            "cerSim1", "xenTro9", "xenTro7", "xenTro3", "xenTro2", "xenTro1",
            "taeGut2", "taeGut1", "danRer11", "danRer10", "danRer7", "danRer7",
            "danRer6", "danRer5",  "danRer4", "danRer3", "ci3", "ci2", "ci1",
            "strPur2", "strPur1", "anoGam1", "apiMel3", "apiMel2", "apiMel1",
            "droAna2", "droAna1", "droEre1", "droGri1", "dm6", "dm3", "dm2", 
            "dm1", "droMoj2", "droMoj1", "droPer1", "dp3", "dp2", "droSec1", "droSim1",
            "droVir2", "droVir1", "droYak2", "droYak1", "caePb2", "caePb1",
            "cb3", "ce10", "ce6", "ce4", "ce2", "caeJap1", "caeRem3", "caeRem2",
            "priPac1", "sacCer3", "sacCer2", "sacCer1", "aplCal1", "eboVir3"]

