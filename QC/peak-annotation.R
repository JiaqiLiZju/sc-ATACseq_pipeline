setwd("~/Desktop/merged-mixed-bam/")

# 然后再打开的命令中设置清华源信息
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 放在前面的话：一般要安装什么包直接搜索包名，对应的包的 manual 说怎么装就怎么装
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("ChIPseeker")
# BiocManager::install("clusterProfiler")

# 加载包
require(ChIPseeker)
require(clusterProfiler) 

require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# 读取当前目录 oldBedFiles 文件夹下的 Peak 文件
bedPeaksFile = './Mouse/peaks/mouse.sorted_peaks.narrowPeak'; 
peak <- readPeakFile( bedPeaksFile) 

# 查看 Peak 文件中染色体信息
seqlevels(peak)

# 过滤掉带有 Het 字眼的染色体
# 请留意 `grepl()` 函数
# keepChr = !grepl('Het',seqlevels(peak)) 
# seqlevels(peak, pruning.mode = "coarse") <- seqlevels(peak)[keepChr]

## ---------------------------------------------------------------------
## 重命名 TxDb 对象的染色体信息
## ---------------------------------------------------------------------
seqlevels(txdb)

# 使用 `sub()` 函数将染色体前缀 `chr` 替换
seqlevels(txdb) <- sub("chr", "MOUSE_", seqlevels(txdb))
seqlevels(txdb)

# 使用 `seqlevels0()` 函数恢复原来的染色体水平
# seqlevels(txdb) <- seqlevels0(txdb)
# seqlevels(txdb)

# 使用 annotatPeak 进行注释，
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), 
                         TxDb = txdb) 

# 转变成 data.frame 格式文件，方便查看与后续操作
peakAnno_df <- as.data.frame(peakAnno)
write.csv(peakAnno_df, file = "Peak.anno.Mouse.csv")

pdf("Peak.anno.Mouse.pdf", width = 10, height = 4)
plotAnnoPie(peakAnno)
dev.off()

plotAnnoBar(peakAnno)

###########################################################################
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# 读取当前目录 oldBedFiles 文件夹下的 Peak 文件
bedPeaksFile = './Human/peaks/human.sorted_peaks.narrowPeak'; 
peak <- readPeakFile( bedPeaksFile) 

# 查看 Peak 文件中染色体信息
seqlevels(peak)

# 过滤掉带有 Het 字眼的染色体
# 请留意 `grepl()` 函数
# keepChr = !grepl('Het',seqlevels(peak)) 
# seqlevels(peak, pruning.mode = "coarse") <- seqlevels(peak)[keepChr]

## ---------------------------------------------------------------------
## 重命名 TxDb 对象的染色体信息
## ---------------------------------------------------------------------
seqlevels(txdb)

# 使用 `sub()` 函数将染色体前缀 `chr` 替换
seqlevels(txdb) <- sub("chr", "HUMAN_", seqlevels(txdb))
seqlevels(txdb)

# 使用 `seqlevels0()` 函数恢复原来的染色体水平
# seqlevels(txdb) <- seqlevels0(txdb)
# seqlevels(txdb)

# 使用 annotatPeak 进行注释，
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), 
                         TxDb = txdb) 

# 转变成 data.frame 格式文件，方便查看与后续操作
peakAnno_df <- as.data.frame(peakAnno)
write.csv(peakAnno_df, file = "Peak.anno.Human.csv")

pdf("Peak.anno.Human.pdf", width = 10, height = 4)
plotAnnoPie(peakAnno)
dev.off()

plotAnnoBar(peakAnno)
