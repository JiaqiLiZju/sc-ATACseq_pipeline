setwd("/media/ggj/Files/scATAC/CHATAC/Axolotl/")
library(SnapATAC);

# Step 1. Barcode selection
x.sp = createSnap(
  file="../Axolotl-batch2-20210604-AllPeak/COL100/bowtie_out/H.snap",
  sample="Axolotl",
  num.cores=10
);

x.sp #5000 cells
barcodes = x.sp@barcode
metadata = x.sp@metaData
sel = sort(metadata$UM, decreasing = T, index.return=T)$ix[1:1000]
metadata.sel = metadata[sel,]
head(metadata.sel)
summary(metadata.sel$UM) #Total number of uniquely mapped
summary(metadata.sel$UQ) #Total number of unique fragments
summary(metadata.sel$TN) #Total number of fragments

# write.csv(metadata.sel, file = "./Axolotl.summary_xin.txt")

x.sp.sel = x.sp[sel,];
# x.sp.sel@metaData = barcodes.sel[x.sp@barcode,];
x.sp <- x.sp.sel
x.sp

x.sp <- readRDS("./x.fater.sp-batch2-BS5k.rds")

# sed 's/"/\t/g' AmexG_v3.0.0.fa.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {if($7=="+") {start=$4-2000; end=$4+2000;} else {if($7=="-") start=$5-2000; end=$5+2000; } if(start<0) start=0; print $1,start,end,$14,$10,$7;}}' >AmexG.promoter.bed
promoters_list = read.csv("./AmexG.promoter.bed", sep='\t', header=F);
promoters_list.gr = GRanges(
  promoters_list[,1],
  IRanges(promoters_list[,2], promoters_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, promoters_list.gr));
if(length(idy) > 0){x.sp.promoter = x.sp.promoter[,-idy, mat="bmat"]};
x.sp.promoter

promoter.covs <- Matrix::rowSums(x.sp.promoter@bmat);
UMI.covs <- Matrix::rowSums(x.sp@bmat);
promoter_ratio <- promoter.covs / UMI.covs 
data = data.frame(UMI=log10(UMI.covs+1), promoter_ratio=promoter_ratio);
# row.names(data) <- x.sp@barcode

# barcodes = read.csv(
#   "atac_v1_adult_brain_fresh_5k_singlecell.csv",
#   head=TRUE
# );
# barcodes = barcodes[2:nrow(barcodes),];
# promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
# UMI = log(barcodes$passed_filters+1, 10);
# data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
# barcodes$promoter_ratio = promoter_ratio;

library(viridisLite);
library(ggplot2);
ggplot(
  data,
  aes(x= UMI, y= promoter_ratio)) +
  geom_point(size=0.1, col="grey") +
  theme_classic() +
  ggtitle("10X Fresh Adult Brain") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio")

barcodes.sel <- which(UMI >= 2 & UMI <= 6 & promoter_ratio >= 0.5 & promoter_ratio <= 0.9);

x.sp.sel = x.sp[barcodes.sel,];
x.sp.sel@metaData$UMI = data$UMI[barcodes.sel]
x.sp.sel@metaData$promoter_ratio = data$promoter_ratio[barcodes.sel]
x.sp.sel

# Step 2. Add cell-by-bin matrix
showBinSizes("../Axolotl-batch2-20210604-AllPeak/COL100/bowtie_out/H.snap");
x.sp = addBmatToSnap(x.sp, bin.size=50000, num.cores=4);
x.sp = addPmatToSnap(x.sp, num.cores=4);
x.sp
# Step 3. Matrix binarization
x.sp = makeBinary(x.sp, mat="bmat");
x.sp@bmat

# Step 4. Bin filtering
# system("wget https://www.axolotl-omics.org/dl/RM_all_repeats.bed.gz");
library(GenomicRanges);
black_list = read.table("mm10.blacklist.bed.gz");
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# bin filter
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 3)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.999);
bin.cutoff 
idy = which(bin.cov <= bin.cutoff & bin.cov > 0.1);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# Step 5. Dimensionality reduction
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);

# Step 6. Determine significant components
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);


# Step 7. Graph-based clustering
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:15,
  k=15
);
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);
x.sp@metaData$cluster = x.sp@cluster;

# Step 8. Visualization
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:15, 
  method="Rtsne",
  seed.use=10
);
# par(mfrow = c(2, 2));

plotViz(
  obj=x.sp,
  method="tsne", 
  main="CH-ATAC(Axolotl)",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);

plotFeatureSingle(
  obj=x.sp,
  feature.value=log(x.sp@metaData$TN+1,10),
  method="tsne", 
  main="Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
); 

#plotFeatureSingle(
#  obj=x.sp,
#  feature.value=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters,
#  method="tsne", 
#  main="FRiP",
#  point.size=0.2, 
#  point.shape=19, 
#  down.sample=10000,
#  quantiles=c(0.01, 0.99) # remove outliers
#);

#plotFeatureSingle(
#  obj=x.sp,
#  feature.value=x.sp@metaData$duplicate / x.sp@metaData$total,
#  method="tsne", 
#  main="Duplicate",
#  point.size=0.2, 
#  point.shape=19, 
#  down.sample=10000,
#  quantiles=c(0.01, 0.99) # remove outliers
#)

# Step 9. Gene based annotation
# system("wget http://renlab.sdsc.edu/r3fang/share/github/Mouse_Brain_10X/gencode.vM16.gene.bed");
genes = read.table("../../Mix.Human.Refseq.bed");
genes.gr = GRanges(genes[,1], 
                   IRanges(genes[,2], genes[,3]), name=genes[,4]
);
marker.genes = c(
  "NM_001101", "NM_002046"
);
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
genes.sel.gr

# re-add the cell-by-bin matrix to the snap object;
# x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=10
);
# normalize the cell-by-gene matrix
x.sp = scaleCountMatrix(
  obj=x.sp, 
  cov=x.sp@metaData$passed_filters + 1,
  mat="gmat",
  method = "RPM"
);
# smooth the cell-by-gene matrix
x.sp = runMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);
par(mfrow = c(3, 3));
for(i in 1:9){
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[, marker.genes[i]],
    method="tsne", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0, 1)
  )}

# Step 10. Heretical clustering
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
plotViz(
  obj=x.sp,
  method="tsne", 
  main="10X Brain Cluster",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
);
plot(hc, hang=-1, xlab="");

# Step 11. Identify peaks
runMACS(
  obj=x.sp[which(x.sp@cluster==1),], 
  output.prefix="atac_v1_adult_brain_fresh_5k.1",
  path.to.snaptools="/home/ggj/anaconda3/envs/scATAC/bin/snaptools",
  path.to.macs="/home/ggj/anaconda3/envs/scATAC/bin/macs2",
  gsize="hs", 
  buffer.size=500, 
  num.cores=5,
  macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
  tmp.folder=tempdir()
);

# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 50)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("atac_col10.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/home/ggj/anaconda3/envs/scATAC/bin/snaptools",
    path.to.macs="/home/ggj/anaconda3/envs/scATAC/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
}, mc.cores=5);
# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));
peak.gr

# Step 12. Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
saveRDS(x.sp, file="atac_v1_adult_brain_fresh_5k.snap.rds");

system("snaptools snap-add-pmat \
       --snap-file ./human.snap \
       --peak-file peaks.combined.bed")

# Step 13. Add cell-by-peak matrix
# x.sp = readRDS("atac_v1_adult_brain_fresh_5k.snap.rds");
x.sp = addPmatToSnap(x.sp);
x.sp = makeBinary(x.sp, mat="pmat");
x.sp

# Step 14. Identify differentially accessible regions
DARs = findDAR(
  obj=x.sp,
  input.mat="pmat",
  cluster.pos=1,
  cluster.neg.method="knn",
  test.method="exactTest",
  bcv=0.1, #0.4 for human, 0.1 for mouse
  seed.use=10
);
DARs$FDR = p.adjust(DARs$PValue, method="BH");
idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
par(mfrow = c(1, 2));
plot(DARs$logCPM, DARs$logFC, 
     pch=19, cex=0.1, col="grey", 
     ylab="logFC", xlab="logCPM",
     main="Cluster 26"
);
points(DARs$logCPM[idy], 
       DARs$logFC[idy], 
       pch=19, 
       cex=0.5, 
       col="red"
);
abline(h = 0, lwd=1, lty=2);
covs = Matrix::rowSums(x.sp@pmat);
vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
vals.zscore = (vals - mean(vals)) / sd(vals);
plotFeatureSingle(
  obj=x.sp,
  feature.value=vals.zscore,
  method="tsne", 
  main="Cluster 26",
  point.size=0.1, 
  point.shape=19, 
  down.sample=5000,
  quantiles=c(0.01, 0.99)
);

idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    bcv=0.1,
    test.method="exactTest",
    seed.use=10
  );
  DARs$FDR = p.adjust(DARs$PValue, method="BH");
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue;
    PValues[DARs$logFC < 0] = 1;
    idy = order(PValues, decreasing=FALSE)[1:2000];
    rm(PValues); # free memory
  }
  idy
})
names(idy.ls) = levels(x.sp@cluster);
par(mfrow = c(3, 3));
for(cluster_i in levels(x.sp@cluster)){
  print(cluster_i)
  idy = idy.ls[[cluster_i]];
  vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
  vals.zscore = (vals - mean(vals)) / sd(vals);
  plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="tsne", 
    main=cluster_i,
    point.size=0.1, 
    point.shape=19, 
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  );
}

# Step 15. Motif analysis identifies master regulators
motifs = runHomer(
  x.sp[,idy.ls[["1"]],"pmat"], 
  mat = "pmat",
  path.to.homer = "/home/ggj/jiaqiLi/general_global_soft/homer/bin/findMotifsGenome.pl",
  result.dir = "./homer/C5",
  num.cores=5,
  genome = 'hg38',
  motif.length = 10,
  scan.size = 300,
  optimize.count = 2,
  background = 'automatic',
  local.background = FALSE,
  only.known = TRUE,
  only.denovo = FALSE,
  fdr.num = 5,
  cache = 100,
  overwrite = TRUE,
  keep.minimal = FALSE
)

# ---- run homer in console ------
system("rm -r ./homer");
result.dir<-"homer"
dir.create(result.dir)

levels(x.sp@cluster)
for (i in levels(x.sp@cluster)) {
  obj <- x.sp[,idy.ls[[i]],"pmat"]
  x <- as.data.frame(obj@peak)[,c(1:3)]
  x$seqnames<-as.character(x$seqnames)
  x$seq='none'
  for (j in 1:2000){  
    x[j,4]<-unlist(strsplit(x$seqnames,"'"))[2*j]
  }
  x$seqnames=x$seq
  x<-x[,-4]
  
  dir.create(paste0(result.dir, '/C', i))
  write.table(x, file=paste0(result.dir, '/C', i, "/target.bed"), 
              row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
}

# # run in terminal
# for i in {1..14};do
#   findMotifsGenome.pl ./C${i}/target.bed /media/ggj/Files/NvWA/PreprocGenome/database_genome/Smed_genome/SmedAsxl_genome_v1.1.nt \
#     ./C${i} -len 10, -size 300 -S 2 -p 5 -cache 100 -fdr 5 &>log.$i &
# done

##merge
motif <- matrix(nrow = 0, ncol = 5, byrow = TRUE)
for (i in levels(x.sp@cluster)) {
  fname <- paste0(result.dir, '/C', i, "/knownResults.txt")
  if(!file.exists(fname)){next}
  x = read.csv(fname, sep="\t", header=TRUE)
  x$Cluster = i
  x <- x[,c(10, 1:4)]
  motif <- rbind(motif,x)
}
motif <- motif[motif['P.value'] < 0.05, ]

WriteXLS::WriteXLS(motif, "Smed_motifs.xls",row.names = F)
write.csv(motif, "Smed_motifs.csv",row.names = F)
