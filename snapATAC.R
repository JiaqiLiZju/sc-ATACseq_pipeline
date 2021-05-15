setwd("/media/ggj/Files/scATAC/CHATAC/HumanMouseMix_20210510/COL10/")

library(SnapATAC);

x.sp = createSnap(
  file="human.snap",
  sample="human",
  num.cores=10
);
barcodes = x.sp@barcode
metadata = x.sp@metaData
sel = sort(metadata$UM, decreasing = T, index.return=T)$ix[1:500]
metadata.sel = metadata[sel,]
head(metadata.sel)
summary(metadata.sel$UM)
summary(metadata.sel$UQ)
summary(metadata.sel$TN)

# write.csv(metadata.sel, file = "./mouse.summary.txt")

x.sp.sel = x.sp[sel,];
# x.sp.sel@metaData = barcodes.sel[x.sp@barcode,];
x.sp <- x.sp.sel
x.sp

barcodes = read.csv(
  "atac_v1_adult_brain_fresh_5k_singlecell.csv",
  head=TRUE
);
barcodes = barcodes[2:nrow(barcodes),];
promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
UMI = log(barcodes$passed_filters+1, 10);
data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
barcodes$promoter_ratio = promoter_ratio;

library(viridisLite);
library(ggplot2);
p1 = ggplot(
  data, 
  aes(x= UMI, y= promoter_ratio)) + 
  geom_point(size=0.1, col="grey") +
  theme_classic() +
  ggtitle("10X Fresh Adult Brain") +
  ylim(0, 1) + xlim(0, 6) +
  labs(x = "log10(UMI)", y="promoter ratio") 
p1 
barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6),];
rownames(barcodes.sel) = barcodes.sel$barcode;

x.sp.sel = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),];
x.sp.sel@metaData = barcodes.sel[x.sp@barcode,];
x.sp.sel

showBinSizes("./tmp/human.demo.snap");
x.sp = addBmatToSnap(x.sp, bin.size=10000, num.cores=1);

x.sp = makeBinary(x.sp, mat="bmat");
x.sp

# system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz");
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


bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 3)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.99);
bin.cutoff = 1
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
);

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

x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
);
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
);
x.sp@metaData$cluster = x.sp@cluster;

x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  seed.use=10
);
par(mfrow = c(2, 2));
plotViz(
  obj=x.sp,
  method="tsne", 
  main="CHATAC Human",
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
  main="10X Brain Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
); 
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters,
  method="tsne", 
  main="10X Brain FRiP",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
);
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$duplicate / x.sp@metaData$total,
  method="tsne", 
  main="10X Brain Duplicate",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99) # remove outliers
);

