# options(future.globals.maxSize=200*1024^3)
# getOption("future.globals.maxSize")

library(SnapATAC)
setwd('/media/ggj/Files/scATAC/CHATAC/Axolotl/')
# setwd("J:/Axolotl/Axolotl")

#Step 0. loading data && Create snap object
file.list = c()
sample.list = c()
# for(i in c(10:81)){
#   print(paste0("Axolotl-batch1-20210524-snATAC_COL",i))
#   file.list <- append(file.list, paste0("../Axolotl-batch1-20210524/COL",i,"/H.snap"))
#   sample.list <- append(sample.list, paste0("Axolotl-batch1-20210524-snATAC_COL",i))
# }
# 
# for(i in c(11:105)){
#   print(paste0("Axolotl-batch1-20210602-snATAC_COL",i))
#   file.list <- append(file.list, paste0("../Axolotl-batch1-20210602/COL",i,"/bowtie_out/H.snap"))
#   sample.list <- append(sample.list, paste0("Axolotl-batch1-20210602-snATAC_COL",i))
# }

for(i in c(10:105)){
  print(paste0("Axolotl-batch2-20210604-snATAC_COL",i))
  file.list <- append(file.list, paste0("../Axolotl-batch2-20210604/COL",i,"/bowtie_out/H.snap"))
  sample.list <- append(sample.list, paste0("Axolotl-batch2-20210604-snATAC_COL",i))
}

# file.list = c("./COL10/H.snap", "./COL11/H.snap",'./COL12/H.snap','./COL13/H.snap','./COL14/H.snap','./COL15/H.snap','./COL16/H.snap','./COL17/H.snap','./COL18/H.snap','./COL19/H.snap','./COL20/H.snap');
# sample.list = c("snATAC_COL10", "snATAC_COL11",'snATAC_COL12','snATAC_COL13','snATAC_COL14','snATAC_COL15','snATAC_COL16','snATAC_COL17','snATAC_COL18','snATAC_COL19','snATAC_COL20');
x.sp.ls = lapply(seq(file.list), function(i){
  x.sp = createSnap(file=file.list[i], sample=sample.list[i]);
  x.sp
})
x.sp.ls
names(x.sp.ls) = sample.list;
x.sp.ls

#Step 1. select barcode 
x.sp.list = lapply(seq(x.sp.ls), function(i){
  x.sp = x.sp.ls[[i]]
  barcodes = x.sp@barcode
  metadata = x.sp@metaData
  sel = sort(metadata$UM, decreasing = T, index.return=T)$ix[1:1000]
  x.sp.sel = x.sp[sel,]
  x.sp <- x.sp.sel
})
names(x.sp.list) = sample.list
x.sp.list

#Step 2. Add cell-by-bin matrix
showBinSizes("../Axolotl-batch2-20210604/COL10/bowtie_out/H.snap");
x.sp.list = lapply(seq(x.sp.list), function(i){
  x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000);
  x.sp
})
x.sp.list

#Step 3. Combine snap objects
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name));
x.sp.list <- lapply(x.sp.list, function(x.sp){
  idy = match(bin.shared, x.sp@feature$name);
  x.sp[,idy, mat="bmat"];
})
x.sp = Reduce(snapRbind, x.sp.list);
x.sp

# rm(x.sp.list)
# gc()
table(x.sp@sample)

# Step 4. Matrix binarization
x.sp = makeBinary(x.sp, mat="bmat");


x.sp <- readRDS("./x.fater.sp-batch2-BS5k.rds")

# Step 5. cutoff low-quality barcodes
barcode.cov = log10(Matrix::rowSums(x.sp@bmat)+1);
hist(barcode.cov, xlab="log10(barcode.cov)", main="log10(barcode.cov)", 
  col="lightblue", xlim=c(0, 8));

barcode.cutoff.left = quantile(barcode.cov, 0.1); # low-quality 
barcode.cutoff.right = quantile(barcode.cov, 0.98); # doublets
barcode.cutoff.left; barcode.cutoff.right

idx = which((barcode.cov >= barcode.cutoff.left) & (barcode.cov <= barcode.cutoff.right));
x.sp = x.sp[idx,];
x.sp

# Step 5. Bin filtering
# system("wget https://www.axolotl-omics.org/dl/RM_all_repeats.bed.gz");
# system("gunzip -c RM_all_repeats.bed.gz |sed -E 's/[[:space:]]+/\t/g' |sed 's/^\t//' |cut -f 5,6,7 >RM_all_repeats.bed");
library(GenomicRanges);
black_list = read.csv("RM_all_repeats.bed", sep='\t', header=F);
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# high variable bins
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(bin.cov, xlab="log10(bin cov)", main="log10(Bin Cov)", 
  col="lightblue", xlim=c(0, 5));
bin.cutoff.left = 0 # quantile(bin.cov, 0.6); # low-quality bins
bin.cutoff.right = 5 # quantile(bin.cov, 0.999); # housekeeping bins
bin.cutoff.left; bin.cutoff.right 

idy = which(bin.cov <= bin.cutoff.right & bin.cov > bin.cutoff.left);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# covs = Matrix::colSums(x.sp@bmat);
# covs.zscore = (covs - mean(covs)) / sd(covs);
# idy <- which(covs.zscore > 1.65)
# 
# x.sp = x.sp[, idy, mat="bmat"];
# x.sp

#Step 6. Reduce dimensionality
row.covs = log10(Matrix::rowSums(x.sp@bmat)+1);
hist(row.covs, xlab="log10(row.covs)", main="log10(row.covs)", 
     col="lightblue", xlim=c(0, 5));
row.covs.dens = density(
  x = row.covs, 
  bw = 'nrd', adjust = 1
);
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(1);
x.sp

# sample landmarks
idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 20000, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];
x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp[Matrix::rowSums(x.landmark.sp@bmat)!=0,],
  input.mat="bmat", 
  num.eigs=100
);
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp[Matrix::rowSums(x.query.sp@bmat)!=0,],
  input.mat="bmat"
);
x.landmark.sp@metaData$landmark = 1;
x.query.sp@metaData$landmark = 0;
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@sample),];
x.sp

rm(x.landmark.sp, x.query.sp);

# Step 5. Dimensionality reduction
x.after.sp = runDiffusionMaps(
  obj=x.after.sp,
  input.mat="bmat", 
  num.eigs=50
);

#Step 7. Determine significant components
plotDimReductPW(
  obj=x.after.sp, 
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
x.after.sp

#Step 8. Remove batch effect
# library(harmony);
# x.after.sp = runHarmony(
#   obj=x.sp,
#   eigs.dim=1:15,
#   meta_data=x.sp@sample # sample index
# );

x.after.sp <- x.sp
# save(x.after.sp,file = './x.after.sp.Rdata')

#Step 9. Graph-based cluster
x.after.sp = runKNN(
  obj= x.after.sp,
  eigs.dim=1:28,
  k=20
);
library(leiden)
x.after.sp = runCluster(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  resolution=.5,
  seed.use=10
);

x.after.sp@metaData$cluster = x.after.sp@cluster;

#Step 10. Visualization
x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:28, 
  method="umap",
  num.cores=10,
  seed.use=10
);

plotViz(
  obj=x.after.sp,
  method="umap", 
  main="After Harmony",
  point.color=x.after.sp@sample, 
  point.size=0.1, 
  text.add=FALSE,
  down.sample=NULL,
  legend.add=TRUE
);

plotViz(
  obj=x.after.sp,
  method="umap", 
  main="Cluster",
  point.color=x.after.sp@cluster, 
  point.size=0.1, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=NULL,
  legend.add=FALSE
)

#Step 10. Visualization
x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  num.cores=10,
  seed.use=10
);

plotViz(
  obj=x.after.sp,
  method="tsne", 
  main="Cluster",
  point.color=x.after.sp@cluster, 
  point.size=0.1, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=NULL,
  legend.add=FALSE
)


# Step 11. Identify peaks
# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.after.sp.sel@cluster))[which(table(x.after.sp.sel@cluster) > 50)];
# clusters.sel = c(1, 3)
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  runMACS(
    obj=x.after.sp.sel[which(x.after.sp.sel@cluster==clusters.sel[i]),], 
    output.prefix=paste0("./Peaks_Cluster/ATAC_Axolotl", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/home/ggj/anaconda3/envs/scATAC/bin/snaptools",
    path.to.macs="/home/ggj/anaconda3/envs/scATAC/bin/macs2",
    gsize="32000000000", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    # macs.options="--nomodel --qval 0.01 -B",
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir(),
    keep.minimal=T
  );
}, mc.cores=2);

# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls ./Peaks_Cluster/ | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(paste0("./Peaks_Cluster/", x))
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));
peak.gr

# Step 12. Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "./Peaks_Cluster/sorted_peaks.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
saveRDS(x.after.sp.sel, file="./Peaks_Cluster/peaks.snap.rds")

# system("/home/ggj/anaconda3/envs/scATAC/bin/snaptools snap-del \
#        --snap-file ../../Axolotl-batch2-20210604/COL88/bowtie_out/H.snap \
#        --session-name PM")

# system("/home/ggj/anaconda3/envs/scATAC/bin/snaptools snap-add-pmat \
#        --snap-file ../Axolotl-batch2-20210604/COL88/bowtie_out/H.snap \
#        --peak-file peaks.COL88.bed")

# Step 13. Add cell-by-peak matrix
x.after.sp = addPmatToSnap(x.after.sp);
x.after.sp = makeBinary(x.after.sp, mat="pmat");
x.after.sp

saveRDS(x.after.sp, file="./x.after.sp-batch2-BS5k-hvb-pmat.rds")
# x.after.sp <- readRDS("./x.fater.sp-batch2-BS5k.rds")
x.after.sp

# save.image("./workspace-batch2-filterbin.Rdata")
load("./workspace-batch2-hvb-pmat.Rdata")