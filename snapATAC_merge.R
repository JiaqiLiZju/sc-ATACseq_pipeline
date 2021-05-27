library(SnapATAC)
setwd('/media/ggj/Files/scATAC/CHATAC/Axolotl/')

#Step 0. loading data && Create snap object
file.list = c()
sample.list = c()
for(i in c(10:81)){
  file.list <- append(file.list, paste0("./COL",i,"/H.snap"))
  sample.list <- append(sample.list, paste0("snATAC_COL",i))
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
  sel = sort(metadata$UM, decreasing = T, index.return=T)$ix[1:500]
  x.sp.sel = x.sp[sel,]
  x.sp <- x.sp.sel
})
names(x.sp.list) = sample.list
x.sp.list

#Step 2. Add cell-by-bin matrix
x.sp.list = lapply(seq(x.sp.list), function(i){
  x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=50000);
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
rm(x.sp.list)
gc()
table(x.sp@sample)

# Step 4. Matrix binarization
x.sp = makeBinary(x.sp, mat="bmat");

# Step 5. Bin filtering
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 10)
);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.999);
bin.cutoff 
idy = which(bin.cov <= bin.cutoff & bin.cov > 0.1);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

#Step 6. Reduce dimensionality
row.covs = log10(Matrix::rowSums(x.sp@bmat)+1);
row.covs.dens = density(
  x = row.covs, 
  bw = 'nrd', adjust = 1
);
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(1);
x.sp
idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 5000, prob = sampling_prob));
x.landmark.sp = x.sp[idx.landmark.ds,];
x.query.sp = x.sp[-idx.landmark.ds,];
x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp,
  input.mat="bmat", 
  num.eigs=50
);
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp, 
  obj2=x.query.sp,
  input.mat="bmat"
);
x.landmark.sp@metaData$landmark = 1;
x.query.sp@metaData$landmark = 0;
x.sp = snapRbind(x.landmark.sp, x.query.sp);
x.sp = x.sp[order(x.sp@sample),];
x.sp
rm(x.landmark.sp, x.query.sp);

#Step 7. Determine significant components
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

#Step 8. Remove batch effect
library(harmony);
x.after.sp = runHarmony(
  obj=x.sp, 
  eigs.dim=1:15, 
  meta_data=x.sp@sample # sample index
);

#Step 9. Graph-based cluster
x.after.sp = runKNN(
  obj= x.after.sp,
  eigs.dim=1:15,
  k=15
);
x.after.sp = runCluster(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  path.to.snaptools=NULL,
  seed.use=10
);
x.after.sp@metaData$cluster = x.after.sp@cluster;

#Step 10. Visualization
x.after.sp = runViz(
  obj=x.after.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:15, 
  method="Rtsne",
  seed.use=10
);

plotViz(
  obj=x.after.sp,
  method="tsne", 
  main="After Harmony",
  point.color=x.sp@sample, 
  point.size=0.1, 
  text.add=FALSE,
  down.sample=10000,
  legend.add=TRUE
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
  down.sample=10000,
  legend.add=FALSE
)

save.image("./workspace.Rdata")
