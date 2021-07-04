setwd('/media/ggj/Files/scATAC/CHATAC/Axolotl/')

x.sp <- readRDS("./x.after.sp-batch2-BS5k-hvb-pmat.rds")

library(SnapATAC);
library(viridisLite);
library(ggplot2);
library(GenomicRanges);

# sample and test
idx <- x.after.sp@sample=="Axolotl-batch2-20210604-snATAC_COL88"

# barcodes <- x.after.sp@barcode
# Digest 13899 cells
# idx <- grep("TCCTACCAGT$|GCGTTGGAGC$|GATCTTACGC$|CTGATGGTCA$|CCGAGAATCC$|GCCGCAACGA$|TGAGTCTGGC$|TGCGGACCTA$|ACGGAGGCGG$|TAGATCTACT$|AATTAAGACT$|TTATTCATTC$|TTGACTTCAG$|AGAGCTATAA$|CTAAGAGAAG$|ACTCAATAGG$|CTTGCGCCGC$|GGTACTGCCT$|TAGAATTAAC$", barcodes)
# kidney 11785
# idx <- grep("GCCATTCTCC$|TGCCGGCAGA$|TTACCGAGGC$|ATCATATTAG$|TGGTCAGCCA$|ACTATGCAAT$|CGACGCGACT$|GATACGGAAC$|TTATCCGGAT$|TAGAGTAATA$|GCAGGTCCGT$|TCGGCCTTAC$|AGAACGTCTC$|CCAGTTCCAA$|ACTTAACCTT$|CAACCGCTAA$|GACCTTGATA$|TCTGATACCA$|GAAGATCGAG$", barcodes)
# Neuron 13906
# idx <- grep("AAGAAGCTAG$|TCCGGCCTCG$|AGAGAAGGTT$|CATACTCCGA$|GCTAACTTGC$|GGCTGAGCTC$|CCGATTCCTG$|ACCGCCAACC$|ATAAGGAGCA$|CGAACGCCGG$|GGTATGCTTG$|AACCTGCGTA$|GGCAGACGCC$|TAGCCGTCAT$|CTAGTAGTCT$|ACGCGAGATT$|GGTATCCGCC$|AACTAGGCGC$|TCGCTAAGCA$", barcodes)
# locomotor 37963
# idx <- grep("TATATACTAA$|ACTTGCTAGA$|AACCATTGGA$|TCGCGGTTGG$|CGTAGTTACC$|TCCAATCATC$|AATCGATAAT$|CCATTATCTA$|TCAACGTAAG$|TCTAATAGTA$|GATCGCTTCT$|CTAACTAGAT$|GCTGGAACTT$|AGGTTAGTTC$|CATTCGACGG$|CATTCAATCA$|CGGATTAGAA$|ATCGGCTATC$|ACGAAGTCAA$", barcodes)
# Lung
# idx <- grep("TTACCTCGAC$|GGAGGATAGC$|GGCTCTCTAT$|CGGTCAAGAA$|CGCTCCTAAC$|ATCCATGACT$|AACCTGGTCT$|GGTACCGGCA$|AAGCCAGTTA$|TCTTGCCGAC$|AAGACCGTTG$|AGGTTAGCAT$|TTCGCCTCCA$|AGAGCCAAGG$|AATACCATCC$|AGCTCTCCTC$|AGCTTATCCG$|CATCTCTGCA$|ACCTGGCCAA$|TAACTGGTTA$", barcodes)
idx[1:5]

x.after.sp.sel <- x.after.sp[idx,]
x.after.sp.sel

# Step 5. cutoff low-quality barcodes
barcode.cov = log10(Matrix::rowSums(x.after.sp.sel@pmat)+1);
hist(barcode.cov, xlab="log10(barcode.cov)", main="log10(barcode.cov)", 
     col="lightblue", xlim=c(0, 5));

barcode.cutoff.left = quantile(barcode.cov, 0.05); # low-quality 
barcode.cutoff.right = quantile(barcode.cov, 0.99); # doublets
barcode.cutoff.left; barcode.cutoff.right

idx = which((barcode.cov >= barcode.cutoff.left) & (barcode.cov <= barcode.cutoff.right));
x.after.sp.sel = x.after.sp.sel[idx,];
x.after.sp.sel

# Step 4. Peak filtering
# system("wget https://www.axolotl-omics.org/dl/RM_all_repeats.bed.gz");
# system("gunzip -c RM_all_repeats.bed.gz |sed -E 's/[[:space:]]+/\t/g' |sed 's/^\t//' |cut -f 5,6,7 >RM_all_repeats.bed");
library(GenomicRanges);
black_list = read.csv("RM_all_repeats.bed", sep='\t', header=F);
black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.after.sp.sel@feature, black_list.gr));
if(length(idy) > 0){
  x.after.sp.sel = x.after.sp.sel[,-idy, mat="bmat"]};
x.after.sp.sel

# chr.exclude = seqlevels(x.after.sp.sel@feature)[grep("random|chrM", seqlevels(x.after.sp.sel@feature))];
# idy = grep(paste(chr.exclude, collapse="|"), x.after.sp.sel@feature);
# if(length(idy) > 0){x.after.sp.sel = x.after.sp.sel[,-idy, mat="pmat"]};
# x.after.sp.sel

# Step 5. high variable features
peak.covs = log10(Matrix::colSums(x.after.sp.sel@pmat)+1);
hist(peak.covs, xlab="log10(peak.covs)", main="log10(peak.covs)", 
     col="lightblue", xlim=c(0, 4));
peak.cutoff.left = 2.5 # quantile(peak.cov, 0.6); # low-quality peaks
peak.cutoff.right = 5 # quantile(peak.cov, 0.999); # housekeeping peaks
peak.cutoff.left; peak.cutoff.right 

# peak.covs.zscore = (peak.covs - mean(peak.covs)) / sd(peak.covs);
# idy <- which(peak.covs.zscore > 2)
# length(idy)

idy = which(peak.covs <= peak.cutoff.right & peak.covs > peak.cutoff.left);
x.after.sp.sel = x.after.sp.sel[, idy, mat="pmat"];
x.after.sp.sel


# as_matrix <- function(mat){
#   
#   tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
#   
#   row_pos <- mat@i+1
#   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
#   val <- mat@x
#   
#   for (i in seq_along(val)){
#     tmp[row_pos[i],col_pos[i]] <- val[i]
#   }
#   
#   row.names(tmp) <- mat@Dimnames[[1]]
#   colnames(tmp) <- mat@Dimnames[[2]]
#   return(tmp)
# }
# 
# get_variable_genes<-function(data) {
#   dat.sub=data
#   med.dat=apply(dat.sub,1,median)
#   var.dat=apply(dat.sub,1,var)
#   quant.med=unique(quantile(med.dat,prob=seq(0,1,length=11),type=5))
#   
#   var.genes1=vector("list")
#   genes.list=vector("list")
#   
#   genes.list=vector("list",length=length(quant.med))
#   for(i in 1:length(quant.med)){
#     if(i==1){
#       filt1=med.dat<=quant.med[i]
#       var.temp=var.dat[filt1]
#       quant.var=quantile(var.temp,na.rm=T)
#       filt2=var.temp > quant.var[4]###### total is 4;TF is3
#       genes.list[[i]]=names(var.temp)[filt2]
#     }
#     else {
#       filt1=med.dat<=quant.med[i]&med.dat>quant.med[i-1]
#       var.temp=var.dat[filt1]
#       quant.var=quantile(var.temp,na.rm=T)
#       filt2=var.temp > quant.var[4]###########
#       genes.list[[i]]=names(var.temp)[filt2]
#     }
#   }
#   temp=length(genes.list)
#   var.genes=unlist(genes.list[1:temp-1])
# 
#   return(var.genes)
# }
# 
# # 
# data <- as.data.frame(as_matrix(x.after.sp.sel@pmat))
# colnames(data) <- c(1:dim(data)[2])
# data <- t(data)
# data[1:5,1:5]
# 
# hvp <- get_variable_genes(data)
# x.after.sp.sel = x.after.sp.sel[, hvp, mat="pmat"];
# x.after.sp.sel

# cutoff low-quality peaks
# idy = which(Matrix::colSums(x.after.sp.sel@pmat) > 1000);
# x.after.sp.sel = x.after.sp.sel[, idy, mat="pmat"];
# x.after.sp.sel

# peak.cov = log10(Matrix::colSums(x.after.sp.sel@pmat)+1);
# hist(
#   peak.cov,
#   xlab="log10(peak cov)",
#   main="log10(peak Cov)",
#   col="lightblue",
#   xlim=c(0, 4)
# );
# # peak.cutoff = quantile(peak.cov[peak.cov > 0], 0.95);
# # peak.cutoff
# # 
# # idy = which(peak.cov <= peak.cutoff & peak.cov > 1);
# # 因为蝾螈太稀疏了，要卡掉低值
# # x.after.sp.sel = x.after.sp.sel[, idy, mat="pmat"];
# # x.after.sp.sel

#Step 6. Reduce dimensionality
row.covs = log10(Matrix::rowSums(x.after.sp.sel@pmat)+1);
row.covs.dens = density(
  x = row.covs,
  bw = 'nrd', adjust = 1
);
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(1);
x.after.sp.sel

# sample landmarks
idx.landmark.ds = sort(sample(x = seq(nrow(x.after.sp.sel)), size = 10000, prob = sampling_prob));
x.landmark.sp = x.after.sp.sel[idx.landmark.ds,];
x.query.sp = x.after.sp.sel[-idx.landmark.ds,];
x.landmark.sp = runDiffusionMaps(
  obj= x.landmark.sp[Matrix::rowSums(x.landmark.sp@pmat)!=0,],
  input.mat="pmat",
  num.eigs=50
);
x.query.sp = runDiffusionMapsExtension(
  obj1=x.landmark.sp,
  obj2=x.query.sp[Matrix::rowSums(x.query.sp@pmat)!=0,],
  input.mat="pmat"
);
x.landmark.sp@metaData$landmark = 1;
x.query.sp@metaData$landmark = 0;
x.after.sp.sel = snapRbind(x.landmark.sp, x.query.sp);
x.after.sp.sel = x.after.sp.sel[order(x.after.sp.sel@sample),];
x.after.sp.sel
rm(x.landmark.sp, x.query.sp);

# Step 5. Dimensionality reduction
# x.after.sp.sel = runDiffusionMaps(
#   obj=x.after.sp.sel[Matrix::rowSums(x.after.sp.sel@pmat)!=0,],
#   input.mat="pmat",
#   num.eigs=50
# );

# Step 6. Determine significant components
plotDimReductPW(
  obj=x.after.sp.sel, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=2000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
);

# Step 7. Graph-based clustering
x.after.sp.sel = runKNN(
  obj=x.after.sp.sel,
  eigs.dims=1:4,
  k=50
);
x.after.sp.sel=runCluster(
  obj=x.after.sp.sel,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  resolution=.5,
  seed.use=10
);
x.after.sp.sel@metaData$cluster = x.after.sp.sel@cluster;

# Step 8. Visualization
x.after.sp.sel = runViz(
  obj=x.after.sp.sel, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:4, 
  method="umap",
  seed.use=10
);
# par(mfrow = c(2, 2));

plotViz(
  obj=x.after.sp.sel,
  method="umap", 
  main="CH-ATAC(Axolotl)",
  point.color=x.after.sp.sel@cluster, 
  point.size=1, 
  # point.shape=10, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=NULL,
  legend.add=FALSE
);

plotFeatureSingle(
  obj=x.after.sp.sel,
  feature.value=log(x.after.sp.sel@metaData$TN+1,10),
  method="umap", 
  main="Read Depth",
  point.size=0.2, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
); 

# Step 8. Visualization
x.after.sp.sel = runViz(
  obj=x.after.sp.sel, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:4, 
  method="Rtsne",
  seed.use=10
);
# par(mfrow = c(2, 2));

plotViz(
  obj=x.after.sp.sel,
  method="tsne", 
  main="CH-ATAC(Axolotl)",
  point.color=x.after.sp.sel@cluster, 
  point.size=1, 
  # point.shape=10, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=NULL,
  legend.add=FALSE
);

saveRDS(x.after.sp.sel, file="./Tissue/Lung-Axolotl-batch2-BS40k-hvb-peak.rds")

# Step 14. Identify differentially accessible regions
DARs = findDAR(
  obj=x.after.sp.sel,
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
     main="Cluster"
);
points(DARs$logCPM[idy], 
       DARs$logFC[idy], 
       pch=19, 
       cex=0.5, 
       col="red"
);
abline(h = 0, lwd=1, lty=2);
covs = Matrix::rowSums(x.after.sp.sel@pmat);
vals = Matrix::rowSums(x.after.sp.sel@pmat[,idy]) / covs;
vals.zscore = (vals - mean(vals)) / sd(vals);
plotFeatureSingle(
  obj=x.after.sp.sel,
  feature.value=vals.zscore,
  method="tsne", 
  main="Cluster",
  point.size=0.1, 
  point.shape=19, 
  down.sample=5000,
  quantiles=c(0.01, 0.99),
  na.rm = T
);

idy.ls = lapply(levels(x.after.sp.sel@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.after.sp.sel,
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
names(idy.ls) = levels(x.after.sp.sel@cluster);
par(mfrow = c(3, 3));
for(cluster_i in levels(x.after.sp.sel@cluster)){
  print(cluster_i)
  idy = idy.ls[[cluster_i]];
  vals = Matrix::rowSums(x.after.sp.sel@pmat[,idy]) / covs;
  vals.zscore = (vals - mean(vals)) / sd(vals);
  plotFeatureSingle(
    obj=x.after.sp.sel,
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
  x.after.sp.sel[,idy.ls[["1"]],"pmat"], 
  mat = "pmat",
  path.to.homer = "/home/ggj/jiaqiLi/general_global_soft/homer/bin/findMotifsGenome.pl",
  result.dir = "./homer/C5",
  num.cores=5,
  genome = 'hb38',
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

