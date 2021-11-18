library(ggplot2)
library(dplyr)
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
library(ggpointdensity) # 绘制密度散点图

metadata <- read.csv("COL10_Aligned.snap.metadata.csv", sep=',', row.names = 1)
head(metadata); dim(metadata)

tss <- read.csv("./COL10_Aligned.snap.bed.barcode_tssEnrich.txt", sep='\t', header=F, row.names = 1)
tss <- tss[rownames(metadata),]
head(tss); dim(tss)

metadata$tss <- tss
metadata$LogUQ <- log(metadata$UQ+1)
head(metadata); dim(metadata)

ggplot(data = metadata, aes(x=LogUQ, y=tss)) +
  geom_pointdensity() + #密度散点图（geom_pointdensity）
  scale_color_viridis() +
  # geom_smooth(method = lm) +  ##省略拟合曲线
  # stat_cor(method = "spearman") +
  xlab("UQ") + ylab("tss") + 
  theme(axis.title.x = element_text(size = 16, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 16, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size=12,  family="serif")) +
  theme(legend.position='none')  ##去除legend
