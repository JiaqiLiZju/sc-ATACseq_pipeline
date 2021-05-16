setwd("/media/ggj/Files/scATAC/CHATAC/HumanMouseMix_20210510/")

all <- list()
for(i in c(10:40)){
  data <- read.table(paste0("COL", i ,"/fragment.length.txt"), header = F)
  all[[i]] = data
}
all %>% reduce(rbind) -> data

# 设置插入片段长度的阈值，过滤掉太长的片段
length_cutoff <- 1000
fragment <- data$V1[data$V1 <= length_cutoff]
# 利用直方图统计频数分布，设置柱子个数
breaks_num <- 500
res <- hist(fragment, breaks = breaks_num, plot = FALSE)
# 添加坐标原点
plot(x = c(0, res$breaks),
     y = c(0, 0, log(res$counts)),
     type = "l", col = "red",
     xlab = "Fragment length(bp)",
     ylab = expression(Log ~ Normalized ~ read ~ density),
     main = "Merged Fragment sizes")

plot(x = c(0, res$breaks),
     y = c(0, 0, res$counts) / 1e2,
     type = "l", col = "red",
     xlab = "Fragment length(bp)",
     ylab = expression(Normalized ~ read ~ density),
     main = "Merged Fragment sizes")


summary(data$V1)
