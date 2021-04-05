library(ggplot2)
library(Biostrings)

movMeanCirc <- function(depths, window = 500, focus = 1){
  linear_end <- length(depths)
  index_left <- focus - window
  index_right <- focus + window
  index_before <- ifelse(index_left >= 1,
                         yes = index_left,
                         no = linear_end + index_left)
  index_after <- ifelse(index_right <= linear_end,
                        yes = index_right,
                        no = index_right - linear_end)
  res <- ifelse(index_before <= index_after,
                yes = mean(depths[index_before:index_after]),
                no = mean(depths[-((index_after + 1):(index_before - 1))]))
  return(res)
}

cov_002 <- read.table("002_depth_short.txt")
cov_004 <- read.table("004_depth_short.txt")
cov_005 <- read.table("005_depth_short.txt")
cov_009 <- read.table("009_depth_short.txt")
cov_012 <- read.table("012_depth_short.txt")
cov_015 <- read.table("015_depth_short.txt")
cov_017 <- read.table("017_depth_short.txt")
cov_018 <- read.table("018_depth_short.txt")
cov_019 <- read.table("019_depth_short.txt")

colnames(cov_002) <- c("Genome", "Locus", "Depth")
colnames(cov_004) <- c("Genome", "Locus", "Depth")
colnames(cov_005) <- c("Genome", "Locus", "Depth")
colnames(cov_009) <- c("Genome", "Locus", "Depth")
colnames(cov_012) <- c("Genome", "Locus", "Depth")
colnames(cov_015) <- c("Genome", "Locus", "Depth")
colnames(cov_017) <- c("Genome", "Locus", "Depth")
colnames(cov_018) <- c("Genome", "Locus", "Depth")
colnames(cov_019) <- c("Genome", "Locus", "Depth")

cov_002$bin <- cut(cov_002$Locus, breaks = ceiling(nrow(cov_002) / 100), labels = FALSE)
cov_004$bin <- cut(cov_004$Locus, breaks = ceiling(nrow(cov_004) / 100), labels = FALSE)
cov_005$bin <- cut(cov_005$Locus, breaks = ceiling(nrow(cov_005) / 100), labels = FALSE)
cov_009$bin <- cut(cov_009$Locus, breaks = ceiling(nrow(cov_009) / 100), labels = FALSE)
cov_012$bin <- cut(cov_012$Locus, breaks = ceiling(nrow(cov_012) / 100), labels = FALSE)
cov_015$bin <- cut(cov_015$Locus, breaks = ceiling(nrow(cov_015) / 100), labels = FALSE)
cov_017$bin <- cut(cov_017$Locus, breaks = ceiling(nrow(cov_017) / 100), labels = FALSE)
cov_018$bin <- cut(cov_018$Locus, breaks = ceiling(nrow(cov_018) / 100), labels = FALSE)
cov_019$bin <- cut(cov_019$Locus, breaks = ceiling(nrow(cov_019) / 100), labels = FALSE)

cov_002_bin <- data.frame("Locus" = tapply(X = cov_002$Locus, INDEX = cov_002$bin, FUN = mean),
                         "Depth" = tapply(X = cov_002$Depth, INDEX = cov_002$bin, FUN = mean))
cov_004_bin <- data.frame("Locus" = tapply(X = cov_004$Locus, INDEX = cov_004$bin, FUN = mean),
                         "Depth" = tapply(X = cov_004$Depth, INDEX = cov_004$bin, FUN = mean))
cov_005_bin <- data.frame("Locus" = tapply(X = cov_005$Locus, INDEX = cov_005$bin, FUN = mean),
                         "Depth" = tapply(X = cov_005$Depth, INDEX = cov_005$bin, FUN = mean))
cov_009_bin <- data.frame("Locus" = tapply(X = cov_009$Locus, INDEX = cov_009$bin, FUN = mean),
                         "Depth" = tapply(X = cov_009$Depth, INDEX = cov_009$bin, FUN = mean))
cov_012_bin <- data.frame("Locus" = tapply(X = cov_012$Locus, INDEX = cov_012$bin, FUN = mean),
                         "Depth" = tapply(X = cov_012$Depth, INDEX = cov_012$bin, FUN = mean))
cov_015_bin <- data.frame("Locus" = tapply(X = cov_015$Locus, INDEX = cov_015$bin, FUN = mean),
                         "Depth" = tapply(X = cov_015$Depth, INDEX = cov_015$bin, FUN = mean))
cov_017_bin <- data.frame("Locus" = tapply(X = cov_017$Locus, INDEX = cov_017$bin, FUN = mean),
                         "Depth" = tapply(X = cov_017$Depth, INDEX = cov_017$bin, FUN = mean))
cov_018_bin <- data.frame("Locus" = tapply(X = cov_018$Locus, INDEX = cov_018$bin, FUN = mean),
                         "Depth" = tapply(X = cov_018$Depth, INDEX = cov_018$bin, FUN = mean))
cov_019_bin <- data.frame("Locus" = tapply(X = cov_019$Locus, INDEX = cov_019$bin, FUN = mean),
                         "Depth" = tapply(X = cov_019$Depth, INDEX = cov_019$bin, FUN = mean))

cov_002_bin$LocusKB <- cov_002_bin$Locus / 1000
cov_004_bin$LocusKB <- cov_004_bin$Locus / 1000
cov_005_bin$LocusKB <- cov_005_bin$Locus / 1000
cov_009_bin$LocusKB <- cov_009_bin$Locus / 1000
cov_012_bin$LocusKB <- cov_012_bin$Locus / 1000
cov_015_bin$LocusKB <- cov_015_bin$Locus / 1000
cov_017_bin$LocusKB <- cov_017_bin$Locus / 1000
cov_018_bin$LocusKB <- cov_018_bin$Locus / 1000
cov_019_bin$LocusKB <- cov_019_bin$Locus / 1000

cov_002_bin$SlidingAverage <- sapply(X = 1:nrow(cov_002_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_002_bin$Depth,
                                    window = 500)
cov_004_bin$SlidingAverage <- sapply(X = 1:nrow(cov_004_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_004_bin$Depth,
                                    window = 500)
cov_005_bin$SlidingAverage <- sapply(X = 1:nrow(cov_005_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_005_bin$Depth,
                                    window = 500)
cov_009_bin$SlidingAverage <- sapply(X = 1:nrow(cov_009_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_009_bin$Depth,
                                    window = 500)
cov_012_bin$SlidingAverage <- sapply(X = 1:nrow(cov_012_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_012_bin$Depth,
                                    window = 500)
cov_015_bin$SlidingAverage <- sapply(X = 1:nrow(cov_015_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_015_bin$Depth,
                                    window = 500)
cov_017_bin$SlidingAverage <- sapply(X = 1:nrow(cov_017_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_017_bin$Depth,
                                    window = 500)
cov_018_bin$SlidingAverage <- sapply(X = 1:nrow(cov_018_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_018_bin$Depth,
                                    window = 500)
cov_019_bin$SlidingAverage <- sapply(X = 1:nrow(cov_019_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_019_bin$Depth,
                                    window = 500)

contig_borders_002 <- cumsum(width(readDNAStringSet("Final_002_startaligned.fasta"))) / 1000
contig_borders_004 <- cumsum(width(readDNAStringSet("Final_004_startaligned.fasta"))) / 1000
contig_borders_005 <- cumsum(width(readDNAStringSet("Final_005_startaligned.fasta"))) / 1000
contig_borders_009 <- cumsum(width(readDNAStringSet("Final_009_startaligned.fasta"))) / 1000
contig_borders_012 <- cumsum(width(readDNAStringSet("Final_012_startaligned.fasta"))) / 1000
contig_borders_015 <- cumsum(width(readDNAStringSet("Final_015_startaligned.fasta"))) / 1000
contig_borders_017 <- cumsum(width(readDNAStringSet("Final_017_startaligned.fasta"))) / 1000
contig_borders_018 <- cumsum(width(readDNAStringSet("Final_018_startaligned.fasta"))) / 1000
contig_borders_019 <- cumsum(width(readDNAStringSet("Final_019_startaligned.fasta"))) / 1000

plot_002 <- ggplot(cov_002_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("002 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_004 <- ggplot(cov_004_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("004 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_005 <- ggplot(cov_005_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("005 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_009 <- ggplot(cov_009_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("009 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_012 <- ggplot(cov_012_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("012 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_015 <- ggplot(cov_015_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("015 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_017 <- ggplot(cov_017_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("017 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_018 <- ggplot(cov_018_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("018 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_019 <- ggplot(cov_019_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("019 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)



plot_002_MovAvg <- ggplot(cov_002_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "002 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_004_MovAvg <- ggplot(cov_004_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "004 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_005_MovAvg <- ggplot(cov_005_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "005 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_009_MovAvg <- ggplot(cov_009_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "009 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_012_MovAvg <- ggplot(cov_012_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "012 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_015_MovAvg <- ggplot(cov_015_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "015 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_017_MovAvg <- ggplot(cov_017_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "017 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_018_MovAvg <- ggplot(cov_018_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "018 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)
plot_019_MovAvg <- ggplot(cov_019_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "019 Read Mapping",
       x = "Locus (kb)") +
  theme_bw(base_size = 16)

plot_002_ConBor <- ggplot(cov_002_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("002 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_002, col = "red", lty = 2, size = 1) +
  labs(title = "002 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_004_ConBor <- ggplot(cov_004_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("004 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_004, col = "red", lty = 2, size = 1) +
  labs(title = "004 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_005_ConBor <- ggplot(cov_005_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("005 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_005, col = "red", lty = 2, size = 1) +
  labs(title = "005 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_009_ConBor <- ggplot(cov_009_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("009 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_009, col = "red", lty = 2, size = 1) +
  labs(title = "009 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_012_ConBor <- ggplot(cov_012_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("012 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_012, col = "red", lty = 2, size = 1) +
  labs(title = "012 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_015_ConBor <- ggplot(cov_015_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("015 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_015, col = "red", lty = 2, size = 1) +
  labs(title = "015 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_017_ConBor <- ggplot(cov_017_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("017 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_017, col = "red", lty = 2, size = 1) +
  labs(title = "017 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_018_ConBor <- ggplot(cov_018_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("018 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_018, col = "red", lty = 2, size = 1) +
  labs(title = "018 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)
plot_019_ConBor <- ggplot(cov_019_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("019 Read Mapping") +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue", size = 2) +
  geom_vline(xintercept = contig_borders_019, col = "red", lty = 2, size = 1) +
  labs(title = "019 Read Mapping",
       x = "Locus (kb)",
       caption = "Grey: mean coverage every 100 bp; \n Blue: moving average of 100bp-means with total sliding window size 1001; \n Red: contig borders.") +
  theme_bw(base_size = 16)

ggsave("002plot.png", plot = plot_002, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("002plot_MovAvg.png", plot = plot_002_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("002plot_ConBor.png", plot = plot_002_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("004plot.png", plot = plot_004, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("004plot_MovAvg.png", plot = plot_004_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("004plot_ConBor.png", plot = plot_004_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("005plot.png", plot = plot_005, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("005plot_MovAvg.png", plot = plot_005_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("005plot_ConBor.png", plot = plot_005_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("009plot.png", plot = plot_009, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("009plot_MovAvg.png", plot = plot_009_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("009plot_ConBor.png", plot = plot_009_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("012plot.png", plot = plot_012, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("012plot_MovAvg.png", plot = plot_012_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("012plot_ConBor.png", plot = plot_012_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("015plot.png", plot = plot_015, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("015plot_MovAvg.png", plot = plot_015_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("015plot_ConBor.png", plot = plot_015_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("017plot.png", plot = plot_017, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("017plot_MovAvg.png", plot = plot_017_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("017plot_ConBor.png", plot = plot_017_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("018plot.png", plot = plot_018, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("018plot_MovAvg.png", plot = plot_018_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("018plot_ConBor.png", plot = plot_018_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("019plot.png", plot = plot_019, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("019plot_MovAvg.png", plot = plot_019_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("019plot_ConBor.png", plot = plot_019_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
