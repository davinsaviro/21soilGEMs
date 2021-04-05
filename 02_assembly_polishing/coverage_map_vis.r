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

#####

cov_005 <- read.table("005_depth_short.txt")

colnames(cov_005) <- c("Genome", "Locus", "Depth")

cov_005$bin <- cut(cov_005$Locus, breaks = ceiling(nrow(cov_005) / 100), labels = FALSE)

cov_005_bin <- data.frame("Locus" = tapply(X = cov_005$Locus, INDEX = cov_005$bin, FUN = mean),
                         "Depth" = tapply(X = cov_005$Depth, INDEX = cov_005$bin, FUN = mean))

cov_005_bin$LocusKB <- cov_005_bin$Locus / 1000

cov_005_bin$SlidingAverage <- sapply(X = 1:nrow(cov_005_bin),
                                    FUN = movMeanCirc,
                                    depths = cov_005_bin$Depth,
                                    window = 500)

contig_borders_005 <- cumsum(width(readDNAStringSet("005_contig.fasta"))) / 1000

plot_005 <- ggplot(cov_005_bin, aes(x = LocusKB, y = Depth)) +
  ggtitle("005 Read Mapping") +
  geom_line() +
  labs(x = "Locus (kb)") +
  theme_bw(base_size = 16)

plot_005_MovAvg <- ggplot(cov_005_bin, aes(x = LocusKB, y = Depth)) +
  geom_line(col = "grey") +
  geom_line(aes(y = SlidingAverage), col = "blue") + # add an extra layer
  labs(title = "005 Read Mapping",
       x = "Locus (kb)") +
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

ggsave("005plot.png", plot = plot_005, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("005plot_MovAvg.png", plot = plot_005_MovAvg, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
ggsave("005plot_ConBor.png", plot = plot_005_ConBor, device = "png", scale = 3, width = 10, height = 6, units = c("cm"), dpi = 300, limitsize = TRUE)
