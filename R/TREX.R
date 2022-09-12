library(FNN)
library(fancycut)
library(ggplot2)
library(RColorBrewer)
library(dbscan)

TREX <- function(embedding.data,                   # dataframe, 2 columns of low dim embedding + file_ID for which dataset a cell belongs to  
                 kvalue = 60,                      # number, for KNN search 
                 bins = 5                          # number of bins to use for TREX 
                 ) {
  
  binned.data = TREX_bin(embedding.data, kvalue, bins)
  
  TREX_plot(binned.data)

  TREX_stats(binned.data)
  
  return(binned.data)

}

TREX_bin <- function(embedding.data, 
                     kvalue = 60,
                     bins = 5) {
  
  # KNN search per cell 
  neighbor.index = knnx.index(embedding.data[, 1:2], embedding.data[, 1:2], k = kvalue)
  neighbor.index[neighbor.index <= nrow(embedding.data)/2] <- 0
  neighbor.index[neighbor.index > nrow(embedding.data)/2] <- 1
  
  # calculate percent change in each KNN region
  percent.change = (rowSums(neighbor.index) / kvalue * 100)
  
  # binning 
  binned.data <- data.frame(x = embedding.data[, 1], y = embedding.data[, 2], percent.change = round(percent.change))
  binned.data$cuts <- wafflecut(binned.data$percent.change, c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]'))
  
  return(binned.data)
}

# create T-REX plot
TREX_plot <- function(binned.data) {

  trex.plot <- binned.data[order(binned.data$cuts), ]
  range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  trex.colors <- c("navyblue","lightskyblue","lightgray","lightcoral","darkred")
  names(trex.colors) <- c("[0,5]","(5,15]","(15,85)","[85,95)","[95,100]")
  
  ggplot(trex.plot) + 
    geom_point(aes(x = x, y = y, colour = cuts), cex = 1) + 
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(
      labels = c("\u2265 95% from 1st dataset",
                 "85-95% from 1st dataset",
                 "from 1st and 2nd dataset",
                 "85-95% from 2nd dataset",
                 "\u2265 95% from 2nd dataset"),
      values = trex.colors) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x = "UMAP1", y = "UMAP2",
         # title = paste(sample_type,"_",sample_id,"_",time_comparison," - Percent Change",sep = ""),
         caption = "Data from Brodin lab in Rodriguez et al., Cell Rep Med. 2020") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank())
  
  # save & export plot as PNG
  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_TREX_plot.png"), 
    width = 8, 
    height = 8
  )
}

# calculate degree of change, and direction of change for pairwise comparison
TREX_stats <- function(binned.data) {

  sample.table <- data.frame(total_cells = nrow(binned.data))
  sums = c(sum(binned.data$cuts == '(15,85)'),
           sum(binned.data$cuts == '(5,15]'),
           sum(binned.data$cuts == '[0,5]'),
           sum(binned.data$cuts == '[85,95)'),
           sum(binned.data$cuts == '[95,100]'))
  
  sample.table[, c(2:6)] <- sums
  colnames(sample.table)[2:6] <- c("(15,85)","(5,15]","[0,5]","[85,95)","[95,100]")
  
  percent = 100*(sums/nrow(binned.data))
  sample.table$degree_of_change = (sum(percent[3] + percent[5]))
  sample.table$direction_of_change = (sums[5] - sums[3]) / (sums[5] + sums[3])
  
  # export stats as CSV
  write.csv(sample.table, paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"), "_TREX_stats.csv"))
  
  return(sample.table) 
}

# use DBSCAN to cluster on regions of great change (5th and 95 percentiles of change)
TREX_cluster <- function(binned.data,                   # output from TREX() or TREX_bin()
                         feature.data = data.frame()    # same nrow as binned.data, not used in any computations
                                                        # optional, used only to export cluster data for downstream use 
                         ) {

  if (nrow(feature.data) > 0) {
    all.data <- cbind(binned.data, feature.data) 
  } else {
    all.data <- binned.data 
  }
  
  regions.of.interest <- all.data %>%
    dplyr::filter(cuts == "[0,5]" | cuts == "[95,100]")

  db.result = dbscan::dbscan(regions.of.interest[, 1:2], eps = 0.3, minPts = 1)
  track.data = cbind(regions.of.interest, cluster = db.result$cluster)
  track.data <- track.data %>%
    filter(cluster != 0)
  
  cluster.data = split(track.data, track.data$cluster)
  median.percent.change = lapply(cluster.data, function(x) median(x[, which(colnames(track.data)=="percent.change")]))
  mean.percent.change = lapply(cluster.data, function(x) mean(x[, which(colnames(track.data)=="percent.change")]))
  write.csv(mean.percent.change, paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"),"_cluster_ave_percent_change.csv"))

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  values = sample(col_vector)

  range <- apply(apply(track.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  ggplot() +
    geom_point(data = binned.data, aes(x = x, y = y), col = "lightgray") + 
    geom_point(data = track.data, aes(x = x, y = y, col = as.factor(cluster)), cex = 1.5) +
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(values = values) +
    labs(
      x = "UMAP 1", 
      y = "UMAP 2", 
      title = "DBSCAN Clusters (5th & 95th percentiles)", 
      color = "DBSCAN Cluster"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5), nrow = 13)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank())

  # save & export plot as PNG
  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_DBSCAN_plot.png"), 
    width = 8, 
    height = 8
  )
  
  return(track.data)
}

