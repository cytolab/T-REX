library(FNN)
library(fancycut)
library(ggplot2)
library(RColorBrewer)
library(dbscan)

source("R/themes_and_palettes.R")

TREX <- function(embedding.data,                   # dataframe, 2 columns of low dim embedding + file_ID for which dataset a cell belongs to
                 kvalue = 60,                      # number, for KNN search 
                 bins = c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]')           
                 ) {
  
  # set up warning about 3 min bins 
  
  binned.data = TREX_bin(embedding.data, kvalue, bins)
  
  # TREX_stats(binned.data)

  return(binned.data)
}

TREX_bin <- function(embedding.data, 
                     kvalue,
                     bins) {
  
  # KNN search per cell 
  neighbor.index = knnx.index(embedding.data[, 1:2], embedding.data[, 1:2], k = kvalue)
  neighbor.index[neighbor.index <= nrow(embedding.data)/2] <- 0
  neighbor.index[neighbor.index > nrow(embedding.data)/2] <- 1
  
  # calculate percent change in each KNN region
  percent.change = (rowSums(neighbor.index) / kvalue * 100)
  
  # binning 
  binned.data <- data.frame(x = embedding.data[, 1], y = embedding.data[, 2], percent.change = round(percent.change))
  binned.data$cuts <- wafflecut(binned.data$percent.change, bins)
  
  return(binned.data)
}

TREX_plot <- function(binned.data,
                      embed.type = "Embedding",
                      caption = NULL) {

  trex.plot = binned.data 
  bins = levels(binned.data$cuts)
  
  # reorder bins so that unchanging areas are on the bottom geom_layer, and regions of interest are on top 
  center.indx = round(length(bins)/2, 0) + 1
  bin.levels = bins[center.indx]
  for (i in 1:(center.indx - 1)) {
    bin.levels <- append(bin.levels, bins[center.indx - i])
    bin.levels <- append(bin.levels, bins[center.indx + i])
  }
  
  # if (length(bins) %% 2 == 0) {} 
  # else {}
  
  trex.plot$cuts <- factor(trex.plot$cuts, levels = bin.levels)
  trex.plot <- trex.plot[order(trex.plot$cuts), ]
  range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  # select a palette based on number of bins 
  if (length(bins) == 5) {
    bin.colors = default_colors
  } else if(length(bins) == )

  names(bin.colors) <- bins
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  ggplot(trex.plot) + 
    geom_point(aes(x = x, y = y, colour = cuts), cex = 1) + 
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(
      labels = c("\u2265 95% from 1st dataset",
                 "85-95% from 1st dataset",
                 "from 1st and 2nd dataset",
                 "85-95% from 2nd dataset",
                 "\u2265 95% from 2nd dataset"),
      values = bin.colors) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x =  embed.x, y = embed.y,
         title = paste0(sample_type, "_", sample_id, "_", time_comparison, " - Percent Change"),
         caption = caption) +
    theme_TREX()
  
  # save & export plot as PNG
  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H%M"), "_TREX_plot.png"), 
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
  write.csv(sample.table, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"), "_TREX_stats.csv"))
  
  return(sample.table) 
}

# use DBSCAN to cluster on regions of great change (5th and 95 percentiles of change)
TREX_cluster <- function(binned.data, 
                         marker.data = NULL) {

  if (!is.null(marker.data)) {
    binned.data <- cbind(binned.data, marker.data)
  }
  
  regions.of.interest <- binned.data %>%
    dplyr::filter(cuts == "[0,5]" | cuts == "[95,100]")
  db.result = dbscan::dbscan(regions.of.interest[, 1:2], eps = 0.3, minPts = 1)
  track.data = cbind(regions.of.interest, cluster = db.result$cluster)
  track.data <- track.data %>%
    filter(cluster != 0)
  
  cluster.data = split(track.data, track.data$cluster)
  median.percent.change = lapply(cluster.data, function(x) median(x[, which(colnames(track.data) == "percent.change")]))
  mean.percent.change = lapply(cluster.data, function(x) mean(x[, which(colnames(track.data) == "percent.change")]))
  write.csv(mean.percent.change, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_ave_percent_change.csv"))

  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  set.seed(1)
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
    theme_TREX()

  # save & export plot as PNG
  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H%M"), "_DBSCAN_plot.png"), 
    width = 8, 
    height = 8
  )
  
  return(track.data)
}

