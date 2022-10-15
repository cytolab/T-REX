library(FNN)
library(fancycut)
library(ggplot2)
library(dbscan)
library(grid)
library(gridExtra)

source("R/themes_and_palettes.R")

TREX <- function(embedding.data,                   
                 kvalue = 60,                       
                 bins = c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]')           
                 ) {
  
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
                      percent.labels = TRUE,
                      title = "dataset 1 vs. dataset 2", 
                      caption = NULL) {

  trex.data = binned.data 
  bins = levels(binned.data$cuts)
  
  # reorder bins so that regions of interest are on top 
  if (length(bins) %% 2 == 0) {
    center.indx = length(bins)/2
    bin.levels = c(bins[center.indx], bins[center.indx + 1]) 
    for (i in 1:(center.indx - 1)) {
      bin.levels <- append(bin.levels, bins[center.indx - i])
      bin.levels <- append(bin.levels, bins[center.indx + i + 1])
    }
  } else {
    center.indx = round(length(bins)/2, 0) + 1
    bin.levels = bins[center.indx]
    for (i in 1:(center.indx - 1)) {
      bin.levels <- append(bin.levels, bins[center.indx - i])
      bin.levels <- append(bin.levels, bins[center.indx + i])
    }
  }
  
  trex.data$cuts <- factor(trex.data$cuts, levels = bin.levels)
  trex.data <- trex.data[order(trex.data$cuts), ]
  range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  if (percent.labels) {
    bin.labels = get_percent_labels(bins, title)
  } else {
    bin.labels = bins
  }
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  trex.plot = ggplot(trex.data) + 
    geom_point(aes(x = x, y = y, colour = cuts), cex = 1) + 
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(
      labels = bin.labels,
      values = get_TREX_colors(bins)) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x =  embed.x, y = embed.y, caption = caption) +
    theme_TREX()
  
  # add red/blue colored title 
  titleGrobs <- grobTree(
    gp = gpar(fontsize = 12, fontface = "bold"),
    textGrob(label = str_split(title, " vs. ")[[1]][1], 
             name = "title1",
             x = unit(5, "lines"), 
             y = unit(-0.5, "lines"), 
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#000080")),
    textGrob(label = " vs. ", name = "title2",
             x = grobWidth("title1") + unit(5, "lines"), 
             y = unit(-0.5, "lines"),
             hjust = 0, vjust = 0),
    textGrob(label = str_split(title, " vs. ")[[1]][2], 
             name = "title3",
             x = grobWidth("title1") + grobWidth("title2") + unit(5, "lines"), 
             y = unit(-0.5, "lines"),
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#8B0000"))
  )
  trex.titled <- arrangeGrob(trex.plot, top = titleGrobs, padding = unit(2.6, "line"))

  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_TREX_plot.png"), 
    plot = trex.titled, 
    width = 8, 
    height = 8
  )
  
  return(trex.titled)
}

TREX_stats <- function(binned.data) {

  sample.table <- data.frame(total_cells = nrow(binned.data))
  
  sums = vector()
  for (i in levels(binned.data$cuts)) {
    sums <- append(sums, sum(binned.data$cuts == i))
  }

  sample.table[, c(2:length(levels(binned.data$cuts)))] <- sums
  colnames(sample.table)[2:length(levels(binned.data$cuts))] <- levels(binned.data$cuts)
  
  percent = 100*(sums/nrow(binned.data))
  sample.table$degree_of_change = (sum(percent[3] + percent[5]))
  sample.table$direction_of_change = (sums[5] - sums[3]) / (sums[5] + sums[3])
  
  write.csv(sample.table, paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"), "_TREX_stats.csv"))
  
  return(sample.table) 
}

# use DBSCAN to cluster on regions of great change (5th and 95 percentiles of change)
TREX_cluster <- function(binned.data, 
                         bins.of.interest = NULL,
                         marker.data = NULL) {           

  if (!is.null(marker.data)) {
    binned.data <- cbind(binned.data, marker.data)
  }
  
  if (!is.null(bins.of.interest)) {
    bin.levels = levels(binned.data$cuts)
    bins.of.interest <- c(bin.levels[1], bin.levels[length(bin.levels)])
  }  

  str(bins.of.interest)
    
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

  return(track.data)
}

TREX_cluster_plot <- function() {
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
  
  ggsave(
    paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_DBSCAN_plot.png"), 
    width = 8, 
    height = 8
  )
  
  
}

