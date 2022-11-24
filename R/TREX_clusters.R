library(ggplot2)
library(dbscan)

TREX_cluster <- function(binned.data, 
                         bins.of.interest = NULL,
                         db.eps = 0.5,
                         marker.data = NULL) {           
  
  if (!is.null(marker.data)) {
    marker.data <- subset(marker.data, select = -c(file_ID))
    binned.data <- cbind(binned.data, marker.data)
  }
  
  if (is.null(bins.of.interest)) {
    bin.levels = levels(binned.data$cuts)
    bins.of.interest <- c(bin.levels[1], bin.levels[length(bin.levels)])
  }  
  
  regions.of.interest = binned.data[binned.data$cuts %in% bins.of.interest, ]
  db.result = dbscan::dbscan(regions.of.interest[, 1:2], eps = db.eps, minPts = 1)
  cluster.data = cbind(regions.of.interest, cluster = db.result$cluster)
  cluster.data <- cluster.data %>%
    filter(cluster != 0)
  
  return(cluster.data)
}

TREX_cluster_plot <- function(cluster.data,
                              binned.data = NULL,
                              embed.type = "Embedding",
                              colors = NULL,
                              export = FALSE) {
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  cluster.plot = ggplot()
  color.values = tatarize_optimized(length(unique(as.factor(cluster.data$cluster))))
  
  if (!is.null(colors)) {
    color.values <- colors
  } 
  
  # use all embedding points to calculate ratio if included 
  if (!is.null(binned.data)) {
    cluster.plot <- cluster.plot + geom_point(data = binned.data, aes(x = x, y = y), col = "lightgray")
    color.values <- append(color.values, "lightgray")
    range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
    graphical.ratio <- (range[1] / range[2])
  } else {
    range <- apply(apply(cluster.data[, 1:2], 2, range), 2, diff)
    graphical.ratio <- (range[1] / range[2])
  }
  
  cluster.plot <- cluster.plot + 
    geom_point(data = cluster.data, aes(x = x, y = y, col = as.factor(cluster)), cex = 1.5) +
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(values = color.values) +
    labs(
      x = embed.x, 
      y = embed.y, 
      title = "DBSCAN Clusters, bins of interest"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_TREX()
  
  if (export) {
    ggsave(
      filename = paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_cluster_plot.png"), 
      plot = cluster.plot,
      width = 8, 
      height = 8
    )
  }
  
  return(cluster.plot)
}

TREX_cluster_results <- function(cluster.data,
                                 export = FALSE) {
  
  results.data = split(cluster.data, cluster.data$cluster)
  median.percent.change = lapply(results.data, function(x) median(x[, which(colnames(x) == "percent.change")]))
  mean.percent.change = lapply(results.data, function(x) mean(x[, which(colnames(x) == "percent.change")]))
  
  if (export) {
    write.csv(mean.percent.change, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_ave_percent_change.csv"))
  }
  
  return(mean.percent.change)
}

TREX_counts <- function(cluster.data,
                        export = TRUE) {
  
  results = data.frame(cluster = unique(cluster.data$cluster))
  
  # count total cells in each cluster
  split.clusters = split(cluster.data, as.factor(cluster.data$cluster))
  results$total_cells <- sapply(split.clusters, NROW)
  
  # find bins of interest
  bins.of.interest = vector()
  for(i in split.clusters) {
    for (z in unique(i$cuts)) {
      if (!z %in% bins.of.interest) {
        bins.of.interest <- append(bins.of.interest, z)
      }
    }
  }
  bins.of.interest <- bins.of.interest[order(match(bins.of.interest, as.vector(levels(cluster.data$cuts))))]
  
  # count how many cells are in each bin 
  results <- cbind(results, data.frame = (matrix(0, ncol = length(bins.of.interest))))
  colnames(results)[3:ncol(results)] <- bins.of.interest
  
  col_indx = 3
  for (w in bins.of.interest) {
    for (i in 1:length(split.clusters)) {
      results[i, col_indx] <- sum(split.clusters[[i]][["cuts"]] == w)
    }
    col_indx <- col_indx + 1
  }
  
  if (export) {
    write.csv(results, paste(strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_cluster_counts.csv"))
  }
  
  return(results)
}