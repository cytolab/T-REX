library(FNN)
library(fancycut)
library(ggplot2)
library(dbscan)
library(grid)
library(gridExtra)
library(Palo)

source("R/themes_and_palettes.R")

TREX <- function(embedding.data,                   
                 kvalue = 60,                       
                 bins = c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]')           
                 ) {
  
  # warnings for improperly formatted data 
  data.name = unique(embedding.data$file_ID)[1]
  if (length(which(embedding.data$file_ID == data.name)) != nrow(embedding.data)/2) {
    stop("Datasets are not equally sampled according to file_ID.")
  }
  
  # create # of bins from integer 
  if (is.numeric(bins)) {
    bin.size = round(100/bins, digits = 1)
    new.bins = vector("character", bins)
    start.b = "("
    end.b = "]"
    new.bins[1] <- paste0("[0," , bin.size, "]")
    for(i in 1:(bins - 1)) {
      if(i == (bins - 1)) {
        new.bins[bins] <- paste0("[", i*bin.size, ",100]")
      } else if(i == floor(bins/2)) {
        new.bins[floor(bins/2) + 1] <- paste0("(", i*bin.size, ",", (i*bin.size + bin.size), ")")
        start.b <- "["
        end.b <- ")"
      } else {
        new.bins[i + 1] <- paste0(start.b, i*bin.size, ",", (i*bin.size + bin.size), end.b)
      }
    } 
    bins <- new.bins
    cat(bins)
  }  
  
  # create bins of a given % 
  
  
  # KNN search per cell 
  neighbor.index = knnx.index(embedding.data[, 1:2], embedding.data[, 1:2], k = kvalue)
  
  # assign dataset belonging by row number
  neighbor.index[neighbor.index <= nrow(embedding.data)/2] <- 0
  neighbor.index[neighbor.index > nrow(embedding.data)/2] <- 1

  # calculate percent change in each KNN region
  percent.change = (rowSums(neighbor.index) / kvalue * 100)
  
  # binning 
  binned.data <- data.frame(
    x = embedding.data[, 1],
    y = embedding.data[, 2],
    file_ID = embedding.data$file_ID,
    percent.change = round(percent.change)
  )
  binned.data$cuts <- wafflecut(binned.data$percent.change, bins)

  return(binned.data)
}

TREX_plot <- function(binned.data,
                      title.height = -3,
                      embed.type = "Embedding",
                      percent.labels = TRUE,
                      caption = NULL) {

  # get dataset names from file_ID column 
  data.names = unique(binned.data$file_ID)
  
  # reorder bins so that regions of interest are on top 
  bins = levels(binned.data$cuts)
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
  
  binned.data$cuts <- factor(binned.data$cuts, levels = bin.levels)
  binned.data <- binned.data[order(binned.data$cuts), ]
  range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  if (percent.labels) {
    bin.labels = get_percent_labels(bins, data.names)
  } else {
    bin.labels = waiver()
  }
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  trex.plot = ggplot(binned.data) + 
    geom_point(aes(x = x, y = y, color = cuts), cex = 1) + 
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(
      labels = bin.labels,
      values = get_TREX_colors(bins),
      breaks = bins
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x =  embed.x, y = embed.y, caption = caption) +
    theme_TREX()
  
  # add red/blue colored title 
  titleGrobs <- grobTree(
    gp = gpar(fontsize = 16, fontface = "bold"),
    textGrob(label = data.names[1], 
             name = "title1",
             x = unit(3, "lines"), 
             y = unit(title.height, "lines"), 
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#000080")),
    textGrob(label = " vs. ", 
             name = "title2",
             x = grobWidth("title1") + unit(3, "lines"), 
             y = unit(title.height, "lines"),
             hjust = 0, vjust = 0),
    textGrob(label = data.names[2], 
             name = "title3",
             x = grobWidth("title1") + grobWidth("title2") + unit(3, "lines"), 
             y = unit(title.height, "lines"),
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#8B0000"))
  )
  trex.titled <- arrangeGrob(trex.plot, top = titleGrobs, padding = unit(2.6, "line"))

  ggsave(
    paste0(
      strftime(Sys.time(), "%Y-%m-%d_%H-%M"),
      data.names[1], " vs ", data.names[2],
      " T-REX plot.png"
    ), 
    plot = trex.titled, 
    width = 8, 
    height = 8
  )
  
  return(ggdraw(trex.titled))
}

TREX_stats <- function(binned.data) {
  sample.table <- data.frame(total_cells = nrow(binned.data))
  sums = vector()
  num_bins = length(levels(binned.data$cuts))
  for (i in levels(binned.data$cuts)) {
    sums <- append(sums, sum(binned.data$cuts == i))
  }
  sample.table[, c(2:(num_bins + 1))] <- sums
  colnames(sample.table)[2:(num_bins + 1)] <- levels(binned.data$cuts)
  percent = 100*(sums/nrow(binned.data))
  sample.table$degree_of_change = (sum(percent[1] + percent[length(percent)]))
  sample.table$direction_of_change = (sums[length(sums)] - sums[1]) / (sums[length(sums)] + sums[1])
  write.csv(sample.table, paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"), "_TREX_stats.csv"), row.names = FALSE)
  return(sample.table) 
}

TREX_cluster <- function(binned.data, 
                         bins.of.interest = NULL,
                         db.eps = 0.5,
                         marker.data = NULL) {           

  if (!is.null(marker.data)) {
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
  
  results.data = split(cluster.data, cluster.data$cluster)
  median.percent.change = lapply(results.data, function(x) median(x[, which(colnames(x) == "percent.change")]))
  mean.percent.change = lapply(results.data, function(x) mean(x[, which(colnames(x) == "percent.change")]))
  write.csv(mean.percent.change, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_ave_percent_change.csv"))
 
  return(cluster.data)
}

TREX_cluster_plot <- function(cluster.data,
                              binned.data = NULL,
                              embed.type = "Embedding",
                              colors = NULL) {
  
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
  
  ggsave(
    filename = paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_DBSCAN_plot.png"), 
    plot = cluster.plot,
    width = 8, 
    height = 8
  )
  
  return(cluster.plot)
}

# TREX_counts <- function(cluster.data,
#                         bins.of.interest = NULL) {
#   
#   if (is.null(bins.of.interest)) {
#     bin.levels = levels(cluster.data$cuts)
#     bins.of.interest <- c(bin.levels[1], bin.levels[length(bin.levels)])
#   } 
#   
#   split.data = split(cluster.data, as.factor(cluster.data$cluster))
#   counts = sapply(split.data, NROW)
#   
#   bin.count = vector()
#   for (i in 1:length(split.data)) {
#     bin.count[i] = sum(split.data[[i]][["cuts"]] == "[0,5]")
#   }
#   counts.95 = counts - bin.count
#   str(counts.95)
  
  # all.clusters = split(cluster.data, as.factor(cluster.data$cluster))
  # counts.total  <- sapply(all.clusters, NROW)
  # counts.5 <- vector()
  # for (i in 1:length(all.clusters)){
  #   counts.5[i] = sum(all.clusters[[i]][["status"]] == '[0,5]')
  # }
  # counts.95 = counts.total-counts.5
  # cluster.data = as.data.frame(counts.total)
  # colnames(cluster.data) <- "total_counts"
  # cluster.data$counts_5 <- counts.5
  # cluster.data$counts_95 <- counts.95
  # write.csv(cluster.data, paste(strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_cluster_counts.csv"))
  # return(cluster.data)
# }

