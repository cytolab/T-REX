library(FNN)
library(fancycut)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

source("R/themes_and_palettes.R")
source("R/TREX_helpers.R")

TREX <- function(embedding.data,                   
                 kvalue = 60,                       
                 bins = c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]')           
                 ) {
  
  # warnings for improperly formatted data 
  if (is.na(match("file_ID", colnames(embedding.data)))) {
    stop("Embedding.data has no column named file_ID.")
  }

  data.name = unique(embedding.data$file_ID)[1]
  if (length(which(embedding.data$file_ID == data.name)) != nrow(embedding.data)/2) {
    stop("Datasets are not equally sampled according to file_ID.")
  }
  
  # create bins vector if given integer or % 
  if (is.numeric(bins)) {
    bins <- get_TREX_bins(round(100/bins, digits = 1))
  } else if(is.character(bins) & length(bins) == 1) {
    if(grepl("%", bins)) { 
      n_bins = as.numeric(substr(bins, start = 1, stop = nchar(bins) - 1))
      if (n_bins > 30) {
        stop("Cannot have less than three bins.")
      }
      bins <- get_TREX_bins(n_bins)
    }
  }
  
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
    percent.change = percent.change
  )
  binned.data$cuts <- wafflecut(round(binned.data$percent.change), bins)
  
  return(binned.data)
}

TREX_plot <- function(binned.data,
                      title.height = -3,
                      embed.type = "Embedding",
                      percent.labels = TRUE,
                      caption = NULL,
                      return.obj = "drawn plot",
                      export = FALSE) {

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
  
  if (export) {
    ggsave(
      paste0(
        strftime(Sys.time(), "%Y-%m-%d_%H_%M "),
        data.names[1], " vs ", data.names[2],
        " T-REX plot.png"
      ), 
      plot = trex.titled, 
      width = 8, 
      height = 8
    )
  }
  
  if (return.obj == "drawn plot") { 
    return(ggdraw(trex.titled))
  } else if (return.obj == "ggplot") {
    return(trex.plot)
  } else if (return.obj == "gtable") {
    return(trex.titled)
  }
  
}

TREX_results <- function(binned.data, 
                         export = FALSE) {
  
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
  
  if (export) {
    write.csv(sample.table, paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"), "_TREX_results.csv"), row.names = FALSE)
  }
  
  return(sample.table) 
}
