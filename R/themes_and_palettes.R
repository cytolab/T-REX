library (tidyverse)
library(RColorBrewer)

get_percent_labels <- function(bin.levels, title) {
  title.one = str_split(title, " vs. ")[[1]][1]
  title.two = str_split(title, " vs. ")[[1]][2]
  
  bin.labels = vector()
  for (bin in bin.levels) {
    num.list = strsplit(bin, split = ",")
    first.num = num.list[[1]][1]
    second.num = num.list[[1]][2]
    
    # check for decimal points 
    if (str_detect(first.num, "\\.")) {
      first.num <- as.numeric(str_extract(first.num, "[0-9]+\\.[0-9]+"))
    } else {
      first.num <- as.numeric(str_extract(first.num, "[0-9]+"))
    }
    if (str_detect(second.num, "\\.")) {
      second.num <- as.numeric(str_extract(second.num, "[0-9]+\\.[0-9]+"))
    } else {
      second.num <- as.numeric(str_extract(second.num, "[0-9]+"))
    }
    
    if (first.num == 0) {
      cluster.label = paste0("\u2265", (100 - second.num), "% from ", title.one)
    } else if (second.num == 100) {
      cluster.label = paste0("\u2265", first.num, "% from ", title.two)
    } else if (second.num <= 50 & first.num <= 50) {
      cluster.label = paste0((100 - second.num), "-", (100 - first.num), "% from ", title.one) 
    } else if (second.num > 50 & first.num > 50) {
      cluster.label = paste0(first.num, "-", second.num, "% from ", title.two)
    } else {
      cluster.label = paste0("from ", title.one, " and ", title.two)
    }
    bin.labels <- append(bin.labels, cluster.label)
  }
  return(bin.labels)
}

default_colors <- c("navyblue","lightskyblue","lightgray","lightcoral","darkred")
thirteen_colors <- c("navyblue","#3E6BCE","lightskyblue","#D6F3FC","#AEC2CC",
                     "#79878E","#5D5D5D","#8E7979","#CCAEAE","#FFDEDE",
                     "lightcoral","#C45555","darkred")

get_TREX_colors <- function(bins) {
  if (length(bins) == 5) {
    bin.colors = default_colors
  } else if(length(bins) > 5 & length(bins) < 13) {
    bin.colors = colorRampPalette(default_colors)(length(bins))
  } else if(length(bins) == 13) {
    bin.colors = thirteen_colors
  } else if (length(bins) > 13) {
    library(RColorBrewer)
    bin.colors = colorRampPalette(thirteen_colors)(length(bins))
  }
  names(bin.colors) <- bins
  return(bin.colors)
}

theme_TREX <- function() {
  theme_bw() %+replace% 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 14),
      legend.title = element_blank()
    )
}