library (tidyverse)
library(RColorBrewer)

get_percent_labels <- function(bin.levels, data.names) {
  
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
      cluster.label = paste0("\u2265", (100 - second.num), "% from ", data.names[1])
    } else if (second.num == 100) {
      cluster.label = paste0("\u2265", first.num, "% from ", data.names[2])
    } else if (second.num <= 50 & first.num <= 50) {
      cluster.label = paste0((100 - second.num), "-", (100 - first.num), "% from ", data.names[1]) 
    } else if (second.num > 50 & first.num > 50) {
      cluster.label = paste0(first.num, "-", second.num, "% from ", data.names[2])
    } else {
      cluster.label = paste0("from ", data.names[1], " and ", data.names[2])
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

tatarize_preferred = c("#3A2465","#7900D7","#671190","#C6005A","#6b002C","#FF1A59","#E83000","#D86A78","#51A058",
                       "#003109","#00D891","#00846F","#1CE6FF","#02525F","#0000A6","#0086ED","#FF6832","#772600",
                       "#FF913F","#372101","#D68E01","#222800","#4C6001","#A4E804")

tatarize_secondary = c("#B88183","#922329","#FF8A9A","#E20027","#A74571","#FF5DA7","#FF90C9","#A30059",
                       "#5B113C","#D157A0","#962B75","#6B3A64","#6F0062","#320033","#BC23FF","#BC65E9",
                       "#72418F","#4A3B53","#9556BD","#A079BF","#9F94F0","#7A7BFF","#3C3E6E","#6367A9",
                       "#3B5DFF","#0045D2","#0060CD","#012C58","#0AA3F7","#29607C","#006FA6","#0AA6D8",
                       "#5EBCD1","#0089A3","#203B3C","#00C6C8","#006A66","#00C2A0","#00FECF","#004B28",
                       "#0CBD66","#1BE177","#006C31","#04F757","#001E09","#5EFF03","#3B9700","#4FC601",
                       "#1B4400","#8BB400","#6B7900","#3A3F00","#FFFF00","#575329","#7E6405","#513A01",
                       "#FFB500","#A77500","#7A4900","#A45B02","#D16100","#5B3213","#CA834E","#BE4700",
                       "#391406","#C86240","#B77B68","#643127","#BF5650","#BA0900","#FF4A46","#000000","#452C2C")

tatarize_avoid = c("#6D80BA","#6367a9","#636375","#bcb1e5","#837393","#D790FF","#FFA0F2","#AA5199",
                   "#402334","#7C6571","#B894A6","#997D87","#FEB2C6","#59738A","#000035","#BCB1E5",
                   "#D83D66","#5B4E51","#797868","#61615A","#6A714A","#658188","#00005F","#F4D749",
                   "#868E7E","#C2ff99","#55813B","#83A485","#A3F3AB","#66796D","#63FFAC","#1e0220",
                   "#8ADBB4","#78AFA1","#809693","#00A6AA","#00A6AA","#C8A1A1","#9C6966","#89412E",
                   "#FAD09F","#CCD87C","#9FA064","#788D66","#15A08A","#E98176","#806C66","#EA8B66",
                   "#7D5A44","#B5F4FF","#6B94AA","#9695C5","#A38469","#CCb87C","#8D8546","#F4ABAA",
                   "#895563","#885578","#FFAA92","#A05837","#CCB87C","#DFFB71","#C2FF99","#953F00",
                   "#88bF4C","#453C23","#B79762","#1E0200","#7B4F4B","#5B4534","#886F4C","#6C8F7D",
                   "#47675D","#9B9700","#83AB58","#98D058","#A97399","#CB7E98","#B05B6F","#943A4D",
                   "#6A3A4C","#5A0007","#E773CE","#CC0744","#374527","#71BB8C","#CCAA35","#64547B",
                   "#7ED379","#4B8160","#3D4F44","#456648","#76912F","#34362D","#300018","#3B000A",
                   "#6B002C","#00B433","#549E79","#02684E","#518A87","#1A3A2A","#1E6E00","#004D43",
                   "#001C1E","#8CD0FF","#353339","#001325","#456D75","#66E1D3","#061203","#008941",
                   "#404E55","#494B5A","#201625","#29201D","#E7AB63","#1D1702","#3E89BE","#324E72",
                   "#8FB0FF","#C895C5","#A76F42","#013349","#00489C") 

tatarize_optimized <- function(num_colors) {
  set.seed(16)
  colors_list = c()
  if (num_colors <= 24) {
    colors_list <- sample(tatarize_preferred, num_colors)
  } else if (num_colors <= 97) {
    colors_list <- sample(tatarize_preferred, 24)
    colors_list <- append(colors_list, sample(tatarize_secondary, (num_colors - 24))) 
  } else {
    colors_list <- sample(tatarize_preferred, 24)
    colors_list <- append(colors_list, sample(tatarize_secondary, 73))
    colors_list <- append(colors_list, sample(tatarize_avoid, (num_colors - 97)))
  }
  return(colors_list)
}



