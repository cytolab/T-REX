library (tidyverse)


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

default_colors <- c("navyblue","lightskyblue","lightgray","lightcoral","darkred")
thirteen_colors <- c("navyblue","#3E6BCE","lightskyblue","#D6F3FC","#AEC2CC",
                     "#79878E","#5D5D5D","#8E7979","#CCAEAE","#FFDEDE",
                     "lightcoral","#C45555","darkred")
