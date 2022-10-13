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
thirteen_colors <- c("navyblue","#3EbBCE")