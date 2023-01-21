get_TREX_bins <- function(bin.size) {
  if (bin.size >= 100 | bin.size <= 0) {
    stop("To generate bins from a percentage or number, the value must be between 1 and 100")
  }
  
  n_bins = round(100/bin.size, digits = 0)
  bin.vector = vector("character", n_bins)
  start.b = "("
  end.b = "]"
  bin.vector[1] <- paste0("[0,",bin.size,"]")
  center.bin = floor(n_bins/2)
  
  # handle special case of only 3 bins 
  if (center.bin == 1) {
    center.bin <- 2 
  }
  
  for(i in 2:n_bins) {
    lower.num = (i - 1)*bin.size
    upper.num = i*bin.size 
    
    # last bin in list
    if(i == n_bins) {              
      bin.vector[i] <- paste0("[",lower.num,",100]") 
      
      # center bin
    } else if(i == center.bin) {   
      bin.vector[i] <- paste0("(",lower.num,",",upper.num,")") 
      start.b <- "["
      end.b <- ")"
      
      # all other bins
    } else {                            
      bin.vector[i] <- paste0(start.b,lower.num,",",upper.num,end.b)         
    }
  } 
  
  return(bin.vector)
}


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