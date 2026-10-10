# clustering

f.consensus.res <- function(temp_res, k, ...){
  
  temp <- temp_res[[k]]
  
  temp_x <- temp$consensusClass
  
  temp_y <- seq_along(temp_x) %>%
    map_dbl(function(i) {
      same_cluster <- which(temp_x == temp_x[i])
      mean(temp$consensusMatrix[i, same_cluster])
    })
  
  res <- tibble(
    Sample = names(temp_x),
    Cluster = as.integer(temp_x),
    Probability = temp_y
  ) %>%
    tibble::column_to_rownames("Sample")
  
  return(res)
}


f.hclust <- function(df1, dist.method="euclidean", hclust.method="ward.D2", ...){
  
  temp_x <- df1 %>%
    dist(
      method = dist.method
    )
  
  temp_res <- temp_x %>%
    hclust(
      method = hclust.method
    )
  
  return(temp_res)
}


f.hclust.res <- function(temp_res, k=NULL, h=NULL){
  
  temp <- dendro_data(temp_res)
  temp$leaf_labels <- NULL
  temp$class <- NULL
  
  temp <- temp %>%
    map(~{
      res <- .x %>% as.data.frame()
      res
    }) %>%
    as.list()
  
  if(!is.null(k) & is.null(h)){
    temp_x <- cutree(temp_res, k = k) %>%
      as.data.frame() %>%
      setNames(nm = "cluster") %>%
      tibble::rownames_to_column("label") %>%
      mutate(
        cluster = factor(cluster, levels = unique(cluster))
      )
    
    temp$labels <- temp$labels %>%
      left_join(
        temp_x,
        by = "label"
      )
    
    return(temp)
  }
  
  if(is.null(k) & !is.null(h)){
    temp_x <- cutree(temp_res, h = 300) %>%
      as.data.frame() %>%
      setNames(nm = "cluster") %>%
      tibble::rownames_to_column("label") %>%
      mutate(
        cluster = factor(cluster, levels = unique(cluster))
      )
    
    temp$labels <- temp$labels %>%
      left_join(
        temp_x,
        by = "label"
      )
    
    return(temp)
  }
  
  
  if(!is.null(k) & !is.null(h)){
    temp_x <- cutree(temp_res, k = k, h = h) %>%
      as.data.frame() %>%
      setNames(nm = "cluster") %>%
      tibble::rownames_to_column("label") %>%
      mutate(
        cluster = factor(cluster, levels = unique(cluster))
      )
    
    temp$labels <- temp$labels %>%
      left_join(
        temp_x,
        by = "label"
      )
    
    return(temp)
  }
  
  return(temp)
}


f.plot.hclust <- function(temp_res, size, ratio, size_tile=5, color, colors, ...){
  # data
  # P.d <- list()
  P.d <- temp_res
  P.d$labels <- P.d$labels %>%
    mutate(
      ymin = -size_tile,
      ymax = 0
    ) %>%
    rename(
      color = !!sym(color)
    )
  
  # plot
  P1 <- ggplot() +
    geom_segment(
      data = P.d$segments, 
      aes(x = x, y = y, xend = xend, yend = yend),
      size = size
    ) +
    geom_tile(
      data = P.d$labels,
      aes(x = x, y = (ymin + ymax)/2, height = ymax - ymin, fill = color),
      width = 1
    ) +
    # geom_text(
    #   data = P.d$labels, 
    #   aes(x = x, y = 0, label = label, color = Group),
    #   angle = 90, 
    #   hjust = 1, 
    #   vjust = 0.5, 
    #   size = 2
    # ) +
    scale_fill_manual(
      values = colors
    ) +
    theme_dendro() +
    theme(
      aspect.ratio = ratio,
      axis.line.y = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      # axis.title.y = element_text(size = 11, color = "black"),
      legend.title = element_blank(),
      # legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    scale_x_continuous(
      expand = c(0, 2, 0, 0)
    ) +
    scale_y_continuous(
      # limits = c(0, NA),
      # breaks = 
      # expand = c(0, 0, 0.05, 0)
    )
  
  # print(P1)
  
  return(P1)
}




























