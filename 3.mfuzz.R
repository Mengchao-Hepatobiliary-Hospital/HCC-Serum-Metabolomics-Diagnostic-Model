## mfuzz ----
f.mfuzz_df <- function(df1, stat="mean", ...){
  group <- names(df1)[sapply(df1, is.factor)]
  
  temp_x <- df1 %>%
    pivot_longer(
      cols = -all_of(group), 
      names_to = "Feature", 
      values_to = "Value"
    ) %>%
    group_by(across(all_of(group)), Feature) %>%
    summarise(
      Value = switch(
        stat,
        "mean" = mean(Value, na.rm = TRUE),
        "median" = median(Value, na.rm = TRUE),
        stop("Invalid choice.")
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = all_of(group),
      values_from = "Value"
    ) %>%
    tibble::column_to_rownames("Feature")
  
  class(temp_x) <- c(class(temp_x), paste0("stat=", stat))
  
  return(temp_x)
}


f.mfuzz_eset <- function(df1, ...){
  eset <- new("ExpressionSet", exprs = as.matrix(df1))
  {
    # eset <- filter.NA(eset, thres = 0.25)
    # eset <- fill.NA(eset, mode = "knn")
    # eset <- filter.std(eset, min.std = 0)
    eset <- standardise(eset)                
  }
  
  return(eset)
}


f.mfuzz <- function(data, centers, ...){
  
  if("ExpressionSet" %in% class(data)){
    eset <- data
  }else{
    eset <- f.mfuzz_eset(data)
  }
  
  set.seed(123)
  temp_res <- eset %>%
    mfuzz(
      eset = .,
      centers = centers, 
      m = mestimate(.)
    )
  
  temp_res$eset <- eset
  
  return(temp_res)
}


f.mfuzz_res <- function(temp_res, cutoff, ...){
  # membership and cluster information
  mc <- temp_res$membership %>%
    data.frame(check.names = F) %>% 
    mutate(
      `Cluster` = names(.)[max.col(., ties.method = "first")],
      `Membership` = apply(., 1, max)
    )
  
  # summary
  sm <- exprs(temp_res$eset) %>% as.data.frame() %>% 
    mutate(
      Cluster = mc$Cluster[match(rownames(.), rownames(mc))],
      Membership = mc$Membership[match(rownames(.), rownames(mc))]
    ) %>%
    tibble::rownames_to_column("Feature") %>%
    pivot_longer(
      cols = all_of(colnames(temp_res$centers)), 
      names_to = "Timer", 
      values_to = "Expression"
    ) %>%
    mutate(
      Cluster = paste0("Cluster ", Cluster),
      Cluster = factor(Cluster, levels = sort(unique(Cluster))),
      Timer = factor(Timer, levels = colnames(temp_res$centers)),
      Feature = factor(Feature, levels = unique(Feature) )
    ) %>%
    arrange(Cluster) %>%
    group_by(Timer, Cluster) %>%
    mutate(
      Count = n()
    ) %>%
    ungroup()
  
  return(sm)
}


f.add_belong <- function(df1, col, cutoff, ...){
  
  temp_res <- df1 %>%
    mutate(
      Belong = if_else(!!sym(col) > cutoff, "T", "F"),
      cutoff = cutoff
    ) %>%
    mutate(
      Belong = factor(Belong, levels = c("T", "F"))
    )
  
  return(temp_res)
}

f.mfuzz_palette <- function(palette="default", ...){
  colors = switch(
    palette,
    "default" = c(
      "#FF8F00", "#FFA700", "#FFBF00",
      "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
      "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00",
      "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28",
      "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F",
      "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF",
      "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF",
      "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF",
      "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
      "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF",
      "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F",
      "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018"
    ),
    "fancy" = {
      fancy.blue  <- c(c(255:0),rep(0,length(c(255:0))),rep(0,length(c(255:150))))
      fancy.green  <- c(c(0:255),c(255:0),rep(0,length(c(255:150))))
      fancy.red  <- c(c(0:255),rep(255,length(c(255:0))),c(255:150))
      rgb(b=fancy.blue/255,g=fancy.green/255,r=fancy.red/255)
    },
    stop("Invalid choice.")
  )

  return(colors)
}


f.plot.mfuzz.hline <- function(temp_res, colors=NULL, palette="default", theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      label = paste0(Cluster, sprintf(" (n = %d)", Count))
    )
  
  if(is.null(colors)){
    colors <- f.mfuzz_palette(palette)
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Timer, y = Expression, color = Membership) ) +
    facet_grid(
      .~ label,
      # labeller = labeller(
      #   Cluster = setNames(P.d$df$label, P.d$df$Cluster)
      # )
    ) +
    geom_line(
      aes(group = Feature)
    )+
    scale_color_gradientn(
      colors = colors,
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1)
    ) +
    stat_summary(
      data = {subset(P.d$df, Belong == "T")},
      aes(x = levels(P.d$df$Timer)[1], y = Expression),
      fun.data = function(y) {
        data.frame(
          y = Inf,
          label = paste0("Membership > ", unique(P.d$df$cutoff), sprintf(" (n = %d)", length(y) ))
        )
      },
      geom = "text",
      hjust = 0, vjust = 1.3,
      size = 3.3,
      color = "firebrick"
    ) +
    f.theme(theme) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      # axis.title.y = element_blank(),
      legend.position = "top",         
      strip.placement = "outside",
      strip.text = element_text(size = 10, color = "black"), 
      # strip.text.x = element_blank()
      # strip.background.x = element_rect(fill = NA, color = NA),
      # strip.background = element_rect(fill = NA, color = NA)
    ) +
    scale_x_discrete(
      expand = c(0.07, 0.03, 0.07, 0.03)
    ) +
    labs(
      y = "Expression Changes"
    ) +
    guides(
      color = guide_colorbar(
        barwidth = 5,
        barheight = 1
      )
    )
    
  # print(P1)
  
  return(P1)
}


f.plot.mfuzz.heat <- function(temp_res, colors=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      label = paste0(Cluster, sprintf(" (n = %d)", Count))
    )
  
  if(is.null(colors)){
    colors <- c("blue", "red")
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Timer, y = Feature, fill = Expression)) +
    facet_wrap(
      . ~ label, 
      ncol = 1, 
      scales = "free_y",
      strip.position = "left",
      # labeller = labeller(
      #   Cluster = setNames(P.d$df$label, P.d$df$Cluster)
      # )
    ) +
    geom_tile() +
    scale_fill_gradient2(
      low = colors[1], high = colors[2], mid = "white",
      # breaks = seq(-3, 3, 0.5),
      limits = c(min(P.d$df$Expression), max(P.d$df$Expression))
    ) +
    scale_x_discrete(
      expand = c(0, 0, 0, 0)
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = 5,
        barheight = 1
      )
    ) +
    f.theme(theme) +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "top",  
      strip.placement = "outside", 
      # strip.text = element_blank(), 
      strip.background = element_rect(fill = NA, color = NA) 
    ) +
    labs(
      fill = "Expression Changes",
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.mfuzz.vlh <- function(temp_res, l.colors=NULL, h.colors=NULL, palette="default", theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- temp %>%
    f.plot.mfuzz.hline(
      palette = palette,
      colors = l.colors,
      theme = theme
    ) +
    facet_wrap(
      .~label, 
      # labeller = labeller(
      #   Cluster = setNames(P.d$df$label, P.d$df$Cluster)
      # ),
      ncol = 1, 
      strip.position = "left"
    ) +
    theme(
      # strip.text = element_blank(), 
      strip.background = element_rect(fill = NA, color = NA)
    )
  
  P2 <- temp %>%
    f.plot.mfuzz.heat(
      colors = h.colors, 
      theme = theme
    ) +
    theme(
      strip.text = element_blank(),
    )
  
  P <- P1 | P2 + 
    plot_layout(guides = 'collect') & theme(legend.position = 'top')
  
  return(P)
}


f.plot.fclust <- function(temp_res, mfrow=NULL, ...) {
  
  # if(is.null(mfrow)){
  #   mfrow <- c(2, round(max(temp_res$cluster)/2 + 0.4, digits = 0))
  # }
  
  mfuzz.plot2(
    eset = temp_res$eset,
    cl = temp_res,
    centre = TRUE,
    time.labels = colnames(temp_res$centers),
    x11 = FALSE,
    xlab = "",
    mfrow = mfrow
  )
  
}


f.plot.fclust.a <- function(temp_res, mfrow=NULL, ...){
  
  for(i in names(temp_res)){
    f.plot.fclust( temp_res[[i]], mfrow )
  }
  
}


f.mfuzz_order <- function(temp, order){
  
  if(length(temp) > 2){
    stop("error: only handle matches of 2 results")
  }
  
  order <- setNames(paste0("Cluster ", order), paste0("Cluster ", names(order)))
  
  lev <- names(order)
  
  temp[[2]] <- temp[[2]] %>%
    mutate(
      Cluster = factor(
        Cluster, 
        labels = as.character(order)
      )
    ) %>%
    # pull(Cluster)
    mutate(
      Cluster = factor( Cluster, levels = lev )
    ) #%>% pull(Cluster)
  
  return(temp)
}





