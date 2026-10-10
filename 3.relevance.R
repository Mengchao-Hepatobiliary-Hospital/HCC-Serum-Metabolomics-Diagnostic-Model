# relevance ----
f.rcorr.res <- function(temp_res, names=NULL, diag=T, ...) {
  
  cormat <- temp_res$r
  pmat <- temp_res$P
  
  ut <- upper.tri(cormat, diag = diag)
  res <- data.frame(
    var1 = rownames(cormat)[row(cormat)[ut]], 
    var2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut], 
    p = pmat[ut],
    method = temp_res$type
  )
  
  return(res)
}


f.rcorr <- function(x, y=NULL, type="pearson", p.adj=T, diag=T, ...){
  
  temp_res <- Hmisc::rcorr(
    x = x, 
    y = y, 
    type = type
  )

  temp <- temp_res %>% 
    f.rcorr.res(
      diag = diag
    )
  
  if(p.adj){
    temp <- temp %>%
      mutate(
        p.signif = f.add_star(p),
        p.adj = p.adjust(p, method = "BH"),
        p.adj.signif = f.add_star(p.adj) 
      )
  }else{
    temp <- temp %>%
      mutate(
        p.signif = f.add_star(p)
      )
  }
  
  return(temp)
}


f.plot.cor.ggheat <- function(temp_res, x, y, point_range=NULL, pos="lt", signif, colors, theme="mini", ...){
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      var1 = !!sym(x),
      var2 = !!sym(y),
      "signif" = !!sym(signif)
    )
  
  if(!is.null(point_range)){
    P.d$df <- P.d$df %>%
      mutate(
        cor_abs = abs(cor),
        size = scales::rescale(cor_abs, to = point_range)
      )
  }
  
  if(grepl("l", pos)){
    lev1 <- P.d$df %>% count(var1) %>%
      arrange(desc(n)) %>% pull(var1)
  } else {
    lev1 <- P.d$df %>% count(var1) %>%
      arrange(n) %>% pull(var1)
  }
  
  if(grepl("t", pos)){
    lev2 <- P.d$df %>% count(var2) %>%
      arrange(n) %>% pull(var2)
  } else {
    lev2 <- P.d$df %>% count(var2) %>%
      arrange(desc(n)) %>% pull(var2)
  }
  
  P.d$df$var1 <- factor(P.d$df$var1, levels = lev1)
  P.d$df$var2 <- factor(P.d$df$var2, levels = lev2)
  
  if(!is.null(point_range)){
    P1 <- ggplot(P.d$df, aes(x = var1, y = var2, fill = cor)) +
      geom_tile(
        # color = "white",
        color = "gray80", fill = "white",
        size = 0.3, 
      ) +
      geom_point(
        aes(fill = cor, size = size), 
        shape = 22, 
        color = "black",
        # show.legend = F
      )
  }else{
    P1 <- ggplot(P.d$df, aes(x = var1, y = var2, fill = cor)) +
      geom_tile(
        # color = "white",
        color = "#f7f7f7", 
        size = 0.3, 
      )
  }
  
  P1 <- P1 +
    geom_text(
      aes(label = signif), 
      color = "black", 
      size = 3
    ) + 
    scale_fill_gradientn(
      name = "Correlation",
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      colors = colors,
      oob = scales::censor
    ) +
    f.theme(theme) +
    theme(
      panel.grid = element_blank(),
      #plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    ) +
    coord_fixed()
  
  if (grepl("r", pos)) {
    P1 <- P1 + scale_y_discrete(position = "right")
  }
  
  if (grepl("t", pos)) {
    P1 <- P1 + scale_x_discrete(position = "top") +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 10, color = "black"),
      )
  }
  
  # print(P1)
  
  return(P1)
}


# temp_res <- temp_rd %>%
#   as.matrix() %>%
#   f.rcorr(
#     y = NULL,
#     type = c("pearson", "spearman")[2],
#     p.adj = T,
#     diag = F
#   )
# 
# # plot
# P1 <- temp_res %>%
#   f.plot.cor.ggheat(
#     x = "var1",
#     y = "var2",
#     signif = c("p.signif", "p.adj.signif")[2],
#     pos = c("lt", "lb", "rt", "rb")[2],
#     colors = c("#1f78b4", "#f7f7f7", "#ff7f00"),
#     theme = "mini",
#     # point_range = c(0, 0.5)
#   )
# 
# print(P1)

f.plot.corrplot <- function(temp_res, text.size=0.6, d.method="ellipse", d.colors, u.size, u.method="number", u.colors, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>% as.matrix()
  
  # plot
  corrplot(
    P.d$df, 
    # method = c("pie", "shade", "color", "square", "ellipse")[5],
    method = d.method,
    type = "lower",
    tl.pos = "n",
    mar = c(1, 1, 1, 1),
    col = d.colors,
    cl.cex = 0.6,         
    cl.ratio = 0.1        
  )
  corrplot(
    P.d$df, 
    method = u.method, 
    # number.cex = 0.2,
    number.cex = u.size,
    type = "upper",      
    tl.pos = "n",        
    diag = FALSE,        
    mar = c(1, 1, 1, 1),
    add = T,             
    cl.pos = "n",
    col = u.colors
  )
  mtext(
    paste0("Range: ", sprintf("%.2f", min(P.d$df[upper.tri(P.d$df, diag = F)])), " ~ ", sprintf("%.2f", max(P.d$df[upper.tri(P.d$df, diag = F)]))), 
    side = 3, line = 2.7, cex = text.size, col = "firebrick", #font = 2
  )
  
}

# temp_res <- temp_x$train %>%
#   cor(
#     method = c("pearson", "kendall", "spearman")[3],
#     use = "complete.obs"
#   )
# 
# temp_func <- function(){
#   temp_res %>%
#     f.plot.corrplot(
#       text.size = 0.7,
#       u.size = 0.6,
#       u.method = "number",
#       u.colors = "black",
#       d.method = c("pie", "shade", "color", "square", "ellipse")[5],
#       d.colors = colorRampPalette(c("dodgerblue", "white", "firebrick1"))(200)
#     )
# }
# 
# pdf(paste0(pf, "/cor=train.pdf"), width = 3.3, height = 3.5)
# temp_func()
# dev.off()


f.plot.cor.ggpair <- function(temp_rd, type="pearson", colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_rd
  
  if(is.null(names(colors))){
    names(colors) <- names(P.d$df)
  }
  
  P.d$lower <- function(data, mapping, ...){
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    
    cor_test <- cor.test(x, y, method = type, na.action = na.omit)
    
    r_squared <- sprintf("%.3f", cor_test$estimate^2)
    p_value <- round(cor_test$p.value, 4)
    
    ggplot(data, mapping) +
      geom_point(color = "grey", size = 0.8) +
      geom_smooth(
        method = "lm", color = "red", se = FALSE,
        linetype = "dashed", linewidth = 0.3
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0(
          "R² = ", r_squared, "\n", 
          ifelse(
            p_value < 0.001, "P < 0.001",
            paste0("P = ", sprintf("%.3f", p_value))
          )
        ),
        size = 2.8,
        hjust = 1.1, vjust = 1.3
      )
  }
  
  P.d$diag <- function(data, mapping, ...){
    var_name <- rlang::as_name(mapping$x)
    current_data <- eval_data_col(data, mapping$x)
    mu <- mean(current_data, na.rm = TRUE)
    sigma <- sd(current_data, na.rm = TRUE)
    
    ggplot(data, mapping) +
      geom_histogram(
        aes(y = after_stat(density)),
        color = "white",
        fill = colors[var_name],
        bins = 15
      ) +
      stat_function(
        fun = dnorm,
        args = list(mean = mu, sd = sigma),
        color = "black", size = 0.5, linewidth = 0.3
      )
  }
  
  P.d$upper <- function(data, mapping, ...){
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    
    cor_test <- cor.test(x, y, method = type, na.action = na.omit)
    
    r_squared <- sprintf("%.3f", cor_test$estimate^2)
    p_value <- round(cor_test$p.value, 4)
    
    ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0(
          "R² = ", r_squared, "\n", 
          ifelse(
            p_value < 0.001, "P < 0.001",
            paste0("P = ", sprintf("%.3f", p_value))
          )
        ),
        size = 2.8
      )
  }
  
  P1 <- ggpairs(
    P.d$df, 
    # columns = 1:4,
    lower = list(continuous = P.d$lower ),
    diag = list(continuous = P.d$diag ),
    # upper = list(continuous = P.d$upper ),
    upper = list(continuous = "blank"),
    switch = "both",
    progress = F
  ) +
    f.theme(theme) +
    theme(
      # aspect.ratio = ratio,
      panel.grid = element_blank(),
      strip.background = element_blank(),  
      strip.placement = "outside",  
      strip.text = element_text(size = 10, color = "black")
    )
  
  P1 <- ggplotify::as.ggplot(ggmatrix_gtable(P1))
  
  # print(P1)
  
  return(P1)
}


f.plot.chord <- function(temp_res, degree=90, gap=30, cex=0.8, label=NULL, col.inner, col.outer, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  P.d$col.inner <- col.inner
  
  if(!is.null(names(col.outer))){
    P.d$col.outer <- c(
      rep(col.outer[1],  length(unique(P.d$df[[names(col.outer)[1]]]))),
      rep(col.outer[2],    length(unique(P.d$df[[names(col.outer)[2]]])))
    ) %>% as.character()
  }else{
    P.d$col.outer <- c(
      rep(col.outer[1],  length(unique(P.d$df[[1]]))),
      rep(col.outer[2],  length(unique(P.d$df[[2]])))
    ) %>% as.character()
  }

  circos.clear()
  circos.par(start.degree = degree)
  # plot
  chordDiagram(
    x = P.d$df,
    big.gap = gap,
    col = P.d$col.inner,
    annotationTrack = c("grid", "axis", "name")[1],
    transparency = 0.4,
    preAllocateTracks = list(
      track.height = 0.2
    ),
    grid.col = P.d$col.outer
  )
  
  # text
  circos.trackPlotRegion(
    track.index = 2, 
    panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      circos.text(
        mean(xlim), max(ylim)+0.5, sector.name,
        facing = "clockwise", niceFacing = TRUE, 
        adj = c(0, 0.5), cex = cex,
        col = if_else(sector.name %in% label, "firebrick", "black")
      )
    }, 
    bg.border = NA
  )
  
  P1 <- Legend(
    title = "Type", 
    at = names(P.d$df)[1:2], 
    legend_gp = gpar(
      fill = unique(P.d$col.outer)
    )
  )
  
  P2 <- Legend(
    col_fun = P.d$col.inner, 
    title = names(P.d$df)[3],
    at = c(-1, -0.5, 0, 0.5, 1),
    direction = c("horizontal", "vertical")[2], 
    title_position = "topleft",
    # title_gp = gpar(
    #   fontsize = 12*0.8,
    #   fontface = c("bold", "plain")[2]
    # ),
    # labels_gp = gpar(
    #   fontsize = 12*0.8,
    #   fontface = c("bold", "plain")[2]
    # )
  )
  
  P <- packLegend(
    P1, 
    P2,
    max_height = unit(0.9, "npc")
  )
  
  draw(
    P, 
    x = unit(0.9, "npc"),
    y = unit(0.5, "npc"),
    just = c("right", "bottom")
  )
  
}

