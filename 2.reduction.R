# reduction ----
f.dr <- function(df1, col, method, ...){
  
  obj <- switch (
    method,
    `PCA` = "f.dr_PCA",
    `tSNE` = "f.dr_tSNE",
    `UMAP` = "f.dr_UMAP",
    `PLS` = "f.dr_PLS",
    `OPLS` = "f.dr_OPLS",
    stop("Invalid choice.")
  )
  
  if( length(levels(df1[[col]])) > 2 & method == "OPLS" ){
    stop("Error: OPLS just need two levels")
  }
  
  temp_res <- get(obj)(df1, col, ...)
  
  class(temp_res) <- unique(c("dr", method, class(temp_res)))
  
  return(temp_res)
}


f.dr.p <- function(temp_rd, p.method, core.method=NULL, col, ...){
  data = temp_rd
  
  temp_res <- f.func.p(
    data = data,
    func = f.dr,
    col = col,
    arg = "method",
    p.arg = p.method,
    core = core.method,
    globals = c( f.find_func("f.") ),
    packages = c("ropls", "umap", "Rtsne", "dplyr")
  )
  
  return(temp_res)
}


f.dr.pa <- function(temp_rd, p.method, core.method=NULL, col, ...){
  
  temp_res <- f.func.a(
    data = temp_rd,
    func = f.dr.p,
    col = col,
    p.method = p.method,
    core.method = core.method,
    ...
  )
  
  return(temp_res)
}


f.dr_PCA <- function(df1, col, ...){
  # X y
  X <- df1 %>% dplyr::select(-all_of(col))
  y <- df1 %>% dplyr::select(all_of(col))
  
  set.seed(001)
  temp_res <- ropls::opls( X )
  
  P.d <- list()
  P.d$df <- temp_res@scoreMN %>%
    data.frame() %>%
    cbind( y ) %>%
    dplyr::mutate(
      method = "PCA"
    )
  
  P.d$labs <- temp_res@modelDF %>%
    dplyr::mutate(
      labs = paste0("P", row_number(), " (", sprintf("%.1f", R2X*100), "%)"),
      method = "PCA"
    )
  
  class(P.d) <- c("dr", "PCA", class(P.d))
  # class(temp_res) <- c("dr", "PCA", class(temp_res))
  
  res <- list(
    res = temp_res,
    P.d = P.d
  )
  
  return(res)
}


f.dr_tSNE <- function(df1, col, ...){
  # X y
  X <- df1 %>% dplyr::select(-all_of(col))
  y <- df1 %>% dplyr::select(all_of(col))
  
  set.seed(002)
  temp_res <- Rtsne::Rtsne( X )
  
  P.d <- list()
  P.d$df <- temp_res$Y %>%
    data.frame() %>%
    dplyr::rename( "p1" = X1, "p2" = X2 ) %>%
    data.frame( row.names = rownames(X) )%>%
    cbind( y ) %>%
    dplyr::mutate(
      method = "tSNE"
    )
  
  P.d$labs <- data.frame(
    labs = c(p1 = "tSNE1", p2 = "tSNE2"),
    method = "tSNE"
  )
  
  class(P.d) <- c("dr", "tSNE", class(P.d))
  # class(temp_res) <- c("dr", "tSNE", class(temp_res))
  
  res <- list(
    res = temp_res,
    P.d = P.d
  )
  
  return(res)
}


f.dr_UMAP <- function(df1, col, ...){
  # X y
  X <- df1 %>% dplyr::select(-all_of(col))
  y <- df1 %>% dplyr::select(all_of(col))
  
  set.seed(003)
  temp_res <- umap::umap(
    X,
    method = 'naive',
    # n_neighbors = 10
  )
  
  P.d <- list()
  P.d$df <- temp_res$layout %>%
    data.frame() %>%
    dplyr::rename( "p1" = X1, "p2" = X2 ) %>%
    cbind( y ) %>%
    dplyr::mutate(
      method = "UMAP"
    )
  
  P.d$labs <- data.frame(
    labs = c(p1 = "UMAP1", p2 = "UMAP2"),
    method = "UMAP"
  )
  
  class(P.d) <- c("dr", "UMAP", class(P.d))
  # class(temp_res) <- c("dr", "UMAP", class(temp_res))
  
  res <- list(
    res = temp_res,
    P.d = P.d
  )
  
  return(res)
}


f.dr_PLS <- function(df1, col, ...){
  # X y
  X <- df1 %>% dplyr::select(-all_of(col))
  y <- df1 %>% dplyr::select(all_of(col))
  
  set.seed(004)
  temp_res <- ropls::opls(
    x = X,
    y = y %>% pull(col),
    orthoI = 0
  )
  
  P.d <- list()
  P.d$df <- temp_res@scoreMN %>%
    data.frame() %>%
    cbind( y ) %>%
    dplyr::mutate(
      method = "PLS"
    )
  
  P.d$labs <- temp_res@modelDF %>%
    dplyr::mutate(
      # Latent Variables
      labs = paste0("t", row_number(), " (", sprintf("%.1f", R2X*100), "%)"),
      method = "PLS"
    )
  
  class(P.d) <- c("dr", "PLS", class(P.d))
  # class(temp_res) <- c("dr", "PLS", class(temp_res))
  
  res <- list(
    res = temp_res,
    P.d = P.d
  )
  
  return(res)
}


f.dr_OPLS <- function(df1, col, ...){
  # X y
  X <- df1 %>% dplyr::select(-all_of(col))
  y <- df1 %>% dplyr::select(all_of(col))
  
  set.seed(005)
  temp_res <- ropls::opls(
    x = X,
    y = y %>% pull(col),
    predI = 1,
    orthoI = NA
  )
  
  P.d <- list()
  P.d$df <- data.frame(
    p1 = temp_res@scoreMN,
    p2 = temp_res@orthoScoreMN[, 1]
  ) %>%
    cbind( y ) %>%
    dplyr::mutate(
      method = "OPLS"
    )
  
  P.d$labs <- data.frame(
    labs = c(
      # Predictive Component
      p1 = paste0("t1 (", sprintf("%.1f", temp_res@modelDF["p1", "R2X"]*100), "%)"),
      # Orthogonal Component
      p2 = paste0("O1 (", sprintf("%.1f", temp_res@modelDF["o1", "R2X"]*100), "%)"),
      method = "OPLS"
    )
  )
  
  class(P.d) <- c("dr", "OPLS", class(P.d))
  # class(temp_res) <- c("dr", "OPLS", class(temp_res))
  
  res <- list(
    res = temp_res,
    P.d = P.d
  )
  
  return(res)
}


f.plot.dr.batch <- function(P.d, ratio=NULL, colors.group, colors.batch, shapes.batch, theme="bw", ...){
  # data
  P.d$df <- P.d$df %>%
    arrange(Batch, Group) %>%
    mutate(
      Combine = paste0(Batch, ": ", Group)
    ) %>%
    mutate(
      Combine = factor(Combine, levels = unique(Combine))
    )
  
  P.d$temp <- P.d$df %>%
    dplyr::select(Batch, Group, Combine) %>%
    filter(!duplicated(across(c(Combine)))) %>%
    mutate(
      colors.group = colors.group[as.character(.$Group)],
      shapes.batch = shapes.batch[as.character(.$Batch)]
    )
  
  P.d$colors.group <- P.d$temp %>% { setNames(.$colors.group, .$Combine) }
  P.d$shapes.batch <- P.d$temp %>% { setNames(.$shapes.batch, .$Combine) }
  
  # plot
  P1 <- ggplot(P.d$df) +
    stat_ellipse(
      aes(x = p1, y = p2, color = Batch),
      type = "norm",
      geom ="polygon",
      alpha = 0.05,
      linewidth = 0.2
    ) +
    geom_point(
      aes(x = p1, y = p2, fill = Combine, shape = Combine),
      size = 1, alpha = 1, stroke = 0,
    ) +
    geom_hline(
      yintercept = 0, color = "gray50",
      linetype = "longdash", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0, color = "gray50",
      linetype = "longdash", linewidth = 0.5
    ) +
    scale_fill_manual(
      values = P.d$colors.group
    ) +
    scale_color_manual(
      values = colors.batch
    ) +
    scale_shape_manual(
      values = P.d$shapes.batch
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(size = 12, hjust = 0),
      # legend.position = "top",
      # legend.direction = c("horizontal", "vertical")[2],
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(size = 2),
        ncol = length(unique(P.d$df$Batch)), byrow = F, #reverse = T
      )
    ) +
    labs(
      x = P.d$labs$labs[1],
      y = P.d$labs$labs[2]
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.dr <- function(P.d, ratio=NULL, colors, shapes, theme="bw", ...){
  # data
  
  # plot
  P1 <- ggplot(P.d$df) +
    facet_wrap(
      .~ method,
      nrow = 1,
      scales = "free_x",
      strip.position = "top"
    ) +
    geom_point(
      aes(x = p1, y = p2, color = Group, shape = Group),
      size = 0.5, alpha = 1
    ) +
    stat_ellipse(
      aes(x = p1, y = p2, color = Group, fill = Group),
      type = "norm",
      geom ="polygon",
      alpha = 0.1,
      linewidth = 0.3
    ) +
    geom_hline(
      yintercept = 0, color = "gray50",
      linetype = "longdash", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0, color = "gray50",
      linetype = "longdash", linewidth = 0.5
    ) +
    scale_color_manual(
      values = colors
    ) +
    scale_fill_manual(
      values = colors
    ) +
    scale_shape_manual(
      values = shapes,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(size = 12, hjust = 0),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    guides(
      fill = guide_none(),
      shape = guide_legend(override.aes = list(size = 2))
    ) +
    labs(
      x = P.d$labs$labs[1],
      y = P.d$labs$labs[2]
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.dr.3D <- function(P.d, color, ...){
  # # data
  # P.d <- list()
  P.d <- temp$PLS
  
  color <- color[match(unique(P.d$df$Group), names(color))]
  
  # library(rgl)
  # open3d()
  # plot3d(
  #   x = P.d$df$p1, 
  #   y = P.d$df$p2, 
  #   z = P.d$df$p3,
  #   xlab = "PC1", 
  #   ylab = "PC2", 
  #   zlab = "PC3",
  #   col = color,
  #   size = 5
  # )
  # legend3d(
  #   "topright", 
  #   legend = unique(P.d$df$Group),
  #   col = color,
  #   pch = 16,
  #   cex = 0.8
  # )
  # 
  # library(plotly)
  # 
  # P1 <- plot_ly(
  #   data = P.d$df,
  #   x = ~p1,
  #   y = ~p2,
  #   z = ~p3,
  #   type = "scatter3d",  # 明确指定3D散点类型
  #   mode = "markers",    # 明确指定显示点
  #   color = ~Group,  # 按分组着色
  #   colors = color,  # 使用你定义的颜色向量
  #   marker = list(size = 1.3)  # 点大小
  # ) %>%
  #   layout(
  #     scene = list(
  #       xaxis = list(title = P.d$labs$labs[1]),
  #       yaxis = list(title = P.d$labs$labs[2]),
  #       zaxis = list(title = P.d$labs$labs[3])
  #     ),
  #     title = unique(P.d$df$method)
  #   )
  # 
  # print(P1)
  # 
  # htmlwidgets::saveWidget(ggplotly(), "pca_3d_plotly.html")
  
  
  # library(plot3D)
  # 
  # plot3D::scatter3D(
  #   x = P.d$df$p1,
  #   y = P.d$df$p2,
  #   z = P.d$df$p3,
  #   xlab = P.d$labs$labs[1],
  #   ylab = P.d$labs$labs[2],
  #   zlab = P.d$labs$labs[3],
  #   pch = 21,
  #   cex = 1,
  #   col = NA,        
  #   bg = color,
  #   ticktype = c("simple", "detailed")[1],
  #   bty = "b2",
  #   box = T,
  #   theta = 45, 
  #   phi = 35,
  #   d = 10,
  #   colkey = F
  # )
  
  
}