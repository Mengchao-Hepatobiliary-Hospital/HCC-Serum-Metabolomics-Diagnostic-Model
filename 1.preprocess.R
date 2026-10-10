# preprocess ----
f.calNA <- function(df1, split=NULL, ...){

  temp <- names(df1)[!sapply(df1, is.numeric)]
  
  if(!is.null(split)){
    if(length(split) == length(temp) && all(sort(split) == sort(temp))){
      
      temp_res <- df1 %>% 
        pivot_longer(
          col = -all_of(split), 
          names_to = "column", 
          values_to = "value"
        ) %>%
        group_by(across(all_of(split)), column) %>%
        reframe(
          n_miss = sum(is.na(value)),
          n_count = sum(!is.na(value)),
          n_total = length(value),
          miss_percent = (n_miss / n_total) * 100,
          miss_paste = paste0(n_miss, "/", n_total, " (", sprintf("%.2f%%", miss_percent), ")"),
          count_percent = (n_count / n_total) * 100,
          count_paste = paste0(n_count, "/", n_total, " (", sprintf("%.2f%%", count_percent), ")") 
        )
      
    }else{
      stop("split value must correspond with the character column.")
    }
  }else{
    
    if(length(temp) != 0){
      stop("when the split value is NULL, character columns are not allowed.")
    }else{
      
      temp_res <- df1 %>% 
        pivot_longer(
          col = everything(), 
          names_to = "column", 
          values_to = "value"
        ) %>%
        group_by(column) %>%
        reframe(
          n_miss = sum(is.na(value)),
          n_count = sum(!is.na(value)),
          n_total = length(value),
          miss_percent = (n_miss / n_total) * 100,
          miss_paste = paste0(n_miss, "/", n_total, " (", sprintf("%.2f%%", miss_percent), ")"),
          count_percent = (n_count / n_total) * 100,
          count_paste = paste0(n_count, "/", n_total, " (", sprintf("%.2f%%", count_percent), ")") 
        )
      
    }
  }
  
  return(temp_res)
}


f.calNA_byCol <- function(df1, split=NULL, name, ...){
  
  temp_res <- df1 %>%
    f.calNA( split ) %>%
    dplyr::rename(all_of(setNames("column", name)))
  
  return(temp_res)
}


f.calNA_byRow <- function(df1, split=NULL, name, ...){
  
  temp_res <- df1 %>%
    dplyr::select(-all_of(split)) %>%
    t() %>% as.data.frame() %>%
    f.calNA( split=NULL ) %>%
    left_join(
      df1 %>% 
        dplyr::select(all_of(split)) %>%
        tibble::rownames_to_column("column"),
      by = "column"
    ) %>%
    dplyr::rename(all_of(setNames("column", name)))
  
  return(temp_res)
}


f.plot.NA <- function(temp_res, ...){
  # data
  P.d <- list()
  P.d$df <- temp
  
  # plot
  ggplot(P.d$df, aes(x = Sample, y = Feature, fill = Value)) +
    facet_row(
      ~ Group,
      scales = "free_x",
      space = "free"
    ) +
    geom_tile(
      color = "white", 
      linewidth = 0.2
    ) +
    scale_fill_manual(
      labels = c("Present", "Missing"),
      values = c(
        "0" = "white",
        "1" = "red"
      )
    ) +
    theme_minimal() +
    theme(
      # aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank(),
    )
  
  
  
}



f.plot.Count_hist <- function(temp_res, cutoff=0.5, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot(data = P.d$df, aes(x = miss_percent)) +
    geom_histogram(
      # bins = 25, 
      binwidth = 5,
      fill = "#69b3a2", 
      color = "white", 
      alpha = 0.9
    ) +
    geom_vline(
      xintercept = cutoff*100,
      color = "firebrick",
      linetype = "twodash", linewidth = 0.4
    ) +
    stat_bin(
      geom = "text",
      aes(label = after_stat(count), y = (max(after_stat(count))+min(after_stat(count)))/12.5),
      breaks = c(0, cutoff, 1)*100,
      # bins = 2, 
      # binwidth = 50,
      vjust = -0.5, 
      size = 3.5,
      color = "firebrick"
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio
    ) +
    labs(
      x = "Missing proportion (%)",
      y = "Count"
    ) +
    scale_x_continuous(
      expand = c(0.05, 0, 0.05, 0),
      # limits = c(0, 100),
      breaks = seq(0, 100, by = 25)
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.05, 0)
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.Count_dot <- function(temp_res, x="Group", color="Group", colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      "x" = !!sym(x),
      "color" = !!sym(color)
    )
  
  # plot
  P1 <- ggplot( P.d$df, aes(x = x, y = n_count, fill = color, color = color )) +
    geom_dotplot(
      binaxis = 'y',
      stackdir = 'center',
      #stackratio = 1,
      dotsize = 0.5,
    ) +
    stat_summary(
      geom = "text",
      fun.data = function(y) {
        data.frame(
          y = max(P.d$df$n_count, na.rm = TRUE) + diff(range(na.omit(P.d$df$n_count)))*0.05,
          label = sprintf("Median: %.0f", median(y, na.rm = TRUE) )
        )
      },
      # fun = median,
      # aes(label = sprintf("Median: %.0f", after_stat(y)) ),
      # position = position_nudge(
      #   y = max(P.d$df$n_count)/2
      # ),
      # y = max(P.d$df$n_count),
      vjust = 0,
      hjust = 0.5,
      size = 3.5,
      # color = "black"
    ) +
    stat_summary(
      fun.data = mean_sdl,
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.15, linewidth = 0.4, color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.5, linewidth = 0.3, color = "black",
    ) +
    scale_fill_manual(
      values = colors,
    ) +
    scale_color_manual(
      values = colors,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none",
    ) +
    labs(
      y = "Count"
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.Int_hist <- function(temp_res, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = log2(Intensity))) +
    geom_histogram(
      aes( y = after_stat(density)),
      bins = 20,
      fill = "#69b3a2", 
      color = "white", 
      alpha = 0.9
    ) +
    geom_density(
      colour = "black",
      fill = NA,
      # alpha = 1,
      linewidth = 0.3,
      linetype = "dashed",
      key_glyph = c("smooth", "rect")[2],
      # show.legend = FALSE
    ) +
    stat_bin(
      bins = 20,
      geom = "text",
      aes(y = after_stat(density), label = after_stat(count)),
      vjust = -0.5, 
      size = 3, 
      color = "firebrick"
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.1, 0)
    ) +
    labs(
      x = parse(text = paste0("Log[2] ~ (", P.d$df$stat[1], " ~ \"of Intensity\")")),
      y = "Density"
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.Int_box <- function(temp_res, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Sample, y = log2(Value), color = Group, fill = Group)) +
    geom_boxplot(
      fill = NA,
      outlier.shape = NA,
      linewidth = 0.15,
      box.linetype = "solid",
      whisker.linetype = "dashed",
      key_glyph = "rect",
      # show.legend = FALSE,
    ) +
    scale_color_manual(
      values = colors
    ) +
    scale_fill_manual(
      values = colors
    ) +
    guides(
      fill = guide_legend(override.aes = list(
        color = colors[names(colors) %in% levels(P.d$df$Group)],
        fill = colors[names(colors) %in% levels(P.d$df$Group)]
      ))
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    labs(
      x = "Samples",
      y = expression(paste(Log[2], " (Intensity)"))
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.CV_vln <- function(temp_res, x="Group", color="Group", ratio=NULL,colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      "x" = !!sym(x),
      "color" = !!sym(color)
    )
  
  # plot
  P1 <- ggplot( P.d$df, aes(x = x, y = CV, fill = color, color = color )) +
    geom_violin(
      trim = F,
      alpha = 0.6,
      # fill = "white",
      linewidth = 0.3
    ) +
    geom_boxplot(
      width = 0.10,
      fill = "white",
      alpha = 1,
      size = 0.3,
      # coef = 0,
      # staplewidth = 0.5,
      outlier.shape = NA
    ) +
    stat_summary(
      geom = "text",
      fun.data = function(y) {
        data.frame(
          y = max(P.d$df$CV, na.rm = TRUE),# + diff(range(na.omit(P.d$df$CV)))*0.05,
          label = sprintf("Median: %.0f %%", median(y, na.rm = TRUE) )
        )
      },
      # fun = median,
      # aes(label = sprintf("Median: %.0f", after_stat(y)) ),
      # position = position_nudge(
      #   y = max(P.d$df$n_count)/2
      # ),
      # y = max(P.d$df$n_count),
      # vjust = 0,
      # hjust = 0.5,
      size = 3.5,
      # color = "black"
    ) +
    scale_fill_manual(
      values = colors,
    ) +
    scale_color_manual(
      values = colors,
    ) +
    # f.theme_bw() + 
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none",
    ) +
    labs(
      y = "CVs (%)"
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.Abund <- function(temp_res, label="Top10", ratio=NULL, colors, theme="bw", p.text=0.2, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # label
  if(length(label) == 1 & any(grepl("Top", label))){  # default Top10
    
    P.d$df <- P.d$df %>%
      group_by(Group) %>%
      mutate(
        label = if_else(Rank <= as.numeric(stringr::str_extract(label, "\\d+")), "Highlight", "Normal")
      ) %>%
      ungroup()
  
  }else{
    
    P.d$df <- P.d$df %>%
      group_by(Group) %>%
      mutate(
        label = if_else(Feature %in% label, "Highlight", "Normal")
      ) %>%
      ungroup()
    
  }
  
  # plot
  P1 <- ggplot(P.d$df) +
    facet_wrap(
      .~ Group,
      nrow = 1,
      #scales = "free_y",
      strip.position = "top"
    ) +
    geom_point(
      data = subset(P.d$df, label == "Normal"),
      aes(x = Rank, y = log10(Abundance)),
      color = "gray", size = 0.5,
    ) +
    geom_point(
      data = subset(P.d$df, label == "Highlight"),
      aes(x = Rank, y = log10(Abundance), color = Group),
      # color = "firebrick",
      size = 0.5,
    ) +
    scale_color_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      legend.position = "none",
    ) +
    labs(
      # y = expression(paste(Log[10], " (Abundance)" )),
      y = parse(text = paste0("Log[10] ~ (", P.d$df$stat[1], " ~ \"of Intensity\")")),
      x = "Rank"
    )
  
  # print(P1)
  
  
  if(length(label) == 1 & any(grepl("Top", label))){ # default Top10
    
    P.d$x.range <- diff(ggplot_build(P1)@layout$panel_params[[1]]$x.range)
    
    P.d$label <- subset(P.d$df, label == "Highlight") %>%
      mutate(
        log_Abund = log10(Abundance),
        global_ymin = min(log_Abund, na.rm = TRUE),
        global_ymax = max(log_Abund, na.rm = TRUE),
        global_diff = diff(range(log_Abund, na.rm = TRUE))
      ) %>%
      group_by(Group) %>%
      reframe(
        Feature = Feature,
        Rank = Rank,
        `Cum.%` = `Cum.%`,
        x = max(Rank) + P.d$x.range*0.005,
        xend = max(Rank) + P.d$x.range*0.05,
        ymin = min(log_Abund, na.rm = TRUE),
        ymax = max(log_Abund, na.rm = TRUE),
        y = seq(unique(ymax), unique(ymin), length.out = n()),
        yend = seq(unique(global_ymax)+unique(global_diff)/2, unique(global_ymin)-unique(global_diff)/2, length.out = n())
      ) %>%
      ungroup() %>%
      dplyr::select(-ymin, -ymax)
    
    P.d$bracket <- P.d$label %>%
      group_by(Group) %>%
      reframe(
        x1 = max(xend) + P.d$x.range*0.01,
        x3 = max(xend) + P.d$x.range*0.02,
        ymin = min(yend),
        ymax = max(yend),
        `Cum.%` = paste0(sprintf("% .1f%%", max(`Cum.%`)), " total  abundance")
      ) %>%
      group_by(Group) %>%
      reframe(
        x = c(x1, x1, x3) + P.d$x.range*p.text,
        xend = c(x3, x3, x3) + P.d$x.range*p.text,
        y = c(ymin, ymax, ymin),
        yend = c(ymin, ymax, ymax),
        y_cum = (ymax + ymin)/2,
        label_cum = `Cum.%`
      ) %>%
      ungroup()
    
    P1 <- P1 + 
      geom_segment(
        data = P.d$label,
        aes(
          x = x, xend = xend,
          y = y, yend = yend,
          color = Group
        ),
        linewidth = 0.35,
        # color = "firebrick",
        # linetype = "dashed",
      ) +
      geom_text(
        data = P.d$label,
        aes(
          x = xend + P.d$x.range*0.01, 
          y = yend, 
          label = Feature,
          color = Group
        ), 
        hjust = 0,
        size = 3.3
      ) +
      geom_segment(
        data = P.d$bracket,
        aes(
          x = x, xend = xend,
          y = y, yend = yend,
          color = Group
        ),
        linewidth = 0.35,
        # color = "firebrick"
      ) +
      geom_text(
        data = P.d$bracket %>% group_by(Group) %>% slice_head(n = 1),
        aes(
          x = xend + P.d$x.range*0.01, 
          y = y_cum, 
          label = label_cum,
          color = Group
        ), 
        hjust = 0,
        size = 3.3
      )
    
  }else{
    
    P1 <- P1 +
      geom_label_repel(
        data = subset(P.d$df, label == "Highlight"),
        aes(label = Feature, x = Rank, y = log10(Abundance), color = Group),
        nudge_x = 2,
        nudge_y = 1,
        alpha = 0.7,
        size = 2,
        min.segment.length = 0,
        direction = "both",
        segment.size = 0.3,
        # segment.linetype = 6,
        segment.color = "black",
      )
    
  }
  
  return(P1)
}


f.plot.linecor <- function(temp_res, x, y, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      x = !!sym(x),
      y = !!sym(y)
    )
  
  # plot
  P1 <- ggplot( P.d$df, aes(x = log2(x), y = log2(y)) ) +
    geom_point(
      size = 1, color = "red"
    ) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = T, 
      linewidth = 0.6,
      color = "blue"
    ) +
    stat_cor(
      method = "pearson",
      r.accuracy = 0.001, p.accuracy = 0.001,
      # label.x = 3, label.y = 0.5
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      panel.grid = element_blank(),
      legend.position = "none",
    )
    # labs(
    #   x = expression(paste(Log[2], " (Urine)")),
    #   y = expression(paste(Log[2], " (Serum)"))
    # )
  
  # print(P1)
  
  return(P1)
}



f.plot.nTier_donut <- function(temp_res, colors, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggdonut(
    data = P.d$df, group_key = "Tier", count_type = c("count", "full")[1],
    fill_color = colors,
    border_color = NA, border_size = 0,
    label_info = "all", label_color = "gray20", label_size = 3, 
    # donut.label.color = "", 
    donut.label.size = 3.5,
    # label_pos = "out", label_type = "horizon",
    label_pos = "in", label_type = "circle",
  ) + 
    theme(
      legend.position = "right",
      legend.title = element_blank(),
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.Count_box <- function(temp_res, width=0.3, x="Group", color="Group", colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      "x" = !!sym(x),
      "color" = !!sym(color)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Group, y = n_count, fill = color, color = color)) +
    stat_boxplot(
      geom = "errorbar", 
      width = 0.15, 
      color = "black", 
      linewidth = 0.55
    ) +
    geom_boxplot(
      outliers = F, 
      width = width, 
      color = "black", 
      fill = "white"
    ) +
    geom_jitter(
      width = 0.2,
      height = 0,
      shape = 21, 
      size = 2, 
      stroke = 0.01, 
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = colors,
    ) +
    scale_color_manual(
      values = colors,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none", 
      plot.title = element_text(size = 12, hjust = 0.5), 
    ) +
    labs(
      y = "Count"
    ) +
    stat_summary(
      geom = "text",
      fun.data = function(y) {
        data.frame(
          y = max(P.d$df$n_count, na.rm = TRUE) + diff(range(na.omit(P.d$df$n_count)))*0.05,
          label = sprintf("Median: %.0f", median(y, na.rm = TRUE) )
        )
      },
      # fun = median,
      # aes(label = sprintf("Median: %.0f", after_stat(y)) ),
      # position = position_nudge(
      #   y = max(P.d$df$n_count)/2
      # ),
      # y = max(P.d$df$n_count),
      vjust = 0,
      hjust = 0.5,
      size = 3.5,
      # color = "black"
    )
  
  # print(P1)
  
  return(P1)
}

f.plot.class_pie <- function(temp_res, main, width = NULL, height = NULL, ...){
  
  P1 <- plot_ly(
    labels = temp_res$category,
    values = temp_res$count,
    type = "pie",
    direction = c("counterclockwise", "clockwise")[2],
    textinfo = "label+percent",
    insidetextorientation = "radial",
    hoverinfo = "text",
    marker = list(
      # colors = f.color()$Pokemon,
      line = list(color = "#FFFFFF", width = 2)
    ),
    # hole = 0.1,  # 创建环形图效果，增强3D视觉感
    rotation = 90,  # 初始旋转角度
  ) %>%
    layout(
      title = list(
        text = main
      ),
      showlegend = FALSE,
      font = list(
        family = "Arial", 
        size = 11,
        color = "black" 
      )
    )
  
  if(!is.null(width)){
    P1 <- P1 %>% layout( width = width, height = height )
  }
  
  return(P1)
}


f.plot.hist <- function(temp_res, x, bins, theme, ratio, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      x = !!sym(x)
    )
  
  P.d$label <- list(
    median = median(P.d$df$x),
    label = paste0("Median: ", sprintf("%.1f", median(P.d$df$x)) )
  )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = x)) +
    geom_density(
      color = "#ea4335",
      linewidth = 0.7,
      alpha = 0.3,
      fill = "#ea4335"
    ) + 
    geom_histogram(
      aes(y = after_stat(density)), 
      fill = "gray80",
      color = NA,
      alpha = 0.3,
      bins = bins,
      # binwidth = 0.9
    ) + 
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      strip.text = element_text(size = 10, color = "black"),
      strip.background.y = element_rect(fill = NA, color = NA),
      legend.position = "none"
    ) +
    labs(
      # x = "Count",
      y = "Density"
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.1, 0),
      limits = c(0, NA)
    ) +
    geom_vline(
      xintercept = P.d$label$median,
      color = "#fbbc05",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = P.d$label$median,
      y = Inf,
      label = P.d$label$label,
      size = 3.3,
      color = "firebrick",
      hjust = 0.5,
      vjust = 1.3
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.ridge <- function(temp_res, x, y, colors, theme="bw", ratio, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      x = !!sym(x),
      y = !!sym(y)
    )
  
  P.d$label <- P.d$df %>%
    group_by(y) %>%
    summarise(
      median = median(x, na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      label = paste0("Median: ", sprintf("%.1f", median) )
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = x, y = y, fill = y, color = y)) +
    geom_density_ridges(
      # stat="binline"
      scale = 4,
      alpha = 0.3,
      # jittered_points = TRUE,
      # quantile_lines = TRUE,
      # quantiles = 2,
    ) + 
    scale_fill_manual(
      values = colors
    ) +
    scale_color_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      
      axis.title.y = element_blank(),
      # axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      
      legend.title = element_blank(),
      # legend.position = "inside",
      # legend.position.inside = c(0.03, 0.97),  
      # legend.justification = c(0, 1), 
      # legend.background = element_rect(colour = NA, fill = NA, linewidth = 0.01),
      # legend.key.size = unit(0.4, "cm")
    ) +
    scale_y_discrete(
      expand = c(0.1, 0, 0.1, 0)
    ) +
    geom_vline(
      data = P.d$label,
      aes(xintercept = median),
      color = "#fbbc05",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    geom_text(
      data = P.d$label,
      aes(x = median, label = label),
      y = Inf,
      size = 3.3,
      color = "firebrick",
      hjust = 0.5,
      vjust = 1.3
    )
  
  # print(P1)
  
  return(P1)
}
