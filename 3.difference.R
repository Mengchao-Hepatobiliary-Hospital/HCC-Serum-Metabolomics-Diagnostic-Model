## difference ----
f.add_star <- function(p_values) {
  case_when(
    p_values >= 0.05 ~ "ns",
    p_values < 0.05 & p_values >= 0.01 ~ "*",
    p_values < 0.01 & p_values >= 0.001 ~ "**",
    p_values < 0.001 ~ "***",
    TRUE ~ NA_character_  
  )
}


f.test_assumptions <- function(data, formula){
  
  value <- f.splitform(formula)$left
  group <- f.splitform(formula)$right
  
  data %>%
    group_modify(function(.x, .y) {
      shapiro_res <- .x %>%
        group_by(!!sym(group)) %>%
        summarise(
          Stat = shapiro.test(!!sym(value))$statistic,
          P_val = shapiro.test(!!sym(value))$p.value,
          .groups = "drop"
        )
      
      shapiro_stars <- f.add_star(shapiro_res$P_val)
      
      lev_test <- leveneTest(formula, data = .x)
      levene_res <- tibble(
        Levene_F = as.numeric(lev_test[1, "F value"]),
        Levene_p = as.numeric(lev_test[1, "Pr(>F)"])
      )
      
      levene_star <- f.add_star(levene_res$Levene_p)
      
      tibble(
        Shapiro_Stats = paste(sprintf("%.4f", shapiro_res$Stat), collapse = " | "),
        Shapiro_Ps = paste(sprintf("%.4f", shapiro_res$P_val), collapse = " | "),
        Shapiro_Sig = paste(shapiro_stars, collapse = " | "),
        Groups = paste(unique(.x[[group]]), collapse = " | "),
        levene_res,
        Levene_Sig = levene_star
      )
      # shapiro_res %>%
      #   mutate(
      #     Levene_F = levene_res$Levene_F,
      #     Levene_p = levene_res$Levene_p,
      #     Groups_Text = paste(unique(.x[[group]]), collapse = " | ")
      #   )
      
    }) %>%
    ungroup()
}


f.test_wilcox <- function(df1, formula, comparisons=NULL, ref.group=NULL, p.adjust.method="holm", paired=F, ...){
  
  if(comparisons == "all"){
    temp_res <- df1 %>% 
      pairwise_wilcox_test(
        formula = formula,
        p.adjust.method = p.adjust.method
      )
  }else{
    temp_res <- df1 %>%
      wilcox_test(
        formula = formula,
        comparisons = comparisons,
        ref.group = ref.group,
        p.adjust.method = p.adjust.method,
        paired = paired
      )
  }
  
  return(temp_res)
}

f.test_t <- function(df1, formula, comparisons=NULL, ref.group=NULL, p.adjust.method="holm", paired=F, ...){
  
  if(comparisons == "all"){
    temp_res <- df1 %>% pairwise_t_test( formula = formula )
  }else{
    temp_res <- df1 %>%
      t_test(
        formula = formula,
        comparisons = comparisons,
        ref.group = ref.group,
        p.adjust.method = p.adjust.method,
        paired = paired
      )
  }
  
  return(temp_res)
}


f.test_tukey <- function(df1, formula, ...){
  
  temp_res <- df1 %>% tukey_hsd(formula)
  
  temp_res <- temp_res %>%
    dplyr::rename(
      "p" = "p.adj",
      "p.signif" = "p.adj.signif"
    )
  
  return(temp_res)
}


f.test_dunn <- function(df1, formula, p.adjust.method="holm", ...){
  
  temp_res <- df1 %>% 
    dunn_test(
      formula = formula,
      p.adjust.method = p.adjust.method
    )
  
  return(temp_res)
}


f.test <- function(df1, method, formula, comparisons=NULL, ref.group=NULL, p.adjust.method="holm", paired=F, ...){
  
  obj <- switch (
    method,
    `wilcox` = list(func="f.test_wilcox", formula=formula, comparisons=comparisons, ref.group=ref.group, p.adjust.method=p.adjust.method, paired=paired),
    `t` = list(func="f.test_t", formula=formula, comparisons=comparisons, ref.group=ref.group, p.adjust.method=p.adjust.method, paired=paired),
    `aov_tukey` = list(func="f.test_tukey", formula=formula),
    `kw_dunn` = list(func="f.test_dunn", formula=formula, p.adjust.method=p.adjust.method),
    `chiq` = "f.test_chiq",
    stop("Invalid choice.")
  )
  
  # temp_res <- get(obj)(df1, ...)
  temp_res <- do.call(
    what = get(obj$func),
    args = c(list(df1), obj[names(obj) != "func"])
  )
  
  temp_res <- temp_res %>%
    dplyr::mutate(
      t.method = method
    )
  
  # class(temp_res) <- unique(c(method, class(temp_res)))
  
  return(temp_res)
}


f.test.p <- function(df1, p.method, core.method=NULL, formula, comparisons=NULL, ref.group=NULL, p.adjust.method="holm", paired=F, ...){
  data <- df1
  
  temp_res <- f.func.p(
    data = data,
    func = f.test,
    formula = formula,
    comparisons = comparisons, 
    ref.group = ref.group,
    p.adjust.method = p.adjust.method, 
    paired = paired,
    arg = "method",
    p.arg = p.method,
    core = core.method,
    globals = c( f.find_func("f.") ),
    packages = c("rstatix", "car", "carData", "dplyr")
  )
  
  return(temp_res)
}


f.fc.o <- function(df1, formula, type="logFC_median", ...) {
  
  temp_res <- rstatix::doo(
    data = df1,
    .f = f.fc,
    formula = formula,
    type = type,
    ... 
  )
  
  # class(temp_res) <- unique(c(class(temp_res), paste0("type=", type)))
  
  return(temp_res)
}


f.fc <- function(df1, formula, type="logFC_median", ...) {
  
  formula <- as.formula(formula)
  var <- all.vars(terms(formula)[[2]])
  group <- all.vars(terms(formula)[[3]])
  
  # var <- as.character(formula[[2]])
  # group <- attr(terms(formula), "term.labels")
  lev <- levels(df1[[group]])
  
  if(length(var) > 2 ){
    stop("Error: must one Feature !")
  }
  
  temp <- combn(lev, 2, simplify = FALSE)
  names(temp) <- map_chr(temp, ~ paste(.x, collapse = "/"))
  
  temp_res <- temp %>%
    purrr::imap_dfr(function(levels_pair, diff_name) {
      group1 <- levels_pair[1]
      group2 <- levels_pair[2]
      
      v1 <- df1[[var]][df1[[group]] == group1]
      v2 <- df1[[var]][df1[[group]] == group2]
      
      logFC <- if (any(v1 < 0 )) {
        if (type == "logFC_median") {
          median(v1, na.rm = TRUE) - median(v2, na.rm = TRUE)
        } else if (type == "logFC_mean") {
          mean(v1, na.rm = TRUE) - mean(v2, na.rm = TRUE)
        } else {
          stop("No such calculation method: ", type)
        }
      } else {
        fc1 <- median(v1, na.rm = TRUE)
        fc2 <- median(v2, na.rm = TRUE)
        if (type == "logFC_median") {
          if (fc2 == 0) fc2 <- 1e-10
          log2(fc1 / fc2)
        } else if (type == "logFC_mean") {
          fc1 <- mean(v1, na.rm = TRUE)
          fc2 <- mean(v2, na.rm = TRUE)
          if (fc2 == 0) fc2 <- 1e-10
          log2(fc1 / fc2)
        } else {
          stop("No such calculation method: ", type)
        }
      }
      
      tibble::tibble(
        `Diff` = diff_name,
        `logFC` = logFC,
        `FC.type` = type
      )
    })
  
  # class(temp_res) <- unique(c(class(temp_res), paste0("FC.type=", type)))
  
  return(temp_res)
}


f.test_fc_adj <- function(temp_rd, param, ...){
  # test
  temp_x <- temp_rd %>%
    group_by(Feature) %>%
    f.test.p(
      # p.method = c("wilcox", "t", "aov_tukey", "kw_dunn")[1:2],
      # core.method = T,
      # formula = Value ~ Group,
      p.method = param$p.method,
      core.method = param$core.method,
      formula = param$formula,
      comparisons = "all", 
      # ref.group = NULL,
      ref.group = param$ref.group,
      p.adjust.method = c("holm", "bonferroni", "BH", "none")[4], 
      paired = param$paired
    )
  
  # logFC
  temp_y <- temp_rd %>%
    group_by(Feature) %>%
    f.fc.o(
      # formula = Value ~ Group,
      # type = c("logFC_mean", "logFC_median")[2],
      formula = param$formula,
      type = param$fc_type
    )
  
  # merge & p.adjust & regulate
  temp_res <- temp_x %>%
    map(~{
      .x %>% 
        mutate( Diff = paste(group1, group2, sep = "/") ) %>%
        group_by(Diff) %>%
        mutate(
          adj.P.Val = p.adjust(
            p, 
            # method = "BH"
            method = param$p.adjust.method
          ),
          # adj.method = "BH"
          adj.method = param$p.adjust.method
        ) %>%
        ungroup() %>%
        left_join( temp_y, by = c("Feature", "Diff") )
    })
  
  return(temp_res)
}


f.add_regulate <- function(temp_res, col, FDR, FC.up, FC.down, ...){
  # logFC
  logFC.up <- log(FC.up, 2)
  logFC.down <- log(FC.down, 2)
  
  res <- temp_res %>%
    mutate(
      `temp_p` = !!sym(col)
    ) %>%
    dplyr::mutate(
      Regulate = if_else(
        temp_p < FDR & logFC > logFC.up, "Up",
        if_else(temp_p < FDR & logFC < logFC.down, "Down", "Normal")
      )
    ) %>%
    dplyr::mutate(
      Regulate = factor(
        Regulate,
        levels = c("Up", "Down", "Normal"),
        exclude = NULL
      )
    ) %>%
    dplyr::select(-temp_p) %>%
    mutate(
      Up.cutoff = round(FC.up, digits = 2),
      Down.cutoff = round(FC.down, digits = 2),
      FDR.cutoff = FDR
    )
  
  return(res)
}


f.plot.vol <- function(temp_res, formula, label=NULL, ratio = NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$rd <- temp_res
  
  P.d$df <- P.d$rd %>%
    as_tibble() %>%
    dplyr::rename(
      "x" = f.splitform(formula)$left,
      "y" = f.splitform(formula)$right
    )
  
  # parameter
  P.d$param <- list(
    adj.method = unique(P.d$rd$adj.method),
    FC.type = unique(P.d$rd$FC.type),
    Up = unique(P.d$rd$Up.cutoff),
    Down = unique(P.d$rd$Down.cutoff),
    FDR = unique(P.d$rd$FDR.cutoff)
  )
  
  P.d$param
  
  # # label
  # P.d$label <- P.d$df %>%
  #   group_by(Diff, Regulate) %>%
  #   summarise(
  #     Count = n(),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     label = case_when(
  #       Regulate == "Up" ~ sprintf("Up (n=%d)", Count),
  #       Regulate == "Down" ~ sprintf("Down (n=%d)", Count),
  #       Regulate == "Normal" ~ sprintf("Normal (n=%d)", Count),
  #       TRUE ~ Regulate
  #     )
  #   )
  # P.d$label
  
  # plot
  P1 <- ggplot( P.d$df, aes(x = x, y = -log10(y), color = Regulate) ) +
    geom_point(
      size = 1, alpha = 1, stroke = 0
    ) +
    geom_hline(
      yintercept = -log10(as.numeric(P.d$param$FDR)), 
      color = "gray", linetype = "twodash", linewidth = 0.4
    ) + 
    geom_vline(
      xintercept = log(as.numeric(P.d$param$Up)), 
      color = "gray", linetype = "twodash", linewidth = 0.4
    ) +
    geom_vline(
      xintercept = log(as.numeric(P.d$param$Down), 2), 
      color = "gray", linetype = "twodash", linewidth = 0.4
    ) +
    facet_wrap(
      .~ Diff,
      nrow = 1,
      #scales = "free_y",
      strip.position = "top"
    ) +
    scale_color_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(
      x = expression(paste(Log[2], " (Fold Change)")),
      y = expression(paste("-", log[10], " (adjust P Value)"))
    ) +
    # geom_text(
    #   data = P.d$label %>% filter(Regulate == "Up"), 
    #   aes(
    #     # x = round(max(P.d$df$x)+0.5, 0), 
    #     # y = round(max(-log10(P.d$df$y))+2, 0), 
    #     label = label
    #   ), 
    #   x = Inf,
    #   y = Inf,
    #   vjust = 1.3, hjust = 1.1,
    #   size = 3.3, 
    #   color = colors["Up_text"],    
    # ) +
    # geom_text(
    #   data = P.d$label %>% filter(Regulate == "Down"),
    #   aes(
    #     # x = round(min(P.d$df$x)-0.5, 0),
    #     # y = round(max(-log10(P.d$df$y))+2, 0),
    #     label = label
    #   ),
    #   x = -Inf,
    #   y = Inf,
    #   vjust = 1.3, hjust = -0.1,
    #   size = 3.3,
    #   color = colors["Down_text"],
    # )
    stat_summary(
      data = {subset(P.d$df, Regulate == "Up")},
      aes(x = max(x), y = y),
      fun.data = function(y) {
        data.frame(
          y = Inf,
          label = sprintf("Up (n=%d)", length(y))
        )
      },
      geom = "text",
      hjust = 1, vjust = 1.3,
      size = 3.3,
      color = colors["Up_text"]
    ) +
    stat_summary(
      data = {subset(P.d$df, Regulate == "Down")},
      aes(x = min(x), y = y),
      fun.data = function(y) {
        data.frame(
          y = Inf,
          label = sprintf("Down (n=%d)", length(y))
        )
      },
      geom = "text",
      hjust = 0, vjust = 1.3,
      size = 3.3,
      color = colors["Down_text"]
    )
  
  # print(P1)
  
  if(!is.null(label)){
    
    P.d$label2 <- P.d$df %>% 
      filter(Feature %in% label)
    
    P1 <- P1 + 
      geom_point(
        data = P.d$label2,
        aes(x = x, y = -log10(y)),
        size = 1, alpha = 1, stroke = 0,
        color = "black",
      ) +
      geom_label_repel(
        data = P.d$label2 ,
        aes(label = Feature),
        nudge_x = ifelse(P.d$label2$x > 0, 1, -1),
        nudge_y = -1,
        alpha = 1,
        size = 2, 
        min.segment.length = 0, 
        direction = "y", 
        segment.size = 0.2, 
        # segment.linetype = 6,
        segment.color = "black",
      )
  }
  
  return(P1)
}


f.plot.venn <- function(l1, seed=46, color=NULL, label=NULL, pf=NULL, ...){
  
  if(!is.null(color)){
    if( length(l1) > length(color) ){
      message("Error: Wrong color length")
      return(l1)
      
    }else{
      if( all(names(l1) %in% names(color)) ){
        color <- color[names(l1)]
      }else{
        color <- color[1:length(l1)]
      }
    }
  }else{
    if(is.null(seed)){
      seed <- sample.int(100, 1)
      set.seed(seed)
    }else{
      if(is.numeric(seed)){
        set.seed(seed)
      }else{
        set.seed(46)
      }
    }
    color <- sample(colors(), length(l1))
    
  }
  
  # label
  if( !is.null(label) ){
    if( (length(l1) != length(color)) ){
      message("Error: Wrong color length")
      return(l1)
      
    }else{
      names(l1) <- paste0(label, "\n(n = ", map_int(l1, length), ")")
    }
  }else{
    names(l1) <- paste0(names(l1), "\n(n = ", map_int(l1, length), ")")
  }
  
  # plot
  P1 <- VennDiagram::venn.diagram(
    l1,
    filename = NULL, 
    # filename = pf,
    # height = 4, width = 4, units = "in",
    euler.d = T,        
    scaled = T,           
    alpha = 0.5,          
    col = NA,             
    fill = color,
    #label.col = "white",    
    margin = 0.1,         
    cex = 0.9,            
    inverted = F,         
    # cat.dist = c(0.05, 0.05, 0.04),   
    # cat.pos = c(-20, 20, 180),        
    #category.names =  ,  
    #cat.default.pos = "outer",         
    cat.cex = 0.8,                      
    disable.logging = T,
    fontfamily = "sans",
    cat.fontfamily = "sans",
    sub.fontfamily = "sans",
    main.fontfamily = "sans"
  )
  
  # grid.newpage()
  # grid.draw(P1)
  
  if(!is.null(pf)){
    dir.create(dirname(pf), recursive = TRUE)
    
    pdf(paste0(pf, "-venn.pdf"), width = 5, height = 5)
    grid.draw(P1)
    dev.off()
  }
  
  class(P1) <- c(class(P1), seed)
  
  return(P1)
}


f.plot.upset <- function(l1, ratio=0.6, nintersects=NA, cutoff=NULL, p.scale=1.05, colors=NULL, ...){
  if(is.null(colors)){
    colors = list(
      `sets_bar` = "#E18727", 
      `main_bar` = "#B09C85", 
      `matrix` = "#a7c2c2", 
      `shade` = "#e0b95e"
    )
  }
  
  mb.ratio <- c(ratio, 1-ratio)
  
  set_size.scale_max <- l1 %>% map(~length(unique(.x))) %>% unlist() %>% max()*p.scale

  # data
  P.d <- list()
  P.d$rd <- l1
  
  P.d$unique <- unique(unlist(P.d$rd))
  
  P.d$df <- P.d$rd %>% 
    map(~{
      as.numeric(P.d$unique %in% .x)
    }) %>%
    bind_cols() %>%
    mutate(
      rowname =  P.d$unique
    ) %>%
    tibble::column_to_rownames("rowname")
  
  P1 <- upset(
    fromList(P.d$rd),
    nsets = length(P.d$rd),
    order.by = "freq", 
    nintersects = nintersects, 
    # empty.intersections = "on",
    mb.ratio = mb.ratio,
    keep.order = T,
    cutoff = cutoff,
    # mainbar.y.label = "Shared feature number",
    # sets.x.label = "Total feature number",
    point.size = 2.5,
    line.size = 0.5,
    sets.bar.color = colors$sets_bar,
    main.bar.color = colors$main_bar,
    matrix.color = colors$matrix,
    shade.color = colors$shade,
    set_size.show = T,
    set_size.scale_max = set_size.scale_max,
    show.numbers = "yes"
  )
  
  # print(P1)
  
  P1$df <- P.d$df
  
  return(P1)
}

