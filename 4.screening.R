# screening ----
## predata ----
f.inTrain <- function(df1, seed=123, y=NULL, p=0.7, ...){
  
  if(p > 1){ p <- p/nrow(df1) }
  
  set.seed(seed)
  # ?stratified sampling
  if(is.null(y)){
    inT <- sample(1:nrow(df1), size = p* nrow(df1))
  }else{
    inT <- caret::createDataPartition(y = df1[[y]], times=1, p = p, list = F, ...)
  }
  
  return(inT)
}


f.part <- function(df1, inT, ...){
  
  temp_res <- list(
    `train` = df1[inT, ],
    `valid` = df1[-inT, ],
    `all` = df1
  )
  
  return(temp_res)
}


f.resampling <- function(df1, method="none", formula, seed=123, over_ratio=NULL, ...){
  
  resp <- f.splitform(formula)$left
  # lev <- levels(df1[[resp]])
  # var <- f.splitform(formula)$right
  y <- f.splitXy(df1, formula)$y
  X <- f.splitXy(df1, formula)$X
  data <- f.splitXy(df1, formula)$data
  
  set.seed(as.integer(seed))
  temp_res <- switch (
    method,
    "none" = df1,
    "up" = upSample(x = X, y = y, yname = resp),
    "down" = downSample(x = X, y = y, yname = resp),
    "rose" = ROSE::ROSE(formula = formula, data  = data)$data,
    "smote" = themis::smote(df = data, var = resp)
  )
  
  temp_res <- temp_res %>%
    dplyr::select(all_of(resp), everything())
  
  class(temp_res) <- c(class(temp_res), "resampling", as.character(method))
  
  return(temp_res)
}

## glmnet ----
### fit ----
f.glmnet.alpha <- function(alpha){
  if (is.numeric(alpha)) {
    if (length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
      stop("alpha must be a single numeric value in [0, 1].")
    }
    return(alpha)
  }
  
  if (is.character(alpha)) {
    alpha <- trimws(alpha)
    if (alpha == "Ridge") {
      return(0)
    } else if (alpha == "Lasso") {
      return(1)
    } else if (grepl("^Enet\\s*\\[", alpha)) {
      match <- sub("^Enet\\s*\\[\\s*([0-9]*\\.?[0-9]+)\\s*\\]$", "\\1", alpha, perl = TRUE)
      if (match == alpha) {
        stop("Invalid Enet format: ", alpha)
      }
      num <- as.numeric(match)
      if (num < 0 || num > 1) {
        stop("Enet alpha must be in [0, 1], got: ", num)
      }
      return(num)
    } else {
      stop("Unknown alpha string: '", alpha, "'. Use 'Ridge', 'Lasso', or 'Enet[x]'.")
    }
  }

  stop("alpha must be numeric or a character string.")
}


f.glmnet <- function(df1, formula, alpha=1, nfold=10, seed=123, ...){
  # data
  traindata <- df1
  
  X <- f.splitXy(traindata, formula)$X %>% as.matrix()
  y <- f.splitXy(traindata, formula)$y
  
  if (inherits(y, "Surv")) {
    family <- c("binomial", "cox")[2]
  }else{
    family <- c("binomial", "cox")[1]
  }
  
  alpha <- f.glmnet.alpha(alpha)
  
  set.seed(seed)
  # modeling
  temp_res <- cv.glmnet(
    x = X,      
    y = y,  
    nfolds = nfold, 
    alpha = alpha,         # 0/Ridge/L2    1/Lasso/L1   (0, 1)/Enet
    family = family,
  )
  # Non-zero variable with coefficient
  var.min <- rownames(coef(temp_res, s = "lambda.min"))[coef(temp_res, s = "lambda.min")[, 1] != 0]
  if (is.element("(Intercept)", var.min)) {
    var.min <- sort(var.min[-1])
  } else {
    var.min <- sort(var.min)
  }
  
  var.1se <- rownames(coef(temp_res, s = "lambda.1se"))[coef(temp_res, s = "lambda.1se")[, 1] != 0]
  if (is.element("(Intercept)", var.1se)) {
    var.1se <- sort(var.1se[-1])
  } else {
    var.1se <- sort(var.1se)
  }
  
  temp_res$var.min <- var.min
  temp_res$var.1se <- var.1se
  
  if(alpha == 1){
    metric <- "Lasso"
  }else if(alpha == 0){
    metric <- "Ridge"
  }else{
    metric <- "Enet"
  }
  
  temp_res$param <- list(
    modeling = "glmnet",
    formula = formula,
    metric = metric,
    alpha = alpha,
    nfold = nfold
  )
  
  res <- list(
    fit = temp_res,
    str = paste(var.min, collapse = " "),
    param = list(
      formula = formula,
      metric = metric,
      alpha = alpha,
      nfold = nfold
    )
  )
  
  return(res)
}


f.glmnet.res <- function(cvfit, ...){
  
  temp <- data.frame(
    lambda = cvfit$lambda,
    nzero = cvfit$nzero,
    mean = cvfit$cvm,
    meanSD = cvfit$cvsd,
    metric = cvfit$name
  ) %>%
    mutate(
      label = case_when(
        lambda == cvfit$lambda.min ~ "min error",
        lambda == cvfit$lambda.1se ~ "error in 1se",
        TRUE ~ NA_character_
      )
    )
  
  return(temp)
}


f.glmnet.p <- function(df1, p.seed, core.seed=NULL, formula, alpha=1, nfold=10, index=c("var.min", "var.1se")[1], ...){
  data <- df1
  formula <- formula
  
  temp_res <- f.func.p(
    data = data,
    func = f.glmnet,
    formula = formula,
    alpha = alpha, 
    nfold = nfold, 
    arg = "seed",
    p.arg = p.seed,
    core = core.seed,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "Matrix", "glmnet", "survival", "pbapply")
  )
  
  temp_res <- temp_res %>%
    f.glmnet.p_res(
      index = index
    )
  
  return(temp_res)
}


f.glmnet.p_res <- function(temp_res, index=c("var.min", "var.1se")[1], ...){
  
  temp <- sapply(temp_res, function(x){
    data.frame(
      `str` = x$str,
      `length` =  length(x$fit[[index]]) #length(unlist(strsplit(x$str, "=")))
    )
  }) %>% 
    t() %>% data.frame() %>% 
    tibble::rownames_to_column("order") %>%
    add_count(str, name = "count") %>%
    distinct(str, .keep_all = TRUE) %>%
    arrange(desc(count)) %>%
    mutate(
      name = paste0("fit_", row_number())
    ) %>%
    mutate(
      name = factor(
        name, 
        levels = paste0("fit_", row_number()), 
        exclude = NULL
      ),
      metric = temp_res[[1]]$param$metric
    )
  # keep only fit
  temp_res <- temp_res[as.numeric(temp$order)]
  names(temp_res) <- paste0("fit_", 1:length(temp_res))
  
  res <- list(
    glmnet = temp_res,
    times = temp,
    fit = temp_res$fit_1$fit,
    Covar = temp_res$fit_1$fit[[index]]
  )
  
  return(res)
}


f.glmnet.mp <- function(df1, m.alpha, p.seed, core.seed, formula, nfold=10, index=c("var.min", "var.1se")[1], ...){
  
  data <- df1
  temp_res <- f.func.m(
    data = data, 
    func = f.glmnet.p, 
    p.seed = p.seed,
    core.seed = core.seed,
    formula = formula,
    nfold = nfold, 
    index = index,
    arg = "alpha", 
    m.arg = m.alpha
  )
  
  return(temp_res)
}


### plot ----
f.plot.glmnet.cv <- function(temp_res, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  P.d$sec_axis <- P.d$df %>%
    mutate(
      log_lambda = log(lambda)
    ) %>%
    arrange(log_lambda) %>%
    mutate(
      # is_key = log_lambda %in% {
      #   # key_x <- log(c(cvfit$lambda.min, cvfit$lambda.1se))
      #   key_x <- log(c(cvfit$lambda.min))
      #   sapply(key_x, function(px) log_lambda[which.min(abs(log_lambda - px))])
      # },
      # is_m = row_number() %in% c(1L, n()),
      is_m = row_number() %in% c(n()),
      is_key = row_number() %in% c(n()),
      is_keep = is_key | is_m
    ) %>%
    {
      key_df <- filter(., is_keep)
      rest_df <- filter(., !is_key)
      n_max <- 10
      n_rest <- max(0, n_max - nrow(key_df))
      if (n_rest > 0 && nrow(rest_df) > 0) {
        step <- max(1, floor(nrow(rest_df) / n_rest))
        rest_sample <- slice(rest_df, seq(1, nrow(rest_df), by = step)[1:n_rest])
      } else {
        rest_sample <- rest_df[0, ]
      }
      bind_rows(key_df, rest_sample) %>% arrange(log_lambda)
    }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = log(lambda), y = mean)) +
    geom_errorbar(
      aes(ymin = mean - meanSD, ymax = mean + meanSD, color = as.factor(nzero)), 
      width = 0.2, linewidth = 0.4, #color = "darkgrey"
      show.legend = F
    ) +
    geom_point(
      size = 0.7, fill = "white", 
      color = "red"
    ) +
    geom_vline(
      xintercept = log(subset(P.d$df, !is.na(label))$lambda), 
      color = "black", linetype = "dotted", linewidth = 0.4
    ) +
    geom_text(
      data = subset(P.d$df, !is.na(label)), 
      aes(
        x = log(lambda), 
        # y = if_else(label == "min error", diff(range(P.d$df$mean))/2, diff(range(P.d$df$mean))/1.5),
        y = Inf,
        label = label
      ), 
      vjust = 1.3, hjust = 0.5,
      size = 3.3, color = "firebrick"    
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank()
    ) +
    labs(
      x = expression(paste(ln, " (Lambda)")) ,
      y = unique(P.d$df$metric)
    ) +
    scale_y_continuous(
      # breaks = seq(from = floor(min(P.d$df$mean)), to = ceiling(max(P.d$df$mean)), by = 0.2)
      breaks = pretty_breaks(n = 10)
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 10),
      expand = c(0.03, 0, 0.03, 0),
      sec.axis = sec_axis(
        ~ .,
        name = "Number of non-zero coefficients",
        breaks = P.d$sec_axis$log_lambda,
        labels = P.d$sec_axis$nzero
      )
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.glmnet.p_times <- function(temp_res, ratio=NULL, colors="skyblue", theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = name, y = count)) +
    geom_bar(
      stat = "identity", width = 0.7, 
      fill = colors
    ) + 
    geom_text(
      aes(label = count), 
      position = position_dodge(width = 0.7), 
      vjust = -0.5, hjust = 0.5,      
      color = "black", size = 3.3  
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
      legend.position = "none", 
    ) + 
    scale_y_continuous(
      expand = c(0, 0, 0.1, 0)
    ) +    
    labs(
      y = "Count"
    )
  
  # print(P1)
  
  return(P1)
}


## boruta ----
f.boruta <- function(df1, formula, ...){
  # data
  traindata <- df1
  
  # boruta
  set.seed(8)
  temp_res <- Boruta(
    formula = formula, 
    data = traindata, 
    # maxRuns = 500,      # iNumber of iterations
    # getImp = getImpRfZ, 
    # default: RF, Z-score(Mean Decrease Accuracy, MDA)
    # (actual_importance - shadow_mean) / shadow_sd   
    doTrace = 1
  )

  res <- list(
    boruta = temp_res,
    Covar = names(temp_res$finalDecision[temp_res$finalDecision == "Confirmed"])
  )
  
  return(res)
}


f.boruta.res <- function(temp_res, ...){
  
  res <- temp_res$boruta
  
  temp <- res$ImpHistory %>%
    as.data.frame() %>%
    pivot_longer(
      cols = everything(), 
      names_to = "Feature",
      values_to = "Imp"
    ) %>%
    dplyr::mutate(
      Decision = case_when(
        grepl("shadow", Feature) ~ "Shadow",
        Feature %in% names(res$finalDecision) ~ res$finalDecision[Feature]
      )
    ) %>%
    dplyr::filter(!is.infinite(Imp)) %>%
    dplyr::group_by(Feature) %>%
    dplyr::mutate(
      Median_Imp = median(Imp, na.rm = TRUE)
    ) %>%
    dplyr::arrange( (Median_Imp) ) %>% 
    # mutate( order = cur_group_id() ) %>%
    dplyr::mutate(
      Decision = factor(
        Decision,
        levels = c("Confirmed", "Tentative", "Rejected", "Shadow"),
        exclude = NULL
      )
    ) %>%
    dplyr::select( -Median_Imp )
  
  temp$Feature <- factor(temp$Feature, levels = unique(temp$Feature))
  
  return(temp)
}


f.plot.boruta <- function(temp_res, label=F, ratio=NULL, colors=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Feature, y = Imp, fill = Decision)) +
    geom_boxplot(
      # fill = "white",
      color = "black",
      linewidth = 0.15,
      box.linetype = "solid",
      whisker.linetype = "dashed",
      # key_glyph = "rect",
      # show.legend = FALSE,
      # outlier.shape = NA,
      outlier.shape = 21,
      outlier.color = "black", 
      outlier.fill = "white",  
      outlier.stroke = 0.15,
      outlier.size = 0.2
    ) +
    # scale_color_manual(
    #   values = colors,
    # ) +
    scale_fill_manual(
      values = colors,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(), 
      axis.text.x = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      legend.title = element_blank(),
      legend.position = "inside",
      # legend.direction = "horizontal",
      legend.position.inside = c(0.03, 0.97),  
      legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA, linewidth = 0.01),
      legend.key.size = unit(0.3, "cm"),
    ) +
    labs(
      x = "Features",
      y = "Importance"
    ) +
    scale_y_continuous(
      expand = c(0.05, 0, 0.2, 0),
      breaks = pretty_breaks(n = 10)
    ) +
    scale_x_discrete(
      expand = c(0.01, 0, 0.01, 0)
    )
  
  if( label ){
    P1 <- P1 + 
      geom_text(
        data = filter(P.d$df, Decision != "Shadow"),
        aes(
          label = Feature,
          x = Feature,
          y = if_else(
            ave(Imp, Feature, FUN = median) > median(P.d$df$Imp),
            ave(Imp, Feature, FUN = min),
            ave(Imp, Feature, FUN = max)
          )
        ),
        vjust = 0.5,
        hjust = 0.5,
        angle = 90,
        size = 2.3,
        color = "gray20"
      )
  }
  
  # print(P1)
  
  return(P1)
}


## rfe ----
### fit ----
f.rfe_auc <- function(df1, formula, seed=005, nfold=10, ...){
  # data
  traindata <- df1
  
  X <- as.matrix(f.splitXy(traindata, formula)$X)
  y <- f.splitXy(traindata, formula)$y
  
  # preprocess
  {
    rfFuncs$summary <- twoClassSummary
    
    set.seed(002)
    rfe_ctrl <- rfeControl(
      functions = rfFuncs,    
      method = "repeatedcv",    
      # number = 10, 
      number = nfold,
      repeats = 1,              
      # rerank = T,
      saveDetails = T,          
      returnResamp = "all",     
      verbose = TRUE            
    )
    
    set.seed(003)
    tr_ctrl <- trainControl(
      classProbs = TRUE,  
      summaryFunction = twoClassSummary 
    )
    
    # rfe
    set.seed(seed)
    temp_res <- rfe(
      X,  
      y,  
      sizes = c(1:ncol(X)), 
      rfeControl = rfe_ctrl, 
      method = "rf",         
      metric = "ROC",        
      trControl = tr_ctrl    
    )
  }
  
  res <- list(
    rfe = temp_res,
    Covar = temp_res$optVariables
  )
  
  return(res)
}


f.rfe_acc <- function(df1, formula, seed=005, nfold=10, ...){
  # data
  traindata <- df1
  
  X <- as.matrix(f.splitXy(traindata, formula)$X)
  y <- f.splitXy(traindata, formula)$y
  
  # preprocess
  {
    set.seed(125)
    rfe_ctrl <- rfeControl(
      functions = rfFuncs,     
      method = "repeatedcv",
      # number = 10, 
      number = nfold,
      repeats = 1, 
      # rerank = T,
      # metric = 
      saveDetails = T,          
      returnResamp = "all",
      verbose = TRUE
    )
    
    # rfe
    set.seed(seed)
    temp_res <- rfe(
      x = X, 
      y = y,
      sizes = c(2:ncol(X)), 
      rfeControl = rfe_ctrl
    )
  }
  
  res <- list(
    rfe = temp_res,
    Covar = temp_res$optVariables
  )
  
  return(res)
}


f.rfe <- function(df1, metric="AUC", formula, seed=005, nfold=10, ...){
  obj <- switch (
    metric,
    `AUC` = "f.rfe_auc",
    `Accuracy` = "f.rfe_acc",
    stop("Invalid choice.") 
  )
  
  temp_res <- get(obj)(df1, formula, seed, nfold, ...)

  return(temp_res)
}


f.rfe.p <- function(df1, p.metric, core.metric=NULL, formula, seed=005, nfold=10, ...){
  data <- df1
  formula <- formula
  
  temp_res <- f.func.p(
    data = data,
    func = f.rfe,
    formula = formula,
    seed = seed,
    nfold = nfold,
    arg = "metric",
    p.arg = p.metric,
    core = core.metric,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "caret")
  )
  
  return(temp_res)
}


f.rfe.res <- function(temp_res, ...){
  res <- temp_res
  
  temp <- cbind(
    res$results %>%
      dplyr::select(Variables, !matches("SD")) %>%
      pivot_longer(
        cols = -c("Variables"), #everything(), 
        names_to = "Metric", 
        values_to = "Mean"  
      ),
    res$results %>%
      dplyr::select(Variables, matches("SD")) %>%
      pivot_longer(
        cols = -c("Variables"), #everything(), 
        names_to = "Metric", 
        values_to = "SD"  
      ) %>%
      dplyr::select(-Metric, -Variables)
  ) %>%
    mutate(
      Metric = case_when(
        Metric == "ROC" ~ "AUC",
        Metric == "Sens" ~ "Sensitivity",
        Metric == "Spec" ~ "Specificity",
        TRUE ~ Metric
      )
    ) %>%
    filter(Metric %in% c("AUC", "Accuracy", "Sensitivity", "Specificity")) %>%
    mutate(
      Metric = factor(
        Metric, 
        levels = c("AUC", "Accuracy", "Sensitivity", "Specificity"), 
        exclude = NULL
      )
    ) %>%
    mutate(
      lev = paste0(levels(res$pred$pred), collapse = "/"),
      optVariables = length(res$optVariables),
      metric = if_else(res$metric == "ROC", "AUC", res$metric)
    ) %>%
    mutate(
      facet = paste0("RFE (", lev, ")")
    )
  
  return(temp)
}


f.rfe.var <- function(temp_res, No, ...){
  
  res <- temp_res$rfe
  
  res$No <- temp_res$rfe$No <- No
  
  temp_res$Covar <- if( length(res$optVariables) >= res$No ){
    head(res$optVariables, res$No)
  }else{
    res$variables %>% filter(Variables == res$No) %>%
      group_by(var) %>% summarise(AverageOverall = mean(Overall, na.rm = TRUE)) %>%
      arrange(desc(AverageOverall)) %>% head(res$No) %>% pull(var)
  }
  
  return(temp_res)
}

f.rfe.var2 <- function(temp_res, No, ...){
  
  if(is.list(temp_res)){
    temp_res %>%
      map(~{
        f.rfe.var(.x, No, ...)
      })
  }else{
    f.rfe.var(temp_res, No, ...)
  }
  
}


### plot ----
f.plot.rfe2 <- function(temp_res, xlim=15, ...){
  
  if(is.list(temp_res)){
    P1 <- temp_res %>%
      map(~{
        plot(
          .x$rfe, type = c("g", "o"), 
          ylim = c(0.45, 1.05), xlim = rev(c(0.1, xlim))
        )
      })
    aplot::plot_list( gglist = P1 )
  }else{
    plot(
      temp_res$rfe, type = c("g", "o"), 
      ylim = c(0.45, 1.05), xlim = rev(c(0.1, xlim))
    )
  }
  
}


f.plot.rfe <- function(temp_res, No=NULL, xlim=NULL, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  # plot
  P1 <- ggplot( P.d$df, aes(x = Variables, y = Mean, color = Metric) ) +
    geom_point(
      size = 1, 
      fill = "white",
    ) +
    geom_ribbon(
      aes(ymin = Mean - SD, ymax = Mean + SD, fill = Metric), 
      color = "lightgray", alpha = 0.1, 
      linewidth = 0.3
    ) +
    geom_line(
      linewidth = 0.3, 
    ) +
    scale_color_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      # panel.grid = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = 10, color = "black"),
      strip.background.y = element_rect(fill = NA, color = NA)
    ) +
    labs(
      x = "Number of Variables"
    ) +
    scale_y_continuous(
      breaks = seq(from = 0, to = 1, by = 0.1),
      limits = c(0.5, NA)
    ) +
    scale_x_reverse()
  
  # print(P1)
  
  if(!is.null(xlim)){
    
    P.d$breaks <- if(max(P.d$df$Variables) < 10){
      rev(seq(from = 1, to = max(P.d$df$Variables), by = 1))
    }else{
      rev(c(seq(from = 1, to = 9, by = 1), seq(from = 11, to = max(P.d$df$Variables), by = 2)))
    }
    
    P1 <- P1 +
      scale_x_reverse(
        breaks = P.d$breaks,
        expand = c(0.05, 0, 0.05, 0)
      )
    
    if( temp_res$optVariables[1] > xlim ){
      P1 <- P1 +
        coord_cartesian(xlim = rev(c(1, xlim)))
    }
    
  }

  if( !is.null(No) ){
    P1 <- P1 + geom_vline(
      aes(xintercept = No),
      # xintercept = No, 
      linetype = "dotted", linewidth = 0.4, 
      color = "black",
    )
  }
  
  return(P1)
}


## auc/cor ----
### calculate ----
f.auc_ind <- function(df1, formula, direction="auto", cutoff="yd", ...){
  resp <- f.splitform(formula)$left
  lev <- levels(df1[[resp]])
  
  # data
  traindata <- df1
  # evaluate
  temp_res <- traindata %>%
    f.eval(
      formula = formula,
      direction = direction, 
      cutoff = cutoff
    )
  # extraction
  temp <- temp_res$metric %>%
    filter(Metric %in% c("AUC", "Sensitivity", "Specificity")) %>%
    tidyr::pivot_wider( names_from = Metric, values_from = Value ) %>%
    mutate(
      lev = paste(lev, collapse = "/")
    )
  
  return(temp)
}


f.auc_ind.o <- function(df1, formula, direction="auto", cutoff="yd", ...) {
  
  temp_res <- rstatix::doo(
    data = df1,
    .f = f.auc_ind,
    formula = formula,
    direction = direction,
    cutoff = cutoff,
    ... 
  )
  
  return(temp_res)
}

# temp_res <- traindata %>%
#   pivot_longer(
#     cols = -Group,
#     names_to = "Feature",
#     values_to = "Value"
#   ) %>%
#   group_by(Feature) %>%
#   f.auc_ind.o(
#     formula = Group ~ Value,
#     direction = c("auto", `down`=">", `up`="<")[1],
#     cutoff = c("yd", 0.5)[1]
#   )

f.auc_ind.p <- function(df1, formula, core.formula=NULL, direction="auto", cutoff="yd", auc.cut=0.7, ...) {
  
  resp <- f.splitform(formula)$left
  lev <- levels(df1[[resp]])
  if(length(lev) > 3){
    stop("must two levels!")
  }
  
  var <- f.splitform(formula)$right
  
  if(all("." %in% var)){
    var <- setdiff(names(df1), resp)
  }
  
  p.formula <- sapply(var, function(v) as.formula(paste(resp, "~", v)))

  data <- df1
  temp_res <- f.func.p(
    data = data,
    func = f.auc_ind,
    direction = direction,
    cutoff = cutoff,
    arg = "formula",
    p.arg = p.formula,
    core = core.formula,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "tidyr", "caret", "pROC")
  )
  
  temp_res <- temp_res %>% bind_rows() %>%
    mutate(
      Feature = var,
      auc.cut = auc.cut
    ) %>%
    dplyr::select("Feature", everything())
  
  
  if(!is.null(auc.cut)){
    Covar <- temp_res %>% filter(AUC > auc.cut) %>% pull(Feature)
  }else{
    Covar <- temp_res %>% pull(Feature)
  }
  
  res <- list(
    auc = temp_res,
    Covar = Covar
  )
  
  return(res)
}

# temp_res <- traindata %>%
#   f.auc_ind.p(
#     formula = Group ~ .,
#     direction = c("auto", `down`=">", `up`="<")[1],
#     cutoff = c("yd", 0.5)[1],
#     auc.cut = 0.7,
#     core.formula = NULL
#   )


f.aucor <- function(df1, formula, core.formula=NULL, direction="auto", cutoff="yd", auc.cut=NULL, cor.cut=NULL, cor.meth=NULL, ...){
  # data
  traindata <- df1
  
  # individual AUC
  temp_res <- traindata %>%
    f.auc_ind.p(
      formula = formula,
      lev = lev,
      direction = direction, 
      cutoff = cutoff,
      auc.cut = auc.cut,
      core.formula = core.formula
    )
  
  if(is.null(auc.cut) | is.null(cor.cut) | is.null(cor.meth)){
    return(temp_res)
  }
  
  if(!is.null(auc.cut)){
    temp_var <- temp_res[["Covar"]]
    
    if(length(temp_var) == 0){
      stop("Error: No feature > ", auc.cut)
    }
  }
  
  # Correlation of stuff (AUC > auc.cut)
  temp_res[["cor"]] <- traindata %>%
    dplyr::select(all_of(temp_var)) %>%  # columns: variables/features   # rows: observations/sample
    cor(method = cor.meth, use = "complete.obs") %>%
    as.data.frame() %>% 
    tibble::rownames_to_column("Covar1") %>% 
    pivot_longer(-Covar1, names_to = "Covar2", values_to = "Correlation") %>%
    filter( Covar1 < Covar2 ) %>%
    mutate(
      cor.cut = cor.cut
    )
  
  # filter: |cor| > cor.cut & keep high auc
  temp <- temp_res[["cor"]] %>%
    filter( abs(Correlation) > cor.cut )
  
  temp_res[["Covar"]] <- pblapply( unique(temp$Covar1), function(i){
    j <- temp %>% filter(Covar1 %in% i) %>% pull(Covar2)
    
    res <- temp_res[["auc"]] %>%
      rename("Covar" = "Feature") %>%
      filter(Covar %in% j) %>%
      filter(AUC == max(AUC)) %>%
      pull(Covar)
    
    return(res)
  }) %>%
    unlist()
  
  temp_res <- temp_res[c("auc", "cor", "Covar")]
  
  return(temp_res)
}


### plot ----
f.plot.auc.hist <- function(temp_res, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  cutoff <- P.d$df$auc.cut[1]
  
  # density
  P.d$density <- P.d$df %>%
    reframe(
      x = density( as.numeric(AUC) )$x,
      y = density( as.numeric(AUC) )$y
    )
  
  # label
  P.d$label <- P.d$density %>%
    summarise(
      down_x = median(x[x < cutoff]),
      up_x = median(x[x >= cutoff]),
      down_y = approx(x[x < cutoff], y[x < cutoff], xout = median(x[x < cutoff]))$y / 2,
      up_y = approx(x[x >= cutoff], y[x >= cutoff], xout = median(x[x >= cutoff]))$y / 2
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("group", ".value"),
      names_sep = "_",
      values_drop_na = TRUE
    ) %>%
    left_join(
      P.d$df %>%
        reframe(
          down = sum(AUC < cutoff),
          up = sum(AUC >= cutoff)
        ) %>%
        pivot_longer(
          cols = everything(),
          names_to = "group",
          values_to = "label"
        ),
      by = "group"
    )
  P.d$label
  
  # plot
  P1 <- ggplot(P.d$df, aes( x = as.numeric(AUC) )) + 
    geom_density(
      colour = "black",
      fill = NA,
      # alpha = 1,
      linewidth = 0.5,
      key_glyph = c("smooth", "rect")[2],
      # show.legend = FALSE
    ) +
    geom_area(
      data = subset(P.d$density, x < cutoff), 
      aes(x, y, fill = "1"), 
      alpha = 0.5
    ) +
    geom_area(
      data = subset(P.d$density, x >= cutoff), 
      aes(x, y, fill = "2"), 
      alpha = 0.5
    ) +
    geom_histogram(
      aes( y = after_stat(density) ),
      # binwidth = 0.05,
      # boundary = 0,
      fill = NA,
      alpha = 0.7,
      color = "gray50"
    ) +
    scale_fill_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      # panel.grid = element_blank(),
      legend.position = "none",
    ) +
    labs(
      x = "Individual AUC",
      y = "Density",
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.05, 0)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.1),
      limits = c(0, 1)
    ) +
    geom_text(
      data = P.d$label,
      aes(label = label, y = y, x = x),
      vjust = 0.5,
      hjust = 0.5,
      #lineheight = 0.7,
      size = 3.3,
      color = "firebrick",
      show.legend = FALSE
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.auc.rank <- function(temp_res, label2, label=NULL, colors, ratio=NULL, theme="bw", ...){

  P.d <- list()
  P.d$df <- temp_res %>%
    pivot_longer(
      cols = c("AUC", "Sensitivity", "Specificity"),  
      names_to = "Metric", 
      values_to = "Value"  
    ) %>%
    group_by(Metric) %>%
    mutate(Rank = rank(dplyr::desc(Value), ties.method = "min")) %>%
    ungroup()
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Rank, y = Value, color = Metric)) +
    # facet_wrap(
    #   .~ Metric,
    #   nrow = 1,
    #   #scales = "free_y",
    #   strip.position = "top"
    # ) +
    geom_point(
      data = subset(P.d$df, !Feature %in% label2),
      color = "gray",
      size = 1,
      alpha = 1
    ) +
    geom_point(
      data = subset(P.d$df, Feature %in% label2),
      # color = "gray",
      size = 1,
      alpha = 1
    ) +
    # geom_segment(
    #   aes(x = Rank, xend = Rank, y = Value, yend = 0, color = Metric),
    #   show.legend = F
    # ) +
    scale_color_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      legend.position = "none",
      # axis.title.y = element_blank(),
    ) +
    labs(
      x = "Rank",
      y = "individual Metric"
    ) +
    scale_y_continuous(
      expand = c(0.05, 0, 0.05, 0),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    )
  
  # print(P1)
  
  if(!is.null(label)){
    P.d$label2 <- P.d$df %>%
      filter(Feature %in% label)
    
    P1 <- P1 + 
      geom_point(
        data = P.d$label2,
        aes(x = Rank, y = Value),
        shape = 21, size = 1.2, alpha = 1, stroke = 0.15,
        color = "black",
        fill = "black"
      ) +
      geom_label_repel(
        data = P.d$label2 ,
        aes(x = Rank, y = Value, label = Feature),
        fill = NA, 
        alpha = 1,
        size = 2.3, 
        nudge_x = 3,
        nudge_y = -0.2,
        min.segment.length = 0, 
        direction = c("both", "y", "x")[1], 
        segment.size = 0.2, 
        # segment.linetype = 6,
        segment.color = "black",
      )
  }
  
  return(P1)
}


## screening ----
f.screening <- function(df1, method, formula, ...){
  if(!is.character(method) || length(method) != 1){
    stop("method must be a single character string")
  }
  
  args <- list(...)
  
  method <- trimws(method)
  
  is_glmnet <- method %in% c("Ridge", "Lasso") || grepl("^Enet\\s*\\[.*\\]$", method)
  if(is_glmnet){
    alpha <- f.glmnet.alpha(method)
    args$alpha <- alpha
  }
  
  is_auc    <- method == "auc"
  is_aucor  <- method == "aucor"
  is_boruta <- method == "Boruta"
  
  is_rfe    <- grepl("^RFE", method)
  if(is_rfe){
    if (grepl("\\[.*\\]", method)) {
      metric <- sub(".*\\[(.*?)\\].*", "\\1", method)
      args$metric <- metric
    }
  }
  
  is_unicox <- method == "uniCox"
  is_multicox <- method == "multiCox"
  is_boostcox <- method == "CoxBoost"
  is.rsf <- method == "rsf"
  
  is_logrank    <- grepl("^logrank", method)
  if(is_logrank){
    if (grepl("\\[.*\\]", method)) {
      cutoff <- sub(".*\\[(.*?)\\].*", "\\1", method)
      args$type <- cutoff
    }
  }
  
  is_stepcox <- grepl("^stepCox", method)
  if(is_stepcox){
    if (grepl("\\[.*\\]", method)){
      step <- sub(".*\\[(.*?)\\].*", "\\1", method)
      args$step <- step
    }
  }

  func <- case_when(
    is_glmnet ~ "f.glmnet.p",
    is_auc    ~ "f.auc_ind.p",
    is_aucor  ~ "f.aucor",
    is_boruta ~ "f.boruta",
    is_rfe    ~ "f.rfe", 
    is_unicox ~ "f.cox_uni",
    is_multicox ~ "f.cox_multi",
    is_boostcox ~ "f.cox_boost",
    is_logrank ~ "f.logrank.p",
    is_stepcox ~ "f.cox_step",
    is.rsf ~ "f.rsf",
    TRUE      ~ NA_character_
  )
  
  if(is.na(func)){
    stop("Unknown method: '", method, "'. Supported: Ridge, Lasso, Enet[...], auc, aucor, Boruta, RFE, RFE[metric] ...")
  }

  temp_res <- do.call(
    what = get(func, mode = "function"),
    args = c(list(df1), list(formula=formula), args)
  )
  
  return(temp_res)
}


f.screening.m <- function(df1, m.method, formula, ...){
  
  data <- df1
  temp_res <- f.func.m(
    data = data, 
    func = f.screening, 
    formula = formula,
    arg = "method", 
    m.arg = m.method,
    ...
  )
  
  return(temp_res)
}


f.screen_seq <- function(df1, method_seq, formula, ...){
  args <- list(...)
  
  resp <- f.splitform(formula)$left
  var <- f.splitform(formula)$right
  
  cat("Total Feature:", length(var), "\n")
  
  temp_res <- list()
  
  for(method in method_seq){
    
    if(length(resp) >= 2){
      formula <- reformulate(termlabels = var, response = paste("Surv(", resp[1], ",", resp[2], ")") )
    }else{
      formula <- reformulate(termlabels = var, response = resp)
    }
    
    temp_res[[method]] <- do.call(
      what = f.screening,
      args = c(list(df1), list(method=method), list(formula=formula), args)
    )
    
    var <- temp_res[[method]]$Covar
    cat(method, "ing refine Feature:", length(var), "\n")
    
    if(length(var) == 0){ break() }
  }

  return(temp_res)
}


f.screen_seq.m <- function(df1, m.method_seq, formula, ...){
  
  data <- df1
  temp_res <- f.func.m(
    data = data, 
    func = f.screen_seq, 
    formula = formula,
    arg = "method_seq", 
    m.arg = m.method_seq,
    ...
  )
  
  return(temp_res)
}


## opt_fea_info ----
f.plot.var.box <- function(temp_res, facet, log=T, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      facet = !!sym(facet)
    )
  
  P.d$median <- P.d$df %>%
    dplyr::group_by(Group, Feature) %>%
    dplyr::summarize(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = Median) %>%
    rowwise() %>%
    mutate(
      max_median = max(c_across(-Feature), na.rm = TRUE),
      max_group = names(pick(-Feature))[which.max(c_across(-Feature))]
    ) %>%
    ungroup() %>%
    mutate(
      max_group = factor(
        max_group,
        levels = levels(P.d$df$Group)
      )
    ) %>%
    arrange(max_group) %>%
    pull(Feature)
  P.d$median
  
  P.d$df$Feature <- factor(P.d$df$Feature, levels = P.d$median )
  
  P.d$df <- P.d$df %>% arrange(Feature)
  
  if(log){
    P.d$df <- P.d$df %>%
      mutate(
        Value = log2(Value)
      )
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Group, y = Value, fill = Group)) +
    facet_grid(
      facet ~ Feature, 
      scales = "free_y",
      # switch = "y"
    ) +
    geom_boxplot(
      # fill = "white",
      color = "black",
      alpha = 0.7,
      width = 0.5,
      linewidth = 0.15,
      box.linetype = "solid",
      whisker.linetype = "solid",
      staplewidth = 0.25,
      # key_glyph = "rect",
      show.legend = FALSE,
      # outlier.shape = NA,
      outlier.shape = 21,
      outlier.color = "black", 
      outlier.fill = "white",  
      outlier.stroke = 0.15,
      outlier.size = 0.3
    ) +
    scale_fill_manual(
      values = colors,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      axis.title.x = element_blank(), 
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    labs(
      y = expression(paste(Log[2], " (Intensity)"))
    ) +
    stat_compare_means(
      method = c("t.test", "wilcox.test")[2], 
      label =  "p.signif",
      paired = F, 
      comparisons = combn(levels(P.d$df$Group), 2, simplify = F)
    ) +
    scale_y_continuous(
      expand = c(0.05, 0, 0.1, 0),
      breaks = pretty_breaks(n = 5)
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.var.box2 <- function(temp_res, facet, log=T, order=T, ratio=NULL, colors, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      facet = !!sym(facet)
    )
  
  if(order){
    P.d$median <- P.d$df %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarize(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(Median)) %>%
      mutate(
        Feature = factor(Feature, levels = unique(Feature))
      )
    
    P.d$df$Feature <- factor(P.d$df$Feature, levels = levels(P.d$median$Feature))
  }else{
    if(class(P.d$df$Feature) != "factor"){
      P.d$df <- P.d$df %>%
        mutate(
          Feature = factor(Feature, levels = unique(Feature))
        )
    }
  }

  if(log){
    P.d$df <- P.d$df %>%
      mutate(
        Value = log2(Value)
      )
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Feature, y = Value, fill = Group)) +
    geom_point(color = "white") +
    guides(
      fill = guide_legend( override.aes = list( color = "black", shape = 22, size = 3) )
    ) + 
    facet_grid(
      facet ~ .,
      scales = "free_y",
      # switch = "both"
    ) +
    geom_boxplot(
      # fill = "white",
      color = "black",
      alpha = 0.7,
      width = 0.5,
      linewidth = 0.15,
      box.linetype = "solid",
      whisker.linetype = "solid",
      staplewidth = 0.25,
      # key_glyph = "rect",
      show.legend = FALSE,
      # outlier.shape = NA,
      outlier.shape = 21,
      outlier.color = "black", 
      outlier.fill = "white",  
      outlier.stroke = 0.15,
      outlier.size = 0.3
    ) +
    scale_fill_manual(
      values = colors,
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 1, color = "black"),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    scale_y_continuous(
      expand = c(0.05, 0, 0.1, 0),
      breaks = pretty_breaks(n = 6)
    ) +
    labs(
      y = expression(paste(Log[2], " (Intensity)"))
    ) +
    stat_compare_means(
      method = c("t.test", "wilcox.test")[2],
      label =  c("p.signif", "p.format")[1],
      paired = F,
      size = 3.3,
    #   comparisons = combn(levels(P.d$df$Group), 2, simplify = F)
    )
  
  # print(P1)
  
  return(P1)
}


## surv ----
### cox ----
f.cox_uni <- function(df1, formula, p.cut=0.05, ph.cut=0.05,...){
  # data
  survdata <- df1
  
  formula <- as.formula(formula)
  
  var <- f.splitform(formula)$right
  resp <- f.splitform(formula)$left
  if(all("." %in% var)){
    var <- setdiff(names(df1), resp)
  }
  
  temp_res <- pblapply(setNames(var, var), function(.x){
    formula <- as.formula(paste("Surv(", resp[1], ",", resp[2], ") ~", .x))
    
    temp_res <- survdata %>% 
      f.coxph(
        formula = formula
      )
    
    temp <- temp_res$sm
    
    temp_x <- cox.zph(temp_res$fit)
    temp$chisq <- temp_x$table[1, "chisq"]
    temp$Ph    <- temp_x$table[1, "p"]
    
    temp
  })
  
  temp_x <- temp_res %>% bind_rows() %>%
    mutate(
      p.cut = p.cut,
      ph.cut = ph.cut
    )

  temp <- list(
    beta = subset(temp_x, P < p.cut) %>% rownames(.),
    ph = subset(temp_x, Ph > ph.cut) %>% rownames(.)
  ) %>%
    c( "all" = list(intersect(.$beta, .$ph)) )
  
  res <- list(
    cox = temp_x,
    filter = temp,
    Covar = temp$all
  )
  
  return(res)
}


f.logrank.p <- function(df1, formula, cutoff="Q2", p.cut=0.05, ...){
  # data
  survdata <- df1
  
  var <- f.splitform(formula)$right
  resp <- f.splitform(formula)$left
  
  temp_res <- pblapply(setNames(var, var), function(.x){
    formula <- as.formula(paste("Surv(", resp[1], ",", resp[2], ") ~", "strata"))
    
    temp_res <- survdata %>% 
      f.strata(
        y = .x,
        cutoff = cutoff,
      ) %>%
      f.logrank(
        formula = formula
      )
    
    temp <- temp_res$sm
    
    temp
  })
  
  temp_x <- temp_res %>% 
    imap_dfr(~{
      .x %>%
        mutate(Feature = .y) %>%
        remove_rownames() %>%
        tibble::column_to_rownames("Feature")
    }) %>%
    mutate(
      p.cut = p.cut
    )
  
  res <- list(
    logrank = temp_x,
    Covar = subset(temp_x, P < p.cut) %>% rownames(.)
  )
  
  return(res)
}


f.logrank.pp <- function(df1, p.cutoff="Q2", core.cutoff=NULL, formula, p.cut=0.05, ...){
  data <- df1
  formula <- formula
  
  temp_res <- f.func.p(
    data = data,
    func = f.logrank.p,
    formula = formula,
    p.cut = p.cut, 
    arg = "cutoff",
    p.arg = p.cutoff,
    core = core.cutoff,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "survival", "pbapply", "purrr")
  )
  
  return(temp_res)
}


f.cox_multi <- function(df1, formula, p.cut=0.05, ...){
  # data
  survdata <- df1
  
  temp_res <- survdata %>% f.coxph(formula = formula)
  
  temp_res$fit$param <- list(
    modeling = "multiCox",
    formula = formula
  )
  
  temp_res$Covar <- subset(temp_res$sm, P < p.cut) %>% rownames(.)
  
  return(temp_res)
}


f.cox_step <- function(df1, formula, step="both", ...){
  # data
  survdata <- df1
  
  set.seed(1234)
  temp <- coxph(
    formula = formula,
    data = survdata,
    x = T, y = T
  )
  
  # temp_fit <- MASS::stepAIC(
  #   object = temp,
  #   direction = step  # c("both", "backward", "forward")
  # )

  temp_fit <- stats::step(
    object = temp,
    direction = step  # c("both", "backward", "forward")
  )
  
  res <- list(
    fit = temp_fit,
    Covar = names(coef(temp_fit))
  )
  
  return(res)
}


f.cox_step.p <- function(df1, core.step=NULL, p.step="both", formula, ...){
  data <- df1
  formula <- formula
  
  temp_res <- f.func.p(
    data = data,
    func = f.cox_step,
    formula = formula,
    arg = "step",
    p.arg = p.step,
    core = core.step,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "purrr", "survival")
  )
  
  return(temp_res)
}


f.cox_boost <- function(df1, formula, nfold=10, ...){
  # data
  survdata <- df1
  
  X <- f.splitXy(survdata, formula)$X %>% as.matrix()
  y <- f.splitXy(survdata, formula)$y
  
  # var <- f.splitform(formula)$right
  resp <- f.splitform(formula)$left
  
  time = survdata[[resp[1]]]
  status = survdata[[resp[2]]]
  
  
  set.seed(1234)
  temp_x <- optimCoxBoostPenalty(
    time = time,
    status = status,
    x = X,
    trace = TRUE, 
    start.penalty = 500, 
    parallel = F
  )
  
  temp_y <- cv.CoxBoost(
    time = time,
    status = status,
    x = X,
    maxstepno = 500,
    K = nfold,
    type = "verweij",
    penalty = temp_x$penalty,
    parallel = F
  )

  temp_fit <- CoxBoost(
    time = time,
    status = status,
    x = X,
    stepno = temp_y$optimal.step,
    penalty = temp_x$penalty
  )

  temp_fit$call <- call(
    "CoxBoost",
    formula = formula,
    penalty = temp_x$penalty,
    stepno = temp_y$optimal.step
  )
  
  temp_fit$param <- list(
    modeling = "CoxBoost",
    formula = formula,
    penalty = temp_x$penalty,
    stepno = temp_y$optimal.step
  )

  Covar <- coef(temp_fit) %>%
    as.data.frame() %>% setNames("coef") %>%
    filter(coef != 0) %>%
    rownames()

  res <- list(
    tune = list(
      penalty = temp_x,
      optimal.step = temp_y
    ),
    fit = temp_fit,
    Covar = Covar
  )
  
  return(res)
}

### other ----
f.gbm <- function(df1, formula, nfold=10, n.cores=1, ...){
  # data
  survdata <- df1
  
  var <- f.splitform(formula)$right
  resp <- f.splitform(formula)$left
  if(all("." %in% var)){
    var <- setdiff(names(df1), resp)
  }
  
  formula = reformulate(termlabels = var, response = paste0("Surv(", resp[1], ",", resp[2], ")") )
  
  set.seed(1234)
  # tuning
  temp_tune <- gbm(
    formula = formula, 
    data = survdata, 
    distribution = 'coxph', 
    n.trees = 10000, 
    interaction.depth = 3,
    n.minobsinnode = 10, 
    shrinkage = 0.001,  
    cv.folds = nfold,  
    n.cores = n.cores
  )
  
  set.seed(1234)
  temp_fit <- gbm(
    formula = formula,
    data = survdata,
    distribution = 'coxph',
    n.trees = which.min(temp_tune$cv.error),
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = nfold,
    n.cores = n.cores
  )
  
  temp_fit$call <- call(
    "gbm",
    formula = formula,
    data = quote(newdata),
    distribution = "coxph",
    n.trees = which.min(temp_tune$cv.error),
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = nfold,
    n.cores = n.cores
  )
  
  temp_fit$param <- list(
    modeling = "gbm",
    formula = formula,
    n.trees = which.min(temp_tune$cv.error),
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = nfold
  )
  
  res <- list(
    tune = temp_tune,
    fit = temp_fit
  )
  
  return(res)
}


f.rsf <- function(df1, formula, var.select=NULL, ...){
  # data
  survdata <- df1
  var <- f.splitform(formula)$right
  
  set.seed(1234)
  temp_fit <- rfsrc(
    formula = formula, 
    data = survdata,
    # ntree = 1000,
    # nodesize = rf_nodesize,
    splitrule = 'logrank',
    importance = T,
    save.memory = T,
    proximity = T,
    forest = T,
    seed = -123
  )
  
  temp_fit$param <- list(
    modeling = "rsf",
    formula = formula
  )
  
  if(!is.null(var.select)){
    temp_res <- var.select(
      object = temp_fit, 
      # method = c("md", "vh", "vh.vimp")[1],
      method = var.select,
      conservative = "high"
    )
    
    res <- list(
      fit = temp_fit,
      var.select = temp_res,
      Covar = temp_res$topvars
    )
    
  }else{
    res <- list(
      fit = temp_fit,
      Covar = var
    )
  }
  
  return(res)
}


f.rfs_select <- function(temp_res, ...){
  
  temp <- temp_res$varselect %>%
    tibble::rownames_to_column("Feature") %>%
    rename(
      "Depth" = depth,
      "Mean Decrease Gini" = vimp
    ) %>%
    pivot_longer(
      cols = -Feature,
      names_to = "Metric",
      values_to = "Importance"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("Depth", "Mean Decrease Gini")
      )
    ) %>%
    arrange(Metric, Importance) %>%
    mutate(
      Status = if_else(Feature %in% temp_res$topvars, "Selected", "Unselected")
    ) %>%
    mutate(
      Status = factor(
        Status,
        levels = c("Selected", "Unselected")
      )
    )
  
  return(temp)
}


f.plot.rsf_select <- function(temp_res, color, colors, ratio, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    group_by(Metric) %>%
    mutate(Rank = rank(dplyr::desc(Importance), ties.method = "min")) %>%
    ungroup() %>%
    filter(
      Metric == "Depth"
    ) %>%
    arrange(desc(Importance)) %>%
    f.grouping("Feature", "z") %>%
    rename(
      color = !!sym(color)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Feature, y = Importance, color = color)) +
    # facet_wrap(
    #   .~ Metric,
    #   nrow = 1,
    #   scales = "free_y",
    #   strip.position = "top"
    # ) +
    # geom_col(
    #   # color = "firebrick",
    #   size = 1,
    #   alpha = 1
    # ) +
    geom_segment(
      aes(xend = Feature, yend = 0),
      linewidth = 1.2,
      alpha = 0.7,
      key_glyph = c("smooth", "rect")[2]
    ) +
    geom_point(
      size = 2.5, 
      alpha = 1,
      fill = "white", 
      stroke = 0.5,
      show.legend = FALSE
    ) +
    coord_flip() +
    # geom_segment(
    #   aes(x = Rank, xend = Rank, y = Value, yend = 0, color = Metric),
    #   show.legend = F
    # ) +
    scale_color_manual(
      values = colors
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.1, 0)
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      # legend.position = "none",
      # axis.title.y = element_blank(),
    ) +
    labs(
      x = "Feature",
      y = "Depth"
    )
    
  # print(P1)
  
  return(P1)
}



# f.imp.rsf <- function(fit, ...){
#   
#   temp_res <- fit$importance %>%
#     as.data.frame() %>%
#     setNames(., "Importance") %>%
#     tibble::rownames_to_column("Feature") %>%
#     mutate(
#       Metric = 
#     )
#   
# }


## repeatedcv ----
f.prog_part <- function(df1, df2=NULL, part.seed, part.p, part.y, diag2prog=NULL, prog, prog2=NULL, timeUnits="d", ...){
  # inT
  if(part.p > 1){ part.p <- part.p/nrow(df1) }
  set.seed(part.seed)  
  inT <- caret::createDataPartition(y = df1[[part.y]], times=1, p = part.p, list = F)
  
  # part
  temp_set <- df1 %>%
    f.part(
      inT = inT
    )
  
  if(!is.null(df2)){
    temp_set <- temp_set %>% setNames(ifelse(names(.) == "all", "c1", names(.)))
    temp_set$c2 <- df2
    lev_set <- c("train", "valid", "c1", "c2")
  }else{
    lev_set <- c("train", "valid", "all")
  }

  # prog
  if(!is.null(diag2prog)){
    left <- f.splitform(diag2prog)$left
    positive <- tail(levels(df1[[left]]), n=1)
    temp <- temp_set %>%
      map(~.x %>% filter(!!sym(left) == positive) )
  }else{
    temp <- temp_set
  }

  surv <- f.splitform(as.formula(paste0(prog, "~.")))$left
  # status distribution
  is_status <- all(unlist(lapply(temp, function(df){
    unique_vals <- unique(df[[surv[2]]])
    if(all(c(0, 1) %in% unique_vals)){ return(TRUE) }else{ return(FALSE) }
  })))
  
  if(!is_status){
    return(temp)
  }
  
  if(is.null(prog2)){
    p.prog <- setNames(prog, prog)
  }else{
    p.prog <- c(prog, prog2) %>% setNames(., .)
  }
  
  # logrank test
  temp_x <- p.prog %>%
    map_dfr(function(pp){
      lev_set %>% f.comb(k = c(2:length(lev_set))) %>%
        imap_dfr(function(set, name){
          temp[set] %>%
            imap_dfr(~{ .x %>% mutate( Set = .y ) }) %>%
            mutate(
              Set = factor(Set, levels = set)
            ) %>%
            f.logrank(
              formula = as.formula(paste(pp, "~", "Set")),
            ) %>%
            {.$sm} %>%
            dplyr::select("P.set"="P") %>%
            mutate(
              comb = name,
              prog = pp,
              seed = part.seed,
              inT = paste0(sort(inT), collapse = " ")
            )
        })
    })
  
  # km
  temp_y <- temp %>%
    imap_dfr(~{ .x %>% mutate( Set = .y ) }) %>%
    mutate(
      Set = factor(Set, levels = lev_set)
    ) %>%
    f.survfit(
      formula = as.formula(paste(prog, "~", "Set"))
    ) %>%
    f.survfit.res(
      timeUnits = timeUnits
    )
  
  if(!is.null(prog2)){
    temp_z <- temp %>%
      imap_dfr(~{ .x %>% mutate( Set = .y ) }) %>%
      mutate(
        Set = factor(Set, levels = lev_set)
      ) %>%
      f.survfit(
        formula = as.formula(paste(prog2, "~", "Set"))
      ) %>%
      f.survfit.res(
        timeUnits = timeUnits
      )
  }
  
  if(!is.null(prog2)){
    res <- list(
      logrank = temp_x,
      km = temp_y,
      km2 = temp_z
    )
  }else{
    res <- list(
      logrank = temp_x,
      km = temp_y
    )
  }

  return(res)
}


f.prog_part.p <- function(p.part.seed, core.part.seed, df1, df2=NULL, part.p, part.y, diag2prog=NULL, prog, prog2=NULL, timeUnits="d", ...){
  data <- df1
  df2 <- df2
  formula <- formula
  
  temp_res <- f.func.p(
    data = data,
    func = f.prog_part,
    df2 = df2,
    part.p = part.p,
    part.y = part.y,
    diag2prog = diag2prog, 
    prog = prog,
    prog2 = prog2,
    timeUnits = timeUnits,
    arg = "part.seed",
    p.arg = p.part.seed,
    core = core.part.seed,
    globals = c(
      f.find_func("f."), "param"
    ),
    packages = c("pbapply", "dplyr", "purrr", "caret", "Matrix", "glmnet", "survival", "survminer", "Boruta", "pROC", "randomForestSRC", "CoxBoost", "timeROC", "riskRegression", "pec")
  )
  
  return(temp_res)
}


f.loocv_inT <- function(
    df1, df2=NULL, var, resp, resampling="none", inT, 
    seq=NULL, No=NULL, comb = T, k1=1, k2=NULL,
    diag=T, diag2prog=NULL, var_prog=NULL, prog=NULL, prog_fit="cox", seqs=NULL, ks1=1, ks2=NULL, timeUnits="d", combs=T, ...
){
  
  # inT
  # part.seed, part.p
  # if(part.p > 1){ part.p <- part.p/nrow(df1) }
  # set.seed(part.seed)  
  # inT <- caret::createDataPartition(y = df1[[resp]], times=1, p = part.p, list = F)
  
  # part
  temp_set <- df1 %>%
    f.part(
      inT = inT
    )
  
  if(!is.null(df2)){
    temp_set <- temp_set %>% setNames(ifelse(names(.) == "all", "c1", names(.)))
    temp_set$c2 <- df2
    lev_set <- c("train", "valid", "c1", "c2")
  }else{
    lev_set <- c("train", "valid", "all")
  }
  
  if(resampling != "none"){
    temp_set$train.rs <- temp_set$train %>%
      f.resampling(
        method = resampling,
        formula = reformulate(termlabels = ".", response = resp), 
        # seed = seed, 
        # over_ratio = NULL
      )
    traindata <- temp_set$train.rs
  }else{
    traindata <- temp_set$train
  }
  
  # temp_set <- temp_set %>%
  #   map(~{
  #     .x %>%
  #       filter(
  #         RFS_status == 0 | RFS_time < 365
  #       )
  #   })
  
  # screen
  if(!is.null(seq)){
    list_temp <- traindata %>%
      f.screen_seq(
        method_seq = seq,
        formula = reformulate(termlabels = var, response = resp),
        # param
        auc.cut = 0.7,    # auc, aucor
        # auc.cut = NULL,
        cor.cut = 0.5,    # aucor
        cor.meth = c("pearson", "kendall", "spearman")[3], # aucor
        seed = 005,       # rfe
        p.seed = 1:200,   # glmnet
        core.seed = 50,    # glmnet
        nfold = 10,       # glmnet, rfe
        index = c("var.min", "var.1se")[1] # glmnet
      )
    temp_var <- list_temp[[length(list_temp)]]$Covar
  }else{
    temp_var <- var
  }
  
  # feature
  if(length(temp_var) == 0){
    message("diag: refining feature is NULL !")
    return(NULL)
  }
  if(!is.null(No)){
    if(length(temp_var) > No){ temp_var <- head(temp_var, No) }
  }
  
  # diag: evaluate
  if(is.null(k2)){
    if(length(temp_var) <= k1){
      k = length(temp_var)
    }else{
      k <- seq(k1, length(temp_var), 1)
    }
  }else{
    k <- seq(k1, k2, 1)
  }
  
  if(diag){
    if(comb){
      temp_x <- temp_var %>% f.comb(k = k) %>%
        map_dfr(function(fea){
          if(length(fea) == 1){
            temp_fit = NULL
            formula = reformulate(termlabels = fea, response = resp)
            
            temp_eval <- f.eval_temp( temp_set, temp_fit, cutoff=c(0.5, "yd")[2], formula=formula, direction="auto" ) %>%
              pivot_longer(
                cols = c("AUC", "Sensitivity", "Specificity"),
                names_to = "Metric",
                values_to = "Value"
              ) %>%
              mutate(
                var = paste0(sort(fea), collapse = " "),
                length_var = length(fea),
                # seed = part.seed,
                inT = paste0(sort(inT), collapse = " ")
              )
          }else{
            temp_fit <- traindata %>%
              f.modeling.rf(
                formula = reformulate(termlabels = fea, response = resp),
                mtry = 1,
                # ntree = 500
              )
            
            temp_eval <- f.eval_temp( temp_set, temp_fit, cutoff=c(0.5, "yd")[1] ) %>%
              pivot_longer(
                cols = c("AUC", "Sensitivity", "Specificity"),
                names_to = "Metric",
                values_to = "Value"
              ) %>%
              mutate(
                var = paste0(sort(fea), collapse = " "),
                length_var = length(fea),
                # seed = part.seed,
                inT = paste0(sort(inT), collapse = " ")
              )
          }
          # reture
          temp_eval
        })
    }else{
      if(length(temp_var) == 1){
        temp_fit = NULL
        formula = reformulate(termlabels = temp_var, response = resp)
        
        temp_x <- f.eval_temp( temp_set, temp_fit, cutoff=c(0.5, "yd")[2], formula=formula, direction="auto" ) %>%
          pivot_longer(
            cols = c("AUC", "Sensitivity", "Specificity"),
            names_to = "Metric",
            values_to = "Value"
          ) %>%
          mutate(
            var = paste0(sort(temp_var), collapse = " "),
            length_var = length(temp_var),
            # seed = part.seed,
            inT = paste0(sort(inT), collapse = " ")
          )
      }else{
        temp_fit <- traindata %>%
          f.modeling.rf(
            formula = reformulate(termlabels = temp_var, response = resp),
            mtry = 1,
            # ntree = 500
          )
        
        temp_x <- f.eval_temp( temp_set, temp_fit, cutoff=c(0.5, "yd")[1] ) %>%
          pivot_longer(
            cols = c("AUC", "Sensitivity", "Specificity"),
            names_to = "Metric",
            values_to = "Value"
          ) %>%
          mutate(
            var = paste0(sort(temp_var), collapse = " "),
            length_var = length(temp_var),
            # seed = part.seed,
            inT = paste0(sort(inT), collapse = " ")
          )
      }
    }
  }else{
    temp_x <- NULL
  }
  
  # prog: 
  if(is.null(prog)){
    res <- list(
      diag = temp_x
    )
    return(res)
  }
  
  # data
  if(!is.null(diag2prog)){
    left <- f.splitform(diag2prog)$left
    positive <- tail(levels(df1[[left]]), n=1)
    temp_set <- temp_set %>%
      map(~.x %>% filter(!!sym(left) == positive) )
  }

  surv <- f.splitform(as.formula(paste0(prog, "~.")))$left
  is_status <- all(unlist(lapply(temp_set, function(df){
    unique_vals <- unique(df[[surv[2]]])
    if(all(c(0, 1) %in% unique_vals)){ return(TRUE) }else{ return(FALSE) }
  }))) %>%
    paste0(., collapse = " ")

  P.set <- temp_set[lev_set] %>%
    imap_dfr(~{ .x %>% mutate( Set = .y ) }) %>%
    mutate(
      Set = factor(Set, levels = lev_set)
    ) %>%
    f.logrank(
      formula = as.formula(paste(prog, "~", "Set")),
    ) %>%
    {.$sm} %>% pull(P)

  if(is.null(var_prog)){
    var_prog <- var
  }

  # screen → feature
  if(!is.null(seqs)){
    list_temp <- temp_set$train %>%
      f.screen_seq(
        method_seq = seqs,
        formula = reformulate(termlabels = var_prog, response = prog),
        # param
        p.cut = 0.05,    # unicox multicox logrank
        ph.cut = 0.05,   # unicox
        var.select = c("md", "vh", "vh.vimp")[1],  # rsf
        p.seed = 1:200,   # glmnet
        core.seed = T,    # glmnet
        nfold = 10,       # glmnet, gbm
        index = c("var.min", "var.1se")[1] # glmnet
      )
    temp_var <- list_temp[[length(list_temp)]]$Covar
  }else{
    temp_var <- var_prog
  }

  # feature
  if(length(temp_var) == 0){
    message("prog: refining feature is NULL !")
    return(NULL)
  }

  # prog: evaluate
  if(is.null(ks2)){
    if(length(temp_var) <= ks1){
      k = length(temp_var)
    }else{
      k <- seq(ks1, length(temp_var), 1)
    }
  }else{
    k <- seq(ks1, ks2, 1)
  }

  if(combs){
    temp_y <- temp_var %>% f.comb(k = k) %>%
      map(function(fea){
        if(prog_fit == "cox"){
          temp_fit <- temp_set$train %>%
            f.cox_multi(
              formula = reformulate(termlabels = fea, response = prog)
            ) %>% {.$fit}
        }else{
          temp_fit <- temp_set$train %>%
            f.rsf(
              formula = reformulate(termlabels = fea, response = prog),
              var.select = NULL
            ) %>% {.$fit}
        }
        temp_km <- temp_set %>%
          imap(function(x, y){
            f.riskscore(temp_fit, newdata = x, center=F) %>%
              f.strata(
                y = "riskscore",
                cutoff = "Q2"
              ) %>%
              f.survfit(
                formula = as.formula(paste(prog, "~", "strata"))
              ) %>%
              f.survfit.res(
                timeUnits = timeUnits
              ) %>%
              map(~{
                .x %>% mutate( Set = y )
              })
          }) %>%
          purrr::transpose() %>%
          map(~{
            .x %>% bind_rows()
          })

        trend <- temp_km$median %>%
          filter(!Set %in% c("all", "c1")) %>%
          dplyr::select(c("median", "strata", "Set")) %>%
          mutate(median = ifelse(is.na(median), 1000, median)) %>%
          pivot_wider(names_from = strata, values_from = median) %>%
          mutate(ratio = High / Low) %>%
          pull(ratio)

        P <- temp_km$df %>% group_by(Set) %>% slice(1) %>% ungroup() %>% pull(P)

        if(all(P < 0.05) & all(trend < 1)){
          temp_eval <- f.evals_r2e(
            l1 = temp_set,
            fit = temp_fit,
            index = "train",
            center = F,
            cutoff = c("Q2", "best")[1],
            weight_l1 = NULL,
            weight_formula = NULL,
            timeUnits = timeUnits,
            timeCut = NULL,
            times = NULL
          )

          temp_metric <- temp_eval %>%
            f.eval.res(index = "metric") %>%
            filter(
              Metric %in% c("P.logrank", "C Index") | grepl("AUC", Metric)
            ) %>%
            # pivot_wider( names_from = Metric, values_from = Value ) %>%
            mutate(
              var = paste0(sort(fea), collapse = " "),
              length_var = length(fea),
              # seed = part.seed,
              inT = paste0(sort(inT), collapse = " "),
              is_status = is_status,
              P.set = P.set
            )

        }else{
          temp_metric <- NULL
          temp_km <- NULL
        }
        # reture
        res <- list(
          sm = temp_metric,
          km = temp_km
        )
      })

    temp_km <- temp_y %>% map(~{.x$km})
    temp_y <- temp_y %>% map_dfr(~{.x$sm}) %>% keep(~ !is.null(.x))
  }else{

    if(prog_fit == "cox"){
      temp_fit <- temp_set$train %>%
        f.cox_multi(
          formula = reformulate(termlabels = temp_var, response = prog)
        ) %>% {.$fit}
    }else{
      temp_fit <- temp_set$train %>%
        f.rsf(
          formula = reformulate(termlabels = temp_var, response = prog),
          var.select = NULL
        ) %>% {.$fit}
    }

    temp_eval <- f.evals_r2e(
      l1 = temp_set,
      fit = temp_fit,
      index = "train",
      center = F,
      cutoff = c("Q2", "best")[1],
      weight_l1 = NULL,
      weight_formula = NULL,
      timeUnits = timeUnits,
      timeCut = NULL,
      times = NULL
    )

    temp_y <- temp_eval %>%
      f.eval.res(index = "metric") %>%
      filter(
        Metric %in% c("P.logrank", "C Index") | grepl("AUC", Metric)
      ) %>%
      # pivot_wider( names_from = Metric, values_from = Value ) %>%
      mutate(
        var = paste0(sort(temp_var), collapse = " "),
        length_var = length(temp_var),
        # seed = part.seed,
        inT = paste0(sort(inT), collapse = " "),
        is_status = is_status,
        P.set = P.set
      )

    temp_km <- temp_eval %>%
      {
        list(
          df = f.eval.res(., index = "km_df"),
          timeCut = f.eval.res(., index = "km_timeCut"),
          median = f.eval.res(., index = "km_median")
        )
      }
  }

  # return
  temp_res <- list(
    diag = temp_x,
    prog = temp_y,
    km = temp_km
  )

  return(temp_res)
}


f.loocv_inT.p <- function(
    p.inT, core.inT,
    df1, df2=NULL, var, resp, resampling="none", 
    seq=NULL, No=NULL, comb = T, k1=1, k2=NULL,
    diag=T, diag2prog=NULL, var_prog=NULL, prog=NULL, prog_fit="cox", seqs=NULL, ks1=3, ks2=NULL, timeUnits="d", combs=T, ...
){
  
  data <- df1
  df2 <- df2
  
  temp_res <- f.func.p(
    data = data,
    func = f.loocv_inT,
    df2 = df2,
    var = var,
    resp = resp,
    resampling = resampling,
    seq = seq,
    No = No,
    comb = comb,
    k1 = k1,
    k2 = k2,
    diag = diag,
    prog = prog, 
    diag2prog = diag2prog,
    var_prog = var_prog,
    prog_fit = prog_fit,
    seqs = seqs,
    ks1 = ks1,
    ks2 = ks2,
    timeUnits = timeUnits,
    combs = combs,
    arg = "inT",
    p.arg = p.inT,
    core = core.inT,
    globals = c( f.find_func("f."), "param" ),
    packages = c("pbapply", "dplyr", "stringr", "tidyr", "purrr", "caret", "Matrix", "glmnet", "survival", "survminer", "Boruta", "pROC", "randomForestSRC", "randomForest", "CoxBoost", "timeROC", "riskRegression", "pec")
  )
  
  return(temp_res)
}


f.cv <- function(
    part.seed, part.y, nflod,
    df1, df2=NULL, var, resp, resampling="none", 
    seq=NULL, No=NULL, comb = T, k1=1, k2=NULL, 
    diag2prog=NULL, var_prog=NULL, prog=NULL, prog_fit="cox", seqs=NULL, ks1=3, ks2=NULL, timeUnits="d", combs=T, ...
  ){
  
  set.seed(part.seed)
  folds <- createFolds(df1[[part.y]], k = nflod, list = T, returnTrain = FALSE)
  
  temp_res <- list()
  
  for(k in 1:nflod){
    # inT <- folds[[k]]
    inT <- setdiff(1:nrow(df1), folds[[k]])
    
    temp_res[[k]] <- f.loocv_inT(
      df1 = df1,
      df2 = df2,
      var = var, 
      resp = resp, 
      resampling = resampling, 
      inT = inT,
      seq = seq, 
      No = No, 
      comb = comb, 
      k1 = k1,
      k2 = k2,
      diag2prog = NULL,  
      var_prog = var_prog, 
      prog = prog, 
      prog_fit = prog_fit,
      seqs = seqs, 
      ks1 = ks1,
      ks2 = ks2,
      timeUnits = timeUnits, 
      combs = combs
    )
    
  }

  # diag
  fold_x <- temp_res %>%
    imap_dfr(~{
      .x$diag %>% filter( Set == "valid" ) %>%
        mutate( folds = paste("fold ", .y) )
    }) %>%
    dplyr::select(-Set) %>%
    dplyr::select(folds, everything())
  
  mean_x <- fold_x %>%
    group_by(var, Metric) %>%
    reframe(
      Value = mean(Value, na.rm = TRUE),
      length_var = dplyr::first(length_var),
    ) %>%
    ungroup()
  # pivot_wider( names_from = Metric, values_from = Value ) %>%
  # dplyr::select(c("folds", "positive", "AUC", "Sensitivity", "Specificity"), everything())
  
  if(is.null(prog)){
    res <- list(
      diag = list(
        folds = fold_x,
        mean = mean_x
      )
    )
    return(res)
  }
  
  # prog
  fold_y <- temp_res %>%
    imap_dfr(~{
      .x$prog %>% filter( Set == "valid" ) %>%
        mutate( folds = paste("fold ", .y) )
    }) %>%
    dplyr::select(-Set) %>%
    dplyr::select(folds, everything())
  
  mean_y <- fold_y %>%
    group_by(var, Metric) %>%
    reframe(
      Value = mean(Value, na.rm = TRUE),
      length_var = dplyr::first(length_var),
      is_status = dplyr::first(is_status),
      P.set = dplyr::first(P.set)
    ) %>%
    ungroup()
  # pivot_wider( names_from = Metric, values_from = Value ) %>%
  # dplyr::select(c("folds", "P.logrank", "C Index"), everything()) %>%
  # dplyr::select(-c("var", "length_var", "inT"), "var", "length_var", "inT")
  
  res <- list(
    diag = list(
      folds = fold_x,
      mean = mean_x
    ),
    prog = list(
      folds = fold_y,
      mean = mean_y
    )
  )
  
  return(res)
}


f.cv.p <- function(
    p.part.seed, core.part.seed, part.y, nflod,
    df1, df2=NULL, var, resp, resampling="none", 
    seq=NULL, No=NULL, comb = T, k1=1, k2=NULL,
    diag2prog=NULL, var_prog=NULL, prog=NULL, prog_fit="cox", seqs=NULL, ks1=3, ks2=NULL, timeUnits="d", combs=T, ...
  ){
  data <- df1
  df2 <- df2
  
  temp_res <- f.func.p(
    data = data,
    func = f.cv,
    df2 = df2,
    var = var,
    resp = resp,
    resampling = resampling,
    seq = seq,
    No = No,
    comb = comb,
    k1 = k1,
    k2 = k2,
    diag2prog = diag2prog,
    var_prog = var_prog,
    prog = prog, 
    prog_fit = prog_fit,
    seqs = seqs,
    ks1 = ks1,
    ks2 = ks2,
    timeUnits = timeUnits,
    combs = combs,
    part.y = part.y, 
    nflod = nflod,
    arg = "part.seed",
    p.arg = p.part.seed,
    core = core.part.seed,
    globals = c( f.find_func("f."), "param" ),
    packages = c("dplyr", "stringr", "tidyr", "purrr", "caret", "Matrix", "glmnet", "survival", "survminer", "Boruta", "pROC", "randomForestSRC", "randomForest", "CoxBoost", "timeROC", "riskRegression", "pec")
  )
  
  return(temp_res)
}


f.plot.comb <- function(temp_res, x, y, colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      "x" = !!sym(x),
      "y" = !!sym(y)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = x, y = y)) +
    facet_wrap(
      ~ Metric, 
      # scales = "free_y"
    ) +
    stat_summary(
      fun = "mean", geom = "bar",
      fill = "gray80",
      color = "black",
      width = 0.8
    ) +
    stat_summary(
      fun.data = "mean_se", 
      geom = "errorbar", 
      width = 0.2
    ) +
    geom_jitter(
      aes(fill = Metric),
      alpha = 0.7,
      color = "black",
      shape = 21,
      width = 0.2,
      size = 1.5
    ) +
    geom_hline(yintercept = 0.8, color = "burlywood4", linetype = "twodash", linewidth = 0.4) + 
    scale_fill_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      panel.grid = element_blank(),
      # axis.text.x = element_text(angle = 45, size = 10, hjust = 1, color = "black"),
      # axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.05, 0),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    scale_x_continuous(
      breaks = seq(floor(min(P.d$df$x)), ceiling(max(P.d$df$x)), by = 1)
    ) +
    labs(
      x = "Order of Combination"
    )
  
  # print(P1)
  
  return(P1)
}



f.plot.prog <- function(temp_res, label, colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    mutate(
      temp = paste0(var, length_var)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = length_var, y = Value, fill = Set)) +
    facet_wrap(
      ~ Metric, 
      # scales = "free_y"
    ) +
    stat_summary(
      fun = "mean", geom = "bar",
      # fill = "gray80",
      color = "black",
      width = 0.6,
      alpha = 0.6,
      position = position_dodge(width = 0.6)
    ) +
    geom_point(
      data = subset(P.d$df, var != label),
      aes(group = Set),
      alpha = 0.7,
      fill = "gray",
      color = "black", 
      shape = 21, size = 1.5,
      position = position_dodge(width = 0.6)
    ) +
    # geom_line(
    #   data = P.d$df,
    #   aes(x = length_var, y = Value, group = temp),
    #   position = position_dodge(width = 0.8),
    #   color = "black", linewidth = 0.5
    # ) +
    # geom_jitter(
    # # geom_point(
    #   aes(group = Set),
    #   alpha = 0.7,
    #   fill = "gray",
    #   color = "black",
    #   shape = 21,
    #   # width = 0.2,
    #   size = 1.5,
    #   position = position_jitterdodge(
    #     dodge.width = 0.5,
    #     jitter.width = 0.15,
    #     jitter.height = 0
    #   )
    # ) +
    stat_summary(
      fun.data = "mean_se", 
      geom = "errorbar", 
      width = 0.2,
      position = position_dodge(width = 0.6)
    ) +
    geom_hline(yintercept = 0.8, color = "burlywood4", linetype = "twodash", linewidth = 0.4) + 
    scale_fill_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      # axis.text.x = element_text(angle = 45, size = 10, hjust = 1, color = "black"),
      # axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 10),
      legend.title = element_blank(),
      # legend.position = "top",
      legend.key.size = unit(0.5, "cm"),
      legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit="pt")
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.05, 0),
      limits = c(0, NA),
      breaks = seq(0, 1, by = 0.2)
    ) +
    labs(
      x = "Order of Combination"
    ) +
    geom_point(
      data = subset(P.d$df, var == label),
      aes(group = Set, fill = Set),
      alpha = 1,
      # fill = "gray",
      color = "black", 
      shape = 21, size = 1.5,
      position = position_dodge(width = 0.6)
    ) +
    stat_summary(
      data = subset(P.d$df, var == label),
      aes(label = var, group = Set), 
      fun = "mean",
      geom = "label_repel", 
      position = position_dodge(width = 0.6), 
      nudge_x = 0.2,
      nudge_y = -0.5,
      fill = NA,
      alpha = 1, size = 2, min.segment.length = 0,
      direction = "both", 
      segment.size = 0.2, segment.color = "black"
    )
    # geom_label_repel(
    #   data = subset(P.d$df, var == label) ,
    #   aes(label = var),
    #   nudge_x = 0.1,
    #   nudge_y = -0.1,
    #   position = position_dodge(width = 0.6),
    #   alpha = 1,
    #   fill = "white",
    #   size = 2,
    #   min.segment.length = 0,
    #   direction = "y",
    #   segment.size = 0.2,
    #   # segment.linetype = 6,
    #   segment.color = "black",
    # )
  
  # print(P1)
  
  return(P1)
}






