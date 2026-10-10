# diagnosis ----
## scores ----
f.score_aMAP <- function(Age, Gender, TBIL, ALB, PLT){
  if(any(is.na(c(Age, Gender, TBIL, ALB, PLT)))) {
    return(NA_real_)
  }
  # Gender (Male：1, Female：0)
  Gender <- ifelse(Gender == "Male", 1, 0)
  
  Age <- as.numeric(Age)
  TBIL <- as.numeric(TBIL)
  ALB <- as.numeric(ALB)
  PLT <- as.numeric(PLT)
  
  # score
  score <-  ((0.06*Age + 0.89*Gender + 0.48*((log10(TBIL)*0.66)+(-0.085*ALB))-0.01*PLT)+7.4)/14.77*100
  
  return(score)
}


f.score_ASAP <- function(Age, Gender, AFP, DCP) {
  if(any(is.na(c(Age, Gender, AFP, DCP)))) {
    return(NA_real_)
  }
  # Gender，0 = Male，1 = Female
  Gender <- ifelse(Gender == "Male", 0, 1)
  
  Age <- as.numeric(Age)
  AFP <- as.numeric(AFP)
  DCP <- as.numeric(DCP)
  
  # ln(AFP)和ln(DCP)
  ln_afp <- log(AFP)
  ln_dcp <- log(DCP)
  # ln(p/(1-p))
  log_odds <- -7.57711770 + 0.04666357 * Age - 0.57611693 * Gender + 0.42243533 * ln_afp + 1.10518910 * ln_dcp
  # probability/score
  p <- exp(log_odds) / (1 + exp(log_odds))
  
  return(p)
}


## tuning ----
f.tune_ctrl <- function(metric){
  
  tune <- list()
  tune$metric <- as.character(metric)
  
  set.seed(1234)
  tune$ctrl <- if(tune$metric %in% c("ROC", "Sens", "Spec")){
    tune$ctrl <- trainControl(
      method = "repeatedcv",
      number = 10,          
      repeats = 1,          
      classProbs = TRUE,    
      # sampling = "up",    
      summaryFunction = twoClassSummary,
    )
  }else{
    tune$ctrl <- trainControl(
      method = "repeatedcv",
      number = 10,          
      repeats = 1,          
      classProbs = TRUE,    
      # sampling = "up",    
    )
  }
  
  return(tune)
}


f.tuning_rf <- function(df1, formula, ctrl, ...){
  resp <- f.splitform(formula)$left
  lev <- levels(df1[[resp]])
  var <- f.splitform(formula)$right
  if(all("." %in% var)){
    var <- setdiff(names(df1), resp)
  }
  # data
  traindata <- df1
  # param
  rf_grid <- expand.grid(
    mtry = 1:(length(var)-1)
  )
  # tuning
  set.seed(12345)
  rf_tune <- train(
    form = formula, #Group ~ .,  
    data = traindata,         
    method = "rf",  
    trControl = ctrl$ctrl, 
    metric = ctrl$metric, 
    tuneGrid = rf_grid,  
    importance = TRUE 
  )
  # plot data
  rf_tune$sm <- cbind(
    rf_tune$results %>%
      dplyr::select(mtry, !matches("SD")) %>%
      pivot_longer(
        cols = -c("mtry"), #everything(), 
        names_to = "Metric", 
        values_to = "Mean"  
      ),
    rf_tune$results %>%
      dplyr::select(mtry, matches("SD")) %>%
      pivot_longer(
        cols = -c("mtry"), #everything(), 
        names_to = "Metric", 
        values_to = "SD"  
      ) %>%
      dplyr::select(-Metric, -mtry)
  ) %>%
    mutate(
      Metric = case_when(
        Metric == "ROC" ~ "AUC",
        Metric == "Sens" ~ "Sensitivity",
        Metric == "Spec" ~ "Specificity",
        TRUE ~ Metric
      )
    )
  
  rf_tune$param <- list(
    modeling = "rf",
    formula = formula,
    lev = lev,
    mtry = rf_tune$bestTune$mtry
  )
  
  # class(rf_tune) <- unique(c("rf", class(rf_tune)))
  
  return(rf_tune)
}


f.tuning_xgb <- function(df1, formula, ctrl, ...){
  resp <- f.splitform(formula)$left
  lev <- levels(df1[[resp]])
  # data
  traindata <- df1
  # auto tuning
  for(y in 1:3){
    # grid range
    obj_xgb <- c("n_e", "g_d", "jz")[y]
    i <- switch(
      obj_xgb, 
      "n_e" = c("2", "0", "0"),  # 0/default，1/best，2/tuning
      "g_d" = c("1", "2", "0"), 
      "jz" = c("1", "1", "2")
    )
    # 1: nrounds, eta，
    # 2: gamma, max_depth
    # 3: others
    xgb_grid <- expand.grid(
      nrounds = switch(i[1], "2" = c(50, 100, 150), "1" = as.numeric(xgb_tune$bestTune[1]), "0" = 100),
      # default/0.3
      eta = switch(i[1], "2" = c(0.01, 0.05, 0.1, 0.3), "1" = as.numeric(xgb_tune$bestTune[3]), "0" = 0.1),
      # default/0
      gamma = switch(i[2], "2" = c(0, 0.5, 1), "1" = as.numeric(xgb_tune$bestTune[4]), "0" = 0), 
      # default/6
      max_depth = switch(i[2], "2" = c(3, 4, 5, 6), "1" = as.numeric(xgb_tune$bestTune[2]), "0" = 6),
      # default/1
      subsample = switch(i[3], "2" = c(0.6, 0.8, 1), "1" = as.numeric(xgb_tune$bestTune[7]), "0" = 1),
      # default/1
      colsample_bytree = switch(i[3], "2" = c(0.6, 0.8, 1), "1" = as.numeric(xgb_tune$bestTune[5]), "0" = 1),
      min_child_weight = 1
    )
    
    rm(i)
    
    set.seed(12346)
    # tuning
    xgb_tune <- train(
      form = formula, #Group ~ ., 
      data = traindata,         
      method = "xgbTree",
      trControl = ctrl$ctrl,
      metric = ctrl$metric,
      tuneGrid = xgb_grid
    )
    
    # plot
    # plot(xgb_tune)
    
    rm(list = c("y", "obj_xgb"))
  }
  
  xgb_tune$param <- list(
    modeling = "xgb",
    formula = formula,
    lev = lev,
    # important
    max_depth = as.numeric(xgb_tune$bestTune[2]),    # tuning
    gamma = as.numeric(xgb_tune$bestTune[4]),        # tuning
    # normal
    eta = as.numeric(xgb_tune$bestTune[3]),
    colsample_bytree = as.numeric(xgb_tune$bestTune[5]),
    min_child_weight = as.numeric(xgb_tune$bestTune[6]),
    subsample = as.numeric(xgb_tune$bestTune[7])
  )
  
  # class(xgb_tune) <- c("xgb", class(xgb_tune))
  
  return(xgb_tune)
}


f.xgb_mat <- function(data, ...){
  
  group <- names(data)[sapply(data, is.factor)]
  if(length(group) !=1){
    stop("must one factor column")
  }
  # data
  X <- data %>% dplyr::select(!all_of(group))
  y <- data[[group]]
  lev <- levels(data[[group]])
  # xgb.DMatrix
  xgb_mat <- xgb.DMatrix(
    data = Matrix(data.matrix(X), sparse = T), 
    label = ifelse(y == lev[2], 1, 0)
  )
  # xgb.DMatrix.save(dtrain, "dtrain.buffer")
  # dtrain2 <- xgb.DMatrix("dtrain.buffer")
  
  return(xgb_mat)
}


f.tuning_ksvm <- function(df1, formula, ctrl, ...){
  resp <- f.splitform(formula)$left
  lev <- levels(df1[[resp]])
  # data
  traindata <- df1
  # param
  svm_grid <- expand.grid(
    C = seq(from = 0.01, to = 1, by = 0.09),
    sigma = seq(from = 0.001, to = 0.5, by = 0.05)
  )
  
  # tuning
  set.seed(12346)
  svm_tune <- train(
    form = formula, #Group ~ .,               
    data = traindata,         
    method = "svmRadial",  
    trControl = ctrl$ctrl,  
    metric = ctrl$metric,   
    tuneGrid = svm_grid,
    preProcess = c("scale")
    # preProcess = c("center", "scale") 
  )
  
  svm_tune$param <- list(
    modeling = "svm",
    formula = formula,
    lev = lev,
    sigma = as.numeric(svm_tune$bestTune[1]),  # tuning
    C = as.numeric(svm_tune$bestTune[2])       # tuning
  )
  
  # class(svm_tune) <- c("svm", class(svm_tune))
  
  return(svm_tune)
}


f.tuning_lr <- function(df1, formula, ctrl, ...){
  
  lr_tune <- list()
  
  lr_tune$param <- list(
    modeling = "lr",
    formula = formula
  )
  
  # class(lr_tune) <- c("lr", class(lr_tune))
  
  return(lr_tune)
}


f.tuning <- function(df1, method, formula, ctrl, ...){
  
  obj <- switch (
    method,
    `rf` = list(
      func = "f.tuning_rf"
    ),
    `xgb` = list(
      func = "f.tuning_xgb"
    ),
    `svm` = list(
      func = "f.tuning_ksvm"
    ),
    `lr` = list(
      func = "f.tuning_lr"
    ),
    stop("Invalid choice.")
  )
  
  temp_tune <- do.call(
    what = get(obj$func),
    args = list(
      df1 = df1,
      formula = formula,
      ctrl = ctrl
    )
  )
  
  return(temp_tune)
}


## modeling ----
f.modeling.rf <- function(df1, formula, mtry, ntree=500, ...){
  # data
  traindata <- df1
  resp <- f.splitform(formula)$left
  lev <- levels(traindata[[resp]])
  
  set.seed(001)
  fit <- randomForest(
    formula = formula, #Group ~ .,
    data = traindata,
    mtry = mtry, #temp_tune$bestTune$mtry,
    ntree = ntree,
    # ntree = 500,
    importance = TRUE,  
    proximity = TRUE
  )
  
  fit$param <- list(
    modeling = "rf",
    formula = formula,
    mtry = mtry,
    ntree = ntree,
    type = "prob",
    direction = "<",
    lev = lev
  )
  
  class(fit) <- unique(c("rf", class(fit)))
  
  return(fit)
}


f.modeling.xgb <- function(df1, formula, watchlist=NULL, max_depth, gamma, eta, colsample_bytree, min_child_weight, subsample, ... ){
  # data
  traindata <- df1
  
  resp <- f.splitform(formula)$left
  lev <- levels(traindata[[resp]])
  var <- f.splitform(formula)$right
  
  traindata <- traindata %>%
    dplyr::select(all_of(c(resp, var)))
  
  if(is.null(watchlist)){
    watchlist = list(
      # valid = temp_set$valid %>% f.xgb_mat(),
      # all = temp_set$c1 %>% f.xgb_mat(),
      train = traindata %>% f.xgb_mat()
    )
  }else{
    watchlist <- watchlist %>%
      map(~{
        .x %>% 
          dplyr::select(all_of(c(resp, var))) %>%
          f.xgb_mat()
      })
  }
  
  set.seed(002)
  fit <- xgb.train(
    # params = as.list(temp_tune$bestTune),
    params = list(
      # eval_metric = "auc",
      # max_depth = as.numeric(temp_tune$bestTune[2]),    # tuning
      max_depth = max_depth,
      # gamma = as.numeric(temp_tune$bestTune[4]),        # tuning
      gamma = gamma,
      
      eta = eta,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      subsample = subsample
      # eta = as.numeric(temp_tune$bestTune[3]),
      # colsample_bytree = as.numeric(temp_tune$bestTune[5]),
      # min_child_weight = as.numeric(temp_tune$bestTune[6]),
      # subsample = as.numeric(temp_tune$bestTune[7])
    ),
    nrounds = 100,
    # nrounds = as.numeric(temp_tune$bestTune[1]),
    data = traindata %>% f.xgb_mat(),
    watchlist = watchlist,
    # eval_metric = "auc",
    verbose = 1, 
    early_stopping_rounds = 10, 
    booster = "gbtree",
    objective = "binary:logistic"
  )
  
  fit$param <- list(
    modeling = "xgb",
    formula = formula,
    max_depth = max_depth, 
    gamma = gamma, 
    eta = eta, 
    colsample_bytree = colsample_bytree, 
    min_child_weight = min_child_weight, 
    subsample = subsample,
    type = "prob",
    direction = "<",
    lev = lev
  )
  
  class(fit) <- unique(c("xgb", class(fit)))
  
  return(fit)
}


f.modeling.ksvm <- function(df1, formula, sigma, C, ...){
  # data
  traindata <- df1
  
  resp <- f.splitform(formula)$left
  lev <- levels(traindata[[resp]])
  
  set.seed(003)
  fit <- ksvm(
    x = formula, #Group ~ .,
    data = traindata,
    kernel = 'rbfdot', 
    kpar = list(
      # sigma = as.numeric(temp_tune$bestTune[1]) # 0.01
      sigma = sigma
    ),
    # C = as.numeric(temp_tune$bestTune[2]),  # 0.1
    C = C,
    # type = "C-svc",
    # scaled = c(), 
    prob.model = TRUE
  )
  
  slot(fit, "param") <- list(
    modeling = "svm",
    formula = formula,
    sigma = sigma,
    C = C,
    type = "probabilities",
    direction = "<",
    lev = lev
  )
  
  # class(fit) <- unique(c("svm", class(fit)))
  # setClass("ksvm", contains = class(fit))
  # fit <- as(fit, "svm")
  
  return(fit)
}


f.modeling.lr <- function(df1, formula, ...){
  # data
  traindata <- df1
  
  resp <- f.splitform(formula)$left
  lev <- levels(traindata[[resp]])
  
  set.seed(004)
  fit <- glm(
    formula = formula, #Group ~ .,
    data = traindata,
    family = binomial()
  )
  
  fit$param <- list(
    modeling = "lr",
    formula = formula,
    type = "response",
    direction = "<",
    lev = lev
  )
  
  class(fit) <- unique(c("lr", class(fit)))
  
  return(fit)
}


f.modeling <- function(df1, tune=NULL, manual=NULL,  ...){
  # data
  traindata <- df1
  
  if(!is.null(tune)){
    param <- tune$param
  }else{
    param <- list()
  }
  
  # manual <- list(
  #   # rf
  #   mtry, ntree,
  #   # xgb
  #   max_depth, gamma, 
  #   # svm
  #   sigma, C
  # )
  if(!is.null(manual)){
    param[names(manual)] <- manual
  }
  
  modeling <- param$modeling
  
  obj <- switch (
    modeling,
    `rf` = list(
      func = "f.modeling.rf", 
      formula = param$formula,
      mtry = param$mtry
      # ntree = param$ntree
    ),
    `xgb` = list(
      func = "f.modeling.xgb",
      max_depth = param$max_depth, 
      gamma = param$gamma, 
      watchlist = param$watchlist, 
      eta = param$eta, 
      colsample_bytree = param$colsample_bytree, 
      min_child_weight = param$min_child_weight, 
      subsample = param$subsample,
      formula = param$formula
    ),
    `svm` = list(
      func = "f.modeling.ksvm", 
      formula = param$formula,
      sigma = param$sigma, 
      C = param$C
    ),
    `lr` = list(
      func = "f.modeling.lr",
      formula = param$formula
    ),
    stop("Invalid choice.")
  )
  
  fit <- do.call(
    what = get(obj$func),
    args = c(list(traindata), obj[names(obj) != "func"])
  )
  
  return(fit)
}


## tuning & modeling ----
f.tundeling <- function(df1, method, formula, ctrl, ...){
  
  temp_tune <- f.tuning(df1, method, formula, ctrl)
  
  fit <- f.modeling(
    df1 = df1, 
    tune = temp_tune,
    manual = list()
  ) 
  
  res <- list(
    tune = temp_tune,
    fit = fit
  )
  
  return(res)
}


f.tundeling.p <- function(df1, p.method, core.method, formula, ctrl, ...){
  data <- df1
  formula <- formula
  ctrl <- ctrl
  
  temp_res <- f.func.p(
    data = data,
    func = f.tundeling,
    formula = formula,
    ctrl = ctrl,
    arg = "method",
    p.arg = p.method,
    core = core.method,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "tidyr", "Matrix", "caret", "kernlab", "xgboost", "randomForest")
  )
  
  return(temp_res)
}


## evaluate ----
f.eval_cal.metric <- function(df1, pos="HCC", ...) {
  positive_class <- pos
  negative_class <- setdiff(unique(df1$Prediction), positive_class)
  
  tp <- sum(df1[df1$Prediction == positive_class & df1$Reference == positive_class, "Freq"])
  fp <- sum(df1[df1$Prediction == positive_class & df1$Reference == negative_class, "Freq"])
  fn <- sum(df1[df1$Prediction == negative_class & df1$Reference == positive_class, "Freq"])
  tn <- sum(df1[df1$Prediction == negative_class & df1$Reference == negative_class, "Freq"])
  total <- tp + fp + fn + tn
  
  sensitivity <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))  # 灵敏度/召回率
  specificity <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))  # 特异度
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))    # 精确率/阳性预测值
  npv <- ifelse((tn + fn) == 0, 0, tn / (tn + fn))          # 阴性预测值
  accuracy <- ifelse((tp + fp + fn + tn) == 0, 0, (tp + tn) / (tp + fp + fn + tn))  # 准确率
  
  # f1 <- ifelse((precision + sensitivity) == 0, 0, 2 * (precision * sensitivity) / (precision + sensitivity))
  # prevalence <- ifelse((tp + fn + tn + fp) == 0, 0, (tp + fn) / (tp + fn + tn + fp))
  # detection_rate <- ifelse((tp + fn + tn + fp) == 0, 0, tp / (tp + fn + tn + fp))
  # detection_prevalence <- ifelse((tp + fn + tn + fp) == 0, 0, (tp + fp) / (tp + fn + tn + fp))
  # balanced_accuracy <- (sensitivity + specificity) / 2
  
  format_percent <- function(value) {
    sprintf("%.1f%%", value * 100)
  }
  
  result <- data.frame(
    Sensitivity = paste0(tp, "/", tp + fn, " (", format_percent(sensitivity), ")"),
    # Recall = paste0(tp, "/", tp + fn, " (", format_percent(sensitivity), ")"),
    Specificity = paste0(tn, "/", tn + fp, " (", format_percent(specificity), ")"),
    `Pos Pred Value` = paste0(tp, "/", tp + fp, " (", format_percent(precision), ")"),
    # Precision = paste0(tp, "/", tp + fp, " (", format_percent(precision), ")"),
    `Neg Pred Value` = paste0(tn, "/", tn + fn, " (", format_percent(npv), ")"),
    Accuracy = paste0((tp + tn), "/", total, " (", format_percent(accuracy), ")"),
    # F1 = paste0(format_percent(f1)),
    # Prevalence = paste0(format_percent(prevalence)),
    # `Detection Rate` = paste0(format_percent(detection_rate)),
    # `Detection Prevalence` = paste0(format_percent(detection_prevalence)),
    # `Balanced Accuracy` = paste0(format_percent(balanced_accuracy)),
    positive = positive_class,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


f.eval <- function(df1, formula, direction="auto", cutoff="yd", lev=NULL, ...){
  formula <- as.formula(formula)
  # x y
  resp <- f.splitform(formula)$left
  if(is.null(lev)){
    lev <- levels(df1[[resp]])
  }
  pred_prob <- f.splitform(formula)$right
  if(length(pred_prob) != 1){
    stop("formula must one right")
  }
  # data
  temp_eval <- df1 %>%
    dplyr::rename(
      "pred_prob" = !!sym(pred_prob),
      "resp" = !!sym(resp)
    )
  
  temp_res <- list()
  
  ROC <- pROC::roc(
    response = temp_eval$resp, predictor = as.numeric(temp_eval$pred_prob),
    ci = T, auc = T, levels = lev,
    direction = direction,
    # smooth = T, smooth.method = "density",
  )
  
  # ROC data
  temp_res$roc <- data.frame(
    `FPR` = 1 - ROC$specificities,
    `TPR` = ROC$sensitivities
  ) %>%
    mutate(
      AUC = ROC$auc
    )
  
  temp_res$roc.s <- data.frame(
    `FPR` = 1 - smooth(ROC, method = "density")$specificities,
    `TPR` = smooth(ROC, method = "density")$sensitivities
  ) %>%
    mutate(
      AUC = ROC$auc
    )
  
  # Maximize Youden's Index
  ROC$yd <- pROC::coords(ROC, "best", "threshold",  best.method = "youden") %>% as.data.frame()
  # maximize sensitivity
  ROC$yd <- ROC$yd %>% 
    filter(sensitivity == max(sensitivity)) %>% 
    dplyr::select(threshold) %>% as.numeric()
  
  # pred_class
  if(any("yd" %in% cutoff)){ cutoff <- ROC$yd }
  
  if(ROC$direction == "<"){  # controls/lev[1] < cases/lev[2]
    temp_eval <- temp_eval %>%
      mutate(
        pred_class = if_else( pred_prob > cutoff, lev[2], lev[1] )
      ) %>%
      mutate(
        cutoff = as.character(round(as.numeric(cutoff), digits = 2)),
        direction = "Up"
      ) %>%
      mutate(
        pred_status = case_when(
          resp == lev[2] & pred_prob > as.numeric(cutoff) ~ "Right",
          resp == lev[2] & pred_prob <= as.numeric(cutoff) ~ "Wrong",
          resp == lev[1] & pred_prob < as.numeric(cutoff) ~ "Right",
          resp == lev[1] & pred_prob >= as.numeric(cutoff) ~ "Wrong",
          TRUE ~ "Undefined"
        )
      )
  }else{
    temp_eval <- temp_eval %>%
      mutate(
        pred_class = if_else( pred_prob < cutoff, lev[2], lev[1] )
      ) %>%
      mutate(
        cutoff = as.character(round(as.numeric(cutoff), digits = 2)),
        direction = "Down"
      ) %>%
      mutate(
        pred_status = case_when(
          resp == lev[2] & pred_prob < as.numeric(cutoff) ~ "Right",
          resp == lev[2] & pred_prob >= as.numeric(cutoff) ~ "Wrong",
          resp == lev[1] & pred_prob > as.numeric(cutoff) ~ "Right",
          resp == lev[1] & pred_prob <= as.numeric(cutoff) ~ "Wrong",
          TRUE ~ "Undefined"
        )
      )
  }
  # factorize
  temp_eval$pred_class <- factor(temp_eval$pred_class, levels = lev)
  
  # Confusion matrix
  Confus <- confusionMatrix(
    data = temp_eval$pred_class,     # Predicted class
    reference = temp_eval$resp,      # Response
    positive = lev[2]
  )
  
  # fuz data
  temp_res$confus <- Confus$table %>% 
    data.frame() %>%
    mutate( Accuracy = as.character(Confus$overall[1]) )
  
  # metric with n()
  temp_res$metric.n <- temp_res$confus %>%
    f.eval_cal.metric( pos = lev[2] )
  
  # metric
  temp_res$metric <- rbind(
    data.frame(
      `Metric` = c("AUC", "AUC_L95CI", "AUC_H95CI"),
      `Value` = as.numeric(c(ROC$auc, ROC$ci[1], ROC$ci[3]))
    ),
    rbind(
      Confus$byClass %>% data.frame() %>%
        setNames("Value") %>%
        tibble::rownames_to_column("Metric"),
      Confus$overall %>% data.frame() %>%
        setNames("Value") %>%
        tibble::rownames_to_column("Metric")
    )
  ) %>%
    mutate(
      positive = lev[2]
    )
  
  # metric t()
  temp_res$metric.t <- temp_res$metric %>%
    pivot_wider( names_from = Metric, values_from = Value ) %>%
    mutate(
      `AUC (95% CI)` = paste0(sprintf("%.3f", AUC), " (", sprintf("%.3f", AUC_L95CI), "-", sprintf("%.3f", AUC_H95CI), ")")
    )
  
  temp_res$metric.n <- temp_res$metric.n %>%
    mutate(
      `AUC (95% CI)` = temp_res$metric.t$`AUC (95% CI)`,
      AUC = temp_res$metric.t$AUC,
      AUC_L95CI = temp_res$metric.t$AUC_L95CI,
      AUC_H95CI = temp_res$metric.t$AUC_H95CI
    )
  
  # add AUC
  temp_res$roc <- temp_res$roc %>%
    mutate(
      AUC = temp_res$metric.t$AUC[1]
    )
  temp_res$roc.s <- temp_res$roc.s %>%
    mutate(
      AUC = temp_res$metric.t$AUC[1]
    )
  
  # predict data
  temp_res$score <- temp_eval #%>% tibble::rownames_to_column("Sample")
  
  temp_res$param <- list(
    formula = formula,
    lev = lev,
    direction = direction,
    cutoff = cutoff
  )
  
  return(temp_res)
}


f.predict <- function(fit, newdata, ...){
  if(isS4(fit)){
    param <- fit@param
  }else{
    param <- fit$param
  }
  
  type <- param$type
  resp <- f.splitform(param$formula)$left
  lev <- levels(newdata[[resp]])
  
  # predict
  if( "xgb.Booster" %in% class(fit) ){
    temp <- predict(
      fit, 
      newdata = newdata %>% 
        dplyr::select(all_of(as.character(unlist(f.splitform(fit$param$formula))))) %>%
        f.xgb_mat(), 
      type = type
    ) %>%
      as.data.frame() %>%
      data.frame( row.names = rownames(newdata) )
  }else{
    temp <- predict(fit, newdata = newdata, type = type) %>% 
      as.data.frame() %>%
      data.frame( row.names = rownames(newdata) )
  }
  
  # setname
  if(ncol(temp) == 1){
    temp <- temp %>% setNames(lev[2])
  }
  
  # class(newdata) <- "data.frame"
  
  # format output
  temp_res <- temp %>%
    dplyr::select(all_of(lev[2])) %>% 
    rename("pred_prob" = lev[2]) %>%
    tibble::rownames_to_column("Sample") %>%
    left_join(
      newdata %>% 
        as_tibble(rownames = "Sample") %>%
        dplyr::select(all_of(c("Sample", "resp"=resp))),# %>% setNames("resp")  # Response
      by = "Sample"
    ) %>%
    # tibble::column_to_rownames("Sample") %>%
    dplyr::select(Sample, resp, pred_prob) %>%
    mutate(
      model = param$modeling
    )
  
  return(temp_res)
}


f.eval_p2e <- function(l1, fit=NULL, index="train", cutoff="yd", formula=NULL, direction="auto", lev_re=NULL, ...){
  data <- l1
  temp_eval <- list()

  for(obj_i in c(index, setdiff(names(data), index)) ){
    # newdata
    newdata <- data[[obj_i]]
    # class(newdata) <- "data.frame"
    
    if(!is.null(fit)){
      
      if(isS4(fit)){ param <- fit@param }else{ param <- fit$param }
      # predict
      temp <- f.predict(fit, newdata)
      # replace lev
      if(!is.null(lev_re)){
        temp <- temp %>%
          mutate(
            resp = factor(resp, labels = lev_re)
          ) %>%
          mutate(
            resp = factor(resp, levels = lev_re)
          )
      }
      # evalute
      for(i in 1:2){
        temp_eval[[obj_i]] <- f.eval(
          df1 = temp,
          formula = resp ~ pred_prob,
          direction = param$direction,
          cutoff = c(cutoff, temp_eval[[index]]$param$cutoff)[i]
        )
      }
      
    }else{
      
      if(is.null(formula)){
        stop("If fit is null, a formula must be provided")
      }
      # replace lev
      if(!is.null(lev_re)){
        resp <- f.splitform(formula)$left
        newdata <- newdata %>%
          mutate(
            !!resp := factor(!!sym(resp), labels = lev_re)
          ) %>%
          mutate(
            !!resp := factor(!!sym(resp), levels = lev_re)
          )
      }
      # evalute
      for(i in 1:2){
        temp_eval[[obj_i]] <- f.eval(
          df1 = newdata,
          formula = formula,
          direction = direction,
          cutoff = c(cutoff, temp_eval[[index]]$param$cutoff)[i]
        )
      }
      
    }
    
  }

  class(temp_eval) <- unique(c("eval", class(temp_eval)))

  return(temp_eval)
}


f.eval.res <- function(temp_eval, index, ...){
  
  if( any( "eval" %in% class(temp_eval[[1]])) ){
    temp_res <- temp_eval %>%
      imap_dfr(~{
        .x %>% 
          map(~{
            .x[[index]]
          }) %>%
          keep(~ !is.null(.x)) %>%
          imap_dfr(~{
            .x %>% mutate( Set = .y )
          }) %>%
          mutate( model = .y )
      })
  }else if( any( "eval" %in% class(temp_eval)) ){
    temp_res <- temp_eval %>%
      map(~{
        .x[[index]]
      }) %>%
      keep(~ !is.null(.x)) %>%
      imap_dfr(~{
        .x %>% mutate( Set = .y )
      })
  }else{
    temp_res <- temp_eval[[index]]
  }
  
  return(temp_res)
}


f.eval_temp <- function(temp_set, temp_fit, cutoff="yd", formula=NULL, direction="<", lev_re=NULL, ...){
  
  temp_eval <- f.eval_p2e(
    l1 = temp_set, 
    fit = temp_fit, 
    index = "train", 
    cutoff = cutoff,
    formula = formula,
    direction = direction,
    lev_re = lev_re
  )
  
  temp_res <- f.eval.res(temp_eval, index = "metric") %>%
    filter(Metric %in% c("AUC", "Sensitivity", "Specificity")) %>%
    pivot_wider( names_from = Metric, values_from = Value ) %>%
    mutate(across(where(is.numeric), ~ signif(., digits = 3)))
  
  return(temp_res)
}


## temp_func ----
f.predata_diag <- function(temp_set, temp_info, param, seed=123, form=F, cut=F, ...){
  
  feature <- switch (
    param$var_bl,
    "AFP" = c("AFP"),
    "DCP" = c("DCP"),
    "BL" =  c("AFP", "DCP"),
    "M_AFP" = c("AFP", param$var),
    "M_DCP" = c("DCP", param$var),
    "M_BL" =  c("AFP", "DCP", param$var),
    "aMAP" = c("aMAP"),
    "ASAP" = c("ASAP"),
    "All" =  c("AFP", "DCP", "aMAP", "ASAP", param$var),
    param$var
  )
  # formula
  if(form){
    formula <- reformulate(termlabels = feature, response = param$resp)
    return(formula)
  }
  
  # filter
  temp_info <- switch(
    param$var_bl,
    "AFP" = temp_info %>% filter(!is.na(AFP)),
    "DCP" = temp_info %>% filter(!is.na(DCP)),
    "BL" =  temp_info %>% filter(!is.na(AFP) & !is.na(DCP)),
    "M_AFP" = temp_info %>% filter(!is.na(AFP)),
    "M_DCP" = temp_info %>% filter(!is.na(DCP)),
    "M_BL" =  temp_info %>% filter(!is.na(AFP) & !is.na(DCP)),
    "aMAP" = temp_info %>% filter(!is.na(aMAP)),
    "ASAP" = temp_info %>% filter(!is.na(ASAP)),
    "All" =  temp_info %>% filter(!is.na(AFP) & !is.na(DCP) & !is.na(aMAP) & !is.na(ASAP)),
    temp_info
  )
  
  cutoff <- switch(
    param$var_bl,
    "AFP" = 20,
    "DCP" = 40,
    # "aMAP" = 50,
    "aMAP" = 60,
    "ASAP" = 0.5256,
    c(0.5, "yd")[1]
  )
  # cutoff
  if(cut){
    return(cutoff)
  }
  
  # subgroup
  sample <- switch(
    param$subgroup,
    "none" = unique(as.character(unlist(lapply(temp_set, row.names)))),
    "AFP" = temp_info %>% filter(AFP < 20) %>% rownames(.),
    "DCP" = temp_info %>% filter(DCP < 40) %>% rownames(.),
    "BL" = temp_info %>%  filter(AFP < 20 & DCP < 40) %>% rownames(.),
    "TNM_1" = temp_info %>% filter(is.na(TNM) | TNM %in% c("I")) %>% rownames(.),
    "TNM_2" = temp_info %>% filter(is.na(TNM) | TNM %in% c("II", "IIB")) %>% rownames(.),
    "TNM_3" = temp_info %>% filter(is.na(TNM) | TNM %in% c("III", "IIIA", "IIIB")) %>% rownames(.),
    "TNM_4" = temp_info %>% filter(is.na(TNM) | TNM %in% c("IV")) %>% rownames(.),
    "BCLC_A" = temp_info %>% filter(is.na(BCLC) | BCLC == "A") %>% rownames(.),
    "BCLC_B" = temp_info %>% filter(is.na(BCLC) | BCLC == "B") %>% rownames(.),
    "BCLC_C" = temp_info %>% filter(is.na(BCLC) | BCLC == "C") %>% rownames(.),
    "CNLC_1" = temp_info %>% filter(is.na(CNLC) | CNLC %in% c("I", "Ia", "Ib")) %>% rownames(.),
    "CNLC_2" = temp_info %>% filter(is.na(CNLC) | CNLC %in% c("II", "IIa", "IIb")) %>% rownames(.),
    "CNLC_3" = temp_info %>% filter(is.na(CNLC) | CNLC %in% c("III", "IIIa", "IIIb", "IIIc")) %>% rownames(.),
    "CNLC_4" = temp_info %>% filter(is.na(CNLC) | CNLC %in% c("IV")) %>% rownames(.),
    stop("Invalid choice.")
  )
  
  # filter
  temp_x <- temp_set %>%
    map(~{
      if("resampling" %in% class(.x)){
        .x
      }else{
        .x %>%
          f.left_join(
            temp_info,
            by = "rownames"
          ) %>%
          filter(rownames(.) %in% sample) %>%
          dplyr::select(all_of(c(param$resp, feature)))
      }
    })
  
  if(param$var_bl != "M"){
    temp_x <- temp_x[!sapply(temp_x, function(x) "resampling" %in% class(x))]
  }
  
  temp_x <- temp_x %>%
    map(~{
      .x %>% filter(complete.cases(.))
    })

  return(temp_x)
}


f.eval_sub <- function(list_temp, str_sub, str_a, temp_set, temp_info, index="train", lev_re=NULL, direction="<", ...){
  
  temp_eval <- pblapply(str_sub, function(sub){
    temp_eval <- list()
    obj_a <- switch(
      sub,
      "AFP" = str_a[!grepl("AFP", names(str_a), ignore.case = TRUE)],
      "DCP" = str_a[!grepl("DCP", names(str_a), ignore.case = TRUE)],
      "BL" = str_a[!grepl("AFP|DCP", names(str_a), ignore.case = TRUE)],
      str_a
    )
    temp_x <- obj_a %>%
      map(~{
        # param
        temp <- if(.x %in% names(list_temp)){list_temp[[.x]]$param}else{list_temp[[1]]$param}
        temp$subgroup <- sub    # subggroup
        temp$var_bl <- .x
        temp$formula <- f.predata_diag(temp_set, temp_info, param=temp, form=T)
        temp$cutoff <- f.predata_diag(temp_set, temp_info, param=temp, cut=T)
        # evaluate
        f.eval_p2e(
          l1 = f.predata_diag(temp_set, temp_info, param=temp),
          formula = temp$formula,
          fit = list_temp[[.x]]$fit,
          index = index, 
          # lev_re = param$lev,
          cutoff = temp$cutoff,
          direction = direction
        )
      })
    
    return(temp_x)
  })
  
  return(temp_eval)
}


## eval_table----
f.table.3line <- function(df1, nrow=NULL, ...){
  
  if(is.null(nrow)){
    nrow <- nrow(df1)
  }
  
  T1 <- tableGrob(df1, rows = NULL, theme = ttheme_minimal(
    base_size = 11, 
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 11)
    ),
    core = list(
      fg_params = list(fontface = "plain")
    )
  ))
  # top and sutitle
  T1 <- gtable_add_grob(
    T1,
    grobs = grobTree(
      segmentsGrob(
        x0 = unit(0, "npc"), y0 = unit(1, "npc"),
        x1 = unit(1, "npc"), y1 = unit(1, "npc"),
        gp = gpar(lwd = 3)),
      segmentsGrob(
        x0 = unit(0, "npc"), y0 = unit(0, "npc"),
        x1 = unit(1, "npc"), y1 = unit(0, "npc"),
        gp = gpar(lwd = 1.5))
    ),
    t = 1, b = 1, l = 1, r = ncol(T1)
  )
  # bottom
  T1 <- gtable_add_grob(
    T1,
    grobs = segmentsGrob(
      x0 = unit(0, "npc"), y0 = unit(0, "npc"),
      x1 = unit(1, "npc"), y1 = unit(0, "npc"),
      gp = gpar(lwd = 3)),
    t = nrow+1, b = nrow+1, l = 1, r = ncol(T1)
  )
  
  # grid.newpage()
  # grid.draw(T1)
  
  return(T1)
}


f.table.add_title <- function(T1, text, ...){
  
  T1 <- gtable_add_rows(T1, height = unit(1, "line"), pos = 0)
  
  T1 <- gtable_add_grob(
    T1, 
    grobs = textGrob(
      text, 
      # x = unit(0, "npc"), hjust = 0,      
      gp = gpar(fontsize = 13, fontface = "plain")
    ),
    t = 1, b = 1, l = 1, r = ncol(T1)
  )
  
  return(T1)
}


f.table.eval.metric <- function(df1, to=NULL, ...){
  # data
  P.d <- list()
  P.d$df <- df1
  
  nrow <- nrow(P.d$df)
  if(is.null(to) || to == nrow){
    to <- nrow
  }else if(to > nrow){
    P.d$df[seq(from = (nrow+1), to = to), ] <- NA
  }else{
    P.d$df[seq(from = (to+1), to = nrow), ] <- NA
    nrow <- to
  }
  
  P.d$df <- P.d$df %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))
  
  T1 <- P.d$df %>% f.table.3line(nrow)
  
  return(T1)
}


## eval_plot ----
f.plot.eval.roc <- function(temp_res, label=F, color, colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      "color" = !!sym(color)
    )
  
  if(label){
    P.d$label <- setNames(unique(P.d$df$label), unique(P.d$df$color))
  }else{
    P.d$label <- levels(P.d$df$color)
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = FPR, y = TPR, fill = color)) +
    geom_point(color = "white") +
    geom_path(
      aes(x = FPR, y = TPR, color = color),
      linewidth = 0.8, show.legend = FALSE,
      # alpha = 0.5
    ) +
    # geom_smooth(
    #   aes(x = FPR, y = TPR, color = Panel), 
    #   linewidth = 0.6, show.legend = FALSE,
    #   # alpha = 0.5
    # ) +
    geom_abline(
      slope = 1, intercept = 0, 
      linewidth = 0.6, linetype = 2,
      color = "grey"
    ) +
    scale_color_manual(
      values = colors,
      label = P.d$label
    ) +
    scale_fill_manual(
      values = colors,
      label = P.d$label
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      legend.position = "inside", 
      legend.position.inside = c(1, 0), 
      legend.justification = c(1, 0), 
      #legend.direction = c("horizontal", "vertical")[1] 
      #legend.byrow =  
      legend.background = element_blank(), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 10),
    ) +
    guides(
      fill = guide_legend( override.aes = list( color = "black", shape = 22, size = 3) )
    ) + 
    labs(
      x = "1 - Specificity", 
      y = "Sensitivity"
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.eval.metric.bar <- function(temp_res, color, colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      "color" = !!sym(color)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Metric, y = Value, fill = color)) +
    geom_bar(
      stat = "identity", 
      position = "dodge", 
      width = 0.7,
      alpha = 0.8
    ) +
    # geom_hline(yintercept = 0.8, color = "burlywood4", linetype = "twodash", linewidth = 0.4) + 
    scale_fill_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      # panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, size = 10, hjust = 1, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 10)
    ) +
    labs(
      fill = color
    ) +
    scale_y_continuous(
      expand = c(0, 0, 0.05, 0),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    )
  
  # print(P1)
  
  return(P1)
}


f.plot.eval.metric.radar <- function(temp_res, color, colors, ncol=NULL, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      "color" = !!sym(color)
    )
  
  P.d$df2 <- P.d$df %>%
    pivot_wider(
      names_from = "Metric",
      values_from = "Value"
    ) %>%
    dplyr::select(-Set, -positive)
  
  P1 <- levels(P.d$df2$color) %>%
    map(~{
      ggradar(
        plot.data = P.d$df2 %>% filter(color == .x),
        # grid.min = 0.5,
        # grid.mid = 0.75,
        # centre.y = 0,
        # grid.max = 1,
        # values.radar = c(0.5, 0.75, 1),
        # group.colours = colors[unique(P.d$df$color)],
        group.colours = colors[.x],
        plot.title = .x,
        fill = T,
        fill.alpha = 0.2,
        base.size = 4,
        axis.label.size = 3.5,
        group.point.size = 2.5,
        group.line.width = 0.7,
        legend.position = "bottom",
        background.circle.colour = c("white", "lightgray")[1],
        background.circle.transparency = 0.1,
        grid.line.width = 0.5,
        grid.label.size = 4,
        legend.text.size = 10,
        gridline.min.colour = "grey",
        gridline.mid.colour = "grey",
        gridline.max.colour = "grey",
        gridline.min.linetype = "longdash",
        gridline.mid.linetype = "longdash",
        gridline.max.linetype = "longdash"
      ) +
        theme(
          # aspect.ratio = 1/1,
          plot.title = element_text(size = 12, hjust = 0.5)
        )
    })
  
  # print(P1[[1]])
  
  if(is.null(ncol)){
    ncol <- length(P1)
  }
  
  P <- wrap_plots(P1, ncol = ncol)
  
  # print(P)
  
  return(P)
  
  # P.d$df2 <- P.d$df %>%
  #   pivot_wider(
  #     names_from = "Metric",
  #     values_from = "Value"
  #   ) %>%
  #   dplyr::select(-Set, -positive) %>%
  #   tibble::column_to_rownames("color") %>%
  #   add_row(
  #     !!!setNames(rep(0.5, ncol(.)), colnames(.)),
  #     .before = 1
  #   ) %>%
  #   add_row(
  #     !!!setNames(rep(1, ncol(.)), colnames(.)),
  #     .before = 1
  #   )
  # 
  # par(mfrow = c(1, length(rownames(P.d$df2)[-c(1:2)])))
  # for(i in rownames(P.d$df2)[-c(1:2)]){
  #   radarchart(
  #     df = P.d$df2 %>%
  #       filter(
  #         rownames(.) %in% c("...1", "...2", i)
  #       ), 
  #     maxmin = T,
  #     axistype = 1, 
  #     seg = 2,
  #     caxislabels = c("50 (%)", "75 (%)"),
  #     axislabcol = "grey", 
  #     pcol = colors[i], 
  #     pfcol = adjustcolor(colors[i], alpha.f = 0.2),
  #     plwd = 2, plty = 1,
  #     cglcol = "grey", cglwd = 1, cglty = 2, vlcex = 1,
  #     title = i
  #   )
  # }
  
}


f.plot.eval.confus <- function(temp_res, colors=NULL, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  P.d$df$Reference <- factor(P.d$df$Reference, levels = rev(levels(P.d$df$Reference)))
  P.d$df$Reference
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile(color = "white") +
    coord_equal() + 
    scale_fill_gradientn(
      colors = colors
    ) +
    geom_text(
      aes(label = Freq), 
      color = c( "white", "black", "black","white"),
      size = 4.5
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio, 
      panel.grid = element_blank(),
      # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.3),
      legend.position = "none",
      panel.border = element_rect(colour = "white", fill = NA),
      strip.placement = "outside",
      strip.text = element_text(size = 10, color = "black"),
      strip.background = element_rect(fill = NA, color = NA),
      plot.title = element_text(size = 11, hjust = 0.5)
    ) + 
    scale_x_discrete(
      name = expression(underline("          Prediction          ")),
      position = "top",
      expand = c(0, 0, 0, 0)
    ) + 
    scale_y_discrete(
      name = expression(underline("          Reference          ")),
      position = "left",
      expand = c(0, 0, 0, 0)
    ) +
    labs(
      title = unique(P.d$df$Combine),
      # parse = TRUE
    )
  
  # print(P1)
  
  P2 <- ggdraw() +
    draw_label(
      paste0("Accuracy = ", sprintf("%.3f", as.numeric(unique(P.d$df$Accuracy)) ) ), 
      x = 0.5, hjust = 0.3, vjust = 0, size = 11
    )
  
  # print(P2)
  
  P <- ggarrange(
    P1, P2,
    ncol = 1, 
    heights = c(1, 0.1)
  )
  
  # print(P)
  
  return(P)
}


f.plot.eval.score <- function(temp_res, width=0.5, label=F, color, colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      color = !!sym(color)
    ) %>%
    arrange(color) %>%
    mutate(
      Combine = paste0(color, " in ", Set),
      x_jittered = as.numeric(factor(resp)) + runif(n(), -0.3/2, 0.3/2)
    )
  
  P.d$df$Combine <- factor(P.d$df$Combine, levels = unique(P.d$df$Combine[order(P.d$df$color)]))
  
  P.d$point <- P.d$df %>%
    mutate(
      resp = if_else(pred_status == "Wrong", "Normal", resp)
    )
  
  P.d$param$lev <- levels(P.d$df$resp)
  P.d$param$cutoff <- P.d$df %>%
    group_by(Combine) %>%
    summarise(
      cutoff_y = round(as.numeric(unique(cutoff_y)), digits = 2),
      cutoff = round(as.numeric(unique(cutoff)), digits = 2),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste("cutoff: ", cutoff)
    )
  
  # plot 
  P1 <- ggplot(P.d$df, aes(x = resp, y = pred_prob)) +
    facet_wrap(
      ~ Combine, 
      #scales = "free_y"
    ) +
    stat_boxplot(
      geom = "errorbar", 
      width = 0.15, 
      color = "black", 
      linewidth = 0.55
    ) +
    geom_boxplot(
      outliers = F, 
      width = width, 
      color = "black", fill = "white"
    ) + 
    geom_point(
      data = P.d$point,
      aes(x = x_jittered, fill = resp),
      shape = 21, size = 2, stroke = 0.01, alpha = 0.8
    ) +
    geom_hline(
      data = P.d$param$cutoff,
      aes(yintercept = cutoff_y),
      color = "burlywood4", 
      linetype = "longdash", linewidth = 0.5
    ) + 
    scale_fill_manual(
      values = c("Normal"="grey", colors),
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
      y = paste0("Probability for ", tail(P.d$param$lev, 1))
    ) +
    geom_text(
      data = P.d$param$cutoff, 
      aes(label = label, x = 0.5, y = cutoff_y), 
      vjust = 0, hjust = 0,          
      size = 3.3, 
      color = "firebrick",
      inherit.aes = FALSE
    ) + 
    stat_compare_means(
      method = c("t.test", "wilcox.test")[2],
      label =  "p.signif",
      paired = F,
      size = 3.3,
      comparisons = combn(P.d$param$lev, 2, simplify = F)
    ) +
    scale_y_continuous(
      expand = c(0.08, 0, 0.08, 0)
    )
  
  if(max(P.d$df$pred_prob) < 1){
    P1 <- P1 +
      scale_y_continuous(
        expand = c(0.08, 0, 0.08, 0),
        breaks = seq(0, 1, by = 0.25)
      )
  }
  
  # print(P1)
  
  if(label){
    P.d$label2 <- P.d$df %>%
      filter(pred_status == "Wrong")
    
    P1 <- P1 + 
      geom_point(
        data = P.d$label2,
        aes(x = x_jittered, y = pred_prob),
        size = 1, alpha = 1, stroke = 0,
        color = "black",
      ) +
      geom_label_repel(
        data = P.d$label2 ,
        aes(x = x_jittered, y = pred_prob, label = Sample),
        alpha = 1, fill = NA,
        size = 2, 
        min.segment.length = 0, 
        direction = c("both", "y", "x")[1], 
        segment.size = 0.2, 
        # segment.linetype = 6,
        segment.color = "black",
      )
  }
  
  return(P1)
}


## importance ----
f.imp <- function(fit, ...){
  UseMethod("f.imp")
}

f.imp.rf <- function(fit, ...){
  
  # varImp(temp_fit, scale = FALSE)
  # varImp(temp_fit)  平均排列重要性 (Permutation Importance)
  # importance(temp_fit)  
  # 特定类别排列重要性 (Permutation Importance)、平均准确度下降（基于排列重要性计算的变量重要性度量）、基于基尼不纯度减少量计算的变量重要性
  # temp_fit$importance
  # temp_fit$importanceSD
  # 1.206764e-03
  # 2.48850676   1.80451333  1.996626590
  # 4.849351e-04
  
  temp_res <- importance(fit) %>%
    data.frame() %>%
    rename(
      "Mean Decrease Accuracy" = MeanDecreaseAccuracy,
      "Mean Decrease Gini" = MeanDecreaseGini
    ) %>%
    dplyr::select(c("Mean Decrease Accuracy", "Mean Decrease Gini")) %>%
    tibble::rownames_to_column("Feature") %>%
    pivot_longer(
      cols = -Feature, 
      names_to = "Metric", 
      values_to = "Importance"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("Mean Decrease Accuracy", "Mean Decrease Gini")
      )
    ) %>%
    arrange(Metric, desc(Importance)) %>%
    mutate(
      model = "Random Forest"
    )
  
  class(temp_res) <- unique(c("imp", "rf", class(temp_res)))
  
  return(temp_res)
}


f.imp.xgb <- function(fit, ...){
  
  temp_res <- xgboost::xgb.importance(model = fit) %>%
    data.frame() %>%
    pivot_longer(
      cols = -Feature,
      names_to = "Metric",
      values_to = "Importance"
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("Gain", "Cover", "Frequency")
      )
    ) %>%
    arrange(Metric, desc(Importance)) %>%
    mutate(
      model = "XGBoost"
    )
  
  class(temp_res) <- unique(c("imp", "xgb", class(temp_res)))
  
  return(temp_res)
}


f.plot.imp.bar <- function(temp_res, color="Metric", colors, ratio=NULL, theme="bw", ...){
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    rename(
      color = !!sym(color)
    )
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = Importance, y = Feature, fill = color)) + 
    geom_bar(
      stat = "identity",
      position = "dodge",
      width = 0.7,
      alpha = 0.8,
    ) +
    scale_fill_manual(
      values = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      # panel.grid = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(size = 10),
      legend.title = element_blank(),
      # legend.position = "none"
    ) +
    scale_x_continuous(
      expand = c(0, 0, 0.1, 0)
    ) +
    facet_wrap(
      ~ color,
      scales = "free_x"
    ) +
    theme(
      aspect.ratio = 1.5/1,
      legend.position = "none"
    )
  
  # print(P1)
  
  return(P1)
}



