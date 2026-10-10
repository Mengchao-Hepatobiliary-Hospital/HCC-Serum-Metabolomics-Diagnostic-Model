# enrich ----
## predata ----
f.add_bitr <- function(df1, fromType=c("Feature"="UNIPROT"), toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Hs.eg.db", ...){
  
  v1 <- df1 %>%
    dplyr::select(all_of(setNames(names(fromType), unname(fromType)))) %>%
    pull(fromType)
  
  temp <- bitr(
    geneID = v1,
    fromType = fromType,
    toType = toType,
    OrgDb = OrgDb,
    drop = F
  )
  
  temp_res <- df1 %>%
    left_join(
      temp,
      by = fromType
    )
  
  return(temp_res)
}


f.enrich.res <- function(temp_res, SYMBOL=T, ...){
  
  # if(SYMBOL){
  #   temp <- temp_res %>% setReadable(OrgDb = "org.Hs.eg.db", keyType=names(temp_res@keytype))
  #   
  #   res <- temp_res@result %>%
  #     as.data.frame() %>%
  #     cbind(
  #       temp@result %>%
  #         as.data.frame() %>%
  #         dplyr::select(any_of(c("geneID", "core_enrichment"))) %>%
  #         dplyr::rename(
  #           "SYMBOL" = any_of(c("geneID", "core_enrichment"))
  #         )
  #     ) %>%
  #     remove_rownames()
  # }else{
  #   res <- temp_res@result %>%
  #     as.data.frame() %>%
  #     remove_rownames()
  # }
  
  tryCatch({
    temp <- temp_res %>% setReadable(OrgDb = "org.Hs.eg.db", keyType = names(temp_res@keytype))
  }, error = function(e) {
    has_error <<- TRUE
    message("Jump: ", e$message)
  })
  
  if (!has_error) {
    res <- temp_res@result %>%
      as.data.frame() %>%
      cbind(
        temp@result %>%
          as.data.frame() %>%
          dplyr::select(any_of(c("geneID", "core_enrichment"))) %>%
          dplyr::rename(
            "SYMBOL" = any_of(c("geneID", "core_enrichment"))
          )
      ) %>%
      remove_rownames()
  } else {
    res <- temp_res@result %>%
      as.data.frame() %>%
      remove_rownames()
  }
  

  if("gseaResult" %in% class(temp_res)){
    temp_x <- res %>%
      mutate(
        organism = temp_res@organism,
        keytype = temp_res@keytype,
        setType = temp_res@setType,
        pvalueCutoff = temp_res@params$pvalueCutoff,
        pAdjustMethod = temp_res@params$pAdjustMethod,
        exponent = temp_res@params$exponent,
        geneList_name = paste(names(temp_res@geneList), collapse = "/"),
        geneList_value = paste(temp_res@geneList, collapse = "/")
      )
    
    temp_y <- setNames(res$ID, res$ID) %>%
      map_dfr(~{
        enrichplot:::gsInfo(temp_res, geneSetID = .x) %>%
          mutate(
            ID = .x
          )
      })
    
    temp_y <- temp_y %>%
      left_join(
        temp_x %>% dplyr::select(-all_of( setdiff(intersect(names(temp_x), names(temp_y)), "ID") )) %>%
          dplyr::select(-c("geneList_name", "geneList_value")),
        by = "ID"
      ) %>%
      mutate(
        Delta = abs(runningScore - enrichmentScore)
      )
    
    res <- list(
      `res` = temp_x,
      `score` = temp_y
    )
    
  }else{
    res <- res %>%
      mutate(
        organism = temp_res@organism,
        keytype = temp_res@keytype,
        pvalueCutoff = temp_res@pvalueCutoff,
        pAdjustMethod = temp_res@pAdjustMethod,
        qvalueCutoff = temp_res@qvalueCutoff
      )
  }
  
  return(res)
}


f.splitRegulate <- function(temp_var, obj=NULL, ...){
  
  temp <- temp_var %>%
    split(.$Regulate) %>%
    c(
      "diff" = list(bind_rows(.[["Up"]], .[["Down"]])),
      "All" = list(bind_rows(.[["Up"]], .[["Down"]], .[["Normal"]])) 
    ) %>%
    map(~{
      .x %>% arrange(desc(logFC))
    })
  
  if(!is.null(obj)){
    temp <- temp[obj]
  }
  
  return(temp)
}


## ora ----
f.ora_uni <- function(v1, USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, universe=NULL, pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.01, minGSSize=10, maxGSSize=500, ...){
  
  if(!is.null(USER_DATA)){
    temp_res <- v1 %>% 
      DOSE:::enricher_internal(
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        USER_DATA = USER_DATA,
        universe = universe
      ) 
  }else{
    temp_res <- v1 %>%
      enricher(
        # gene,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        gson  = gson,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        universe = universe
      )
  }
  
  temp_res@result <- temp_res@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = "ORA"
    )
  
  return(temp_res)
}


f.ora_kegg <- function(v1, keyType="UNIPROT", organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.01, minGSSize=10, maxGSSize=500, ...){
  
  temp_x <- enrichKEGG(
    gene = v1,
    keyType = c(`KEGG`="kegg", `UNIPROT`="uniprot", `SYMBOL`="SYMBOL")[keyType],
    organism = organism,
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  temp_x@result <- temp_x@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = "KEGG"
    )
  
  return(temp_x)
}


f.ora_go <- function(v1, ont="ALL", keyType="UNIPROT", OrgDb="org.Hs.eg.db", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.01, minGSSize=10, maxGSSize=500, ...){
  
  temp_y <- enrichGO(
    gene = v1, 
    keyType = c(`ENTREZID`="ENTREZID", `UNIPROT`="UNIPROT", `SYMBOL`="SYMBOL")[keyType], 
    OrgDb = OrgDb,
    ont = ont,
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  temp_y@result <- temp_y@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = paste0("GO: ", ONTOLOGY)
    )
  
  return(temp_y)
}


f.ora <- function(v1, enrich, ont="ALL", USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, universe=NULL, keyType="UNIPROT", organism="hsa", OrgDb="org.Hs.eg.db", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.01, minGSSize=10, maxGSSize=500, ...){
  
  obj <- switch (
    enrich,
    `KEGG` = list(
      func = "f.ora_kegg", 
      v1 = v1, 
      keyType = keyType,
      organism = organism,
      pAdjustMethod = pAdjustMethod, 
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    ),
    `GO` = list(
      func = "f.ora_go", 
      ont = ont, 
      v1 = v1, 
      keyType = keyType, 
      OrgDb = OrgDb,
      pAdjustMethod = pAdjustMethod, 
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    ),
    `ORA` = list(
      func = "f.ora_uni", 
      v1 = v1, 
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      USER_DATA = USER_DATA,
      gson  = gson,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      universe = universe
    ),
    stop("Invalid choice.")
  )
  
  temp_res <- do.call(
    what = get(obj$func),
    args = c(list(v1), obj[names(obj) != "func"])
  )
  
  return(temp_res)
}


f.ora.p <- function(v1, p.enrich, core.enrich=F, ont="ALL", USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, universe=NULL, keyType="UNIPROT", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.01, minGSSize=10, maxGSSize=500, ...){
  data <- v1
  
  temp_res <- f.func.p(
    data = data,
    func = f.ora,
    ont = ont,
    keyType = keyType,
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    USER_DATA = USER_DATA,
    gson  = gson,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    universe = universe,
    arg = "enrich",
    p.arg = p.enrich,
    core = core.enrich,
    globals = c( f.find_func("f.") ),
    packages = c("dplyr", "clusterProfiler", "org.Hs.eg.db")
  )
  
  return(temp_res)
}


## gsea ----
f.gsea_uni <- function(v1, USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, pAdjustMethod="BH", pvalueCutoff=0.05, exponent=1, seed=T, minGSSize=10, maxGSSize=500, by="fgsea", ...){
  if(!is.null(USER_DATA)){
    temp_res <- v1 %>%
      DOSE:::GSEA_internal(
        # geneList,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        exponent = exponent,
        eps  = 1e-10,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        USER_DATA = USER_DATA,
        verbose = TRUE,
        seed = seed,
        # nPerm = 10000,
        # by = c("fgsea", "DOSE")[1]
        by = by
      )
  }else{
    temp_res <- v1 %>%
      GSEA(
        # geneList,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        exponent = exponent,
        eps  = 1e-10,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        gson  = gson,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        verbose = TRUE,
        seed = seed,
        # nPerm = 10000,
        by = c("fgsea", "DOSE")[1],
      )
  }
  
  temp_res@result <- temp_res@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = "GSEA"
    )
  
  return(temp_res)
}


f.gsea_kegg <- function(v1, keyType="UNIPROT", organism = "hsa", pAdjustMethod="BH", pvalueCutoff=0.05, exponent=1, seed=T, minGSSize=10, maxGSSize=500, by="fgsea", ...){
  
  temp_x <- gseKEGG(
    gene = v1,
    keyType = c(`KEGG`="kegg", `UNIPROT`="uniprot")[keyType],
    organism = organism,
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    exponent = exponent,
    seed = seed,
    # by = c("fgsea", "DOSE")[1]
    by = by
  )
  
  temp_x@result <- temp_x@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = "KEGG"
    )
  
  return(temp_x)
}


f.gsea_go <- function(v1, ont="ALL", keyType="UNIPROT", OrgDb="org.Hs.eg.db", pAdjustMethod="BH", pvalueCutoff=0.05, exponent=1, seed=T, minGSSize=10, maxGSSize=500, by="fgsea", ...){
  
  temp_y <- gseGO(
    gene = v1, 
    ont = ont,
    keyType = c(`ENTREZID`="ENTREZID", `UNIPROT`="UNIPROT")[keyType], 
    OrgDb = OrgDb, 
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    exponent = exponent,
    seed = seed,
    # by = c("fgsea", "DOSE")[1]
    by = by
  )
  
  temp_y@result <- temp_y@result %>%
    as.data.frame() %>%
    mutate(
      enrichType = paste0("GO: ", ONTOLOGY)
    )
  
  return(temp_y)
}


f.gsea <- function(v1, enrich, ont="ALL", by="fgsea", USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, keyType="UNIPROT", organism="hsa", OrgDb="org.Hs.eg.db", pAdjustMethod="BH", pvalueCutoff=0.05, exponent=1, seed=T, minGSSize=10, maxGSSize=500, ...){
  
  obj <- switch (
    enrich,
    `KEGG` = list(
      func = "f.gsea_kegg", 
      v1 = v1, 
      keyType = keyType,
      organism = organism,
      pAdjustMethod = pAdjustMethod, 
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent, 
      seed = seed,
      by = by
    ),
    `GO` = list(
      func = "f.gsea_go", 
      ont = ont, 
      v1 = v1, 
      keyType = keyType, 
      OrgDb = OrgDb,
      pAdjustMethod = pAdjustMethod, 
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent, 
      seed = seed,
      by = by
    ),
    `GSEA` = list(
      func = "f.gsea_uni", 
      v1 = v1, 
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent, 
      seed = seed,
      USER_DATA = USER_DATA,
      gson  = gson,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      by = by
    ),
    stop("Invalid choice.")
  )
  
  temp_res <- do.call(
    what = get(obj$func),
    args = c(list(v1), obj[names(obj) != "func"])
  )
  
  return(temp_res)
}


f.gsea.p <- function(v1, p.enrich, core.enrich=F, by="fgsea", ont="ALL", USER_DATA=NULL, gson=NULL, TERM2GENE=NULL, TERM2NAME=NULL, keyType="UNIPROT", pAdjustMethod="BH", pvalueCutoff=0.05, exponent=1, seed=T, minGSSize=10, maxGSSize=500, ...){
  data <- v1
  
  temp_res <- f.func.p(
    data = data,
    func = f.gsea,
    ont = ont,
    keyType = keyType,
    pAdjustMethod = pAdjustMethod, 
    pvalueCutoff = pvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    exponent = exponent, 
    seed = seed,
    USER_DATA = USER_DATA,
    gson  = gson,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    by = by,
    arg = "enrich",
    p.arg = p.enrich,
    core = core.enrich,
    globals = c(f.find_func("f."), ".Random.seed"),
    packages = c("dplyr", "stats", "fgsea", "clusterProfiler", "org.Hs.eg.db")
  )
  
  return(temp_res)
}


f.df2gsea <- function(temp_res, geneSet=NULL, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  P.d$geneList <- data.frame(
    name = str_split(P.d$df$geneList_name[1], "/") %>% unlist(),
    value = str_split(P.d$df$geneList_value[1], "/") %>% unlist() %>% as.numeric()
  ) %>%
    tibble::deframe()
  
  if(any(c("GO: BP", "GO: CC", "GO: MF", "KEGG") %in% P.d$df$enrichType)){
    
    if( any(unique(P.d$df$enrichType) %in% "KEGG") ){
      keytype <- c(ENTREZID="kegg", UNIPROT="uniprot")[names(P.d$df$keytype[1])] %>% as.character()
      P.d$geneSet <- clusterProfiler:::prepare_KEGG(P.d$df$organism[1], "KEGG", keytype )
    }else{
      keytype <- c(ENTREZID="kegg", UNIPROT="UNIPROT")[names(P.d$df$keytype[1])]
      P.d$geneSet <- clusterProfiler:::get_GO_data(OrgDb = "org.Hs.eg.db", ont = "ALL", keytype)
    }
    
    P.d$geneSet <- GseaVis:::getGeneSet(P.d$geneSet)
    
  }else{
    
    if(is.null(geneSet)){
      stop("not KEGG or GO, please provide geneSet")
    }else{
      P.d$geneSet <- geneSet
    }
    
  }
  
  gsea <- methods::new(
    "gseaResult",
    result = P.d$df %>%
      dplyr::select(-any_of(c("enrichType", "SYMBOL", "organism", "keytype", "setType", "pvalueCutoff", "pAdjustMethod", "geneList_name", "geneList_value"))),
    organism = P.d$df$organism[1],
    setType = P.d$df$setType[1],
    geneSets = P.d$geneSet,
    geneList = P.d$geneList,
    keytype = P.d$df$keytype[1],
    permScores = matrix(),
    params = list(
      pvalueCutoff = P.d$df$pvalueCutoff[1],
      eps = 1e-10,
      pAdjustMethod = P.d$df$pAdjustMethod[1],
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500
    ),
    gene2Symbol = character(),
    readable = FALSE,
    termsim = matrix(),
    method = character(),
    dr = list()
  )
  
  return(gsea)
}


## plot ----
f.ora.pd <- function(model, showCategory=5, by = "Count", order=FALSE, drop=T, split=NULL, ...) {
  
  if("data.frame" %in% class(model)){
    res <- model %>%
      filter(
        p.adjust < pvalueCutoff
      )
    
    if("qvalueCutoff" %in% names(model)){
      res <- res %>%
        filter(
          qvalue < qvalueCutoff
        )
    }
    
  }else if("enrichResult" %in% class(model)){
    res <- as.data.frame(model)
  }else{
    stop("error input")
  }
  
  res <- res[!is.na(res$Description), ]
  
  if("core_enrichment" %in% names(res)){
    res$Count <- str_count(res$core_enrichment, "/")
    res$.sign <- "activated"
    res$.sign[res$NES < 0] <- "suppressed"
  }
  
  if (drop) {
    res <- res[res$Count != 0, ]
  }
  
  if("GeneRatio" %in% names(res)){
    res$GeneRatio <- parse_ratio(res$GeneRatio)
    if ("BgRatio" %in% colnames(res)) {
      ## groupGO output doesn't have this column
      res$BgRatio <- parse_ratio(res$BgRatio)
    }
  }else if("core_enrichment" %in% names(res)){
    res$GeneRatio <- res$Count / res$setSize
  }else{
    stop("error no GeneRatio")
  }
  
  # showCategory
  if (is.null(split)) {
    
    if( is.numeric(showCategory) ){
      if( showCategory <= nrow(res) ){
        res <- res[1:showCategory, ]
      }
    }else{
      res <- res[res$Description %in% showCategory, ]
    }
    
  }else{
    if( is.numeric(showCategory) ){
      res <- res %>%
        mutate(
          "temp" = !!sym(split)
        ) %>%
        split(.$temp) %>%
        map_dfr(~{
          if(nrow(.x) < showCategory){
            .x %>% dplyr::select(-temp)
          }else{
            .x[1:showCategory, ] %>% dplyr::select(-temp)
          }
        })
    }else{
      res <- res %>%
        mutate(
          "temp" = !!sym(split)
        ) %>%
        split(.$temp) %>%
        map_dfr(~{
          .x[.x$Description %in% showCategory, ] %>% dplyr::select(-temp)
        })
    }
    
  }
  
  if (order) {
    if(!is.null(split)){
      res <- res %>%
        mutate(
          "temp" = !!sym(split),
          "by" = !!sym(by)
        ) %>%
        group_by(temp) %>%
        arrange(by) %>%
        ungroup() %>% 
        dplyr::select(-"temp", -"by") %>%
        mutate(
          order = row_number()
        ) %>%
        mutate(
          order = factor(order)
        )
    }else{
      res <- res %>%
        mutate(
          "by" = !!sym(by)
        ) %>%
        arrange(by) %>% 
        dplyr::select(-"by") %>%
        mutate(
          order = row_number()
        ) %>%
        mutate(
          order = factor(order)
        )
    }
  }
  
  if("GeneRatio" %in% names(res)){
    if(!is.numeric(res$GeneRatio[1])){
      res <- res %>%
        separate(GeneRatio, into = c("Count", "setSize"), sep = "/", convert = TRUE) %>%
        mutate(GeneRatio = Count / setSize)
    }
  }
  
  return(res)
}


f.plot.enrich <- function(temp_res, param, type, ratio=NULL, pf=NULL, autoTimer=T, ...){
  
  res <- list()
  if("bar" %in% type){
    P1 <- temp_res %>%
      f.plot.enrich.bar(
        colors = param$colors,
        ratio = ratio,
        theme = param$theme,
        text.inside = param$text.inside,
        x = param$x, 
        color = param$color,
        showCategory = param$showCategory,
        label_format = param$label_format,
        order = param$order,
        split = param$split
      )
    
    if(!is.null(pf)){
      f.save_pdf( P1, width = 7.3, height = 10, pf = paste0(pf, "-bar" ), autoTimer = autoTimer )
    }
    
    res[["bar"]] <- P1
  }
  
  if("dot" %in% type){
    P1 <- temp_res %>%
      f.plot.enrich.dot(
        colors = param$colors,
        ratio = ratio,
        theme = param$theme,
        x = param$x, 
        color = param$color,
        showCategory = param$showCategory,
        label_format = param$label_format,
        order = param$order,
        split = param$split
      ) +
      scale_size_continuous(
        breaks = c(0.05, 0.10, 0.15, 0.20),
        limits = c(0, 0.3),
        range = c(2, 6)
      )
    
    if(!is.null(pf)){
      f.save_pdf( P1, width = 7.3, height = 10, pf = paste0(pf, "-dot" ), autoTimer = autoTimer )
    }
    
    res[["dot"]] <- P1
  }
  
  if("net" %in% type){
    P1 <- temp_res %>%
      f.plot.enrich.net(
        colors = param$colors,
        foldChange = param$foldChange,
        x = param$x, 
        color = param$color,
        showCategory = param$showCategory,
        label_format = param$label_format,
        order = param$order,
        split = param$split
      )
    
    if(!is.null(pf)){
      f.save_pdf( P1, width = 8, height = 8, pf = paste0(pf, "-net" ), autoTimer = autoTimer )
    }
    
    res[["net"]] <- P1
  }
  
  return(res)
}


f.plot.enrich.bar <- function(temp_res, ratio=NULL, colors=NULL, theme="bw", text.inside=F, split=NULL, x="GeneRatio", color="p.adjust", showCategory=10, label_format=30, order=F, ...){
  if(is.null(colors)){
    colors <- rev(c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00"))
  }else{
    colors <- rev(colors)
  }
  # data
  P.d <- list()
  P.d$df <- temp_res %>%
    f.ora.pd(
      showCategory = showCategory, 
      by = x, 
      split = split,
      order = order
    )  %>%
    dplyr::rename(
      "x" = !!sym(x),
      "color" = !!sym(color)
    )
  
  if(!is.null(split)){
    P.d$df <- P.d$df %>% 
      dplyr::rename(
        "split" = !!sym(split)
      )
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes( x = x, y = order )) +
    geom_col(
      aes(fill = color),
      width = 0.7
    ) +
    # coord_flip() +
    # geom_bar(
    #   stat = "identity", 
    #   position = "dodge", 
    #   width = 0.7
    # ) +
    f.theme(theme) + 
    theme(
      aspect.ratio = ratio
    ) +
    scale_x_continuous(
      expand = c(0, 0, 0.1, 0)
    ) +
    scale_y_discrete(
      breaks = P.d$df$order,
      labels = P.d$df$Description %>% stringr::str_wrap(width = label_format)
    ) +
    scale_fill_gradientn(
      # limits = c(0, 0.05),
      # breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
      colors = colors,
      oob = scales::censor
    ) + 
    labs(
      fill = color,
      x = x
    )
  
  if(text.inside){
    P1 <- P1 +
      geom_text(
        aes(label = Description, x = max(P.d$df$x)*0.05), 
        vjust = 0.5,    
        hjust = 0, 
        #lineheight = 0.7,    
        size = 3.3,           
        color = "gray20",
      ) +
      theme(
        # aspect.ratio = 1/1,
        # panel.grid = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()
      )
  }else{
    P1 <- P1 +
      theme(
        # aspect.ratio = 1/1,
        # panel.grid = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
      )
    
  }
  
  if(!is.null(split)){
    P1 <- P1 +
      facet_grid(
        split~.,
        switch = "x",
        scale = "free_y"
      )
  }
  
  return(P1)
}


f.plot.enrich.dot <- function(temp_res, ratio=NULL, colors=NULL, theme="bw", split=NULL, x="GeneRatio", color="p.adjust", showCategory=10, label_format=30, order=F, ...){
  if(is.null(colors)){
    colors <- rev(c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00"))
  }else{
    colors <- rev(colors)
  }
  # data
  x2 <- setdiff(c("Count", "GeneRatio"), x)
  P.d <- list()
  P.d$df <- temp_res %>%
    f.ora.pd(
      showCategory = showCategory, 
      by = x, 
      split = split,
      order = order
    )  %>%
    dplyr::rename(
      "x" = !!sym(x),
      "x2" = !!sym(x2),
      "color" = !!sym(color)
    )
  
  if(!is.null(split)){
    P.d$df <- P.d$df %>% 
      dplyr::rename(
        "split" = !!sym(split)
      )
  }
  
  # plot
  P1 <- ggplot(P.d$df, aes(x = x, y = order)) +
    geom_point(
      aes(size = x2, fill = color),
      shape = 21,
      color = "black"
    ) +
    f.theme(theme) + 
    theme(
      aspect.ratio = ratio,
      axis.title.y = element_blank()
    ) +
    scale_x_continuous(
      expand = c(0.05, 0, 0.1, 0)
    ) +
    scale_y_discrete(
      breaks = P.d$df$order,
      labels = P.d$df$Description %>% stringr::str_wrap(width = label_format)
    ) +
    scale_fill_gradientn(
      # limits = c(0, 0.05),
      # breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
      colors = colors,
      oob = scales::censor
    ) + 
    labs(
      fill = color,
      size = x2,
      x = x
    )
  
  if(!is.null(split)){
    P1 <- P1 +
      facet_grid(
        split~.,
        switch = "x",
        scale = "free_y"
      )
  }
  
  return(P1)
}


f.plot.enrich.net <- function(temp_res, colors=NULL, foldChange, split=NULL, x="GeneRatio", color="p.adjust", showCategory=10, label_format=30, order=F, ...){
  if(is.null(colors)){
    colors <- c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00")
  }
  # data
  P.d <- list()
  if("data.frame" %in% class(temp_res)){
    P.d$df <- temp_res %>%
      f.ora.pd(
        showCategory = showCategory, 
        by = x, 
        split = split,
        order = order
      )
    
    P.d$geneSet <- P.d$df %>%
      split(.$Description) %>%
      map(~{
        .x %>%
          pull(!!sym(foldChange$index)) %>%
          str_split("/") %>%
          unlist() %>%
          purrr::discard(~ .x == "")
      })
    
    if( any(unlist(P.d$geneSet) %in% names(foldChange$value)) ){
      
      P.d$fc <- foldChange$value
      
    }else{
      
      P.d$mapping <- P.d$df %>%
        dplyr::rename(
          "index" = !!sym(foldChange$index)
        ) %>%
        rename_with(~ "temp", matches("geneID|core_enrichment")) %>%
        separate_rows(temp, index, sep = "/") %>%
        distinct(temp, index) %>%
        {setNames(.$index, .$temp)}
      
      P.d$fc <- setNames(
        foldChange$value[names(foldChange$value) %in% names(P.d$mapping)], 
        P.d$mapping[names(foldChange$value)[names(foldChange$value) %in% names(P.d$mapping)]] 
      )
      
    }
    
  }else{
    
    P.d$geneSet <- temp_res %>% enrichplot:::extract_geneSets(showCategory)
    P.d$fc <- temp_res %>% enrichplot:::fc_readable(foldChange$value)
    
  }
  
  
  P1 <- P.d$geneSet %>%
    enrichplot::cnetplot(
      showCategory = showCategory,
      layout = c("star", "circle", "fr", "kk")[4],
      shadowtext = c("none", "all", "category", "gene")[1],
      color.params = list(
        foldChange = P.d$fc,
        # edge = T,
        category = c("#E5C494", "#dc9ee5")[2],
        gene = "#B3B3B3"
      ),
      cex.params = list(
        category_node = 1, 
        gene_node = 1, 
        category_label = 0.7, 
        gene_label = 0.5
      )
    ) +
    # ggtangle::cnetplot(
    #   foldChange = P.d$fc, 
    #   layout = c("star", "circle", "fr", "kk")[4], 
    #   showCategory = 5,
    #   color_category = c("#E5C494", "#dc9ee5")[2],
    #   size_category = 1,
    #   color_item = "#B3B3B3",
    #   size_item = 0.5,
    #   color_edge = "grey",
    #   size_edge = 0.5,
    #   node_label = "all",
    #   hilight = "none",
    #   hilight_alpha = 0.3,
    # ) +
    # enrichplot::set_enrichplot_color(
    scale_color_gradientn(
      name = expression(paste(Log[2], " (Fold Change)")),
      # colors = c("#1f78b4", "#a6cee3", "#fdbf6f", "#ff7f00")
      colors = colors
    )
  
  # print(P1)
  
  return(P1)
}


f.gsea.pd <- function(temp_res, typeID="ID", geneSetID=NULL, ...){
  # data
  P.d <- list()
  P.d$df <- temp_res
  
  P.d$runScore <- P.d$df
  
  if(is.null(geneSetID)){
    geneSetID <- unique(P.d$runScore[[typeID]])
  }else{
    
    if(any(geneSetID %in% "all")){
      geneSetID <- unique(P.d$runScore[[typeID]])
    }
    
  }
  
  P.d$runScore <- P.d$runScore %>%
    filter(
      .data[[typeID]] %in% geneSetID
    )
  
  P.d$min <- P.d$runScore %>%
    group_by(Description) %>%
    slice_min(Delta) %>%
    reframe(
      x = x,
      y = runningScore
    )
  
  if(length(geneSetID)==1 & !any(geneSetID %in% "all")){
    P.d$pvalue <- P.d$runScore %>%
      group_by(Description) %>%
      reframe(
        quantile_80 = quantile(runningScore, probs = 0.80, na.rm = TRUE),
        max_x = max(x),
        min_running = min(runningScore, na.rm = TRUE),
        max_running = max(runningScore, na.rm = TRUE),
        NES = dplyr::first(NES),
        pvalue = dplyr::first(pvalue),
        p.adjust = dplyr::first(p.adjust)
      ) %>%
      mutate(
        x = if_else(quantile_80 <= 0, max_x * 0.02, max_x * 0.98),
        y = if_else(quantile_80 <= 0, min_running * 0.98, max_running * 0.98),
        hjust = if_else(quantile_80 <= 0, 0, 1),
        vjust = if_else(quantile_80 <= 0, 0, 1),
        label = paste0(
          "NES: ",
          sprintf("%.3f", NES),
          "\n",
          "Pvalue: ",
          ifelse(pvalue < 0.001, "< 0.001", sprintf("%.3f", pvalue)),
          "\n",
          "Adjusted Pvalue: ",
          ifelse(p.adjust < 0.001, "< 0.001", sprintf("%.3f", p.adjust))
        )
      )
    
    has_ontology <- "ONTOLOGY" %in% colnames(P.d$runScore)
    
    P.d$title <- P.d$runScore %>%
      mutate(
        title = if(has_ontology){ paste0(ONTOLOGY, ": ", Description) }else{Description}
      ) %>%
      group_by(Description) %>% 
      slice_head(n = 1) %>%
      ungroup()
    
  }
  
  P.d$geneList <- P.d$df %>%
    filter(
      ID == .$ID[1]
    ) %>%
    dplyr::select(gene, geneList) %>%
    tibble::deframe()
  
  P.d$heat <- P.d$geneList %>%
    enframe() %>%
    dplyr::rename(
      "gene" = "name",
      "geneList" = "value"
    ) %>%
    mutate(
      break_intervals = cut(geneList, breaks = 50, include.lowest = TRUE),
      interval_str = str_remove_all(as.character(break_intervals), "[()\\[\\]]")
    ) %>%
    separate(interval_str, into = c("Start", "End"), sep = ",", convert = TRUE) %>%
    group_by(break_intervals, Start, End) %>%
    summarise(
      geneList = median(geneList),
      .groups = "drop"
    ) %>%
    mutate(
      Interval = row_number(),
      xmax = sapply(Start, function(b) which.min(abs(P.d$geneList - b))),
      xmin = sapply(End, function(b) which.min(abs(P.d$geneList - b))),
      xmax = ifelse(Interval == 1, length(P.d$geneList), xmax)
    )
  
  return(P.d)
}


f.plot.enrich.gsea.runScore <- function(temp_res, typeID="ID", geneSetID=NULL, addGene=NULL, addSeg=F, addPval=F, label_format=30, ratio=NULL, theme="bw", ...){
  # data
  if(is.null(geneSetID)){
    geneSetID <- unique(temp_res$ID)
  }else{
    
    if(any(geneSetID %in% "all")){
      geneSetID <- unique(temp_res$ID)
    }
    
  }
  
  P.d <- f.gsea.pd(temp_res, typeID, geneSetID)
  
  # plot
  P1 <- ggplot(P.d$runScore, aes(x = x, y = runningScore, color = Description)) +
    geom_hline(
      yintercept = 0, color = "black",
      linetype = c("solid", "dashed")[2],
      show.legend = if_else(length(geneSetID)==1 & !any(geneSetID %in% "all"), F, T)
    ) +
    geom_line(
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = rainbow(length(unique(P.d$runScore$Description))),
      labels = P.d$df$Description %>% stringr::str_wrap(width = label_format)
    ) + 
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      # axis.title.x = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
    ) +
    scale_x_continuous(
      limits = c(0, max(P.d$runScore$x)),
      expand = c(0, 0, 0, 0)
    ) + 
    labs(
      y = "Running Enrichment Score",
      x = "Rank in Ordered Dataset"
    )
  
  # print(P1)
  
  # add Seg
  if(addSeg){
    P1 <- P1 + 
      geom_segment(
        data = P.d$min,
        aes(x = x, xend = x, y = 0, yend = y, color = Description),
        # color = "black",
        linetype = "dashed",
        linewidth = 0.6,
        show.legend = F
      )
  }
  
  if(length(geneSetID)==1 & !any(geneSetID %in% "all")){
    # add PValue
    if(addPval){
      P1 <- P1 +
        annotate(
          geom = "text",
          x = P.d$pvalue$x,
          y = P.d$pvalue$y,
          hjust = P.d$pvalue$hjust,
          vjust = P.d$pvalue$vjust,
          label = P.d$pvalue$label,
          size = 3.3,
          color = "black"
        )
      # geom_text(
      #   data = P.d$pvalue,
      #   aes(label = label, x = max(P.d$df$x)*0.8), 
      #   vjust = 0.5,    
      #   hjust = 0, 
      #   #lineheight = 0.7,    
      #   size = 3.3,           
      #   color = "gray20",
      # )
    }
    
    # add Title
    P1 <- P1 + 
      labs(
        title = P.d$title$title %>% stringr::str_wrap(width = label_format)
      )
  }
  
  # add Gene
  if(!is.null(addGene)){
    
    TopMetric <- c("geneList", "runningScore", "NES")[1]
    
    if(length(addGene) == 1 & grepl("Top", addGene)){
      P.d$addGene <- P.d$runScore %>%
        filter(
          position == 1
        ) %>%
        group_by(Description) %>%
        mutate(
          quantile_80 = quantile(.data[[TopMetric]], 0.8, na.rm = TRUE),
          sort_key = ifelse(quantile_80 > 0, -.data[[TopMetric]], .data[[TopMetric]])
        ) %>%
        arrange(sort_key) %>% 
        slice_head(n = as.numeric(stringr::str_extract(addGene, "\\d+"))) %>%  
        ungroup()
    }else{
      P.d$addGene <- P.d$runScore %>%
        filter(
          position == 1,
          gene %in% addGene
        )
    }
    
    P1 <- P1 +
      geom_text_repel(
        data = P.d$addGene,
        aes(x = x, y = runningScore, label = gene),
        force = 10,
        max.overlaps = 20,
        size = 3.3,
        color = "black",
        min.segment.length = 0, 
        # fontface = "italic",
        # arrow = ggplot2::arrow(
        #   angle = 20,
        #   length = ggplot2::unit(0.2, "cm"),
        #   ends = "first",
        #   type = c("closed", "open")[2]
        # )
      )
    
  }
  
  return(P1)
}


f.plot.enrich.gsea.heat <- function(temp_res, typeID="ID", geneSetID=NULL, colors=NULL, ratio=NULL, theme="bw", htHeight=0.3, ...){
  if(is.null(colors)){
    colors <- c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00")
  }
  # data
  P.d <- f.gsea.pd(temp_res, typeID, geneSetID)
  
  # plot
  P2 <- ggplot() +
    facet_wrap(
      ~ Description, 
      ncol = 1
    ) +
    geom_segment(
      data = subset(P.d$runScore, position == 1),
      aes(x = x, xend = x, y = 0, yend = htHeight, color = Description),
      # color = "black",
      show.legend = F
    ) +
    # geom_hline(
    #   yintercept = 0
    # ) +
    scale_color_manual(
      values = rainbow(length(unique(P.d$df$Description)))
    ) + 
    geom_rect(
      data = P.d$heat,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = htHeight/3, fill = geneList),
      alpha = 0.75,
      inherit.aes = FALSE,
      show.legend = F
    ) +
    scale_fill_gradientn(
      name = expression(paste(Log[2], " (Fold Change)")),
      values = scales::rescale(c(min(P.d$df$geneList), 0, max(P.d$df$geneList))),
      # colors = c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00")
      colors = colors
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # axis.title.x = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
    ) +
    scale_x_continuous(
      limits = c(0, max(P.d$runScore$x)),
      expand = c(0, 0, 0, 0)
    ) + 
    scale_y_continuous(
      expand = c(0.05, 0, 0.05, 0)
    ) +
    labs(
      x = "Rank in Ordered Dataset"
    )
  
  # print(P2)
  
  return(P2)
}


f.plot.enrich.gsea.rank <- function(temp_res, typeID="ID", geneSetID=NULL, colors=NULL, ratio=NULL, theme="bw", ...){
  if(is.null(colors)){
    colors <- c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00")
  }
  # data
  P.d <- f.gsea.pd(temp_res, typeID, geneSetID)
  
  P3 <- ggplot() +
    geom_segment(
      data = P.d$df, 
      aes(x = x, xend = x, y = geneList, yend = 0, color = geneList),
      show.legend = F
    ) +
    scale_color_gradientn(
      name = expression(paste(Log[2], " (Fold Change)")),
      values = scales::rescale(c(min(P.d$df$geneList), 0, max(P.d$df$geneList))),
      # colors = c("#1f78b4", "#a6cee3", "#f7f7f7", "#fdbf6f", "#ff7f00")
      colors = colors
    ) +
    # scale_color_gradient2(
    #   low = "#1f78b4", mid = "white", high = "#ff7f00", midpoint = 0
    # ) +
    geom_hline(
      yintercept = 0, color = "black",
      linetype = c("solid", "dashed", "longdash")[2]
    ) +
    f.theme(theme) +
    theme(
      aspect.ratio = ratio,
      panel.grid = element_blank(),
      # legend.position = "none"
    ) +
    scale_x_continuous(
      limits = c(0, max(P.d$runScore$x)),
      expand = c(0, 0, 0, 0)
    ) + 
    scale_y_continuous(
      limits = c(-max(abs(max(P.d$df$geneList)), abs(min(P.d$df$geneList))), max(abs(max(P.d$runScore$geneList)), abs(min(P.d$runScore$geneList))) ),
      expand = c(0, 0, 0, 0)
    ) + 
    labs(
      y = expression(paste(Log[2], " (Fold Change)")),
      x = "Rank in Ordered Dataset"
    )
  
  # print(P3)
  
  return(P3)
}



f.plot.enrich.gsea <- function(temp_res, ratio=NULL, typeID="ID", geneSetID=NULL, by="all", addGene=NULL, addSeg=F, addPval=F, label_format=30, colors=NULL, theme="bw", htHeight=0.3, ...){
  # theme
  de_x <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
  
  de_mar <- theme(
    plot.margin = ggplot2::margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="pt")
  )
  
  # plot
  P1 <- temp_res %>% 
    f.plot.enrich.gsea.runScore(
      typeID = typeID,
      geneSetID = geneSetID,
      addGene = addGene, 
      addSeg = addSeg,
      addPval = addPval,
      label_format = label_format,
      ratio = ratio$runScore, 
      theme = theme
    )
  
  P2 <- temp_res %>% 
    f.plot.enrich.gsea.heat(
      typeID = typeID,
      geneSetID = geneSetID,
      colors = colors, 
      ratio = ratio$heat,
      theme = theme,
      htHeight = htHeight
    )
  
  P3 <- temp_res %>% 
    f.plot.enrich.gsea.rank(
      typeID = typeID,
      geneSetID = geneSetID,
      colors = colors,
      ratio = ratio$rank,
      theme = theme
    )
  
  P <- switch (
    by,
    "runScore" = list(P1),
    "heat" = list(P2),
    "rank" = list(P3),
    "runScore + heat" = list(
      runScore = P1 + de_x + de_mar,
      heat = P2 + de_mar
    ),
    "runScore + rank" = list(
      runScore = P1 + de_x + de_mar,
      rank = P3 + de_mar
    ),
    "all" = list(
      runScore = P1 + de_x + de_mar,
      heat = P2 + de_x + de_mar,
      rank = P3 + de_mar
    ),
    stop("Invalid choice.")
  )
  
  P <- wrap_plots( gglist = P, ncol = 1) + 
    plot_layout(guides = 'collect') & theme(legend.position = 'right')
  
  return(P)
}


## pathway -----
f.pathStuff <- function(ID, ...){
  
  temp <- keggGet(dbentries = ID) %>% {.[[1]]}
  
  # GENE
  temp$GENE <- temp$GENE %>%
    tibble(raw = .) %>%
    mutate(
      group_id = ceiling(row_number() / 2),
      type = ifelse(row_number() %% 2 == 1, "ENTREZID", "annotation")
    ) %>%
    pivot_wider(
      id_cols = group_id,
      names_from = type,
      values_from = raw
    ) %>%
    mutate(
      SYMBOL = str_split(annotation, "; ", simplify = TRUE)[, 1],
      Description = str_extract(annotation, "(?<=; ).*(?= \\[KO:)"),
      Orthology = str_extract(annotation, "KO:K\\d{5}"),
      Enzyme = str_extract(annotation, "EC:[\\d.]+(?: [\\d.]+)*"), 
    ) %>%
    dplyr::select(-c("group_id", "annotation"))
  
  # COMPOUND
  temp$COMPOUND <- temp$COMPOUND %>%
    enframe( 
      name = "Compound",
      value = "Name"
    )
  
  temp$DRUG <- temp$DRUG %>%
    tibble(raw = .) %>%
    mutate(
      group_id = ceiling(row_number() / 2),
      type = ifelse(row_number() %% 2 == 1, "Drug", "Name")
    ) %>%
    pivot_wider(
      id_cols = group_id,
      names_from = type,
      values_from = raw
    )
  
  return(temp)
}


f.pathway <- function(ID, stuff, limit, bins, both.dirs, colors, pf, ...){
  
  colors <- colors %>%
    {
      list(
        low = map(., "low"),
        mid = map(., "mid"),
        high = map(., "high")
      )
    }
  
  pathview(
    gene.data = stuff$gene,
    cpd.data = stuff$cpd,
    pathway.id = ID,
    gene.idtype = c("UNIPROT", "SYMBOL", "ENSEMBL")[2],
    gene.annotpkg = "org.Hs.eg.db",
    species = "hsa",
    out.suffix = basename(pf),
    # kegg.dir = paste0(pf, "pathway/"),
    na.col = c("transparent", "gray")[1],
    # both.dirs = FALSE,
    # discrete = T,
    limit = limit,
    bins = bins,
    both.dirs = both.dirs,
    low = colors$low,
    mid = colors$mid,
    high = colors$high,
    # kegg.native = F,
    same.layer = F,
    map.null = F,
    map.symbol = T,
    map.cpdname = T,
    pdf.size = c(15, 9),
    rankdir = "LR", 
    # cex = 0.9,
    # text.width = 20,   
  )
  
  # dir(getwd())
  dir.create(dirname(pf), recursive = TRUE)
  
  list.files(
    path = getwd(),
    pattern = paste0(ID),
    recursive = F
  ) %>%
    walk(~{
      file.copy(
        from = file.path(getwd(), .x), 
        to = file.path(dirname(pf), basename(.x)),
        overwrite = TRUE
      )
      file.remove(.x)
    })
  
}
