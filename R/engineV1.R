# fixSEM v1
#' fixsem
#' @import future
#' @import listenv
#' @import kaefa
#' @import semTools
#' @import lavaan
#' @import progress
#' @import toOrdinal
#'
#' @param model a lavaan-style syntax
#' @param data data for calibration
#' @param group specify the group if you want
#' @param growth logical; your model is growth model? default is FALSE
#' @param remote see future::plan()
#' @param fastrun logical; default is FALSE. If you turn on to TRUE, it will be run with ML family.
#' @param method calibration method what can lavaan run
#'
#' @return calibrated lavaan models
#' @export
#'
#' @examples
#' \dontrun{
#'geocode("3817 Spruce St, Philadelphia, PA 19104")
#'geocode("Philadelphia, PA")
#'dat <- data.frame(value=runif(3),address=c("3817 Spruce St, Philadelphia, PA 19104","Philadelphia, PA","Neverneverland"))
#'geocode(dat)
#'}
  fixsem <- function(model, data, group = NULL, growth = F,
                     remote = getOption("kaefaServers"), fastrun = F,
                     method = if(fastrun) {((c("ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR")))}
                     else {((c("WLS", "DWLS", "WLSM", "WLSMV", "WLSMVS",
                             "ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR",
                             "ULS", "ULSM", "ULSMV", "ULSMVS")))}){
    fitPre <- lavaan::sem(model, data)
    CatOn <- T
    for(i in 1:length(data[fitPre@Data@ov$name])){
      if(length(levels(as.factor(data[fitPre@Data@ov$name][,i]))) >= 30){
        CatOn <- F
      }
    }
    
    group <- c("NULL", group)
    TickIterCount <- 0
    # ticktock counter
    for(groups in group){
      for(optim.method in c("nlminb", "BFGS","L-BFGS-B")){
        for(mimic in c("Mplus", "EQS", "lavaan")){
          for(stdlv in c(T, F)){
            for(meanstructure in c(T, F)){
              for(calibMethod in method){
                if(calibMethod %in% c("WLSM", "WLSMV", "WLSMVS", "ULSM", "ULSMV", "ULSMVS") && CatOn){
                  for(parameterisation in c("detla", "theta")){
                    TickIterCount <- TickIterCount + 1
                  }
                }
                TickIterCount <- TickIterCount + 1
              }
            }
          }
        }
      }
    }

    kaefa::aefaInit(remote)
    iterate <- T
    grandCount <- 0
    while(iterate){
      invisible(gc())
      pb <- progress_bar$new(
        format = "  (:spin) estimating :type SEM model [:bar] :percent (:current/:total) in :elapsed, eta: :eta, in :grandCount fix trials / Estimator: :calibMethod, MIMIC: :mimic, Optimiser: :optim.method, Group: :group",
        total = TickIterCount, clear = T, width= 170)
      grandCount <- grandCount + 1
      iterCount <- 0
      models <- listenv::listenv()
      invisible(gc())

      for(groups in group){
        if(groups == "NULL"){
          iterGroup <- NULL
          groupEqual <- NULL
        } else {
          iterGroup <- groups
          groupEqual <- c("loadings", "intercepts", "residuals", "lv.variances",
                          "lv.covariances")
        }
        for(optim.method in c("nlminb", "BFGS","L-BFGS-B")){
          for(mimic in c("Mplus", "EQS", "lavaan")){
            for(stdlv in c(T, F)){
              for(meanstructure in c(T, F)){
                for(calibMethod in method){
                  if(calibMethod %in% c("WLSM", "WLSMV", "WLSMVS", "ULSM", "ULSMV", "ULSMVS") && CatOn){
                    for(parameterisation in c("detla", "theta")){
                      invisible(gc())
                      iterCount <- iterCount + 1
                      pb$tick(tokens = list(grandCount = toOrdinal::toOrdinal(grandCount), optim.method = optim.method, mimic = mimic, calibMethod = calibMethod, type = 'Discrete', group = groups))
                      if(growth){
                        models[[iterCount]] %<-% tryCatch(lavaan::growth(model = model, data = data, method = calibMethod, ordered = names(data)[!names(data) %in% iterGroup], meanstructure=meanstructure, parameterization = parameterisation, std.lv = stdlv, mimic = mimic, optim.method = optim.method, group = iterGroup, group.equal = groupEqual), error = function(e) {})
                      } else {
                        models[[iterCount]] %<-% tryCatch(lavaan::sem(model = model, data = data, method = calibMethod, ordered = names(data)[!names(data) %in% iterGroup], meanstructure=meanstructure, parameterization = parameterisation, std.lv = stdlv, mimic = mimic, optim.method = optim.method, group = iterGroup, group.equal = groupEqual), error = function(e) {})
                      }
                    }
                  }
                  iterCount <- iterCount + 1
                  pb$tick(tokens = list(grandCount = toOrdinal::toOrdinal(grandCount), optim.method = optim.method, mimic = mimic, calibMethod = calibMethod, type = 'Continuous', group = groups))

                  invisible(gc())
                  if(growth){
                    models[[iterCount]] %<-% tryCatch(lavaan::growth(model = model, data = data, method = calibMethod, meanstructure=meanstructure, std.lv = stdlv, mimic = mimic, optim.method = optim.method, group = iterGroup, group.equal = groupEqual), error = function(e) {})
                  } else {
                    models[[iterCount]] %<-% tryCatch(lavaan::sem(model = model, data = data, method = calibMethod, meanstructure=meanstructure, std.lv = stdlv, mimic = mimic, optim.method = optim.method, group = iterGroup, group.equal = groupEqual), error = function(e) {})
                  }
                }
              }
            }
          }
        }
      }


      models <- .convergenceChecker(as.list(models))

      MIpower %<-% semTools::miPowerFit(models)
      # MIpower[MIpower$op != "=~",]
      if(NROW(MIpower) != 0){
        MIpower <- MIpower[MIpower$decision.pow == "NM",]
      } else {
        # STOP <- T
      }
      if(NROW(MIpower) != 0){
        fixedModel <- paste(model, '\n\t', paste(MIpower[which(max(MIpower$mi) == MIpower$mi),1:3], collapse = ""))
        rm(MIpower)
      } else {
        # STOP <- T
      }
      if(fixedModel == model){
        iterate <- F
      } else {
        model <- fixedModel
      }
    }
    models
  }


  #' @export
  .convergenceChecker <- function(fitlist){
    good <- list()
    count <- 0
    for(i in fitlist){
      if(!is.null(i)){
        if(i@optim$converged){
          count <- count + 1
          good[[count]] <- i
        } else {

        }
      }
    }
    
    if(NROW(good) == 0){
      stop('all model was non-convergence.')
    }
    
    value_RMSEA_lower <- listenv::listenv()
    valueBIC <- listenv::listenv()
    count <- 0
    pb <- progress_bar$new(
      format = "  (:spin) getting model fit indices [:bar] :percent (:current/:total) in :elapsed, eta: :eta",
      total = NROW(good), clear = T, width= 150)
    for(i in good){
      count <- count + 1
      pb$tick()
      value_RMSEA_lower[[count]] %<-% lavaan::fitmeasures(i, "rmsea.ci.lower")
    }

    value_RMSEA_lower <- unlist(as.list(value_RMSEA_lower))

    final_RMSEA_lower1 <- which(min(value_RMSEA_lower) == value_RMSEA_lower)[1]
    if(length(which(min(value_RMSEA_lower) == value_RMSEA_lower)) > 1){
      final_RMSEA_lower2 <- which(min(value_RMSEA_lower) == value_RMSEA_lower)[2]
    } else {
      final_RMSEA_lower2 <- final_RMSEA_lower1
    }
    if(final_RMSEA_lower1 != final_RMSEA_lower2){
      try(scaledRMSEA1 <- lavaan::fitMeasures(good[[final_RMSEA_lower1]], "rmsea.scaled"))
      try(scaledRMSEA2 <- lavaan::fitMeasures(good[[final_RMSEA_lower2]], "rmsea.scaled"))

      if(length(scaledRMSEA1) == 0 | length(scaledRMSEA2) == 0){

        try(scaledRMSEA1 <- lavaan::fitMeasures(good[[final_RMSEA_lower1]], "rmsea.robust"))
        try(scaledRMSEA2 <- lavaan::fitMeasures(good[[final_RMSEA_lower2]], "rmsea.robust"))

        if(length(scaledRMSEA1) == 0 | length(scaledRMSEA2) == 0){

          try(scaledRMSEA1 <- lavaan::fitMeasures(good[[final_RMSEA_lower1]], "rmsea"))
          try(scaledRMSEA2 <- lavaan::fitMeasures(good[[final_RMSEA_lower2]], "rmsea"))
        }
      }

      if(scaledRMSEA1 > scaledRMSEA2){
        good[[final_RMSEA_lower2]]
      } else {
        good[[final_RMSEA_lower1]]
      }
    } else {
      good[[final_RMSEA_lower1]]
    }

  }
