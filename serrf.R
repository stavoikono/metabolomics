normalize = function(file_location = file.choose(new = FALSE), methods=c("SERRF","nomis"), scatter_plot = T, detectcores_ratio = 0.1){
  
  
  # file_location = file.choose(new = FALSE);
  # methods="SERRF";
  # scatter_plot = F
  
  if(identical(methods,"all")){
    methods = c("mTIC","sum","median","PQN","contrast","quantile","linear",'liwong','cubic','batch_ratio','batch_loess','SERRF','svm','nomis','bmis')
  }
  
  if(!"pacman" %in% rownames(installed.packages())){
    cat("Installing packages...\n")
    install.packages('pacman')
  }
  pacman::p_load(affy, parallel,ranger,caret,pcaMethods,ggplot2,tidyr,graphics,grDevices,Hmisc,gtools,cowplot,RColorBrewer,readr,plotly,stringr,GGally,dplyr,e1071,officer,bootstrap,pdist,metabolomics)
  remove_outlier = function(v){
    out = boxplot.stats(v)$out
    return(list(value = v[!v%in%out],index = which(v%in%out)))
  }
  loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.25, 1, by = 0.05), folds = 5){
    # Do model selection using mean absolute error, which is more robust than squared error.
    mean.abs.error <- numeric(length(span.vals))
    
    # Quantify error for each span, using CV
    loess.model <- function(x, y, span){
      loess(y ~ x, span = span, control=loess.control(surface="interpolate",statistics='exact'),family = "gaussian")
    }
    
    loess.predict <- function(fit, newdata) {
      predict(fit, newdata = newdata)
    }
    
    span.index <- 0
    
    for (each.span in span.vals) {
      span.index <- span.index + 1
      mean.abs.error[span.index] = tryCatch({
        y.hat.cv <- bootstrap::crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
        non.empty.indices <- !is.na(y.hat.cv)
        diff = (y[non.empty.indices] / y.hat.cv[non.empty.indices]) * mean(y[non.empty.indices])
        sd(diff)/mean(diff)
      },error = function(er){
        NA
      })
    }
    best.span <- span.vals[which.min(mean.abs.error)]
    if(length(best.span)==0){
      best.span = 0.75
    }
    
    best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="interpolate",statistics='exact'),family = "gaussian")
    
    return(list(best.model, min(mean.abs.error, na.rm = TRUE),best.span))
  }
  shiftData = function(ori,norm){
    ori.min = apply(ori,1,min,na.rm=T)
    norm.min = apply(norm,1,min,na.rm=T)
    return(norm - c(norm.min - ori.min))
  }
  RSD = function(data){
    return(apply(data,1,function(x){
      x = remove_outlier(x)[[1]]
      return(sd(x,na.rm=T)/mean(x,na.rm=T))
    }))
  }
  source("https://raw.githubusercontent.com/slfan2013/rcodes/master/read_data.R")
  data = read_data(file_location)
  e = data$e_matrix
  f = data$f
  p = data$p
  
  
  
  
  
  
  # check if sampleType is in the dataset
  if(!"sampleType" %in% colnames(p)){
    stop("Your data must have 'sampleType'. Please see example data for more information.")
  }else{
    # check if qc, sample are in the data.
    if(any(!c('qc','sample') %in% p$sampleType)){
      stop("The 'sampleType' must contain at least 'qc' and 'sample'. Please see example data for more information.")
    }
  }
  
  # check if time is in the datasheet
  if(!"time" %in% colnames(p)){
    stop("Your data must have 'time'. Please see example data for more information.")
  }
  # check if time has duplicated value.
  if(any(duplicated(p$time))){
    time_duplicated = duplicated(p$time)
    stop("Your dataset has ",sum(time_duplicated)," duplicated 'time' value. 'time' of each sample should be unique. The duplicated 'time' values are ", paste0(p$time[time_duplicated],collapse = ', '),'.')
  }
  
  # check if batch is in the datasheet
  if(!"batch" %in% colnames(p)){
    stop("Your data must have 'batch'. Please see example data for more information.")
  }
  
  # check with missing value.
  num_miss = sum(is.na(e))
  if(num_miss>0){
    cat(paste0("Your dataset has ",num_miss," missing values (i.e. empty cell). These missing values will be replaced by half of the minimum of non-missing value for each compound.\n"))
    for(i in 1:nrow(e)){
      
      e[i, is.na(e[i,])] = 0.5 * min(e[i, !is.na(e[i,])])
      
    }
  }
  
  # check with zero values.
  num_zero = sum(e == 0)
  if(num_zero>0){
    cat(paste0("Your dataset has ",num_zero," zeros. These zeros will be kept zeros in the final normalized data. \n"))
  }
  
  cat(paste0("In your data, \n"))
  cat(paste0("Number of compounds: ",nrow(f)," \n"))
  cat(paste0("Number of samples: ",sum(p$sampleType=='sample', na.rm = TRUE)," \n"))
  cat(paste0("Number of QCs: ",sum(p$sampleType=='qc', na.rm = TRUE)," \n"))
  
  
  
  
  p$sample_index = paste0('p',1:nrow(p))
  
  original_colnames_e = colnames(e)
  colnames(e) = p$sample_index
  
  # Empty sampleType will not be used for normalization, but will be put back in the final sheet.
  if(any(is.na(p$sampleType))){
    cat(paste0("There are ",sum(is.na(p$sampleType))," empty cells in the 'sampleType'. They will not be used for normalization, but will be put back in the final sheet. \n"))
  }
  
  
  # !!! when there is no na, this needs modification.
  p_empty_sampleType = p[is.na(p$sampleType), ]
  e_empty_sampleType = e[,is.na(p$sampleType)]
  if(class(e_empty_sampleType) == 'numeric'){
    e_empty_sampleType = matrix(e_empty_sampleType, ncol = 1)
  }
  colnames(e_empty_sampleType) = p_empty_sampleType$sample_index
  e = e[, !is.na(p$sampleType)]
  p = p[!is.na(p$sampleType), ]
  
  
  
  infinite_index = which(apply(e, 1, function(x){
    sum(is.infinite(x)) == length(x)
  }))
  if(!length(infinite_index)==0){
    e_infinite = e[infinite_index,]
    f_infinite = f[infinite_index,]
    
    e = e[-infinite_index,]
    f = f[-infinite_index,]
  }
  
  
  
  
  
  
  
  
  qc_RSDs = list()
  normalized_dataset = list()
  calculation_times = list()
  with_validate = any(!p$sampleType %in% c('qc','sample'))
  
  # split e, and p to different sample type.
  e_qc = e[, p$sampleType == 'qc']
  e_sample = e[, p$sampleType == 'sample']
  
  p_qc = p[p$sampleType == 'qc',]
  p_sample = p[p$sampleType == 'sample',]
  
  
  e_validates = list()
  
  p_validates = list()
  if(with_validate){
    
    
    val_RSDs = list()
    
    validate_types = unique(p$sampleType[!p$sampleType %in% c('qc','sample')])
    
    for(validate_type in validate_types){
      e_validates[[validate_type]] = e[, p$sampleType %in% validate_type]
      p_validates[[validate_type]] = p[p$sampleType %in% validate_type, ]
      
      val_RSDs[[validate_type]] = list()
    }
  }else{
    validate_types = NULL
  }
  
  
  aggregate_e = function(e_qc,e_sample,e_validates){
    e = do.call('cbind',c(list(e_qc, e_sample), e_validates))
    e = e[,order(as.numeric(gsub("p","",colnames(e))))]
    return(e)
  }
  
  # start
  start = Sys.time()
  normalized_dataset[['none']] = aggregate_e(e_qc,e_sample,e_validates)
  qc_RSDs[['none']] = RSD(e_qc)
  calculation_times[['none']] = Sys.time() - start
  cat("<!--------- raw data --------->\n")
  cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['none']],na.rm = TRUE),4)*100,"%.\n"))
  cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['none']]<0.2,na.rm = TRUE),".\n"))
  
  
  
  # pca = prcomp(t(e), scale = TRUE)
  # 
  # plot(pca$x[,1],pca$x[,2], col = factor(p$sampleType))
  
  
  if(with_validate){
    # for each type of validate.
    for(validate_type in validate_types){
      val_RSDs[[validate_type]][['none']] = RSD(e_validates[[validate_type]])
      cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['none']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['none']]<0.2,na.rm = TRUE),".\n"))
    }
  }
  
  if('mTIC' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['mTIC']] = tryCatch({
      cat("<!--------- mTIC --------->\n")
      mTIC_skip = FALSE
      if(!'compoundType' %in% colnames(f)){
        cat(paste0("warning: 'compountType' is not in the dataset. mTIC skipped.\n"))
        mTIC_skip = TRUE
      }else{
        mTIC_column = f[['compoundType']]
        if(!'known' %in% unique(f[['compoundType']])){
          cat(paste0("'known' (case-sensitive) is not found in the 'compoundType'. mTIC skipped.\n"))
          mTIC_skip = TRUE
        }
      }
      if(!mTIC_skip){
        start = Sys.time()
        index = mTIC_column %in% "known"
        sums = apply(e[index,], 2, sum, na.rm=T)
        mean_sums = mean(sums, na.rm = TRUE)
        e_norm = t(t(e)/(sums/mean_sums))
        
        qc_RSDs[['mTIC']] = RSD(e_norm[,p$sampleType=='qc'])
        calculation_times[['mTIC']] = Sys.time() - start
        cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['mTIC']],na.rm = TRUE),4)*100,"%.\n"))
        cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['mTIC']]<0.2,na.rm = TRUE),".\n"))
        
        
        
        
        if(with_validate){
          
          for(validate_type in validate_types){
            
            val_RSDs[[validate_type]][['mTIC']] = RSD(e_norm[,p$sampleType %in% validate_type])
            cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['mTIC']],na.rm = TRUE),4)*100,"%.\n"))
            cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['mTIC']]<0.2,na.rm = TRUE),".\n"))
          }
        }
      }else{
        e_norm = NA
      }
      
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('sum' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['sum']] = tryCatch({
      cat("<!--------- sum normalizaion --------->\n")
      start = Sys.time()
      sums = colSums(e,na.rm = T)
      e_norm = t(t(e)/sums*mean(sums,na.rm = T))
      qc_RSDs[['sum']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['sum']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['sum']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['sum']]<0.2,na.rm = TRUE),".\n"))
      
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['sum']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['sum']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['sum']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      
      
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('median' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['median']] = tryCatch({
      cat("<!--------- median normalizaion --------->\n")
      start = Sys.time()
      medians = apply(e, 2, median, na.rm = T)
      e_norm = t(t(e)/medians*median(medians,na.rm = T))
      qc_RSDs[['median']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['median']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['median']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['median']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        for(validate_type in validate_types){
          val_RSDs[[validate_type]][['median']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average Validate Sample RSD:",signif(median(val_RSDs[[validate_type]][['median']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% Validate Sample RSD:",sum(val_RSDs[[validate_type]][['median']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('PQN' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['PQN']] = tryCatch({
      cat("<!--------- PQN --------->\n")
      start = Sys.time()
      reference <- apply(e, 1, median, na.rm = TRUE)
      reference[reference==0] = 1
      quotient <- e/reference
      quotient.median <- apply(quotient, 2, median, na.rm = TRUE)
      e_norm <- t(t(e)/quotient.median)
      qc_RSDs[['PQN']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['PQN']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['PQN']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['PQN']]<0.2,na.rm = TRUE),".\n"))
      
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['PQN']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['PQN']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['PQN']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
  }
  
  if("contrast" %in% methods){
    Sys.sleep(2)
    normalized_dataset[['contrast']] = tryCatch({
      cat("<!--------- contrast --------->\n")
      start = Sys.time()
      threshold=1e-11
      e[e <= 0] <- threshold
      maffy.data <- maffy.normalize(data.matrix(e),
                                    subset=1:nrow(e),
                                    span=0.75,
                                    verbose=FALSE,
                                    family="gaussian",
                                    log.it=FALSE)
      subtract <- function(x){
        t(t(x)-apply(x,2,quantile,0.1,na.rm = TRUE))
      }
      e_norm <- subtract(maffy.data)
      rownames(e_norm) = rownames(e)
      colnames(e_norm) = colnames(e)
      qc_RSDs[['contrast']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['contrast']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['contrast']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['contrast']]<0.2,na.rm = TRUE),".\n"))
      
      
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['contrast']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['contrast']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['contrast']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      
      
      e_norm
    },error = function(error){
      error
    })
    
  }
  
  if("quantile" %in% methods){
    Sys.sleep(2)
    normalized_dataset[['quantile']] = tryCatch({
      cat("<!--------- quantile --------->\n")
      start = Sys.time()
      normalize.quantile <- get("normalize.quantiles",en=asNamespace("affy"))
      e_norm <- normalize.quantile(data.matrix(e))
      rownames(e_norm) <- rownames(e)
      colnames(e_norm) <- colnames(e)
      qc_RSDs[['quantile']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['quantile']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['quantile']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['quantile']]<0.2,na.rm = TRUE),".\n"))
      
      
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['quantile']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['quantile']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['quantile']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      
      
      
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('linear' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['linear']] = tryCatch({
      cat("<!--------- linear --------->\n")
      start = Sys.time()
      linear.baseline <- apply(e, 1, median,na.rm = TRUE)
      baseline.mean <- mean(linear.baseline,na.rm = TRUE)
      sample.means <- apply(e, 2, mean,na.rm = TRUE)
      linear.scaling <- baseline.mean/sample.means
      e_norm <- t(t(e) * linear.scaling)
      qc_RSDs[['linear']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['linear']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['linear']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['linear']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['linear']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['linear']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['linear']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('liwong' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['liwong']] = tryCatch({
      cat("<!--------- liwong --------->\n")
      start = Sys.time()
      average.intensity <- apply(e,2,mean,na.rm = TRUE)
      median.number <- round(ncol(e)/2 + 0.1)
      ordering <- order(average.intensity)
      median.sample.number <- ordering[median.number]
      median.sample <- e[,median.sample.number]
      e_norm <- vector()
      for(i in 1:ncol(e)){
        tryCatch({
          liwong.model <- normalize.invariantset(data=e[,i],
                                                 ref=median.sample,
                                                 prd.td=c(0.003,0.007))
          liwong.sample <- predict(liwong.model$n.curve$fit, e[,i])
        },error = function(er){
          liwong.sample = list();liwong.sample$y = e[,i]
        })
        
        e_norm <- cbind(e_norm,liwong.sample$y)
      }
      colnames(e_norm) = colnames(e)
      qc_RSDs[['liwong']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['liwong']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['liwong']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['liwong']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['liwong']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['liwong']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['liwong']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
  }
  
  if("cubic" %in% methods){
    Sys.sleep(2)
    normalized_dataset[['cubic']] = tryCatch({
      cat("<!--------- cubic --------->\n")
      start = Sys.time()
      e_norm <- normalize.qspline(e,samples=0.02,target=apply(e,1,mean,na.rm = TRUE),verbose = FALSE)
      rownames(e_norm) <- rownames(e)
      colnames(e_norm) <- colnames(e)
      qc_RSDs[['cubic']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['cubic']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['cubic']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['cubic']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['cubic']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['cubic']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['cubic']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('batch_ratio' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['batch_ratio']] = tryCatch({
      cat("<!--------- batch ratio --------->\n")
      start = Sys.time()
      e_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
      QC.index = p[["sampleType"]]
      batch = p[["batch"]]
      for(i in 1:nrow(f)){
        means = by(as.numeric(e[i,QC.index=='qc']),batch[QC.index=='qc'], mean, na.rm=T)
        mean_means = mean(means, na.rm = TRUE)
        for(b in 1:length(unique(batch))){
          e_norm[i,batch%in%unique(batch)[b]] = as.numeric(e[i,batch%in%unique(batch)[b]])/(means[[unique(batch)[b]]]/mean_means)
        }
      }
      rownames(e_norm) = rownames(e)
      colnames(e_norm) = colnames(e)
      qc_RSDs[['batch_ratio']] = RSD(e_norm[,p$sampleType=='qc'])
      calculation_times[['batch_ratio']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['batch_ratio']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['batch_ratio']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['batch_ratio']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['batch_ratio']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['batch_ratio']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      
      
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('batch_loess' %in% methods){
    Sys.sleep(2)
    normalized_dataset[['batch_loess']] = tryCatch({
      cat("<!--------- batch loess --------->\n(This may take a while...)\n")
      start = Sys.time()
      e_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
      QC.index = p[["sampleType"]]
      batch = p[["batch"]]
      time = as.numeric(p[["time"]])
      n_CV = 5
      QC.index.train = QC.index.test = list()
      n_CV = 5
      seed = 1
      set.seed(seed)
      QC.index.train = QC.index.test = list()
      e_qc_only = e[,p$sampleType=='qc']
      if(any(table(p$batch[p$sampleType=='qc']))<7){
        ratio = 0.7
      }else{
        ratio = 0.8
      }
      for(j in 1:n_CV){
        QC.index.train.temp = sample(1L:ncol(e_qc_only),round(ncol(e_qc_only)*ratio))
        QC.index.test.temp = c(1L:ncol(e_qc_only))[!c(1L:ncol(e_qc_only))%in%QC.index.train.temp]
        QC.index.train. = rep('sample',ncol(e_qc_only))
        QC.index.test. = rep('sample',ncol(e_qc_only))
        QC.index.train.[QC.index.train.temp] = 'qc'
        QC.index.test.[QC.index.test.temp] = 'qc'
        QC.index.train[[j]] = QC.index.train.
        QC.index.test[[j]] = QC.index.test.
      }
      
      
      loess_fun_cv = function(e,train.index = QC.index,test.index=NULL,batch,time){
        cl = makeCluster(detectCores() * detectcores_ratio)
        e_norm = parSapply(cl, X=1:nrow(e), function(i,e,train.index,batch,time,remove_outlier,loess_wrapper_extrapolate){
          
          e_norm = tryCatch({
            # for(i in 1:nrow(e)){
            line = e[i,]
            for(b in 1:length(unique(batch))){
              outlier_remove = remove_outlier(e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')])
              if(length(outlier_remove$index) == 0){
                lm = loess_wrapper_extrapolate(x=time[(batch %in% unique(batch)[b]) & (train.index=='qc')], y = e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')])[[1]]
              }else{
                lm = loess_wrapper_extrapolate(x=time[(batch %in% unique(batch)[b]) & (train.index=='qc')][-outlier_remove$index], y = e[i,(batch %in% unique(batch)[b]) & (train.index=='qc')][-outlier_remove$index])[[1]]
              }
              line[batch %in% unique(batch)[b]] = predict(lm,newdata  = time[batch %in% unique(batch)[b]])
              
              
              if(length(which(is.na(line[batch %in% unique(batch)[b]])))>0){
                for(j in which(is.na(line[batch %in% unique(batch)[b]]))){
                  time_notNA = time[batch %in% unique(batch)[b]][-which(is.na(line[batch %in% unique(batch)[b]]))]
                  closest_time = time_notNA[which.min(abs(time_notNA - time[batch %in% unique(batch)[b]][j]))]
                  line[batch %in% unique(batch)[b]][j] = line[batch %in% unique(batch)[b]][which(time[batch %in% unique(batch)[b]]==closest_time)]
                }
              }
              
            }
            
            if(sum(line<0)>(length(line)/5)){
              stop("too many negative value. LOESS failed.")
            }else{
              line[line<0] = runif(sum(line<0), min = max(c(median(e[i,], na.rm = TRUE) -  0.1 * sd(e[i,], na.rm = TRUE),0)), max = max(c(median(e[i,], na.rm = TRUE) +  0.1 * sd(e[i,], na.rm = TRUE),1)))
            }
            
            # if(length(which(is.na(line)))>0){
            #   for(j in which(is.na(line))){
            #     time_notNA = time[-which(is.na(line))]
            #     closest_time = time_notNA[which.min(abs(time_notNA - time[j]))]
            #     line[j] = line[which(time==closest_time)]
            #   }
            # }
            e[i,] / (line / median(e[i,], na.rm = TRUE))
            # if(sum((e[i,] / (line / median(e[i,], na.rm = TRUE)))<0)){
            #   stop(i)
            # }
            # }
            
            
          },error = function(er){
            e[i,]
          })
          return(e_norm)
        },e,train.index,batch,time,remove_outlier,loess_wrapper_extrapolate)
        stopCluster(cl)
        e_norm = t(e_norm)
        if(!is.null(test.index)){
          rsd = RSD(e_norm[,test.index=='qc'])
        }else(
          rsd = NULL
        )
        return(list(data = e_norm,rsd = rsd))
      }
      
      
      rsds = matrix(,nrow = nrow(e), ncol = n_CV)
      for(i in 1:n_CV){
        rsds[,i] = loess_fun_cv(e=e_qc_only,train.index=QC.index.train[[i]],test.index = QC.index.test[[i]],batch = batch[p$sampleType=='qc'],time=time[p$sampleType=='qc'])[[2]]
      }
      qc_RSD = apply(rsds,1,mean,na.rm = TRUE)
      e_norm = loess_fun_cv(e,train.index = QC.index,test.index=NULL,batch=batch,time=time)[[1]]
      rownames(e_norm) = rownames(e)
      colnames(e_norm) = colnames(e)
      qc_RSDs[['batch_loess']] = qc_RSD
      calculation_times[['batch_loess']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['batch_loess']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['batch_loess']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['batch_loess']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['batch_loess']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['batch_loess']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
    
    # p_val = c()
    # for(i in 1:nrow(normalized_dataset[['batch_loess']])){
    #
    #   x = normalized_dataset[['batch_loess']][i,]
    #
    #   p_val[i] = tryCatch({t.test(x ~ p$Gender)$p.value}, error = function(e){return(1)})
    #
    # }
    
    
    
  }
  
  if('SERRF' %in% methods){
    
    normalized_dataset[['SERRF']] = tryCatch({
      cat("<!--------- SERRF --------->\n(This may take a while...)\n")
      
      e_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
      QC.index = p[["sampleType"]]
      batch = p[["batch"]]
      time = p[["time"]]
      batch = factor(batch)
      num = 10
      start = Sys.time();
      
      cl = makeCluster(detectCores() * detectcores_ratio)
      
      # e_train = e_qc
      # e_target = e_validates[["Biorec"]]
      # num = 10
      # p_train = p_qc
      # p_target = p_validates[["Biorec"]]
      
      serrfR = function(train = e[,p$sampleType == 'qc'],
                        target = e[,p$sampleType == 'sample'],
                        num = 10,
                        batch. = factor(c(batch[p$sampleType=='qc'],batch[p$sampleType=='sample'])),
                        time. = c(time[p$sampleType=='qc'],time[p$sampleType=='sample']),
                        sampleType. = c(p$sampleType[p$sampleType=='qc'],p$sampleType[p$sampleType=='sample']),cl){
        
        
        all = cbind(train, target)
        normalized = rep(0, ncol(all))
        for(j in 1:nrow(all)){
          for(b in 1:length(unique(batch.))){
            current_batch = levels(batch.)[b]
            all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0] = rnorm(length(all[j,batch.%in%current_batch][all[j,batch.%in%current_batch] == 0]),mean = min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
            all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])] = rnorm(length(all[j,batch.%in%current_batch][is.na(all[j,batch.%in%current_batch])]),mean = 0.5*min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+1,sd = 0.1*(min(all[j,batch.%in%current_batch][!is.na(all[j,batch.%in%current_batch])])+.1))
          }
        }
        
        corrs_train = list()
        corrs_target = list()
        for(b in 1:length(unique(batch.))){
          
          current_batch = levels(batch.)[b]
          
          train_scale = t(apply(train[,batch.[sampleType.=='qc']%in%current_batch],1,scale))
          if(is.null(target[,batch.[!sampleType.=='qc']%in%current_batch])){
            target_scale = t(apply(target[,batch.[!sampleType.=='qc']%in%current_batch],1,scale))
          }else{
            target_scale = scale(target[,batch.[!sampleType.=='qc']%in%current_batch])
          }
          
          # all_scale = cbind(train_scale, target_scale)
          
          # e_current_batch = all_scale
          corrs_train[[current_batch]] = cor(t(train_scale), method = "spearman")
          corrs_target[[current_batch]] = cor(t(target_scale), method = "spearman")
          # corrs[[current_batch]][is.na(corrs[[current_batch]])] = 0
        }
        
        
        
        
        pred = parSapply(cl, X = 1:nrow(all), function(j,all,batch.,ranger, sampleType., time., num,corrs_train,corrs_target){
          # for(j in 1:nrow(all)){
          # j = j+1
          print(j)
          normalized  = rep(0, ncol(all))
          qc_train_value = list()
          qc_predict_value = list()
          sample_value = list()
          sample_predict_value = list()
          
          for(b in 1:length(levels(batch.))){
            current_batch = levels(batch.)[b]
            e_current_batch = all[,batch.%in%current_batch]
            corr_train = corrs_train[[current_batch]]
            corr_target = corrs_target[[current_batch]]
            
            
            corr_train_order = order(abs(corr_train[,j]),decreasing = TRUE)
            corr_target_order = order(abs(corr_target[,j]),decreasing = TRUE)
            
            sel_var = c()
            l = num
            while(length(sel_var)<(num)){
              sel_var = intersect(corr_train_order[1:l], corr_target_order[1:l])
              sel_var = sel_var[!sel_var == j]
              l = l+1
            }
            
            
            
            train.index_current_batch = sampleType.[batch.%in%current_batch]
            train_data_y = scale(e_current_batch[j, train.index_current_batch=='qc'],scale=F)
            train_data_x = apply(e_current_batch[sel_var, train.index_current_batch=='qc'],1,scale)
            
            if(is.null(dim(e_current_batch[sel_var, !train.index_current_batch=='qc']))){
              test_data_x = t(scale(e_current_batch[sel_var, !train.index_current_batch=='qc']))
            }else{
              test_data_x = apply(e_current_batch[sel_var, !train.index_current_batch=='qc'],1,scale)
            }
            
            train_NA_index  = apply(train_data_x,2,function(x){
              sum(is.na(x))>0
            })
            
            train_data_x = train_data_x[,!train_NA_index]
            test_data_x = test_data_x[,!train_NA_index]
            
            if(!class(test_data_x)=="matrix"){
              test_data_x = t(test_data_x)
            }
            
            good_column = apply(train_data_x,2,function(x){sum(is.na(x))==0}) & apply(test_data_x,2,function(x){sum(is.na(x))==0})
            train_data_x = train_data_x[,good_column]
            test_data_x = test_data_x[,good_column]
            if(!class(test_data_x)=="matrix"){
              test_data_x = t(test_data_x)
            }
            train_data = data.frame(y = train_data_y,train_data_x )
            
            if(ncol(train_data)==1){# some samples have all QC constent.
              norm = e_current_batch[j,]
              normalized[batch.%in%current_batch] = norm
            }else{
              colnames(train_data) = c("y", paste0("V",1:(ncol(train_data)-1)))
              model = ranger(y~., data = train_data)
              
              test_data = data.frame(test_data_x)
              colnames(test_data) = colnames(train_data)[-1]
              
              norm = e_current_batch[j,]
              
              
              
              norm[train.index_current_batch=='qc'] = e_current_batch[j, train.index_current_batch=='qc']/((predict(model, data = train_data)$prediction+mean(e_current_batch[j,train.index_current_batch=='qc'],na.rm=TRUE))/mean(all[j,sampleType.=='qc'],na.rm=TRUE))
              # norm[!train.index_current_batch=='qc'] =(e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model, data = test_data)$prediction + mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))/mean(e_current_batch[j,!train.index_current_batch=='qc'],na.rm=TRUE))
              
              norm[!train.index_current_batch=='qc'] =(e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))/(median(all[j,!sampleType.=='qc'],na.rm = TRUE)))
              norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]
              
              
              # plot(p$time[batch.%in%b][!train.index_current_batch=='qc'], (e_current_batch[j,!train.index_current_batch=='qc'])/((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, train.index_current_batch=='qc'],na.rm=TRUE))/(median(e_current_batch[j,!train.index_current_batch=='qc'],na.rm = TRUE))))
              
              
              norm[train.index_current_batch=='qc'] = norm[train.index_current_batch=='qc']/(median(norm[train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,sampleType.=='qc'],na.rm=TRUE))
              norm[!train.index_current_batch=='qc'] = norm[!train.index_current_batch=='qc']/(median(norm[!train.index_current_batch=='qc'],na.rm=TRUE)/median(all[j,!sampleType.=='qc'],na.rm=TRUE))
              norm[!is.finite(norm)] = rnorm(length(norm[!is.finite(norm)]),sd = sd(norm[is.finite(norm)],na.rm=TRUE)*0.01)
              
              
              
              
              out = boxplot.stats(norm, coef = 3)$out
              norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']%in%out] = ((e_current_batch[j,!train.index_current_batch=='qc'])-((predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'],na.rm=TRUE))-(median(all[j,!sampleType.=='qc'],na.rm = TRUE))))[norm[!train.index_current_batch=='qc']%in%out];
              norm[!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]=e_current_batch[j,!train.index_current_batch=='qc'][norm[!train.index_current_batch=='qc']<0]
              normalized[batch.%in%current_batch] = norm
              # points(current_time, norm, pch = (as.numeric(factor(train.index_current_batch))-1)*19, col = "blue", cex = 0.7)
              
              # qc_train_value[[b]] = train_data_y + mean(e_current_batch[j, train.index_current_batch=='qc'])
              # qc_predict_value[[b]] = predict(model,data = train_data)$predictions + mean(e_current_batch[j, train.index_current_batch=='qc'])
              # sample_value[[b]] = e_current_batch[j,!train.index_current_batch=='qc']
              # sample_predict_value[[b]] = predict(model,data = test_data)$predictions  + mean(e_current_batch[j, !train.index_current_batch=='qc'])
            }
            
            
            
            
            
            
          }
          
          
          # par(mfrow=c(1,2))
          # ylim = c(min(e[j,],norm), max(e[j,],norm))
          # plot(time.[sampleType.=='qc'], unlist(qc_train_value),col = "red",ylim = ylim,main=j)
          # points(time.[sampleType.=='qc'],unlist(qc_predict_value),col = "yellow")
          #
          # points(time.[!sampleType.=='qc'],unlist(sample_value),col = "blue")
          # points(time.[!sampleType.=='qc'],unlist(sample_predict_value),col = "green")
          #
          # plot(time.,normalized, col = factor(sampleType.), ylim = ylim,main=f$label[j])
          #
          # j = j + 1
          
          #
          
          
          
          # }
          
          
          
          
          return(normalized)
        },all,batch.,ranger, sampleType., time., num,corrs_train,corrs_target)
        
        
        
        
        normed = t(pred)
        
        normed_target = normed[,!sampleType.=='qc']
        
        
        for(i in 1:nrow(normed_target)){
          normed_target[i,is.na(normed_target[i,])] = rnorm(sum(is.na(normed_target[i,])), mean = min(normed_target[i,!is.na(normed_target[i,])], na.rm = TRUE), sd = sd(normed_target[i,!is.na(normed_target[i,])])*0.1)
        }
        for(i in 1:nrow(normed_target)){
          normed_target[i,normed_target[i,]<0] = runif(1) * min(normed_target[i,normed_target[i,]>0], na.rm = TRUE)
        }
        
        
        normed_train = normed[,sampleType.=='qc']
        
        
        for(i in 1:nrow(normed_train)){
          normed_train[i,is.na(normed_train[i,])] = rnorm(sum(is.na(normed_train[i,])), mean = min(normed_train[i,!is.na(normed_train[i,])], na.rm = TRUE), sd = sd(normed_train[i,!is.na(normed_train[i,])])*0.1)
        }
        for(i in 1:nrow(normed_train)){
          normed_train[i,normed_train[i,]<0] = runif(1) * min(normed_train[i,normed_train[i,]>0], na.rm = TRUE)
        }
        return(list(normed_train=normed_train,normed_target=normed_target))
      }
      
      
      serrf_normalized = e
      serrf_normalized_modeled = serrfR(train = e_qc, target = e_sample, num = num,batch. = factor(c(p_qc$batch, p_sample$batch)),time. = c(p_qc$time, p_sample$time),sampleType. = c(p_qc$sampleType, p_sample$sampleType),cl)
      
      serrf_qc = serrf_normalized_modeled$normed_train
      colnames(serrf_qc) = colnames(e_qc)
      serrf_sample = serrf_normalized_modeled$normed_target
      colnames(serrf_sample) = colnames(e_sample)
      
      serrf_cross_validated_qc = e_qc
      
      cv = 5
      RSDs = list()
      if(any(table(p_qc$batch))<7){
        ratio = 0.7
      }else{
        ratio = 0.8
      }
      
      test_indexes = split(1L:nrow(p_qc), c(1L:nrow(p_qc))%%cv)
      
      for(k in 1:cv){
        
        test_index = test_indexes[[k]]
        train_index = c(1L:nrow(p_qc))[-test_index]
        
        train_index = sample(1L:sum(p$sampleType=='qc'),round(sum(p$sampleType=='qc')*ratio))
        test_index = c(1L:sum(p$sampleType=='qc'))[!(c(1L:sum(p$sampleType=='qc'))%in%train_index)]
        
        
        while(length(unique(p_qc$batch[test_index]))<length(unique(batch))){
          train_index = sample(1L:nrow(p_qc),round(nrow(p_qc)*ratio))
          test_index = c(1L:nrow(p_qc))[!(c(1L:nrow(p_qc))%in%train_index)]
        }
        serrf_normalized_on_cross_validate = serrfR(train = e_qc[,train_index], target = e_qc[,test_index], num = num,batch. = factor(c(p_qc$batch[train_index],p_qc$batch[test_index])),time. = c(p_qc$time[train_index],p_qc$time[test_index]),sampleType. = rep(c("qc","sample"),c(length(train_index),length(test_index))),cl)
        
        serrf_cross_validated_qc[,test_index] = serrf_normalized_on_cross_validate$normed_target
        
        RSDs[[k]] = RSD(serrf_normalized_on_cross_validate$normed_target)
      }
      
      
      
      qc_RSD = apply(do.call("cbind",RSDs),1,mean)
      qc_RSDs[['SERRF']] = qc_RSD
      calculation_times[['SERRF']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['SERRF']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['SERRF']]<0.2,na.rm = TRUE),".\n"))
      
      
      serrf_validates = list()
      if(with_validate){
        
        
        for(validate_type in validate_types){
          
          serrf_validates[[validate_type]] = serrfR(train = e_qc, target = e_validates[[validate_type]], num = num,batch. = factor(c(p_qc$batch, p_validates[[validate_type]]$batch)),time. = c(p_qc$time, p_validates[[validate_type]]$time),sampleType. = rep(c("qc","sample"),c(nrow(p_qc),nrow(p_validates[[validate_type]]))),cl)$normed_target
          
          colnames(serrf_validates[[validate_type]]) = colnames(e_validates[[validate_type]])
          
          val_RSDs[[validate_type]][['SERRF']] = RSD(serrf_validates[[validate_type]])
          cat(paste0("Average ",validate_type," RSD:",signif(median( val_RSDs[[validate_type]][['SERRF']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum( val_RSDs[[validate_type]][['SERRF']]<0.2,na.rm = TRUE),".\n"))
        }
        aggregate_e(serrf_qc,serrf_sample,serrf_validates)
      }else{
        aggregate_e(serrf_qc,serrf_sample,NULL)
      }
      
      
      
    }, error = function(error_message){
      error_message
    })
    
  }
  
  if('svm' %in% methods){
    normalized_dataset[['svm']] = tryCatch({
      cat("<!--------- svm --------->\n(This may take a while...)\n")
      start = Sys.time()
      # https://pubs.rsc.org/en/content/articlepdf/2015/an/c5an01638j
      e_norm = matrix(,nrow=nrow(e),ncol=ncol(e))
      QC.index = p[["sampleType"]]
      batch = p[["batch"]]
      time = as.numeric(p[["time"]])
      n_CV = 5
      QC.index.train = QC.index.test = list()
      n_CV = 5
      seed = 1
      set.seed(seed)
      QC.index.train = QC.index.test = list()
      e_qc_only = e[,p$sampleType=='qc']
      if(any(table(p$batch[p$sampleType=='qc']))<7){
        ratio = 0.7
      }else{
        ratio = 0.8
      }
      for(j in 1:n_CV){
        QC.index.train.temp = sample(1L:ncol(e_qc_only),round(ncol(e_qc_only)*ratio))
        QC.index.test.temp = c(1L:ncol(e_qc_only))[!c(1L:ncol(e_qc_only))%in%QC.index.train.temp]
        QC.index.train. = rep('sample',ncol(e_qc_only))
        QC.index.test. = rep('sample',ncol(e_qc_only))
        QC.index.train.[QC.index.train.temp] = 'qc'
        QC.index.test.[QC.index.test.temp] = 'qc'
        QC.index.train[[j]] = QC.index.train.
        QC.index.test[[j]] = QC.index.test.
      }
      svm_fun_cv = function(e,train.index = QC.index,test.index=NULL,batch,time){
        cl = makeCluster(detectCores() * detectcores_ratio)
        e_norm = parSapply(cl, X=1:nrow(e), function(i,e,train.index,batch,time,remove_outlier,trainControl,train){
          # for(i in 1:nrow(e)){
          tryCatch({
            line = e[i,]
            for(b in 1:length(unique(batch))){
              outlier_remove = remove_outlier(e[i,(batch %in% unique(batch)[b]) & train.index=='qc'])
              if(length(outlier_remove$index) == 0){
                dta = data.frame(x = time[(batch %in% unique(batch)[b]) & train.index=='qc'], y = e[i,(batch %in% unique(batch)[b]) & train.index=='qc'])
                lm = train(y~., data=dta, method = "svmLinear", trControl = trainControl(method = "cv", savePred=T))
              }else{
                dta = data.frame(x = time[(batch %in% unique(batch)[b]) & train.index=='qc'][-outlier_remove$index], y = e[i,(batch %in% unique(batch)[b]) & train.index=='qc'][-outlier_remove$index])
                lm = train(y~., data=dta, method = "svmLinear", trControl = trainControl(method = "cv", savePred=T))
              }
              line[batch %in% unique(batch)[b]] = predict(lm,newdata  = data.frame(x = time[batch %in% unique(batch)[b]]))
            }
            if(length(which(is.na(line)))>0){
              for(j in which(is.na(line))){
                time_notNA = time[-which(is.na(line))]
                closest_time = time_notNA[which.min(abs(time_notNA - time[j]))]
                line[j] = line[which(time==closest_time)]
              }
            }
            e_norm = e[i,] / (line / median(e[i,], na.rm = TRUE))
          },error = function(er){
            e[i,]
          })
          
          # }
          
        },e,train.index,batch,time,remove_outlier,trainControl,train)
        stopCluster(cl)
        e_norm = t(e_norm)
        if(!is.null(test.index)){
          rsd = RSD(e_norm[,test.index=='qc'])
        }else(
          rsd = NULL
        )
        return(list(data = e_norm,rsd = rsd))
      }
      rsds = matrix(,nrow = nrow(e), ncol = n_CV)
      for(i in 1:n_CV){
        rsds[,i] = svm_fun_cv(e=e_qc_only,train.index=QC.index.train[[i]],test.index = QC.index.test[[i]],batch = batch[p$sampleType=='qc'],time=time[p$sampleType=='qc'])[[2]]
      }
      qc_RSD = apply(rsds,1,mean,na.rm = TRUE)
      e_norm = svm_fun_cv(e,train.index = QC.index,test.index=NULL,batch=batch,time=time)[[1]]
      rownames(e_norm) = rownames(e)
      colnames(e_norm) = colnames(e)
      qc_RSDs[['svm']] = qc_RSD
      calculation_times[['svm']] = Sys.time() - start
      cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['svm']],na.rm = TRUE),4)*100,"%.\n"))
      cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['svm']]<0.2,na.rm = TRUE),".\n"))
      if(with_validate){
        
        for(validate_type in validate_types){
          
          val_RSDs[[validate_type]][['svm']] = RSD(e_norm[,p$sampleType %in% validate_type])
          cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['svm']],na.rm = TRUE),4)*100,"%.\n"))
          cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['svm']]<0.2,na.rm = TRUE),".\n"))
        }
      }
      e_norm
    },error = function(error){
      error
    })
  }
  
  if('nomis' %in% methods){
    
    Sys.sleep(2)
    cat("<!--------- nomis --------->\n")
    normalized_dataset[['nomis']] = tryCatch({cat("nomis normalization.\n")
      start = Sys.time()
      skip_NOMIS = FALSE
      if(!'compoundType' %in% colnames(f)){
        cat(paste0("warning: 'compoundType' is not found. NOMIS is skipped.\n"))
        skip_NOMIS = TRUE
      }else{
        IS_column = f[['compoundType']]
        if(!'istd' %in% unique(f[['compoundType']])){
          cat(paste0("'istd' (case-sensitive) is not found in the 'compoundType'. NOMIS is skipped.\n"))
          skip_NOMIS = TRUE
        }
        
      }
      if(!skip_NOMIS){
        IS_column = f[['compoundType']]
        istd_index =  which(IS_column%in%'istd')
        
        inputdata = data.frame(Group = "A", t(log((e + sqrt(e^2 + 4)) * 0.5, base  = exp(1))))
        colnames(inputdata) = c("Group",f$label)
        rownames(inputdata) = paste0("S",1:nrow(inputdata))
        normed = exp(t(Normalise(inputdata,method = 'nomis',nc = istd_index)$output[,-1]))
        qc_RSDs[['nomis']] = RSD(normed[,p$sampleType=='qc'])
        e_norm = rbind(e[istd_index,],normed)
        calculation_times[['nomis']] = Sys.time() - start
        cat(paste0("Average QC RSD:",signif(median(qc_RSDs[['nomis']],na.rm = TRUE),4)*100,"%.\n"))
        cat(paste0("Number of compounds less than 20% QC RSD:",sum(qc_RSDs[['nomis']]<0.2,na.rm = TRUE),".\n"))
        if(with_validate){
          
          for(validate_type in validate_types){
            
            val_RSDs[[validate_type]][['nomis']] = RSD(e_norm[,p$sampleType %in% validate_type])
            cat(paste0("Average ",validate_type," RSD:",signif(median(val_RSDs[[validate_type]][['nomis']],na.rm = TRUE),4)*100,"%.\n"))
            cat(paste0("Number of compounds less than 20% ",validate_type," RSD:",sum(val_RSDs[[validate_type]][['nomis']]<0.2,na.rm = TRUE),".\n"))
          }
        }
        
        
      }else{
        e_norm = NA
      }
      e_norm
    },error = function(error){
      error
    })
    
    
    
  }
  
  
  
  if(!length(infinite_index)==0){
    e_full = rbind(e_infinite, e)
    f_full = rbind(f_infinite, f)
    
    e_full[infinite_index,] = NA
    e_full[-infinite_index,] = e
    e = e_full
    
    f_full[infinite_index,] = f_infinite
    f_full[-infinite_index,] = f
    f = f_full
    
  }
  
  
  
  
  comb_p = rbind(p_empty_sampleType, p)
  comb_e = cbind(e_empty_sampleType, e)
  
  comb_p$sample_index = as.numeric(substring(comb_p$sample_index,2))
  
  comb_e = comb_e[,order(comb_p$sample_index)]
  
  if(!length(infinite_index)==0){
    for(i in 1:length(normalized_dataset)){
      
      e_full = rbind(e_infinite, normalized_dataset[[i]])
      f_full = rbind(f_infinite, f)
      
      e_full[infinite_index,] = NA
      e_full[-infinite_index,] = normalized_dataset[[i]]
      normalized_dataset[[i]] = e_full
      
      
    }
  }
  
  
  
  normalized_dataset = normalized_dataset[!is.na(normalized_dataset)]  
  
  
  for(i in 1:length(normalized_dataset)){
    
    normalized_dataset[[i]] = cbind(e_empty_sampleType, normalized_dataset[[i]])
    
    normalized_dataset[[i]] = normalized_dataset[[i]][,order(comb_p$sample_index)]
    
  }
  
  comb_p = comb_p[order(sample_index),]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(grepl("\\\\",file_location)){
    comp = strsplit(file_location,"\\\\")[[1]]
  }else{
    comp = strsplit(file_location,"/")[[1]]
  }
  
  filename = gsub("\\.csv|\\.xlsx","",comp[length(comp)])
  root = paste0(paste0(comp[-length(comp)],collapse = "\\"),"\\")
  dir = paste0(root,filename," - normalization result")
  dir.create(dir)
  
  
  png(filename=paste0(dir,"\\","Bar Plot and PCA plot.png"), width = 2000, height = 1000 * ifelse(with_validate,3,2))
  qc_RSD_performance = sapply(qc_RSDs,median, na.rm = TRUE)
  qc_RSD_performance = sort(qc_RSD_performance,decreasing = TRUE)
  qc_RSD_performance_color = rep("grey",length(qc_RSD_performance))
  qc_RSD_performance_color[length(qc_RSD_performance_color)-1] = "red"
  qc_RSD_performance_color[length(qc_RSD_performance_color)] = "#ffbf00"
  qc_RSD_performance_color[names(qc_RSD_performance)=='none'] = 'black'
  
  
  val_RSD_performances = list()
  val_RSD_performances_color = list()
  if(with_validate){
    for(validate_type in validate_types){
      val_RSD_performances[[validate_type]] = sapply(val_RSDs[[validate_type]],median, na.rm = TRUE)
      val_RSD_performances[[validate_type]] = sort(val_RSD_performances[[validate_type]],decreasing = TRUE)
      val_RSD_performances_color[[validate_type]] = rep("grey",length(val_RSD_performances[[validate_type]]))
      val_RSD_performances_color[[validate_type]][length(val_RSD_performances_color[[validate_type]])-1] = "red"
      val_RSD_performances_color[[validate_type]][length(val_RSD_performances_color[[validate_type]])] = "#ffbf00"
      val_RSD_performances_color[[validate_type]][names(val_RSD_performances[[validate_type]])=='none'] = 'black'
    }
    layout(matrix(c(rep(1:(length(validate_types)+1), each = 2), length(validate_types)+2, length(validate_types)+3), ncol = 2, byrow = TRUE))
  }else{
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  }
  
  par(lty = 0)
  par(mar=c(5,4,4,2)*3)
  bp = barplot(qc_RSD_performance*100, main="QC RSD", xlab="", ylab="RSD (%)",col = qc_RSD_performance_color,width = 1,las=2,cex.axis =5, cex.names=5,cex.main = 5)
  text(bp[which(names(qc_RSD_performance)=='none'),1], qc_RSD_performance['none']*100, paste0(signif(qc_RSD_performance['none'],4)*100,"%"), cex = 5, pos = 3)
  text(bp[nrow(bp),1], qc_RSD_performance[length(qc_RSD_performance)]*100, paste0(signif(qc_RSD_performance[length(qc_RSD_performance)],4)*100,"%"), cex = 5, pos = 3)
  
  
  if(with_validate){
    for(validate_type in validate_types){
      
      par(lty = 0)
      par(mar=c(5,4,4,2)*3)
      bp = barplot(val_RSD_performances[[validate_type]]*100, main=paste0(validate_type," Sample RSD"), xlab="", ylab="RSD (%)",col = val_RSD_performances_color[[validate_type]],width = 1,las=2,cex.axis =5, cex.names=5,cex.main = 5)
      text(bp[which(names(val_RSD_performances[[validate_type]])=='none'),1], val_RSD_performances[[validate_type]]['none']*100, paste0(signif(val_RSD_performances[[validate_type]]['none'],4)*100,"%"), cex = 5, pos = 3)
      text(bp[nrow(bp),1], val_RSD_performances[[validate_type]][length(val_RSD_performances[[validate_type]])]*100, paste0(signif(val_RSD_performances[[validate_type]][length(val_RSD_performances[[validate_type]])],4)*100,"%"), cex = 5, pos = 3)
      
    }
  }
  
  comb_p$sampleType[is.na(comb_p$sampleType)] = "NA"
  pca_color = factor(comb_p$sampleType, levels = c('sample','qc',validate_types,"NA"))
  dots = c(1,16,rep(16, length(validate_types)),1)[as.numeric(pca_color)]
  if(!length(infinite_index)==0){
    
    sds = apply(comb_e[-infinite_index,],1,sd)
    pca_before = prcomp(t(comb_e[-infinite_index,][!sds==0,]),scale. = TRUE)
  }else{
    sds = apply(comb_e,1,sd)
    pca_before = prcomp(t(comb_e[!sds==0,]),scale. = TRUE)
  }
  
  
  
  par(mar=c(4,2,4,2)*3)
  plot(pca_before$x[,1],pca_before$x[,2], col = factor(comb_p$sampleType, levels = c('sample','qc',validate_types,"NA")),main = 'Before',xlab='raw data',cex.lab=5,yaxt='n', cex.axis=5, cex.main=5, cex.sub=5,ylab="", xaxt='n',cex = 5,pch = dots)
  
  if(!length(infinite_index)==0){
    sds = apply(normalized_dataset[[names(qc_RSD_performance)[length(qc_RSD_performance)]]][-infinite_index,],1,sd)
    
    pca_after = prcomp(t(normalized_dataset[[names(qc_RSD_performance)[length(qc_RSD_performance)]]][-infinite_index,][!sds==0,]),scale. = TRUE) 
  }else{
    sds = apply(normalized_dataset[[names(qc_RSD_performance)[length(qc_RSD_performance)]]],1,sd)
    
    pca_after = prcomp(t(normalized_dataset[[names(qc_RSD_performance)[length(qc_RSD_performance)]]][!sds==0,]),scale. = TRUE)
  }
  
  plot(pca_after$x[,1],pca_after$x[,2], col = pca_color,main = 'After',xlab = names(qc_RSD_performance)[length(qc_RSD_performance)],cex.lab=5, cex.axis=5, cex.main=5, cex.sub=5,ylab="",yaxt='n', xaxt='n',cex = 5,pch = dots)
  
  dev.off()
  
  
  dir.create(paste0(dir,"\\normalized datasets"))
  for(i in 1:length(normalized_dataset)){
    method = names(normalized_dataset[i])
    if(identical(class(normalized_dataset[[i]]),"matrix")){
      
      colnames(normalized_dataset[[i]]) = comb_p$label
      
      
      normalized_dataset[[i]][data$e_matrix==0] = 0
      normalized_dataset[[i]][is.na(data$e_matrix)] = NA
      
      
      
      fwrite(data.table(label = f$label,normalized_dataset[[i]]),paste0(dir,"\\normalized datasets\\normalized by - ",method,'.csv'))
    }else{
      cat("Error in ",method,": ", normalized_dataset[[i]])
    }
  }
  if(!length(infinite_index)==0){
    for(i in 1:length(qc_RSDs)){
      qc_RSD_full = rep(0,nrow(f))
      qc_RSD_full[-infinite_index] = qc_RSDs[[i]]
    }
  }
  
  
  
  
  
  fwrite(data.table(label = f$label, do.call('cbind',qc_RSDs)),paste0(dir,"//","QC - RSD.csv"))
  if(with_validate){
    
    for(validate_type in validate_types){
      
      
      
      if(!length(infinite_index)==0){
        
        for(i in 1:length(val_RSDs[[validate_type]])){
          qc_RSD_full = rep(0,nrow(f))
          val_RSDs[[validate_type]][-infinite_index] = val_RSDs[[validate_type]][[i]]
        }
      }
      
      
      
      
      fwrite(data.table(label = f$label, do.call('cbind',val_RSDs[[validate_type]])),paste0(dir,"//","Validate Samples - ",validate_type," - RSD.csv"))
      
    }
    
    
  }
  
  
  
  
  
  
  
  if(scatter_plot){ # time consuming.
    comb_p$sampleType[is.na(comb_p$sampleType)] = "NA"
    cat("Generating scatter plots for each compounds. Please be patient...\n")
    dir.create(paste0(dir,"\\scatter plots"))
    for(i in 1:length(normalized_dataset)){
      
      method = names(normalized_dataset)[i]
      
      if(identical(class(normalized_dataset[[i]]),"matrix")){
        normalized = normalized_dataset[[i]]
        dir.create(paste0(dir,"\\scatter plots\\",method))
        e =     normalized_dataset[['none']] 
        for(j in 1:nrow(e)){
          
          if(!j %in% infinite_index){
            
            
            png(paste0(dir,"\\scatter plots\\",method,"\\",j,"th.png"), width = 480*2, height = 480)
            par(mfrow=c(1,2))
            ylim = c(min(e[j,],normalized[j,]), max(e[j,],normalized[j,]))
            if(Inf %in% ylim){
              ylim[2] = max(e[j,!is.infinite(e[j,])],normalized[j,!is.infinite(normalized[j,])])*1.1
            }
            if(sum(is.na(ylim))<1){
              plot(comb_p$time,normalized_dataset[['none']][j,], col = factor(comb_p$sampleType, levels = c('sample','qc',validate_types, "NA")), ylim = ylim, main = f$label[j])
              plot(comb_p$time,normalized[j,], col = factor(comb_p$sampleType, levels = c('sample','qc',validate_types, "NA")), ylim = ylim)
            }
            dev.off()
            
            
            
          }
          
          
        }
      }}
  }
  
  if(!exists("serrf_cross_validated_qc")){
    serrf_cross_validated_qc = FALSE
  }
  
  
  cat(paste0("Good Job! All the normalizations are finished!\nPlease check your folder: '",dir,"'.\n"))
  return(list(normalized_dataset = normalized_dataset,qc_RSDs = qc_RSDs,calculation_times=calculation_times,serrf_cross_validated_qc = serrf_cross_validated_qc))
}

o = normalize(methods = "batch_loess",scatter_plot = F,detectcores_ratio = 1)
# methods can be either "all" a c() which can include mTIC,sum,median,PQN,contrast,quantile,linear,liwong,cubic,batch_ratio,batch_loess,SERRF,svm,nomis or bmis.

df_loess <- o$normalized_dataset$batch_loess

df_loess <- data.frame(t(df_loess))
names(df_loess) <- non_zero_met

df_serrf <- read.csv("./data_files/df_serrf2.csv")

df_loess2 <- cbind(df_serrf[,1:4],df_loess)

write.csv(df_loess2, "./data_files/df_loess_bo_sample_070421.csv", row.names = F)
