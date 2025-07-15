

###########################################################  
###########################################################  4. Ecological Niche Model function

# Required packages
library(maxnet)
library(ranger)
library(e1071)
library(ParBayesianOptimization) 
library(dismo)
library(MLmetrics)
library(doParallel)
library(reshape2)

# Ecological Niche Modeling and evaluation function
ENM.develope <- function(DB.train, DB.syn, 
                         DB.valid, # Essential DB; For ENM validation
                         DB.test, # Essential DB; For ENM test
                         ENM.type = "Max", # Select one among Max, RFC, RFR, SVC, and SVR.
                         r.n = 1234 # Random seed setting
                         ){

###  Initial setting   
RNGkind("Mersenne-Twister")
res.var <- "OC" # response variable; dependent variable

#
if (!is.null(DB.syn)) {
  D.train <- rbind(DB.train, DB.syn)
} else {
  D.train <- DB.train
}

#  Hyperparameter range for each model
{
  MAX.para.range <- list(regmult = c(0.5, 3) # Default: 1
                         )  
  
  RFR.para.range <- list(mtry = c(2L,5L), #  default: sqrt(ncol)
                         ntree = c(1000L,2000L), # default: 500
                         nodesize = c(5L,15L) # default: 5 for regression
                         )
  
  RFC.para.range <- list(mtry = c(2L,5L), #  default: sqrt(ncol)
                         ntree = c(1000L,2000L), # default: 500
                         nodesize = c(5L,15L) # Default 1 for classification, 10 for probability
                         )
  
  SVR.para.range <- list(gamma =c(0.1, 2),
                         cost= c(0.1, 3.0),
                         epsilon =c(0.05,0.5)
                         )
  
  SVC.para.range <- list(gamma =c(0.1, 2),
                         cost= c(0.1, 3.0)
                         )
}

##### Bayesian optimization function (BOF) for Each ENM
{
  ev.score.f <- function(y_true, y_pred){ # By TSS
    ev.score <- MLmetrics::Sensitivity(y_true=y_true, y_pred=y_pred, positive = "1") 
    + MLmetrics::Specificity(y_true=y_true, y_pred=y_pred, positive = "1") -1
    if( is.na(ev.score) ){ ev.score <- -999 }
    return(ev.score + rnorm(1)/1000)
  }
  
  # BOF for Maxent
  MAX.BOF <- function(regmult){
    p=D.train[,res.var]
    m.data=D.train[!names(D.train) %in% res.var]
    target.model <- maxnet(p, m.data, f = maxnet.formula(p, m.data, classes = "default"), regmult=regmult, addsamplestobackground=T)
    e.model <- dismo::evaluate(p = predict(target.model, DB.valid, type="cloglog")[DB.valid$OC == 1], 
                               a = predict(target.model, DB.valid, type="cloglog")[DB.valid$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    pred.raw <- dismo::predict(target.model, DB.valid, type="cloglog")
    pred.OC <- ifelse(pred.raw >= t.model, 1, 0)
    
    ev.score <- ev.score.f(y_true=DB.valid$OC, y_pred=pred.OC) 
    return(list(Score =ev.score))
  } 
  
  ### BOF for RF regression
  RFR.BOF <- function(mtry, ntree, nodesize){
    form.in <- paste0(res.var, " ~ .")
    target.model <- ranger(formula(form.in), data = D.train, mtry= mtry, num.trees= ntree, min.node.size= nodesize)
    e.model <- dismo::evaluate(p = predict(target.model, DB.valid)$predictions[DB.valid$OC == 1], 
                               a = predict(target.model, DB.valid)$predictions[DB.valid$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    pred.raw <- dismo::predict(target.model, DB.valid)$predictions 
    pred.OC <- ifelse(pred.raw  >=  t.model, 1, 0)

    ev.score <- ev.score.f(y_true=DB.valid$OC, y_pred=pred.OC) 
    return(list(Score =ev.score))
  } 
  
  ### BOF for RF classification
  RFC.BOF <- function(mtry, ntree, nodesize){
    form.in <- paste0("factor(",res.var, ") ~ .")
    target.model <- ranger(formula(form.in), data = D.train, probability = T, mtry= mtry, num.trees= ntree, min.node.size= nodesize)
    e.model <- dismo::evaluate(p = predict(target.model, DB.valid)$predictions[,"1"][DB.valid$OC == 1], 
                               a = predict(target.model, DB.valid)$predictions[,"0"][DB.valid$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    pred.raw <- dismo::predict(target.model, DB.valid)$predictions[,"1"] 
    pred.OC <- ifelse(pred.raw  >=  t.model, 1, 0)
    
    ev.score <- ev.score.f(y_true=DB.valid$OC, y_pred=pred.OC) 
    return(list(Score =ev.score))
  } 
  
  ### BOF for SVM regression
  SVR.BOF <- function(gamma, cost, epsilon){
    form.in <- paste0(res.var, " ~ .")
    target.model <- svm(formula(form.in), data = D.train, gamma= gamma, cost= cost, epsilon= epsilon, 
                        kernel = "radial", scale = T, type="eps-regression")
    e.model <- dismo::evaluate(p = predict(target.model, DB.valid)[DB.valid$OC == 1], 
                               a = predict(target.model, DB.valid)[DB.valid$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    pred.raw <- dismo::predict(target.model, DB.valid)
    pred.OC <- ifelse(pred.raw  >= t.model, 1, 0)
    
    ev.score <- ev.score.f(y_true=DB.valid$OC, y_pred=pred.OC) 
    return(list(Score =ev.score))
  } 
  
  ### BOF for SVM classification
  SVC.BOF <- function(gamma, cost){
    form.in <- paste0("factor(",res.var, ") ~ .")
    target.model <- svm(formula(form.in), data = D.train, gamma= gamma, cost= cost, kernel = "radial", 
                        scale = T, type="C-classification", probability =T)
    e.model <- dismo::evaluate(p = attr(predict(target.model, DB.valid, probability =T), "probabilities")[,"1"][DB.valid$OC == 1], 
                               a = attr(predict(target.model, DB.valid, probability =T), "probabilities")[,"0"][DB.valid$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    
    pred.raw <- attr(dismo::predict(target.model, DB.valid, probability =T), "probabilities")[,"1"]# for SVM # 1(P) 기준 확률값
    pred.OC <- ifelse(pred.raw  >= t.model, 1, 0)
    ev.score <- ev.score.f(y_true=DB.valid$OC, y_pred=pred.OC) 
    return(list(Score =ev.score))
  } # BOF
}


##### Bayesian optimization performed in parallel 
{
  
  assign("D.train", D.train, envir = .GlobalEnv)
  assign("DB.valid", DB.valid, envir = .GlobalEnv)
  assign("ev.score.f", ev.score.f, envir = .GlobalEnv)  
  assign("res.var", res.var, envir = .GlobalEnv)  
  
  
  if( detectCores()-2 > 5 ){ use.clu.n <- 5 } else { use.clu.n<- detectCores()-2 }
  cl <- makeCluster(use.clu.n)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c('D.train', 'DB.valid', 'ev.score.f', 'res.var'))
  set.seed(r.n*2)
  if( ENM.type == "Max"){
    clusterEvalQ(cl,expr= { library(glmnet)})
    clusterEvalQ(cl,expr= { library(maxnet)})
    optObj <- bayesOpt(FUN = MAX.BOF, bounds = MAX.para.range, initPoints = 5, iters.n = 10, iters.k = 5, parallel = TRUE)
  }
  
  if( ENM.type == "RFR"){
    clusterEvalQ(cl,expr= { library(ranger)})  
    optObj <- bayesOpt(FUN = RFR.BOF, bounds = RFR.para.range, initPoints = 5, iters.n = 10, iters.k = 5, parallel = TRUE)
  }
  
  if( ENM.type == "RFC"){
    clusterEvalQ(cl,expr= { library(ranger)})  
    optObj <- bayesOpt(FUN = RFC.BOF, bounds = RFC.para.range, initPoints = 5, iters.n = 10, iters.k = 5, parallel = TRUE)
  }
  
  if( ENM.type == "SVR"){
    clusterEvalQ(cl,expr= { library(e1071)})  
    optObj <- bayesOpt(FUN = SVR.BOF, bounds = SVR.para.range, initPoints = 5, iters.n = 10, iters.k = 5, parallel = TRUE)
  }
  
  if( ENM.type == "SVC"){
    clusterEvalQ(cl,expr= { library(e1071)})  
    optObj <- bayesOpt(FUN = SVC.BOF, bounds = SVC.para.range, initPoints = 5, iters.n = 10, iters.k = 5, parallel = TRUE)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  rm(D.train, DB.valid, ev.score.f, res.var, envir = .GlobalEnv)
  
}

# Optimization result: hyperparameter selection (to data.frame)
para.selected <-  do.call(cbind, lapply(getBestPars(optObj), as.data.frame))
colnames(para.selected) <- names(getBestPars(optObj))
print("Complete : Bayesian optimization")

##### Development of the optimal model
{
  set.seed(r.n*2)
  #
  if( ENM.type == "Max"){
    p=D.train[,res.var]
    m.data=D.train[!names(D.train) %in% res.var]
    now.model <- maxnet(p, m.data, f = maxnet.formula(p, m.data, classes = "default"), 
                        regmult= para.selected$regmult, addsamplestobackground=T)
    e.model <- dismo::evaluate(p = predict(now.model, DB.test, type="cloglog")[DB.test$OC == 1], 
                               a = predict(now.model, DB.test, type="cloglog")[DB.test$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    #
    train.pred.raw <- predict(now.model, D.train, type="cloglog")
    valid.pred.raw <- predict(now.model, DB.valid, type="cloglog")
    test.pred.raw  <- predict(now.model, DB.test, type="cloglog")
    #
    rm(p, m.data)
  }
  
  if( ENM.type == "RFR"){
    form.in <- paste0(res.var, " ~ .")
    now.model <- ranger(formula(form.in), data = D.train, 
                        mtry= para.selected$mtry, num.trees= para.selected$ntree, min.node.size= para.selected$nodesize)
    e.model <- dismo::evaluate(p = predict(now.model, DB.test)$predictions[DB.test$OC == 1], 
                               a = predict(now.model, DB.test)$predictions[DB.test$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    #
    train.pred.raw <- dismo::predict(now.model, D.train)$predictions
    valid.pred.raw <- dismo::predict(now.model, DB.valid)$predictions
    test.pred.raw  <- dismo::predict(now.model, DB.test)$predictions
  }
  
  if( ENM.type == "RFC"){
    form.in <- paste0("factor(",res.var, ") ~ .")
    now.model <- ranger(formula(form.in), data = D.train, probability = T, 
                        mtry= para.selected$mtry, num.trees= para.selected$ntree, min.node.size= para.selected$nodesize)
    e.model <- dismo::evaluate(p = predict(now.model, DB.test)$predictions[,"1"][DB.test$OC == 1], 
                               a = predict(now.model, DB.test)$predictions[,"0"][DB.test$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    #
    train.pred.raw <- dismo::predict(now.model, D.train)$predictions[,"1"]
    valid.pred.raw <- dismo::predict(now.model, DB.valid)$predictions[,"1"]
    test.pred.raw  <- dismo::predict(now.model, DB.test)$predictions[,"1"] 
  }
  
  if( ENM.type == "SVR"){
    form.in <- paste0(res.var, " ~ .")
    now.model <- svm(formula(form.in), data = D.train, kernel = "radial", scale = T, type="eps-regression",
                     gamma= para.selected$gamma, cost= para.selected$cost, epsilon= para.selected$epsilon)
    e.model <- dismo::evaluate(p = predict(now.model, DB.test)[DB.test$OC == 1], 
                               a = predict(now.model, DB.test)[DB.test$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    #
    train.pred.raw <- dismo::predict(now.model, D.train)
    valid.pred.raw <- dismo::predict(now.model, DB.valid)
    test.pred.raw  <- dismo::predict(now.model, DB.test)
  }
  
  if( ENM.type == "SVC"){
    form.in <- paste0("factor(",res.var, ") ~ .")
    now.model <- svm(formula(form.in), data = D.train, kernel = "radial", scale = T, type="C-classification", probability =T,
                     gamma= para.selected$gamma, cost= para.selected$cost)
    e.model <- dismo::evaluate(p = attr(predict(now.model, DB.test, probability =T), "probabilities")[,"1"][DB.test$OC == 1], 
                               a = attr(predict(now.model, DB.test, probability =T), "probabilities")[,"0"][DB.test$OC == 0])
    t.model <- dismo::threshold(e.model, "spec_sens")
    #
    train.pred.raw <- attr(dismo::predict(now.model, D.train, probability =T), "probabilities")[,"1"]
    valid.pred.raw <- attr(dismo::predict(now.model, DB.valid, probability =T), "probabilities")[,"1"]
    test.pred.raw  <- attr(dismo::predict(now.model, DB.test, probability =T), "probabilities")[,"1"]
  }   
  
  ENM.D <- now.model
}
print("Complete : Development of the optimal model")

##### Calculation of evaluation metrics
{
  train.OC <- D.train$OC
  valid.OC <- DB.valid$OC
  test.OC <- DB.test$OC
  #
  train.pred.OC <- ifelse(train.pred.raw >= t.model, 1, 0)
  valid.pred.OC <- ifelse(valid.pred.raw >= t.model, 1, 0)
  test.pred.OC  <- ifelse(test.pred.raw  >= t.model, 1, 0)
  
  DB.ev.metrics <- function(v.true, v.pred.raw, v.pred.label, label.txt = NULL){
    
    m.Sensitivity <- MLmetrics::Sensitivity(y_true = v.true, y_pred =v.pred.label, positive = "1") 
    m.Specificity <- MLmetrics::Specificity(y_true = v.true, y_pred =v.pred.label, positive = "1") 
    
    DF <- data.frame("m.ACC" = MLmetrics::Accuracy(y_pred = v.pred.label, y_true= v.true), 
                     "m.AUC" = MLmetrics::AUC(y_pred = v.pred.raw, y_true= v.true), 
                     "m.PRAUC" = MLmetrics::PRAUC(y_pred = v.pred.raw, y_true= v.true), 
                     "m.Sens.1"= m.Sensitivity,
                     "m.Spec.1"= m.Specificity,
                     "m.Prec.1" = MLmetrics::Precision(y_true = v.true, y_pred =v.pred.label, positive = "1"),
                     "m.Recall.1" = MLmetrics::Recall(y_true = v.true, y_pred =v.pred.label, positive = "1"),
                     "m.TSS" = m.Sensitivity+m.Specificity-1,  # Specificity + Sensitivity - 1
                     "m.F1.1" = MLmetrics::F1_Score(y_true= v.true, y_pred= v.pred.label, positive = "1")
    )
    if( !is.null(label.txt)){
      names(DF) <- gsub("m.", paste0(label.txt,"."), names(DF))  
    } 
    return(DF)
  }
  
  tDB.ev.train <- DB.ev.metrics(v.true = train.OC, 
                                v.pred.raw = train.pred.raw,
                                v.pred.label = train.pred.OC, 
                                label.txt = "train")
  
  tDB.ev.valid <- DB.ev.metrics(v.true = valid.OC, 
                                v.pred.raw = valid.pred.raw,
                                v.pred.label = valid.pred.OC, 
                                label.txt = "valid")
  
  tDB.ev.test <- DB.ev.metrics(v.true = test.OC, 
                               v.pred.raw = test.pred.raw,
                               v.pred.label = test.pred.OC, 
                               label.txt = "test")
  
  
  tDB.info <- data.frame("ENM.method" = ENM.type, 
                         "N.train.OC0.ori" = table(D.train$OC)["0"],
                         "N.train.OC1.ori" = table(D.train$OC)["1"],
                         "N.syn.OC1.ori" = ifelse(is.null(nrow(DB.syn)), NA, nrow(DB.syn)),
                         "N.valid.OC0" = table(DB.valid$OC)["0"],
                         "N.valid.OC1" = table(DB.valid$OC)["1"],
                         "N.test.OC0" = table(DB.test$OC)["0"],
                         "N.test.OC1" = table(DB.test$OC)["1"]
  )
  
  ENM.evaluation <- cbind(tDB.info, tDB.ev.train, tDB.ev.valid, tDB.ev.test)
  
  ENM.opt.para <- cbind(data.frame("ENM.method" = ENM.type),
                        reshape2::melt(para.selected)
  )
  
}
print("Complete : Calculation of evaluation metrics")

##### Result save
list.ENM.results <- list(DB.train = DB.train, 
                         DB.syn = DB.syn,
                         DB.valid = DB.valid,
                         DB.test = DB.test,
                         ENM.D = ENM.D,
                         ENM.evaluation = ENM.evaluation,
                         ENM.opt.para = ENM.opt.para
                         )
save(list.ENM.results, file=paste0("Sample data.C5_ENM result_",ENM.type,".RData")) # Modeling result
write.csv(ENM.evaluation, file = paste0("ENM_evaluation metrics_",ENM.type,".csv"))
write.csv(ENM.opt.para, file = paste0("ENM_selected parameter_",ENM.type,".csv"))
          
print("Complete : Result save")
return(list.ENM.results)
}


# Input data description
# (1) DB.train: DB for ENM training, only partitioned data (before oversampling process)
# (2) DB.syn: DB for ENM training, only synthetic data generated through oversampling
# (3) DB.valid: DB for ENM validation
# (4) DB.test: DB for ENM test
# (5) ENM.type: Selected one of the five models; 1. Max (maximum entropy), 2. RFC (random forest classification), 
# 3. RFR (random forest regression), 4. SVC (support vector machine for classification), 5. SVR (support vector machine for regression)
# (6) r.n: Initial seed for random number generation

# Output description (one R file and two CSV files) 
# - (R) Sample data.C5_ENM result_[Model type]: Contains the list.ENM.results object, which includes the four input datasets, 
# the generated model (ENM.D), evaluation results (ENM.evaluation), and model parameters (ENM.opt.para).
# - (CSV) ENM_evaluation metrics_[ENM.type]: model evaluation results
# - (CSV) ENM_evaluationENM_selected parameter_[ENM.type]: model parameters used for ENM development

### Sample code 
load("Sample data.C3-4_sp1_p0.20_after oversampling.RData") 

# For maximum entropy model
ENM.max <- ENM.develope(DB.train =  DB.list.k1$train,
                        DB.syn = DB.list.k1$syn.SMOTE,
                        DB.valid = DB.list.k1$valid,
                        DB.test = DB.list.k1$test,
                        ENM.type = "Max",
                        r.n = 1234 
                        )

# For Random forest regression
ENM.RFR <- ENM.develope(DB.train =  DB.list.k1$train,
                        DB.syn = DB.list.k1$syn.SMOTE,
                        DB.valid = DB.list.k1$valid,
                        DB.test = DB.list.k1$test,
                        ENM.type = "RFR",
                        r.n = 1234
                        ) 
