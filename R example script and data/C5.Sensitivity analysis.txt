

###########################################################  
###########################################################  5.1. Variable importance 

### (Model-agnostic) Variable importance function
library(maxnet) 
library(ranger) 
library(e1071) 
library(DALEX) 

ENM.varimp.f <- function(ENM, ENM.type, DB, res.var = "OC", r.n = 1234){

if(ENM.type == "Max"){ custom_predict <-function(model, ...) { predict(model, ..., type="cloglog") } } 
if(ENM.type == "RFR"){ custom_predict <-function(model, ...) { predict(model, ...)$prediction } }
if(ENM.type == "RFC"){ custom_predict <-function(model, ...) { predict(model, ...)$prediction[,"1"] } }
if(ENM.type == "SVR"){ custom_predict <-function(model, ...) { predict(model, ...) } }
if(ENM.type == "SVC"){ custom_predict <-function(model, ...) { attr(predict(model, ..., probability =T), "probabilities")[,"1"] } }

# model explainer
m.explainer <- DALEX::explain(ENM, 
                              data = DB[!names(DB) %in% res.var], 
                              y = DB[,res.var], 
                              label = ENM.type,
                              predict_function = custom_predict)

if( ENM.type =="Max" ){ m.explainer$model_info$type <- "regression" } 

# Permutation importance
set.seed(r.n)
t.var.imp <- model_parts(m.explainer, N=NULL)
VarImp <- data.frame(t.var.imp[!t.var.imp$variable %in% c("_baseline_", "_full_model_"), ]) 
VarImp <- VarImp[VarImp$permutation == 0, ]

return(VarImp)
}

# Input data description
# (1) ENM: The generated ecological niche model. Output model of C4
# (2) ENM.type: ENM.type: Selected one of the five models; 1. Max (maximum entropy), 2. RFC (random forest classification), 
# 3. RFR (random forest regression), 4. SVC (support vector machine for classification), 5. SVR (support vector machine for regression)
# (3) DB: The dataset used for modeling; used for calculating variable importance
# (4) res.var: Name of the response variable in the DB (format: character)
# (5) r.n: Initial seed for random number generation

# Output description (format: data.frame)
# - variable: DB 내 변수(모델의 설명변수)
# - permutation: 0 = dropout_loss average value
# - dropout_loss: 각 변수별로 구해진 Permutation-based variable importance value값
# - label: 모델 종류

### Sample code 
load("Sample data.C5_ENM result_Max.RData") # object name: list.ENM.results

DB.ENM.varimp <- ENM.varimp.f(ENM = list.ENM.results$ENM.D, 
                              ENM.type = list.ENM.results$ENM.evaluation$ENM.method, # Automatically selected (="Max")
                              DB = rbind(list.ENM.results$DB.train, list.ENM.results$DB.syn),
                              res.var = "OC",
                              r.n = 1234)

DB.ENM.varimp



###########################################################  
###########################################################  5.2. Variable importance-based distance 

### Variable importance-based distance function
library(reshape2)

varimp.dist <- function(DB, n.var = 5){ 
  ### Mean of variable importance over 10 repetitions (K)
  mean.DB <- aggregate(DB, dropout_loss ~ variable + label + overampling, mean) 
  
  ##### Distance calculation
  # Top N variables based on the Baseline model (in descending order of importance)
  enm.BM <- mean.DB[mean.DB$overampling == "Baseline", ]
  Bm.TN <- enm.BM$variable[order(enm.BM$dropout_loss, decreasing = T)[1:n.var] ]
  
  # Row-wise descending rank DB (rank 1 = largest value)
  td1 <- dcast(mean.DB, label + overampling ~ variable, value.var = "dropout_loss")
  DF.rank <- as.data.frame(t(apply(td1[unique(DB$variable)], 1, function(x) rank(-x, ties.method = "first")))) 
  rownames(DF.rank) <- td1$overampling
  
  # Weighted Euclidean distance function (using the 'Basline' row as reference)
  weighted.Euclidean.dist <- function(df, weights) {
    if(length(weights) != ncol(df)) {
      stop("Weight vector length must equal number of dataframe columns")
    }
    
    DF.ref <- as.numeric(df[rownames(df) == "Baseline", ])
    DF.other <-df[rownames(df) != "Baseline", ]
    
    # Compute weighted Euclidean distance between ref and each row from the second onward
    WED <- apply(DF.other, 1, function(row) {
      sqrt(sum(weights * (as.numeric(row) - DF.ref)^2))
    })
    return(WED)
  }
  
  # Weights based on variable ranks (simple inverse rank)
  weights <- c(1/(1:n.var)) 
  
  # Distance matrix
  dist.matrix <- weighted.Euclidean.dist(DF.rank[Bm.TN], weights)
  
  Result <- data.frame("ENM"= unique(DB$label), "oversampling"= names(dist.matrix), "dist"= dist.matrix)
  return(Result)
}

# Input data description
# (1) DB: Variable importance database for the models, which use the same ENM and presence proportion (format: data.frame)
# (2) n.var: Number of top variables to use (e.g., n.var = 5: use top 5 variables)

# Output description (format: data.frame)
# - ENM: Name of ecological niche model
# - oversampling: Name of used oversampling method
# - dist: Distance between each oversampled ENM and the baseline ENM

### Sample code 
DB <- read.csv("Sample data.C5_sp1_0.20_Max.csv") # sp 1 (Baetis fuscatus) / presence proportion: 20% / ENM type: Max (maximum entropy) 
head(DB)

varimp.dist(DB=DB, n.var = 5)
