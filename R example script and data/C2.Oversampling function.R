###########################################################  
###########################################################  2. Oversampling function
library(smotefamily)

### Oversampling function

oversampling.f <- function(DB.train, Oversampling = "SWR", r.n = 1){
  ### Tuning the oversampling parameter
  n.min <- min(table(DB.train$OC))
  n.max <- max(table(DB.train$OC))
  n.length <- n.max-n.min # Number of samples required via oversampling
  if(n.min <= 5 ){ Kn = n.min -1 } else{ Kn = 5 }
  
  ##### Oversampling 
  o.result <- c() # For synthesized data 
  
  ### Oversampling 1: SWR
  if( Oversampling == "SWR" ){
    print("Applying SWR")  
    
    # Temporary DB
    tDB.minor <- DB.train[DB.train$OC == 1, ] # Minority class (P) #Synthesis target
    tDB.major <- DB.train[DB.train$OC == 0, ] # Majority class (A)
    
    # SWR 
    set.seed(919+r.n)
    temp.num <- sample(1:nrow(tDB.minor), size = n.length, replace = T) 
    o.result <- tDB.minor[temp.num,]
  }
  
  target.X <- (DB.train)[!names(DB.train) %in% "OC"]
  
  ### Oversampling 2: SMOTE
  if( Oversampling == "SMOTE" ){
    print("Applying SMOTE")  
    o.result <- try(SMOTE(X = target.X, target = DB.train$OC, K = Kn))
  }
  
  ### Oversampling 3: BLSMOTE
  if( Oversampling == "BLSMOTE" ){
    print("Applying BLSMOTE")  
    o.result <-  try(BLSMOTE(X = target.X, target = DB.train$OC, dupSize = 0, K = Kn))
  }
  
  ### Oversampling 4: ADASYN
  if( Oversampling == "ADASYN" ){
    print("Applying ADASYN")  
    o.result <- try(ADAS(X = target.X, target = DB.train$OC, K = Kn))
  }
  
  ### Oversampling 5: DBSMOTE
  if( Oversampling == "DBSMOTE" ){
    print("Applying DBSMOTE")  
    o.result <- try(DBSMOTE(X=target.X, target = DB.train$OC))
  }
  
  # If the data was not synthesized
  if(class(o.result) == "try-error"){
    o.result <- c()
  } 
  
  ### Final output: synthesized data
  new.data <- o.result$syn_data
  if(is.null(new.data)){
    new.data <- cbind(target.X[0, ], data.frame("OC" = numeric(0)))
  }
  names(new.data)[ncol(new.data)] <- "OC" 
  new.data$OC <- as.numeric(new.data$OC) # Column attribute change: Chr. -> Num.
  
  return(new.data)
}

# Input data description
# (1) DB.train: The partitioned training data (format: data.frame).
# (2) Oversampling: Select one among SWR, SMOTE, BLSMOTE, ADASYN, and DBSMOTE.
# (3) r.n: Initial seed for random number generation

# Output description
# -synthesized data via oversampling Oversampling technique (format: data.frame).
# -If no data was synthesized, the data.frame has zero rows.


### Sample code 
load("Sample data.C2_sp1_p0.20-after partition.RData") # Data generated through Data preprocessing code (C1) 
DB.train <- DB.list.k1$train
head(DB.train)

oversampling.f(DB.train, Oversampling = "SMOTE") # Example of successful data synthesis
oversampling.f(DB.train, Oversampling = "BLSMOTE") # Example of failed data synthesis

# Example data
DB.list.k1$syn.SMOTE <- oversampling.f(DB.train, Oversampling = "SMOTE", r.n = 1) 
DB.list.k2$syn.SMOTE <- oversampling.f(DB.train, Oversampling = "SMOTE", r.n = 10) 
save(DB.list.k1, DB.list.k2, 
     file= "Sample data.C3_sp1_p0.20_after oversampling.RData") 
