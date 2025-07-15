###########################################################  
###########################################################  1.1 Outlier detection using Hampel filter

### Hampel filter function  
outlier.Hampel <- function(DB){  
  for(i in 1:ncol(DB)){
    lower.limit <- median(DB[,i]) - 3 * mad(DB[,i], constant = 1) 
    upper.limit <- median(DB[,i]) + 3 * mad(DB[,i], constant = 1)
    outlier.row <- which(DB[,i] < lower.limit | DB[,i] > upper.limit)
    
    if(i == 1){
      outlier.row2 <- outlier.row 
    } else {
      outlier.row2 <- c(outlier.row2, outlier.row) 
    }
  } 
  remove.row <- sort(unique(outlier.row2)) # 
  DB2 <- DB[-remove.row,] 
  return(DB2)
}

# Input data description
# (1) DB: environmental DB (format: data.frame)

# Output description
# - DB with outlier instances removed (format: data.frame)

### Sample code 
DB.sp.env <- read.csv("Sample data.C1_D1.raw data.csv", fileEncoding = "EUC-KR") # 99 rows
head(DB.sp.env)

env.DB <- DB.sp.env[7:14] # Only for Env. data

outlier.Hampel(env.DB) # 51 rows


###########################################################  
###########################################################  1.2. variable selection 

### Variable selection function
library(usdm) # Package for variable selection

variable.selection <- function(DB, th.r = 0.7, th.VIF = 5){
  vs.cor <- vifcor(DB, th= th.r) # by Pearson's coefficient
  DB.cor <- exclude(DB, vs.cor)
  vs.vif <- vifstep(DB.cor, th= th.VIF) # by Variance Inflation Factor
  DB.vif <- exclude(DB.cor, vs.vif)
  return(DB.vif)
} 

# Input data description
# (1) DB: environmental DB (format: data.frame)

# Output description
# - DB excluding variables with multicollinearity (format: data.frame)

### Sample code 
DB.sp.env <- read.csv("Sample data.C1_D1.raw data.csv", fileEncoding = "EUC-KR") # 99 rows
env.DB <- DB.sp.env[7:14] # Only for Env. data

variable.selection(env.DB) 

###########################################################  
###########################################################  1.3. Environmentally filter 

# This function is a modified version of the built-in function ‘VarelaSample’ 
# from the ‘OccurrenceManagement’ function in the ‘megaSDM’ package.

# The 'dplyr' package is required.
library(dplyr) 

### Environmental filtering function  
VarelaSample.modified <- function (EnvOccur, no_bins = 25) {
  
  # (Modifications) Generate dummy column and setting to use PCA 
  {
    EnvOccur <- cbind(data.frame(x = seq(nrow(EnvOccur)), 
                                 y = seq(nrow(EnvOccur))), 
                      EnvOccur)
    PCA ="Y"
    PCAxes = "Y"
  }
  
  # The code below is the same as the original.
  ClimOccur <- EnvOccur
  ClimOccur <- ClimOccur[stats::complete.cases(ClimOccur), ]
  if (PCA == "Y") {
    PCAEnv <- stats::prcomp(ClimOccur[, 3:ncol(ClimOccur)], scale = TRUE)
    PCAImp <- summary(PCAEnv)$importance
    #Determine the number of PC axes to use for subsampling
    if (is.numeric(PCAxes)) {
      NumberAxes <- PCAxes
    } else {
      NumberAxes <- max(2, min(which(PCAImp[3,] > 0.95)))
    }
    
    #Add PCA values to the unsubsampled data frame
    EnvOccur <- data.frame(cbind(ClimOccur[, 1:2], PCAEnv$x[, 1:NumberAxes]))
  }
  
  #make a landing spot for the bin membership vectors
  
  #cycle through all of the environmental variables (columns 3 to end)
  nsamples <- c()
  for (j in 1:length(no_bins)) {
    out_ptz <- EnvOccur[, 1:2]
    for(i in 3:length(names(EnvOccur))) { 
      #make a data frame that is this variable with no NA values
      k <- EnvOccur[!is.na(EnvOccur[, i]), i]
      #calculate the observed range of this variable
      rg <- range(k)
      #figure out the resolution from the number of bins
      res <- (rg[2] - rg[1]) / no_bins[j] 
      #rescale the axis by the range and bin size, so the value is just a
      #number from 1 to no_bins for its bin membership
      d <- (EnvOccur[, i] - rg[1]) / res
      #d is now a vector of values ranging from 0 to no_bins
      f <- ceiling(d)
      #f is a vector of bin membership
      f[f == 0] <- 1 #move the zeros into the 1 bin
      #correct the name of the vector, so it will carry over to the output
      names(f) <- names(EnvOccur)[i]
      #add the bin membership vector to the output df for this section
      out_ptz <- cbind(out_ptz, f)
      #get the names correct
      names(out_ptz)[length(names(out_ptz))] <- names(EnvOccur)[i]
    }
    
    #subsample the bin membership df to come up with the filled bins
    sub_ptz <- dplyr::distinct(out_ptz[, -1:-2])
    
    #count the number of filled bins
    no_grps <- nrow(sub_ptz)
    #add a column with the group membership number; this number is arbitrary
    sub_ptz$grp <- c(1:no_grps)
    
    #join the out_ptz with the subsample to capture the group membership info
    #note: join() will automatically match the variable names from these two dfs
    out_ptz <- suppressMessages(dplyr::left_join(out_ptz, sub_ptz))
    #out_ptz now has a group membership  for each input point
    
    #select a random point for each group -- this is an additional improvement on the
    #Varela et al. function, because it does not pick the same points each time.
    
    #make a landing spot for the data
    final_out <- data.frame(x = numeric(), y = numeric())
    
    #cycle through each group
    for(i in 1:no_grps) {
      
      #subset to the members of the ith group, keep only the Latitude and Longitude
      grp_mbrs <- out_ptz[out_ptz$grp == i, c(1, 2)]
      
      #pick one of these group members to output
      grp_out <- grp_mbrs[sample(1:nrow(grp_mbrs), 1), ]
      #bind this sampled row to the output df
      final_out <- rbind(final_out, grp_out)
    }
    
    #return the subsampled points as a df of Latitude and Longitude values
    final_out <- data.frame(x = final_out[, 1], y = final_out[, 2])
    final_out <- merge(final_out, ClimOccur, by = c("x", "y"), all.x = TRUE)
    nsamples <- c(nsamples, nrow(final_out))
  }
  
  if (length(no_bins) == 1) {
    
    
    # (Modifications) Generate dummy column
    {
      final_out <- final_out[-1]
      names(final_out)[1] <- "Row.ID"
    }
    
    return(final_out)
  } else {
    return(data.frame(NumberofSamples = nsamples, NumberOfClimateBins = no_bins))
  }
}

# Input data description
# (1) EnvOccur: environmental DB (format: data.frame)
# (2) no_bins: The bin count per environmental combination for Varela subsampling

# Output description
# - DB after Environmental filtering. A column (i.e., Row.ID) for row IDs has been added (format: data.frame)

### Sample code 
DB.sp.env <- read.csv("Sample data.C1_D1.raw data.csv", fileEncoding = "EUC-KR")
env.DB <- DB.sp.env[7:14] # Only for Env. data

VarelaSample.modified(env.DB)

###########################################################  
###########################################################  1.4. Data partition 

### Data partition function to Training & Validation & Test datasets  

data.partition <- function(DB.o, test.p = 0.2, train.vp = 0.2, DPP = 0.20, r.n = 1){
  # Dataset division
  DB.P <- DB.o[DB.o$OC == 1, ]
  DB.A <- DB.o[DB.o$OC == 0, ]
  
  ##### Test dataset
  set.seed(6808+r.n)
  ind.test.P <- sort(sample(1:nrow(DB.P), size = test.p*nrow(DB.P)))
  ind.test.A <- sort(sample(1:nrow(DB.A), size = test.p*nrow(DB.A)))
  
  DB.test <- rbind(DB.P[ind.test.P, ], DB.A[ind.test.A, ]) 
  
  # Extra dataset (Training & validation set)
  DB.P2 <- DB.P[-ind.test.P, ]
  DB.A2 <- DB.A[-ind.test.A, ]
  
  ##### Training dataset composition by DPP
  # DPP = P/(P+A) 
  size.P <- round( nrow(DB.A2)*DPP/(1-DPP), 0 ) # Number of P relative to the number of A samples
  
  # Verifying size.P 
  print(paste0("P data ratio: ", round(size.P/(size.P+nrow(DB.A2)), 5))) # == DPP
  
  ##### Training & validation data partitioning by DPP
  set.seed(910520+r.n)
  ind.DPP <- sort(sample(1:nrow(DB.P2), size= size.P))
  DB.P3 <- DB.P2[ind.DPP, ]
  
  set.seed(6808+910520+r.n)
  ind.valid.P <- sample(1:nrow(DB.P3), size= train.vp*nrow(DB.P3))
  ind.valid.A <- sample(1:nrow(DB.A2), size= train.vp*nrow(DB.A2))
  
  DB.valid <- rbind(DB.P3[ind.valid.P, ], DB.A2[ind.valid.A, ])
  DB.train <- rbind(DB.P3[-ind.valid.P, ], DB.A2[-ind.valid.A, ]) 
  
  list.DP <- list("test" = DB.test,
                  "valid" = DB.valid,
                  "train" = DB.train,
                  "original" = DB.o
  )
  return(list.DP)
}

# Input data description
# (1) DB.o: Target DB (format: data.frame)
# (2) test.p: Test set ratio: (e.g.) 20% with the remaining 80% used as the extra dataset
# (3) train.vp: Validation set ratio; Extra dataset = Validation dataset + training dataset
# (4) DPP: Desired presence proportion applied to training data only (e.g., 0.20 == 20%) 
# (5) r.n: Initial seed for random number generation

# Output description
# A list labeled with 'original', 'train', 'valid', and 'test' (format: list)
# Each item is structured as a data frame
# -'original': Same as the target DB
# -'train': DB for model training
# -'valid': DB for model validation
# -'test': DB for model testing

### Sample code 
DB.o <- read.csv("Sample data.C1_D2.sp1_Preprocessed data.csv", fileEncoding = "EUC-KR")
head(DB.o)

DB.o <- DB.o[-1]
data.partition(DB.o, DPP = 0.20, r.n = 1) 

# Repetition (e.g., k=2)
DB.list.k1 <- data.partition(DB.o, DPP = 0.20, r.n = 1) 
DB.list.k2 <- data.partition(DB.o, DPP = 0.20, r.n = 5)

save(DB.list.k1, DB.list.k2, 
     file = "Sample data.C2_sp1_p0.20-after partition.RData")


