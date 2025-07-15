
### (Common function) Data arrangement function for statistical test 
# CD: Combined data / OD: Original data / PD: Partitioned data
pre.CD.OD.f <- function(DB.list, over.m = "SMOTE"){
  DB.original <- DB.list$original[DB.list$original$OC == 1, ] # only for minority class (P)
  DB.train.syn <- rbind(DB.list$train, DB.list[[paste0("syn.", over.m)]]) # combined training and synthetic data
  DB.train.syn[DB.train.syn$OC == 1, ] # # only for minority class (P)
  DB.original$group <- "OD" 
  DB.train.syn$group <- "CD"
  t.DB <- rbind(DB.original, DB.train.syn)
  return(t.DB)
}

pre.CD.PD.f <- function(DB.list, over.m = "SMOTE"){
  DB.train <- DB.list$train[DB.list$train$OC == 1, ] # only for minority class (P)
  DB.train.syn <- rbind(DB.list$train, DB.list$syn) # combined training and synthetic data
  DB.train.syn[DB.train.syn$OC == 1, ] # # only for minority class (P)
  
  DB.train$group <- "PD" 
  DB.train.syn$group <- "CD"
  t.DB <- rbind(DB.train, DB.train.syn)
  return(t.DB)
}

###########################################################  
###########################################################  3.1. MRPP tests with restricted permutation (multivariate test)

### Multi-Response Permutation Procedure (MRPP) tests with restricted permutation function
library(vegan)
library(permute)

MRPP.restricted <- function(DB, groups, blocks, nperm = 1000, r.n = 2025){
  ctrl <- how(blocks = blocks, nperm = nperm)
  set.seed(r.n) 
  DB.mrpp <- mrpp(DB, grouping = groups, permutations = ctrl)
  
  # See the mrpp function documentation in the vegan package for details on each component
  tcd <- data.frame(t(DB.mrpp$classdelta)) 
  names(tcd) <- paste0("delta.", names(tcd))
  
  tn <- data.frame(t(DB.mrpp$n)) 
  names(tn) <- paste0("n.", names(tn))
  
  tdf <- data.frame("p.value" = DB.mrpp$Pvalue, 
                    "expected.delta" = DB.mrpp$E.delta, 
                    "overall.delta" = DB.mrpp$delta, 
                    "mrpp.A" = DB.mrpp$A 
  )
  result <- cbind(tdf, tcd, tn)
  return(result)
}

# Input data description:
# (1) DB: Comparison target database. Must consist only of numeric values to allow distance calculation (format: data.frame)
# (2) groups: Comparison target groups (format: vector)
# (3) blocks: Homogeneous clusters (=blocks) (format: vector)
# (4) nperm: Number of permutation
# (5) r.n: Initial seed for random number generation

# Output description (See the mrpp function documentation in the vegan package for details on each component)
# - p.value: Significance of the test, under the null hypothesis of no group structure
# - expected.delta: expected delta, under the null hypothesis of no group structure.
# - overall.delta: The overall weighted mean of group mean distances
# - mrpp.A: A chance-corrected estimate of the proportion of the distances explained by group identity
# - delta.[group]: Mean dissimilarities within classes (i.e., groups)
# - n.[group] Number of observations in each class.


### Sample code 
load("Sample data.C3-4_sp1_p0.20_after oversampling.RData") # Including oversampled Data 

tDB1 <- pre.CD.OD.f(DB.list.k1); tDB2 <- pre.CD.OD.f(DB.list.k2) # comparison between CD and OD
# tDB1 <- pre.CD.PD.f(DB.list.k1); tDB2 <- pre.CD.PD.f(DB.list.k2) # comparison between CD and PD
tDB1$k <-1; tDB2$k <-2

total.DB <- rbind(tDB1, tDB2)
head(total.DB)

# See the mrpp function documentation in the vegan package for details on each component
MRPP.restricted(DB = total.DB[2:9], groups = total.DB$group, blocks = total.DB$k, r.n=2025)


###########################################################  
###########################################################  3.2. Restricted permutation test (univariate test)

### Restricted permutation test function (Calculates for all variables in the DB)
library(permute)

Perm.restricted <- function(DB, groups, blocks, nperm = 1000, r.n = 2025){
  
  result <- c()
  for(j in seq(names(DB))){
  values <- DB[,names(DB)[j]]
  
  ## mean difference function
  meanDif <- function(i, x, grp) {
    grp <- grp[i]
    mean(x[grp == unique(groups)[1]]) - mean(x[grp == unique(groups)[2]])
  }
  
  # permutation control option
  n_total <- length(values)
  ctrl <- how(
    within = Within(type = "free"), 
    blocks = blocks
  )
  
  # permutation index 
  set.seed(r.n) 
  perm_matrix <- shuffleSet(n_total, nset = nperm, control = ctrl)
  
  ## iterate over the set of permutations applying meanDif
  perm_diffs  <- apply(perm_matrix, 1, meanDif, x = values, grp = groups)
  ## add on the observed mean difference
  obs_diff  <- meanDif(seq_len(n_total), x = values, grp = groups)
  
  ## compute & return the p-value
  pval <- (sum(perm_diffs >= obs_diff) + 1) / (length(perm_diffs) + 1)
  
  tdf <- data.frame("variable" = names(DB)[j], 
                    "p.value" = pval)
  result <- rbind(result, tdf)
  }# for j
  
  return(result)
}

# Input data description:
# (1) DB: Comparison target database. Must consist only of numeric values (format: data.frame)
# (2) groups: Comparison target groups (format: vector)
# (3) blocks: Homogeneous clusters (=blocks) (format: vector)
# (4) nperm: Number of permutation
# (5) r.n: Initial seed for random number generation

# Output description (See the mrpp function documentation in the vegan package for details on each component)
# - variable: Target variable for comparison (variable names ian the DB)
# - p.value: Significance of the test

### Sample code 
load("Sample data.C3-4_sp1_p0.20_after oversampling.RData") # Including oversampled Data 

tDB1 <- pre.CD.OD.f(DB.list.k1); tDB2 <- pre.CD.OD.f(DB.list.k2) # comparison between CD and OD
# tDB1 <- pre.CD.PD.f(DB.list.k1); tDB2 <- pre.CD.PD.f(DB.list.k2) # comparison between CD and PD
tDB1$k <-1; tDB2$k <-2

total.DB <- rbind(tDB1, tDB2)

Perm.restricted(DB = total.DB[2:9], groups = total.DB$group, blocks = total.DB$k, r.n=2025) # for all variables in the DB 


