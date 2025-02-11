# FRQN simulation 
library(dplyr)
library(e1071)
library(caTools)
library(caret)
library(scutr)
library(janitor)
library(broom)
library(xgboost)
library(vioplot)

#### Hello I've changed

constrained_fy_shuffle <- function(vec, m) { 
  n <- length(vec)
  shuffled <- sample(vec) # First shuffle without constraints
  
  for (i in 1:n) {
    if (abs(shuffled[i] - vec[i]) > m) {
      for (j in 1:n) {
        if (i != j && abs(shuffled[j] - vec[i]) <= m && abs(shuffled[i] - vec[j]) <= m) {
          # Swap values
          temp <- shuffled[i]
          shuffled[i] <- shuffled[j]
          shuffled[j] <- temp
          break
        }
      }
    }
  }
  
  return(shuffled)
}

############# Loading microarray data ######

gse4475 <- read.csv("validation_data/final_gse4475.csv")
rownames(gse4475) <- gse4475[,1]
gse4475 <- gse4475[,-1]
gse4475$dataset_id <- "gse4475"
gse4475$panel_id <- "gpl96"

gse7788 <- read.csv("validation_data/final_gse7788.csv")
rownames(gse7788) <- gse7788[,1]
gse7788 <- gse7788[,-1]
gse7788 <- gse7788[,-1]
gse7788$dataset_id <- "gse7788"
gse7788$panel_id <- "gpl570"

gse10138 <- read.csv("validation_data/final_gse10138.csv")
rownames(gse10138) <- gse10138[,1]
gse10138 <- gse10138[,-1]
gse10138 <- gse10138[,-1]
gse10138$dataset_id <- "gse10138"
gse10138$panel_id <- "gpl570"

#gse12453 <- read.csv("validation_data/final_gse12453.csv")
#rownames(gse12453) <- gse12453[,1]
#gse12453 <- gse12453[,-1]
#gse12453 <- gse12453[,-1]
#gse12453$dataset_id <- "gse12453"
#gse12453$panel_id <- "gpl570"

gse13996 <- read.csv("validation_data/final_gse13996.csv")
rownames(gse13996) <- gse13996[,1]
gse13996 <- gse13996[,-1]
gse13996$dataset_id <- "gse13996"
gse13996$panel_id <- "gpl571"

gse16746 <- read.csv("validation_data/final_gse16746.csv")
rownames(gse16746) <- gse16746[,1]
gse16746 <- gse16746[,-1]
gse16746$dataset_id <- "gse16746"
gse16746$panel_id <- "gpl96"

gse29605 <- read.csv("validation_data/final_gse29605.csv")
rownames(gse29605) <- gse29605[,1]
gse29605 <- gse29605[,-1]
gse29605 <- gse29605[,-1]
gse29605$dataset_id <- "gse29605"
gse29605$panel_id <- "gpl570"

gse35082 <- read.csv("validation_data/final_gse35082.csv")
rownames(gse35082) <- gse35082[,1]
gse35082 <- gse35082[,-1]
gse35082 <- gse35082[,-1]
gse35082$dataset_id <- "gse35082"
gse35082$panel_id <- "gpl570"
gse35082$subtype <- "MZL"

gse53820 <- read.csv("validation_data/final_gse53820.csv")
rownames(gse53820) <- gse53820[,1]
gse53820 <- gse53820[,-1]
gse53820$dataset_id <- "gse53820"
gse53820$panel_id <- "gpl570"

gse87371 <- read.csv("validation_data/final_gse87371.csv")
rownames(gse87371) <- gse87371[,1]
gse87371 <- gse87371[,-1]
gse87371 <- gse87371[,-1]
gse87371$dataset_id <- "gse87371"
gse87371$panel_id <- "gpl570"

gse93291 <- read.csv("validation_data/final_gse93291.csv")
rownames(gse93291) <- gse93291[,1]
gse93291 <- gse93291[,-1]
gse93291$dataset_id <- "gse93291"
gse93291$panel_id <- "gpl570"


gse4475_shared <- gse4475[,intersect(colnames(gse93291), colnames(gse4475))]
gse7788_shared <- gse7788[,intersect(colnames(gse93291), colnames(gse4475))]
gse10138_shared <- gse10138[,intersect(colnames(gse93291), colnames(gse4475))]
#gse12453_shared <- gse12453[,intersect(colnames(gse93291), colnames(gse4475))]
gse13996_shared <- gse13996[,intersect(colnames(gse93291), colnames(gse4475))]
gse16746_shared <- gse16746[,intersect(colnames(gse93291), colnames(gse4475))]
gse29605_shared <- gse29605[,intersect(colnames(gse93291), colnames(gse4475))]
gse35082_shared <- gse35082[,intersect(colnames(gse93291), colnames(gse4475))]
gse53820_shared <- gse53820[,intersect(colnames(gse93291), colnames(gse4475))]
gse87371_shared <- gse87371[,intersect(colnames(gse93291), colnames(gse4475))]
gse93291_shared <- gse93291[,intersect(colnames(gse93291), colnames(gse4475))]

microarray <- do.call("rbind", list(gse4475_shared, gse7788_shared, gse10138_shared, 
                                    gse13996_shared, gse16746_shared,
                                    gse29605_shared, gse35082_shared, gse53820_shared,
                                    gse87371_shared, gse93291_shared))


microarray$subtype <- as.factor(microarray$subtype)
levels(microarray$subtype) <- list(
  reactive = "Adenite",                                                  
  burkitt = "Burkitt",                         
  HL = "CHL",                       
  CLL =  "CLL",
  LBCL =  "DLBCL",                                      
  remove = "DLBCL 10% and FL 90%",                                      
  remove = "DLBCL 30% and FL 70%",                             
  remove = "DLBCL 50% and FL 50%",                                  
  FL = "FL",                           
  FL = "follicular lymphoma",                                 
  HL = "Hodgkin's Lymphoma",
  HL = "lymph node affected by NLPHL",
  remove = "Lymphoma Cell line",           
  MCL = "Mantle Cell Lymphoma (MCL)",   
  MZL = "MZL",
  remove = "Morphologic.diagnosis : aggressive B-NHL.unclassifiable",                                        
  burkitt = "Morphologic.diagnosis : atypical BL", 
  burkitt = "Morphologic.diagnosis : BL",
  remove = "Morphologic.diagnosis : BL-leukemia",
  LBCL = "Morphologic.diagnosis : DLBCL",
  HL = "NLPHL",
  LBCL = "PMBL",
  remove = "T-Cell/Histiocyte-Rich Large B-Cell Lymphomas"
)

microarray <- microarray[which(microarray$subtype!="remove"),]
microarray$subtype <- droplevels(microarray$subtype)

table(microarray$subtype)
table(microarray$dataset_id)
table(microarray$panel_id)

microarray_filt <- microarray[,colSums(is.na(microarray))==0]

############# Loading PanB data ####


setwd("R:/HMRN/Analysis/Projects/Lymphoid Gene Expression/Analysis/Data")

genes <- read.csv("gene_exp_processed_wo_validation_290124.csv")
ann <- read.csv("new_annotations_220124.csv")

colnames(genes)[1] <- "HMDS.number"

ann_sub <- ann[,c("HMDS.number", "Diagnosis")]
ann_sub$HMDS.ref <- gsub("/", "_", ann_sub$HMDS.number)
df <- merge(ann_sub[,-3], genes, "HMDS.number")

df_broad <- df

df_broad$Diagnosis <- as.factor(df_broad$Diagnosis)
levels(df_broad$Diagnosis) <- list(
  edge = "?",                                                  
  edge = "?11q (gain of 11q no loss)",                         
  edge = "?11q (gain of chr11 no loss)",                       
  edge =  "?Aggressive B cell Lymphoma.",
  edge =  "?Burkitt 11q?",                                      
  edge = "?Burkitt 11Q?",                                      
  edge = "?Burkitt, atypical MYC",                             
  edge = "?CCND1 neg mantle",                                  
  edge = "?CD5+ FL, ?CCND1- mantle",                           
  edge = "?DHL, atypical MYC",                                 
  edge = "?DLBCL (Maple early progressor study) Atypical BCL2",
  edge = "?DLBCL (Maple early progressor study), atypical DHL",
  edge = "?DLBCL (Maple early progressor study)DHL",           
  edge = "?DLBCL, atypical MYC",                               
  edge = "?DLBCL, DHL",                                        
  edge = "?DLBCL, leg type, atypical MYC",                     
  edge = "?High grade B-cell lymphoma-BCL2 negative",          
  edge = "?Monomorphic PTLD (B&T/NK types)",    
  edge = "?Granulomatous inflammation (MaPLe)",
  edge = "?Suspected 11Q",
  LBCL = "11q",                                                
  burkitt = "Burkitt",                                            
  burkitt = "Burkitt. MYC gene rearrangement. ",                  
  MCL = "CCND1 neg MCL", 
  HL = "CHL",
  CLL = "CLL",                                                
  LBCL = "DHL",                                                
  LBCL = "DHL ",                                               
  LBCL = "DLBCL",                                              
  LBCL =  "DLBCL ",                                             
  "border/transformed DLBCL/FL" = "DLBCL underlying FL",                                
  "border/transformed DLBCL/FL" = "evidence of DLBCL and FL",                           
  FL = "FL",                                                 
  FL = "FL 3B",                                              
  FL = "FL BCL2-neg",                                        
  FL = "FLCL",                                               
  FL = "Follicular lymphoma grade 3",
  MCL = "MCL",                                                
  MCL =  "MCL blastic",                                        
  LBCL = "MHG",                                                
  LBCL = "MHG ",                                               
  LBCL = "MHG/DHL",   
  HL = "NLPHL",
  LBCL = "MYC-R, known partner",                               
  MZL = "MZL",                                                
  LBCL = "Plasmablastic",                                      
  LBCL = "PMBL",      
  plasmacytoma = "Plasmacytoma bone",
  plasmacytoma = "Plasmacytoma Bone",
  plasmacytoma = "Plasmacytoma bone ",
  plasmacytoma = "soft tissue plasmacytoma",
  plasmacytoma = "Soft tissue plasmacytoma",
  plasmacytoma = "Soft tissue plasmacytoma ",
  reactive = "Reactive"
)




df_broad <- df_broad %>% 
  filter(Diagnosis != "edge") %>%
  filter(Diagnosis != "border/transformed DLBCL/FL")

df_broad$Diagnosis <- droplevels(df_broad$Diagnosis)

df_broad <- df_broad[-c(829,920),]
rownames(df_broad) <- df_broad$HMDS.number
df_broad <- df_broad[,-1]

df_broad_reduced <- df_broad[,intersect(colnames(microarray_filt), colnames(df_broad))]
df_broad_reduced$subtype <- df_broad$Diagnosis

# 

############# Identifying reference distance distributions #####


rank_lbcl <- microarray_filt[which(microarray_filt$subtype=="LBCL"),-c(254,255,256)]
rank_lbcl <- data.frame(t(apply(rank_lbcl,1,rank,ties.method="min")))

df_broad_reduced_lbcl <- df_broad_reduced[which(df_broad_reduced$subtype=="LBCL"),-254]
df_broad_reduced_lbcl <- data.frame(t(apply(df_broad_reduced_lbcl,1,rank,ties.method="min")))
df_broad_reduced_lbcl <- unlist((lapply(df_broad_reduced_lbcl, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_lbcl - rep(df_broad_reduced_lbcl, each=nrow(rank_lbcl)))
test_d1_lbcl <- apply(test, 1, max)
test_d2_lbcl <- rowSums(test)
test_d3_lbcl <- rowSums(test != 0)

rank_burk <- microarray_filt[which(microarray_filt$subtype=="burkitt"),-c(254,255,256)]
rank_burk <- data.frame(t(apply(rank_burk,1,rank,ties.method="min")))

df_broad_reduced_burk <- df_broad_reduced[which(df_broad_reduced$subtype=="burkitt"),-254]
df_broad_reduced_burk <- data.frame(t(apply(df_broad_reduced_burk,1,rank,ties.method="min")))
df_broad_reduced_burk <- unlist((lapply(df_broad_reduced_burk, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_burk - rep(df_broad_reduced_burk, each=nrow(rank_burk)))
test_d1_burk <- apply(test, 1, max)
test_d2_burk <- rowSums(test)
test_d3_burk <- rowSums(test != 0)


rank_hl <- microarray_filt[which(microarray_filt$subtype=="HL"),-c(254,255,256)]
rank_hl <- data.frame(t(apply(rank_hl,1,rank,ties.method="min")))

df_broad_reduced_hl <- df_broad_reduced[which(df_broad_reduced$subtype=="HL"),-254]
df_broad_reduced_hl <- data.frame(t(apply(df_broad_reduced_hl,1,rank,ties.method="min")))
df_broad_reduced_hl <- unlist((lapply(df_broad_reduced_hl, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_hl - rep(df_broad_reduced_hl, each=nrow(rank_hl)))
test_d1_hl <- apply(test, 1, max)
test_d2_hl <- rowSums(test)
test_d3_hl <- rowSums(test != 0)


rank_cll <- microarray_filt[which(microarray_filt$subtype=="CLL"),-c(254,255,256)]
rank_cll <- data.frame(t(apply(rank_cll,1,rank,ties.method="min")))

df_broad_reduced_cll <- df_broad_reduced[which(df_broad_reduced$subtype=="CLL"),-254]
df_broad_reduced_cll <- data.frame(t(apply(df_broad_reduced_cll,1,rank,ties.method="min")))
df_broad_reduced_cll <- unlist((lapply(df_broad_reduced_cll, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_cll - rep(df_broad_reduced_cll, each=nrow(rank_cll)))
test_d1_cll <- apply(test, 1, max)
test_d2_cll <- rowSums(test)
test_d3_cll <- rowSums(test != 0)




rank_mcl <- microarray_filt[which(microarray_filt$subtype=="MCL"),-c(254,255,256)]
rank_mcl <- data.frame(t(apply(rank_mcl,1,rank,ties.method="min")))

df_broad_reduced_mcl <- df_broad_reduced[which(df_broad_reduced$subtype=="MCL"),-254]
df_broad_reduced_mcl <- data.frame(t(apply(df_broad_reduced_mcl,1,rank,ties.method="min")))
df_broad_reduced_mcl <- unlist((lapply(df_broad_reduced_mcl, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_mcl - rep(df_broad_reduced_mcl, each=nrow(rank_mcl)))
test_d1_mcl <- apply(test, 1, max)
test_d2_mcl <- rowSums(test)
test_d3_mcl <- rowSums(test != 0)


rank_fl <- microarray_filt[which(microarray_filt$subtype=="FL"),-c(254,255,256)]
rank_fl <- data.frame(t(apply(rank_fl,1,rank,ties.method="min")))

df_broad_reduced_fl <- df_broad_reduced[which(df_broad_reduced$subtype=="FL"),-254]
df_broad_reduced_fl <- data.frame(t(apply(df_broad_reduced_fl,1,rank,ties.method="min")))
df_broad_reduced_fl <- unlist((lapply(df_broad_reduced_fl, function (x) quantile(x,0.5,type=1))))

test <- abs(rank_fl - rep(df_broad_reduced_fl, each=nrow(rank_fl)))
test_d1_fl <- apply(test, 1, max)
test_d2_fl <- rowSums(test)
test_d3_fl <- rowSums(test != 0)


test_d1_lbcl <- data.frame(value=test_d1_lbcl, var="LBCL")
test_d1_cll <- data.frame(value=test_d1_cll, var="CLL")
test_d1_burk <- data.frame(value=test_d1_burk, var="burk")
test_d1_fl <- data.frame(value=test_d1_fl, var="FL")
test_d1_mcl <- data.frame(value=test_d1_mcl, var="MCL")
test_d1_hl <- data.frame(value=test_d1_hl, var="HL")

d1_combo <- rbind(test_d1_lbcl, test_d1_cll, test_d1_burk, test_d1_fl, test_d1_mcl, test_d1_hl)

vioplot(value ~ var, data=d1_combo, main="Comparison of observed D1 by subtype", xlab="Subtype", ylab="Max distance")


test_d2_lbcl <- data.frame(value=test_d2_lbcl, var="LBCL")
test_d2_cll <- data.frame(value=test_d2_cll, var="CLL")
test_d2_burk <- data.frame(value=test_d2_burk, var="burk")
test_d2_fl <- data.frame(value=test_d2_fl, var="FL")
test_d2_mcl <- data.frame(value=test_d2_mcl, var="MCL")
test_d2_hl <- data.frame(value=test_d2_hl, var="HL")

d2_combo <- rbind(test_d2_lbcl, test_d2_cll, test_d2_burk, test_d2_fl, test_d2_mcl, test_d2_hl)

vioplot(value ~ var, data=d2_combo, main="Comparison of observed D2 by subtype", xlab="Subtype", ylab="Sum of distances")

test_d3_lbcl <- data.frame(value=test_d3_lbcl, var="LBCL")
test_d3_cll <- data.frame(value=test_d3_cll, var="CLL")
test_d3_burk <- data.frame(value=test_d3_burk, var="burk")
test_d3_fl <- data.frame(value=test_d3_fl, var="FL")
test_d3_mcl <- data.frame(value=test_d3_mcl, var="MCL")
test_d3_hl <- data.frame(value=test_d3_hl, var="HL")

d3_combo <- rbind(test_d3_lbcl, test_d3_cll, test_d3_burk, test_d3_fl, test_d3_mcl, test_d3_hl)

vioplot(value ~ var, data=d3_combo, main="Comparison of observed D3 by subtype", xlab="Subtype", ylab="Number of positions changed")


repo <- repository("~/frqn_simulation")
pull(repo)
write_vc(d1_combo, file = "rel_path/filename", root = repo, stage = TRUE)
commit(repo, "My message")
push(repo)


############# Identifying optimal constraint parameter ####

M <- 253
N <- 1000
x <- df_broad_reduced_burk
d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 150)
  t <- x-y
  t_abs <- abs(t)
  d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d1_10_burk <- data.frame(value=d1_vals, var=10)
d2_10_burk <- data.frame(value=d2_vals, var=10)
d3_10_burk <- data.frame(value=d3_vals, var=10)

d1_25_burk <- data.frame(value=d1_vals, var=25)
d2_25_burk <- data.frame(value=d2_vals, var=25)
d3_25_burk <- data.frame(value=d3_vals, var=25)

d1_50_burk <- data.frame(value=d1_vals, var=50)
d2_50_burk <- data.frame(value=d2_vals, var=50)
d3_50_burk <- data.frame(value=d3_vals, var=50)

d1_75_burk <- data.frame(value=d1_vals, var=75)
d2_75_burk <- data.frame(value=d2_vals, var=75)
d3_75_burk <- data.frame(value=d3_vals, var=75)

d1_100_burk <- data.frame(value=d1_vals, var=100)
d2_100_burk <- data.frame(value=d2_vals, var=100)
d3_100_burk <- data.frame(value=d3_vals, var=100)

d1_150_burk <- data.frame(value=d1_vals, var=150)
d2_150_burk <- data.frame(value=d2_vals, var=150)
d3_150_burk <- data.frame(value=d3_vals, var=150)

burkitt_d1 <- rbind(test_d1_burk, d1_10_burk, d1_25_burk, d1_50_burk, d1_75_burk, d1_100_burk, d1_150_burk)
burkitt_d2 <- rbind(test_d2_burk, d2_10_burk, d2_25_burk, d2_50_burk, d2_75_burk, d2_100_burk, d2_150_burk)
burkitt_d3 <- rbind(test_d3_burk, d3_10_burk, d3_25_burk, d3_50_burk, d3_75_burk, d3_100_burk, d3_150_burk)

par(mfrow=c(2,1),mai = c(1, 1, 0.6,0.6))
vioplot(value ~ var, data=burkitt_d2, main="Identifying optimal constraint parameter - Burkitt", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=burkitt_d3, xlab="m (constraint parameter)", ylab="D3 - # Positions changed")


M <- 253
N <- 1000
x <- df_broad_reduced_cll
#d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 150)
  t <- x-y
  t_abs <- abs(t)
  #d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d2_10_cll <- data.frame(value=d2_vals, var=10)
d3_10_cll <- data.frame(value=d3_vals, var=10)

d2_25_cll <- data.frame(value=d2_vals, var=25)
d3_25_cll <- data.frame(value=d3_vals, var=25)

d2_50_cll <- data.frame(value=d2_vals, var=50)
d3_50_cll <- data.frame(value=d3_vals, var=50)

d2_75_cll <- data.frame(value=d2_vals, var=75)
d3_75_cll <- data.frame(value=d3_vals, var=75)

d2_100_cll <- data.frame(value=d2_vals, var=100)
d3_100_cll <- data.frame(value=d3_vals, var=100)

d2_150_cll <- data.frame(value=d2_vals, var=150)
d3_150_cll <- data.frame(value=d3_vals, var=150)

cll_d2 <- rbind(test_d2_cll, d2_10_cll, d2_25_cll, d2_50_cll, d2_75_cll, d2_100_cll, d2_150_cll)
cll_d3 <- rbind(test_d3_cll, d3_10_cll, d3_25_cll, d3_50_cll, d3_75_cll, d3_100_cll, d3_150_cll)


vioplot(value ~ var, data=cll_d2,main="Identifying optimal constraint parameter - CLL", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=cll_d3, xlab="K (constraint parameter)", ylab="D3 - # Positions changed")



M <- 253
N <- 1000
x <- df_broad_reduced_fl
#d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 150)
  t <- x-y
  t_abs <- abs(t)
  #d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d2_10_fl <- data.frame(value=d2_vals, var=10)
d3_10_fl <- data.frame(value=d3_vals, var=10)

d2_25_fl <- data.frame(value=d2_vals, var=25)
d3_25_fl <- data.frame(value=d3_vals, var=25)

d2_50_fl <- data.frame(value=d2_vals, var=50)
d3_50_fl <- data.frame(value=d3_vals, var=50)

d2_75_fl <- data.frame(value=d2_vals, var=75)
d3_75_fl <- data.frame(value=d3_vals, var=75)

d2_100_fl <- data.frame(value=d2_vals, var=100)
d3_100_fl <- data.frame(value=d3_vals, var=100)

d2_150_fl <- data.frame(value=d2_vals, var=150)
d3_150_fl <- data.frame(value=d3_vals, var=150)

fl_d2 <- rbind(test_d2_fl, d2_10_fl, d2_25_fl, d2_50_fl, d2_75_fl, d2_100_fl, d2_150_fl)
fl_d3 <- rbind(test_d3_fl, d3_10_fl, d3_25_fl, d3_50_fl, d3_75_fl, d3_100_fl, d3_150_fl)

vioplot(value ~ var, data=fl_d2, main="Identifying optimal constraint parameter - FL", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=fl_d3, xlab="K (constraint parameter)", ylab="D3 - # Positions changed")



M <- 253
N <- 1000
x <- df_broad_reduced_hl
#d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 100)
  t <- x-y
  t_abs <- abs(t)
  #d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d2_10_hl <- data.frame(value=d2_vals, var=10)
d3_10_hl <- data.frame(value=d3_vals, var=10)

d2_25_hl <- data.frame(value=d2_vals, var=25)
d3_25_hl <- data.frame(value=d3_vals, var=25)

d2_50_hl <- data.frame(value=d2_vals, var=50)
d3_50_hl <- data.frame(value=d3_vals, var=50)

d2_75_hl <- data.frame(value=d2_vals, var=75)
d3_75_hl <- data.frame(value=d3_vals, var=75)

d2_100_hl <- data.frame(value=d2_vals, var=100)
d3_100_hl <- data.frame(value=d3_vals, var=100)

d2_150_hl <- data.frame(value=d2_vals, var=150)
d3_150_hl <- data.frame(value=d3_vals, var=150)

hl_d2 <- rbind(test_d2_hl, d2_10_hl, d2_25_hl, d2_50_hl, d2_75_hl, d2_100_hl, d2_150_hl)
hl_d3 <- rbind(test_d3_hl, d3_10_hl, d3_25_hl, d3_50_hl, d3_75_hl, d3_100_hl, d3_150_hl)

vioplot(value ~ var, data=hl_d2, main="Identifying optimal constraint parameter - HL", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=hl_d3, xlab="K (constraint parameter)", ylab="D3 - # Positions changed")




M <- 253
N <- 1000
x <- df_broad_reduced_lbcl
#d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 150)
  t <- x-y
  t_abs <- abs(t)
  #d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d2_10_lbcl <- data.frame(value=d2_vals, var=10)
d3_10_lbcl <- data.frame(value=d3_vals, var=10)

d2_25_lbcl <- data.frame(value=d2_vals, var=25)
d3_25_lbcl <- data.frame(value=d3_vals, var=25)

d2_50_lbcl <- data.frame(value=d2_vals, var=50)
d3_50_lbcl <- data.frame(value=d3_vals, var=50)

d2_75_lbcl <- data.frame(value=d2_vals, var=75)
d3_75_lbcl <- data.frame(value=d3_vals, var=75)

d2_100_lbcl <- data.frame(value=d2_vals, var=100)
d3_100_lbcl <- data.frame(value=d3_vals, var=100)

d2_150_lbcl <- data.frame(value=d2_vals, var=150)
d3_150_lbcl <- data.frame(value=d3_vals, var=150)

lbcl_d2 <- rbind(test_d2_lbcl, d2_10_lbcl, d2_25_lbcl, d2_50_lbcl, d2_75_lbcl, d2_100_lbcl, d2_150_lbcl)
lbcl_d3 <- rbind(test_d3_lbcl, d3_10_lbcl, d3_25_lbcl, d3_50_lbcl, d3_75_lbcl, d3_100_lbcl, d3_150_lbcl)

vioplot(value ~ var, data=lbcl_d2, main="Identifying optimal constraint parameter - LBCL", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=lbcl_d3, xlab="K (constraint parameter)", ylab="D3 - # Positions changed")



M <- 253
N <- 1000
x <- df_broad_reduced_mcl
#d1_vals <- rep(0, N)
d2_vals <- rep(0, N)
d3_vals <- rep(0, N)
#distance <- vector("list", length = N)

for (i in 1:N){
  set.seed(i)
  y <- constrained_fy_shuffle(x, 150)
  t <- x-y
  t_abs <- abs(t)
  #d1_vals[i] <- max(t_abs)
  d2_vals[i] <- sum(t_abs)
  d3_vals[i] <- sum(x!=y)
  #distance[[i]] <- t
}

d2_10_mcl <- data.frame(value=d2_vals, var=10)
d3_10_mcl <- data.frame(value=d3_vals, var=10)

d2_25_mcl <- data.frame(value=d2_vals, var=25)
d3_25_mcl <- data.frame(value=d3_vals, var=25)

d2_50_mcl <- data.frame(value=d2_vals, var=50)
d3_50_mcl <- data.frame(value=d3_vals, var=50)

d2_75_mcl <- data.frame(value=d2_vals, var=75)
d3_75_mcl <- data.frame(value=d3_vals, var=75)

d2_100_mcl <- data.frame(value=d2_vals, var=100)
d3_100_mcl <- data.frame(value=d3_vals, var=100)

d2_150_mcl <- data.frame(value=d2_vals, var=150)
d3_150_mcl <- data.frame(value=d3_vals, var=150)

mcl_d2 <- rbind(test_d2_mcl, d2_10_mcl, d2_25_mcl, d2_50_mcl, d2_75_mcl, d2_100_mcl, d2_150_mcl)
mcl_d3 <- rbind(test_d3_mcl, d3_10_mcl, d3_25_mcl, d3_50_mcl, d3_75_mcl, d3_100_mcl, d3_150_mcl)

vioplot(value ~ var, data=mcl_d2, main="Identifying optimal constraint parameter - MCL", xlab="", ylab="D2 - Sum of distances")
vioplot(value ~ var, data=mcl_d3, xlab="K (constraint parameter)", ylab="D3 - # Positions changed")

