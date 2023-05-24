## code to prepare `DATASET` dataset goes here
# Generate data toy example
rm(list=ls())
library(data.table)
library(truncnorm)
library(nlme)
library(devtools)
library(mice)

#0. Load original dataset ----

# data <- as.data.table(read.csv('../ObesityDataSet_raw_and_data_sinthetic.csv', stringsAsFactors=FALSE, fileEncoding="latin1"))
# save(data, file = "obesity_raw.rda")
load("../obesity_raw.rda")

#1. Modify dataset in order to get a variable weight with MNAR missing mechanism ----
# Generate countries where survey was taken
data <- data[,c("Weight","Gender","Age","Height","FAVC")] # Select some of the variables of the original dataset
data[,Cluster:=c(rep(1,500),rep(2,400),rep(3,450),rep(4,400),rep(5,361))] # in the imputation method is necessary cluster is included as numeric
N<-length(unique(data$Cluster))
#  Generation of the exclusion restriction (ERV): The time variable describes the time it takes for a person to answer the weight question.
set.seed(12345)
data[,Time:=rtruncnorm(n=nrow(data), a = 1, b = 10, mean = 5, sd = 3)]
data[,Gender:=(ifelse(Gender=="Female",1,0))]
data[,FAVC:=(ifelse(FAVC=="yes",1,0))]
data[,Heights:=scale(Height)]
data[,Ages:=scale(Age)]
data[,Times:=scale(Time)]

MAR_model <- with(data, lme4::lmer(Weight~Gender+Ages+Heights+FAVC+(1+Heights+Gender+Ages|Cluster)))
summary(MAR_model)
# Generate outcome variable
data[,XOBO:=predict(MAR_model)+15]

# Generate random slope for Gender Age on selection equation
alpha_gender <- rnorm(n = N, mean=0, sd = 0.01)
alpha_age <- rnorm(n = N, mean=0, sd = 0.01)
nobs<-c(500,400,450,400,361)
data[,Genderi:=unlist(mapply(rep,alpha_gender,each=nobs))]
data[,Agei:=unlist(mapply(rep,alpha_age,each=nobs))]
data[,XSBS:=1.5-(0.7+Genderi)*Gender-(0.5+Agei)*Ages-1.2*Times] # Simulated selection model

# Simulate bivariate normal correlated errors
v <- VarCorr(MAR_model) #Get varcov matrix
eps <- NULL
rho=-0.8 # We assume that non-responders are more likely to have a high weight.
set.seed(12345)
sigmae <- exp( rnorm( n=N, mean=log(attr(v, "sc")), sd = 0.02))
for( i in 1: N){
  d <- diag(2)
  d[2,2] <- sigmae[i]^2
  d[2,1] <- d[1,2] <- sqrt(d[1,1])*sqrt(d[2,2])*rho
  eps_i <- mvtnorm::rmvnorm(n = nrow(data[Cluster==i]), mean = rep(0,2), sigma = d)
  eps <- rbind(eps,eps_i)
}

# Generate latent outcome and selection variables
data[,Weight.star:=XOBO+eps[,2]]
data[,ry.star:=XSBS+eps[,1]]
data[,ry:=ifelse(ry.star>0,1,0)]

# Generate observed variables
data[,ry:=ifelse(ry.star>0&Cluster!=5,1,0)] #Include systematic missingness on cluster 4
data[,Weight :=ifelse(ry==1,Weight.star,NA)]

#Generate missing in other variables
data_other<-data[,c("Gender","Age","Height","FAVC","Weight.star")]
# ampute the complete data once for every mechanism
ampdata0 <- ampute(data_other, patterns = c(1,1,0,0,1), prop = 0.2, mech = "MAR")$amp
ampdata1 <- ampute(data_other, patterns = c(1,0,0,0,1), prop = 0.1, mech = "MAR")$amp
# create a random allocation vector
# use the prob argument to specify how much of each mechanism should be created
# here, 0.5 of the missingness should be MAR and 0.5 should be MCAR
indices <- sample(x = c(1, 2), size = nrow(ampdata0), replace = TRUE, prob = c(0.7, 0.3))
ampdata <- matrix(NA, nrow = nrow(data_other), ncol = ncol(data_other))
ampdata[indices == 1, ] <- as.matrix(ampdata0[indices == 1, ])
ampdata[indices == 2, ] <- as.matrix(ampdata1[indices == 2, ])
colnames(ampdata)<-colnames(data_other)
ampdata<-data.table(ampdata)
ampdata[,Weight:=data$Weight]
ampdata[,Time:=data$Time]
ampdata[,Cluster:=data$Cluster]
ampdata[,Weight.star:=NULL]
ampdata[,FAVC :=as.factor(FAVC)]
ampdata[,Gender :=as.factor(Gender)]

# Save dataObs
dataObs<-as.data.frame(ampdata)
# save(dataObs,file="dataObs.Rdata")

# save obesity data for heckman package
obesity <- janitor::clean_names(dataObs)
obesity <- cbind(cluster = obesity$cluster, obesity[, 1:6])
# save(obesity, file = "obesity.rda")


usethis::use_data(obesity, overwrite = TRUE)
