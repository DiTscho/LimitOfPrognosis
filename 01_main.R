library(plyr)
library(tibble)
library(breastCancerMAINZ)
##### files #####
#'Overall 3288 patients.'
files = c(
  "../data/GSE11121/data/GSE11121.RData",
  "../data/GSE19615/data/GSE19615.RData",
  "../data/GSE22226/data/GSE22226_gpl1708.RData",
 # "../data/GSE23988/data/GSE23988.RData",
  "../data/GSE26971_aka_E_TABM_158/data/E_TABM_158.RData",
  "../data/GSE96058/data/GSE96058.Rdata",
  "../data/GSE1456/data/GSE1456A.RData",
  "../data/GSE42568/data/GSE42568.RData",
  #"../data/GSE43365/data/GSE43365.RData",
  "../data/GSE45255/data/GSE45255.RData",
  "../data/GSE4922/data/GSE4922A.RData",
  "../data/GSE7390/data/GSE7390.RData",
  "../data/GSE7849/data/GSE7849.RData",
  "../data/GSE9195/data/GSE9195.RData",
  "../data/GSE9893/data/GSE9893.RData",
  "../data/GSE2603/data/GSE2603.RData",
  "../data/GSE_XXXX_aka_NKI/data/NKI.RData",
 #"../data/GSE6532/data/GSE6532A.RData",
 "../data/TCGA/data/tcga.Rdata"
)

###### parameters #####
d_time = 30.4167
node_status = "node_status"
grade = "grade"
tumor_size = "tumor_size"
survival_time = "survival_time" # Distant Metastasis-Free Survival
event = "event"
age = "age"
chemo = "chemo" # any kind of Neoadjuvant treatment
hormonal = "hormonal"
radiation = "radiation"
stage = "stage"
ER = "ER"
PR = "PR"
HER2 = "HER2"

##### helper functions #####
isevent = function(x){ifelse(x==0|is.na(x),F,T)}
isnoevent = function(x){ifelse(x==1|is.na(x),F,T)}

sb = function(df){return(names(df)[1:100])}
se = function(df){
  n = names(df)[(ncol(df)-50):ncol(df)]
  return(n)}

reorder_df = function(i, j, df){
  df1 = df[,1:i-1]
  df2 = df[,i:j-1]
  df3 = df[,j:ncol(df)]
  return(cbind(df3,df1,df2))
}

get_clin = function(i, e, df){return(cbind(df[,1:i], df[,e:ncol(df)]))}
get_expr = function(i, e, df){return(df[,(i+1):(e-1)])}
empty_df = function(df_length = 100){
  na_vec = rep(NA,df_length)
  df = data.frame(node_status = na_vec,
                  grade = na_vec,
                  tumor_size = na_vec,
                  survival_time = na_vec,
                  event = na_vec,
                  age = na_vec,
                  chemo = na_vec,
                  hormonal = na_vec,
                  radiation = na_vec,
                  stage = na_vec,
                  ER = na_vec,
                  PR = na_vec,
                  HER2 = na_vec)
  return(df)
}

##### GSE11121 #####
X1 = get(load(files[1]))

c = 50
names(X1)[1:10]
X1[1:5,1:20]

X1.c = empty_df(df_length = dim(X1)[1])
dim(X1.c)
X1.c$node_status =   as.factor(X1$`node:ch1`)
X1.c$grade =         as.factor(X1$`grade:ch1`)
X1.c$tumor_size =    as.numeric(X1$`size_in_cm:ch1`)
X1.c$survival_time = as.numeric(X1$`t.dmfs:ch1`)
X1.c$event =         as.factor(X1$`e.dmfs:ch1`)
X1.c$stage =         as.factor(NA)
X1.c$chemo = as.factor(0)
X1.c$PR = as.factor(NA)
X1.c$HER2 = as.factor(NA)
X1.c$hormonal = as.factor(0)
X1.c$age = as.numeric(NA)
X1.c$radiation = as.factor(NA)


data("mainz")
X1.cM = pData(mainz)
sum(X1.cM$id == X1$title)==dim(X1)[1] #check order
head(X1.cM)
X1.c$age = X1.cM$age
X1.c$ER = X1.cM$er

X1 = cbind(X1.c, X1)

##### GSE19615 #####
X2 = get(load(files[2]))

c = 100
names(X2)[1:c]
X2[1:5,1:c]

X2.c = empty_df(df_length = dim(X2)[1])

X2.c$node_status = X2$"lymph nodes:ch1"
table(X2.c$node_status)
X2.c$node_status = ifelse(X2.c$node_status=="negative", 0, 1)
X2.c$node_status = as.factor(X2.c$node_status)

X2.c$grade = X2$"grade (modified, bloom, richardson):ch1"
table(X2.c$grade)
class(X2.c$grade)
X2.c$grade[X2.c$grade=="I"] = 1
X2.c$grade[X2.c$grade=="II"] = 2
X2.c$grade[X2.c$grade=="III"] = 3
X2.c$grade = as.factor(X2.c$grade)

X2.c$stage = as.factor(NA)

X2.c$tumor_size = X2$"tumor size (cm):ch1"
X2.c$tumor_size
X2.c$tumor_size = as.numeric(X2.c$tumor_size)

X2.c$event = X2$"distant recur (yn):ch1"
X2.c$event
X2.c$event = as.factor(ifelse(X2.c$event=="Y", 1, 0))

X2.c$survival_time = NA
X2.c$survival_time[X2.c$event==1] = X2$`distant recurrence free survival (mo):ch1`[X2.c$event==1]
X2.c$survival_time[X2.c$event==0] = X2$`time of followup (mo):ch1`[X2.c$event==0]
X2.c$survival_time = as.numeric(X2.c$survival_time)

X2.c$age = X2$"age:ch1"
X2.c$age = as.numeric(X2.c$age)
X2.c$ER = X2$"er:ch1"
X2.c$ER= as.factor(ifelse(X2.c$ER=="neg", 0, 1))
X2.c$PR = X2$"pr:ch1"
X2.c$PR= as.factor(ifelse(X2.c$PR=="neg", 0, 1))
X2.c$HER2 = X2$"her.2:ch1"
X2.c$HER2 = as.factor(ifelse(X2.c$HER2=="neg", 0, 1))
X2.c$chemo = as.factor(ifelse(X2$"chemo class:ch1"=="none", 0, 1))
X2.c$hormonal = as.factor(ifelse(X2$"hormonal rx:ch1"=="none", 0, 1))
X2.c$radiation = as.factor(NA)

X2 = cbind(X2.c, X2)
X2[1:5,1:100]

##### GSE22226 #####

X3A = get(load(files[3])) #"../data/GSE22226/data/GSE22226_gpl1708.RData"

c = 20
names(X3A)[1:c]
X3A[1:5,1:c]

X3A.c = empty_df(df_length = dim(X3A)[1])

X3A.c$node_status = as.factor(X3A.c$node_status)

X3A.c$grade = X3A$"histologic.grade"
table(X3A.c$grade)
class(X3A.c$grade)
X3A.c$grade = as.factor(X3A.c$grade)
X3A.c$stage = as.factor(X3A$clinical.tumor.stage)

X3A.c$tumor_size = X3A$"clinical.tumor.size"
X3A.c$tumor_size
X3A.c$tumor_size = as.numeric(X3A.c$tumor_size)
X3A.c$tumor_size[X3A.c$tumor_size==99.0 | X3A.c$tumor_size<=0.0] = NA # there is evidently a mistake

X3A.c$event = as.factor(X3A$"rfs.e")

X3A.c$survival_time = as.numeric(X3A$rfs.t)/d_time #days 2 months 

X3A.c$age = NA
X3A.c$age = as.numeric(X3A.c$age)
X3A.c$ER= as.factor(X3A$"er")
X3A.c$PR= as.factor(X3A$"pr")
X3A.c$HER2= as.factor(X3A$"her2")
X3A.c$chemo = as.factor(1)
X3A.c$hormonal = as.factor(NA)
X3A.c$radiation = as.factor(NA)

X3A = cbind(X3A.c, X3A)
X3A[1:5,1:100]


X3B = get(load("../data/GSE22226/data/GSE22226_gpl4133.RData"))
c = 20
names(X3B)[1:c]
X3B[1:5,1:c]

X3B.c = empty_df(df_length = dim(X3B)[1])

X3B.c$node_status = as.factor(X3B.c$node_status)

X3B.c$grade = X3B$"histologic.grade"
table(X3B.c$grade)
class(X3B.c$grade)
X3B.c$grade = as.factor(X3B.c$grade)
X3B.c$stage = as.factor(X3B$clinical.tumor.stage)

X3B.c$tumor_size = X3B$"clinical.tumor.size"
X3B.c$tumor_size
X3B.c$tumor_size = as.numeric(X3B.c$tumor_size)
X3B.c$tumor_size[X3B.c$tumor_size==99.0 | X3B.c$tumor_size<=0.0] = NA # there is evidently a mistake

X3B.c$event = as.factor(X3B$"rfs.e")

X3B.c$survival_time = as.numeric(X3B$rfs.t)/d_time #days 2 months 

X3B.c$age = NA
X3B.c$age = as.numeric(X3B.c$age)
X3B.c$ER= as.factor(X3B$"er")
X3B.c$PR= as.factor(X3B$"pr")
X3B.c$HER2= as.factor(X3B$"her2")
X3B.c$chemo = as.factor(1)
X3B.c$hormonal = as.factor(NA)
X3B.c$radiation = as.factor(NA)

X3B = cbind(X3B.c, X3B)
X3B[1:5,1:100]

##### GSE26971 #####
X4 = get(load(files[4]))

c = 20
names(X4)[1:c]
X4[1:5,1:c]

X4.c = empty_df(df_length = dim(X4)[1])

X4.c$node_status = X4$"Characteristics..node.positive."
table(X4.c$node_status)
X4.c$node_status = ifelse(X4.c$node_status=="no", 0, 1)
X4.c$node_status = as.factor(X4.c$node_status)

X4.c$grade = X4$"Characteristics..TumorGrading."
table(X4.c$grade)
class(X4.c$grade)
X4.c$grade = as.factor(X4.c$grade)

X4.c$stage = X4$"Characteristics..TumorStaging."
table(X4.c$stage)
class(X4.c$stage)
X4.c$stage = as.factor(X4.c$stage)

X4.c$tumor_size = X4$"Characteristics..tumor.size..mm.."
X4.c$tumor_size
class(X4.c$tumor_size)
X4.c$tumor_size = as.numeric(X4.c$tumor_size)

X4.c$event = X4$"Characteristics..dead.of.disease."
X4.c$event
X4.c$event[X4.c$event=="yes"] = 1
X4.c$event[X4.c$event=="no"] = 0
X4.c$event[X4.c$event=="NA"] = NA
X4.c$event[X4.c$event=="n/a"] = NA
X4.c$event = as.factor(X4.c$event)

X4.c$survival_time = NA
X4.c$survival_time[isevent(X4.c$event)] = X4$Characteristics..distal.recurrence.time.[isevent(X4.c$event)]
X4.c$survival_time[isnoevent(X4.c$event)] = X4$Characteristics..followup.time.[isnoevent(X4.c$event)]
X4.c$survival_time = as.numeric(X4.c$survival_time)*12



X4.c$age = X4$"Characteristics..age.at.diagnosis."
X4.c$age = as.numeric(X4.c$age)

X4.c$ER = X4$"Characteristics..EstrogenReceptorStatus."
X4.c$ER[X4.c$ER=="negative" & !is.na(X4.c$ER)] = 0
X4.c$ER[X4.c$ER=="positive" & !is.na(X4.c$ER)] = 1
X4.c$ER = as.factor(X4.c$ER)

X4.c$PR = X4$"Characteristics..Progesterone.Receptor.status."
X4.c$PR[X4.c$PR=="negative" & !is.na(X4.c$PR)] = 0
X4.c$PR[X4.c$PR=="positive" & !is.na(X4.c$PR)] = 1
X4.c$PR[X4.c$PR=="n/a"] = NA
X4.c$PR = as.factor(X4.c$PR)

X4.c$HER2 = X4$"Characteristics..ErbB2.positive..IHC.."
X4.c$HER2[X4.c$HER2=="n/a"] = NA
X4.c$HER2[X4.c$HER2=="no" & !is.na(X4.c$HER2)] = 0
X4.c$HER2[X4.c$HER2=="yes" & !is.na(X4.c$HER2)] = 1
X4.c$HER2 = as.factor(X4.c$HER2)

X4.c$chemo = X4$"Characteristics..chemotherapy.treatment."
X4.c$chemo[X4.c$chemo=="n/a"] = NA
X4.c$chemo[X4.c$chemo=="no" &  !is.na(X4.c$chemo)] = 0
X4.c$chemo[X4.c$chemo=="yes" & !is.na(X4.c$chemo)] = 1
X4.c$chemo = as.factor(X4.c$chemo)

X4.c$hormonal = X4$"Characteristics..hormonal.therapy."
X4.c$hormonal[X4.c$hormonal=="n/a"] = NA
X4.c$hormonal[X4.c$hormonal=="no" &  !is.na(X4.c$hormonal)] = 0
X4.c$hormonal[X4.c$hormonal=="yes" & !is.na(X4.c$hormonal)] = 1
X4.c$hormonal = as.factor(X4.c$hormonal)

X4.c$radiation = X4$"Characteristics..radiation.treatment."
X4.c$radiation[X4.c$radiation=="n/a"] = NA
X4.c$radiation[X4.c$radiation=="no" &  !is.na(X4.c$radiation)] = 0
X4.c$radiation[X4.c$radiation=="yes" & !is.na(X4.c$radiation)] = 1
X4.c$radiation = as.factor(X4.c$radiation)

X4 = cbind(X4.c, X4)
X4[1:5,1:100]


##### GSE96058 #####

X5 = get(load(files[5]))
X5$tumor_size = X5$tumor_size/10
X5$node_status = as.character(X5$node_status)
X5$node_status[X5$node_status=="NodePositive"] = 1
X5$node_status[X5$node_status=="NodeNegative"] = 0
X5$node_status = as.factor(X5$node_status)
table(X5$node_status)


##### GSE1456 #####
### X6A
X6A = get(load(files[6]))

names(X6A)[1:15]
X6A[1:5,1:50]

X6A.c = empty_df(df_length = dim(X6A)[1])
X6A.c$node_status = as.factor(NA)

X6A.c$grade =X6A$`ELSTON:ch1`
table(X6A.c$grade)
class(X6A.c$grade)

X6A.c$stage = NA
X6A.c$tumor_size = NA

X6A$`RELAPSE:ch1`
X6A$`DEATH_BC:ch1`

X6A.c$event = as.integer(as.integer(X6A$`RELAPSE:ch1`) | as.integer(X6A$`DEATH_BC:ch1`))
X6A.c$survival_time = as.numeric(X6A$`SURV_RELAPSE:ch1`)*12
X6A.c$age = NA

X6A.c$ER = 1
X6A.c$ER[X6A$`SUBTYPE:ch1`=="ERBB2"] = 0
X6A.c$ER[X6A$`SUBTYPE:ch1`=="Basal"] = 0
X6A.c$ER[X6A$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6A.c$ER)

X6A.c$PR = 1
X6A.c$PR[X6A$`SUBTYPE:ch1`=="ERBB2"] = 0
X6A.c$PR[X6A$`SUBTYPE:ch1`=="Basal"] = 0
X6A.c$PR[X6A$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6A.c$PR)

X6A.c$HER2 = NA
X6A.c$HER2[X6A$`SUBTYPE:ch1`=="ERBB2"] = 1
X6A.c$HER2[X6A$`SUBTYPE:ch1`=="Basal"] = 0
X6A.c$HER2[X6A$`SUBTYPE:ch1`=="Normal Like"] = 0
X6A.c$HER2[X6A$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6A.c$HER2)
X6A.c$chemo =     NA
X6A.c$hormonal =  NA
X6A.c$radiation = NA

### X6B

X6B = get(load("../data/GSE1456/data/GSE1456B.RData"))

X6B.c = empty_df(df_length = dim(X6B)[1])
X6B.c$node_status = as.factor(NA)

X6B.c$grade =X6B$`ELSTON:ch1`
table(X6B.c$grade)
class(X6B.c$grade)

X6B.c$stage = NA
X6B.c$tumor_size = NA

X6B$`RELAPSE:ch1`
X6B$`DEATH_BC:ch1`

X6B.c$event = as.integer(as.integer(X6B$`RELAPSE:ch1`) | as.integer(X6B$`DEATH_BC:ch1`))
X6B.c$survival_time = as.numeric(X6B$`SURV_RELAPSE:ch1`)*12
X6B.c$age = NA

X6B.c$ER = 1
X6B.c$ER[X6B$`SUBTYPE:ch1`=="ERBB2"] = 0
X6B.c$ER[X6B$`SUBTYPE:ch1`=="Basal"] = 0
X6B.c$ER[X6B$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6B.c$ER)

X6B.c$PR = 1
X6B.c$PR[X6B$`SUBTYPE:ch1`=="ERBB2"] = 0
X6B.c$PR[X6B$`SUBTYPE:ch1`=="Basal"] = 0
X6B.c$PR[X6B$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6B.c$PR)

X6B.c$HER2 = NA
X6B.c$HER2[X6B$`SUBTYPE:ch1`=="ERBB2"] = 1
X6B.c$HER2[X6B$`SUBTYPE:ch1`=="Basal"] = 0
X6B.c$HER2[X6B$`SUBTYPE:ch1`=="Normal Like"] = 0
X6B.c$HER2[X6B$`SUBTYPE:ch1`=="No Subtype"] = NA
table(X6B.c$HER2)
X6B.c$chemo =     NA
X6B.c$hormonal =  NA
X6B.c$radiation = NA
# clinical info is the same for both platforms and match row order:
X6A = cbind(X6A.c, X6A)
X6B = cbind(X6A.c, X6B)

 
##### GSE42568 #####
X7 = get(load(files[7]))
X7 = subset(X7, X7$characteristics_ch1=="tissue: breast cancer")

X7.c = empty_df(df_length = dim(X7)[1])
X7.c$node_status = X7$"lymph node status:ch1"
table(X7.c$node_status)
X7.c$node_status = as.factor(X7.c$node_status)

X7.c$grade = X7$"grade:ch1"
table(X7.c$grade)
X7.c$grade = as.factor(X7.c$grade)

X7.c$stage = as.factor(NA)

X7.c$tumor_size =X7$"size:ch1"
X7.c$tumor_size
X7.c$tumor_size = as.numeric(X7.c$tumor_size)

X7.c$event = X7$"relapse free survival event:ch1"
X7.c$event = as.factor(as.numeric(X7.c$event))

X7.c$survival_time = as.numeric(X7$"relapse free survival time_days:ch1")/d_time

X7.c$age = X7$"age:ch1"
X7.c$age = as.numeric(X7.c$age)

X7.c$ER = X7$"er_status:ch1"
X7.c$ER = as.factor(as.numeric(X7.c$ER))

X7.c$PR =        as.factor(NA)
X7.c$HER2 =      as.factor(NA)
X7.c$chemo =     as.factor(0)
X7.c$hormonal =  as.factor(0)
X7.c$radiation = as.factor(NA)


X7 = cbind(X7.c, X7)

  
##### GSE45255 #####
X8 = get(load(files[8]))
dim(X8)
names(X8)[53:71]
X8 = subset(X8, !grepl("NUH", X8[,1],))


X8.c = empty_df(df_length = dim(X8)[1])
X8.c$node_status = X8$"ln status:ch1"
table(X8.c$node_status)
X8.c$node_status = ifelse(X8.c$node_status=="LN+", 1, 0)
X8.c$node_status = as.factor(X8.c$node_status)

X8.c$grade =X8$"histological grade:ch1"
table(X8.c$grade)
X8.c$grade[X8.c$grade=="G1"] = 1
X8.c$grade[X8.c$grade=="G2"] = 2
X8.c$grade[X8.c$grade=="G3"] = 3
X8.c$grade[X8.c$grade=="NA"] = NA
X8.c$grade = as.factor(X8.c$grade)

X8.c$stage = as.factor(NA)

X8.c$tumor_size =X8$"size (mm):ch1" 
X8.c$tumor_size
X8.c$tumor_size = as.numeric(X8.c$tumor_size)/10.0

X8.c$event = X8$"dmfs event (defined as distant metastasis or death from breast cancer):ch1"
X8.c$event = as.factor(as.numeric(X8.c$event))

X8.c$survival_time = as.numeric(X8$"dmfs time:ch1")*12

X8.c$age = X8$"patient age:ch1"
X8.c$age = as.numeric(X8.c$age)

X8.c$ER =X8$"er status:ch1"
X8.c$ER[X8.c$ER=="ER+"] = 1
X8.c$ER[X8.c$ER=="ER-"] = 0
X8.c$ER[X8.c$ER=="NA"] = NA
X8.c$ER = as.factor(X8.c$ER)

X8.c$PR = X8$"pgr status:ch1"
X8.c$PR[X8.c$PR=="PgR+"] = 1
X8.c$PR[X8.c$PR=="PgR-"] = 0
X8.c$PR[X8.c$PR=="NA"] = NA
X8.c$PR = as.factor(X8.c$PR)

X8.c$HER2 = X8$"her2 status:ch1"
X8.c$HER2[X8.c$HER2=="He+"] = 1
X8.c$HER2[X8.c$HER2=="He-"] = 0
X8.c$HER2[X8.c$HER2=="NA"] = NA
X8.c$HER2 = as.factor(X8.c$HER2)

X8.c$chemo =X8$"chemo? (0=no, 1=yes):ch1"
X8.c$chemo = as.factor(as.numeric(X8.c$chemo))

X8.c$hormonal = X8$"characteristics:ch1"
X8.c$hormonal = as.factor(as.numeric(gsub("(.*): ", "", X8.c$hormonal)))
X8.c$radiation = as.factor(NA)

X8 = cbind(X8.c, X8)

##### GSE4922 #####

###GSE4922A

X9A = get(load(files[9]))
names(X9A)[51:67]


X9A.c = empty_df(df_length = dim(X9A)[1])
X9A.c$node_status = X9A$"Lymph node status:ch1"
table(X9A.c$node_status)
X9A.c$node_status[X9A.c$node_status=="LN+"] = 1
X9A.c$node_status[X9A.c$node_status=="LN-"] = 0
X9A.c$node_status[X9A.c$node_status=="LN?"] = NA
X9A.c$node_status = as.factor(X9A.c$node_status)

X9A.c$grade = X9A$"Elston (NGS) histologic grade:ch1"
table(X9A.c$grade)
X9A.c$grade = as.factor(X9A.c$grade)

X9A.c$stage = as.factor(NA)

X9A.c$tumor_size = X9A$"tumor size (mm):ch1" 
X9A.c$tumor_size
X9A.c$tumor_size = as.numeric(X9A.c$tumor_size)/10.0

X9A.c$event =X9A$"DFS EVENT (0=censored; 1=event defined as any type of recurrence (local, regional or distant) or death from breast cancer:ch1"
X9A.c$event = as.factor(as.numeric(X9A.c$event))

X9A.c$survival_time = as.numeric(X9A$"DFS TIME (yrs):ch1")*12.0

X9A.c$age = X9A$"age at diagnosis:ch1"
X9A.c$age = as.numeric(X9A.c$age)

X9A.c$ER = X9A$"ER status:ch1"
X9A.c$ER[X9A.c$ER=="ER+"] = 1
X9A.c$ER[X9A.c$ER=="ER-"] = 0
X9A.c$ER[X9A.c$ER=="ER?"] = NA
X9A.c$ER[X9A.c$ER=="NA"] = NA
X9A.c$ER = as.factor(X9A.c$ER)

X9A.c$PR = as.factor(NA)
X9A.c$HER2 = as.factor(NA)

X9A.c$chemo = as.factor(NA)

X9A.c$hormonal = X9A$"ER+, endocrine therapy only (1=included in survival analysis):ch1" 
X9A.c$hormonal[X9A.c$hormonal=="1"]=1
X9A.c$hormonal[X9A.c$hormonal==""]=0
X9A.c$hormonal[X9A.c$hormonal=="NA"]=NA
X9A.c$hormonal= as.factor(X9A.c$hormonal)

X9A.c$radiation = as.factor(NA)

X9A = cbind(X9A.c, X9A)


###GSE4922B

X9B = get(load("../data/GSE4922/data/GSE4922B.RData"))
X9B = cbind(X9A.c, X9B)





##### GSE7390 #####

X10 = get(load(files[10]))
names(X10)[63:89]


X10.c = empty_df(df_length = dim(X10)[1])
X10.c$node_status =X10$"node:ch1"
table(X10.c$node_status)
X10.c$node_status = as.factor(X10.c$node_status)

X10.c$grade = X10$"grade:ch1"
table(X10.c$grade)
X10.c$grade[X10.c$grade=="NA"] = NA
X10.c$grade = as.factor(X10.c$grade)

X10.c$stage = as.factor(NA)

X10.c$tumor_size = X10$"size:ch1" 
X10.c$tumor_size
X10.c$tumor_size = as.numeric(X10.c$tumor_size)

X10.c$event = X10$"e.dmfs:ch1"
X10.c$event = as.factor(as.numeric(X10.c$event))

X10.c$survival_time = as.numeric(X10$"t.dmfs:ch1")
X10.c$survival_time = X10.c$survival_time/d_time

X10.c$age =X10$"age:ch1"
X10.c$age = as.numeric(X10.c$age)

X10.c$ER = X10$"er:ch1"
X10.c$ER = as.factor(as.numeric(X10.c$ER))

X10.c$PR = as.factor(NA)
X10.c$HER2 = as.factor(NA)

X10.c$chemo = as.factor(0)
X10.c$hormonal= as.factor(0)

X10.c$radiation = as.factor(NA)

X10 = cbind(X10.c, X10)




##### GSE7849 #####

X11 = get(load(files[11]))
X11[1:5,1:68]

X11.c = empty_df(df_length = dim(X11)[1])

pat = "DFS ="
tmp=NA
col = X11$"characteristics_ch1.15"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.16"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.17"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.18"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.19"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.20"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.21"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.22"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.23"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.24"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$survival_time = as.numeric(gsub("([^-]*)= ", "", x = tmp ))

X11.c$grade =X11$"characteristics_ch1.5"
table(X11.c$grade)
X11.c$grade[X11.c$grade=="Histo_Grade = Moderately Diff"] = 2 
X11.c$grade[X11.c$grade=="Histo_Grade = Not Graded / Unknown"] = NA 
X11.c$grade[X11.c$grade=="Histo_Grade = Poorly Diff"] = 3
X11.c$grade[X11.c$grade=="Histo_Grade = Undifferentiated"] = 3 
X11.c$grade[X11.c$grade=="Histo_Grade = Well Diff"] =  1
X11.c$grade = as.factor(X11.c$grade)

pat = "^[Tumor_size]"
tmp=NA
col = X11$"characteristics_ch1.8"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.9"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.10"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$tumor_size  = as.numeric(gsub("([^-]*)= ", "", x = tmp ))

pat = "^ER = "
tmp=NA
col = X11$"characteristics_ch1.10"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.11"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.12"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.13"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$ER = tmp
table(X11.c$ER)
X11.c$ER[X11.c$ER=="ER = Borderline,Undetermined"] = NA 
X11.c$ER[X11.c$ER=="ER = Negative/Normal"] = 0 
X11.c$ER[X11.c$ER=="ER = Positive/elevated"] = 1
X11.c$ER[X11.c$ER=="ER = Unknown if Done"] = NA 
X11.c$ER = as.factor(X11.c$ER)

pat = "^PR = "
tmp=NA
col = X11$"characteristics_ch1.10"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.11"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.12"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.13"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$PR = tmp
table(X11.c$PR)
X11.c$PR[X11.c$PR=="PR = Borderline,Undetermined"] = NA 
X11.c$PR[X11.c$PR=="PR = Negative/Normal"] = 0 
X11.c$PR[X11.c$PR=="PR = Positive/elevated"] = 1
X11.c$PR[X11.c$PR=="PR = Unknown if Done"] = NA 
X11.c$PR = as.factor(X11.c$PR)

pat = "^Age"
tmp=NA
col = X11$"characteristics_ch1.15"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.16"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.17"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.18"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.19"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.20"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.21"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$age  = as.numeric(gsub("([^-]*)= ", "", x = tmp ))


pat = "^LVI "
tmp=NA
col = X11$"characteristics_ch1.14"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.15"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.16"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.17"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.18"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.19"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.20"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
tmp[tmp=="LVI = present"] = 1
tmp[tmp=="LVI = absent"] = 0
X11.c$node_status = as.factor(as.numeric(tmp))

X11.c$stage = as.factor(NA)

pat = "1st_Distant_Site_Recur"
tmp=NA
col = X11$"characteristics_ch1.18"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.19"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.20"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.21"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.22"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.23"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.24"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$event = ifelse(tmp=="1st_Distant_Site_Recur = None", 1, 0)
X11.c$event = as.factor(X11.c$event)

X11$HER2 = as.factor(NA)

pat = "XRT"
tmp=NA
col = X11$"characteristics_ch1.22"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.23"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.24"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.25"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.26"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.27"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.28"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$radiation = ifelse(tmp=="Type_1st_XRT = None", 0, 1)
X11.c$radiation = as.factor(X11.c$radiation)

pat = "Hormone"
tmp=NA
col = X11$"characteristics_ch1.25"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.26"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.27"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.28"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.29"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.30"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.31"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.32"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$hormonal = ifelse(tmp=="Hormone_Tx = Yes", 1, 0)
X11.c$hormonal = as.factor(X11.c$hormonal)

pat = "Chemotherapy"
tmp=NA
col = X11$"characteristics_ch1.25"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.26"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.27"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.28"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.29"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.30"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.31"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
col = X11$"characteristics_ch1.32"
tmp[grep(pattern = pat, col)] = col[grep(pattern = pat, col)]
X11.c$chemo = ifelse(tmp=="Chemotherapy = Yes", 1, 0)
X11.c$chemo = as.factor(X11.c$chemo)

X11 = cbind(X11.c, X11)


##### GSE9195 #####

X12 = get(load(files[12]))
names(X12)[49:65]

X12.c = empty_df(df_length = dim(X12)[1])
X12.c$node_status = X12$"node:ch1"
table(X12.c$node_status)
X12.c$node_status = as.factor(X12.c$node_status)

X12.c$grade =X12$"grade:ch1"
table(X12.c$grade)
X12.c$grade[X12.c$grade=="NA"] = NA
X12.c$grade = as.factor(X12.c$grade)

X12.c$stage = as.factor(NA)

X12.c$tumor_size = X12$"size:ch1" 
X12.c$tumor_size
X12.c$tumor_size = as.numeric(X12.c$tumor_size)

X12.c$event = X12$"e.dmfs:ch1"
X12.c$event = as.factor(as.numeric(X12.c$event))

X12.c$survival_time = as.numeric(X12$"t.dmfs:ch1")
X12.c$survival_time = X12.c$survival_time/(d_time)

X12.c$age = X12$"age:ch1"
X12.c$age = as.numeric(X12.c$age)

X12.c$ER = X12$"er:ch1"
X12.c$ER = as.factor(as.numeric(X12.c$ER))

X12.c$PR = X12$"pgr:ch1"
X12.c$PR = as.factor(as.numeric(X12.c$PR))
X12.c$HER2 = as.factor(NA)

X12.c$chemo = as.factor(0)
X12.c$hormonal= as.factor(1)

X12.c$radiation = as.factor(NA)

X12 = cbind(X12.c, X12)


##### GSE9893 #####

X13 = get(load(files[13]))
names(X13)[53:71]

X13.c = empty_df(df_length = dim(X13)[1])
X13.c$node_status = X13$"n dich+:ch1"
table(X13.c$node_status)
X13.c$node_status[X13.c$node_status=="NA"] = NA
X13.c$node_status = as.factor(X13.c$node_status)

X13.c$grade =X13$"sbr grade:ch1"
table(X13.c$grade)
X13.c$grade[X13.c$grade=="NA"] = NA
X13.c$grade = as.factor(X13.c$grade)

X13.c$stage = as.factor(NA)

X13.c$tumor_size = X13$"tumor size (mm):ch1" 
X13.c$tumor_size
X13.c$tumor_size = as.numeric(X13.c$tumor_size)/10

X13.c$event =X13$"state of health:ch1"
X13.c$survival_time = NA
# decease
X13.c$survival_time[X13.c$event=="deceased"] = X13$`decease delay after surgery (months):ch1`[X13.c$event=="deceased"]
# recurrence
X13.c$survival_time[X13$`time before recurrence (months):ch1`!=""]=X13$`time before recurrence (months):ch1`[X13$`time before recurrence (months):ch1`!=""]
# alive
X13.c$survival_time[X13.c$event=="alive"] = X13$`follow-up period (months):ch1`[X13.c$event=="alive"]
X13.c$survival_time[X13.c$event=="alive with meta"] = X13$`follow-up period (months):ch1`[X13.c$event=="alive with meta"]
X13.c$survival_time[X13.c$survival_time==""] = X13$`follow-up period (months):ch1`[X13.c$survival_time==""]
X13.c$event = ifelse(X13.c$event=="alive", 0, 1)
X13.c$event = as.factor(X13.c$event)
X13.c$survival_time = as.numeric(X13.c$survival_time)



X13.c$age = X13$"age:ch1"
X13.c$age = as.numeric(X13.c$age)


# see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2912334/ for threshlding ER and PR
X13.c$ER = X13$"characteristics_ch1.11"
X13.c$ER = gsub("ER ", "", X13.c$ER)
X13.c$ER = gsub("<", "", X13.c$ER)
X13.c$ER = as.numeric(X13.c$ER)
X13.c$ER[is.na(X13.c$ER)] = 1
X13.c$ER = as.factor(ifelse(X13.c$ER<20, 0, 1))

X13.c$PR = X13$"pr:ch1"
X13.c$PR = ifelse(X13.c$PR=="<10", 0, 1)
X13.c$PR = as.factor(as.numeric(X13.c$PR))

X13.c$HER2 = as.factor(NA)

X13.c$chemo = as.factor(0)
X13.c$hormonal= as.factor(1)

X13.c$radiation = as.factor(ifelse(X13$`adjuvant therapy:ch1`=="X-ray+Tam", 1, 0))

X13 = cbind(X13.c, X13)


##### GSE2603 #####

X14 = get(load(files[14]))
names(X14)[56:71]

X14.c = empty_df(df_length = dim(X14)[1])
X14.c$node_status = X14$"pos lymph nodes:ch1"
table(X14.c$node_status)
X14.c$node_status = as.factor(ifelse(X14.c$node_status==0, 0, 1))

X14.c$grade = as.factor(NA)
X14.c$stage = as.factor(NA)

X14.c$tumor_size = X14$"tumor size (cm):ch1" 
X14.c$tumor_size
X14.c$tumor_size = as.numeric(X14.c$tumor_size)

X14.c$event =X14$"met event:ch1"
X14.c$event[X14.c$event=="1"] = 1
X14.c$event[X14.c$event=="0"] = 0
X14.c$event[X14.c$event=="--"] = NA
X14.c$event = as.factor(X14.c$event)


X14.c$survival_time = X14$"mfs (yr):ch1"
X14.c$survival_time = as.numeric(X14.c$survival_time)*12



X14.c$age = X14$"age at dx:ch1"
X14.c$age = as.numeric(X14.c$age)

X14.c$ER = as.factor(ifelse(X14$"path er status:ch1"=="N", 0, 1))
X14.c$PR = as.factor(ifelse(X14$"path pr status:ch1"=="N", 0, 1))
X14.c$HER2 = as.factor(ifelse(X14$"her2 status:ch1"=="N",0,1))

X14.c$chemo = as.factor(0)
X14.c$hormonal= as.factor(1)

X14.c$radiation = as.factor(NA)

X14 = cbind(X14.c, X14)

##### NKI #####

X15 = get(load(files[15]))
names(X15)[1:22]

X15.c = empty_df(df_length = dim(X15)[1])
X15.c$node_status = X15$"node"
table(X15.c$node_status)
X15.c$node_status = as.factor(X15.c$node_status)

X15.c$grade = as.factor(X15$grade)

X15.c$stage = as.factor(NA)

X15.c$tumor_size =X15$"size" 
X15.c$tumor_size
X15.c$tumor_size = as.numeric(X15.c$tumor_size)

X15.c$event = X15$"e.rfs"
X15.c$event = as.factor(X15.c$event)

X15.c$survival_time = X15$"t.rfs"
X15.c$survival_time = as.numeric(X15.c$survival_time)/d_time

X15.c$age = as.numeric(X15$age)

X15.c$ER = as.factor(X15$"er")
X15.c$PR = as.factor(X15$"pgr")
X15.c$HER2 = as.factor(X15$her2)

#N Engl J Med 2002; 347:1999-2009
#DOI: 10.1056/NEJMoa021967
#Ten of the 151 patients who had lymph-node–negative disease and 120 of the 144 who had lymph-node–positive 
#disease had received adjuvant systemic therapy consisting of chemotherapy (90 patients), 
#hormonal therapy (20), or both (20). 

X15.c$chemo = as.factor(ifelse(X15$treatment==1, 1, 0))
X15.c$hormonal= as.factor(ifelse(X15$treatment==2, 1, 0))

X15.c$radiation = as.factor(0)

X15 = cbind(X15.c, X15)

##### TCGA ####
X16 = get(load(files[16]))
dim(X16)
X16[1:5,1:100]

X16.c = empty_df(df_length = dim(X16)[1])

X16.c$node_status = 1
X16.c$node_status[grep("N0", X16$ajcc_pathologic_n)] = 0
X16.c$node_status[grep("NX", X16$ajcc_pathologic_n)] = NA

X16.c$grade = NA
X16.c$tumor_size = NA
X16.c$survival_time[X16$vital_status=="Dead" 
                    & !is.na(X16$vital_status)] = X16$days_to_death[X16$vital_status=="Dead"
                                                                    & !is.na(X16$vital_status)]/d_time
X16.c$survival_time[X16$vital_status=="Alive" 
                    & !is.na(X16$vital_status)] = X16$days_to_last_follow_up[X16$vital_status=="Alive"
                                                                             & !is.na(X16$vital_status)]/d_time
X16.c$event = ifelse(X16$vital_status=="Dead" & !is.na(X16$vital_status), 1, 0)
X16.c$event[is.na(X16$vital_status)] =  NA

X16.c$age = X16$age_at_diagnosis/(d_time*12)
X16.c$chemo = 0
X16.c$chemo[is.na(X16$therapy)] = NA
X16.c$chemo[X16$therapy=="Chemotherapy" | X16$therapy=="ChemotherapyHormone Therapy"] =  1
X16.c$hormonal = 0
X16.c$hormonal[is.na(X16$therapy)] = NA
X16.c$hormonal[X16$therapy=="Hormone Therapy" | X16$therapy=="ChemotherapyHormone Therapy"] =  1
X16.c$stage[X16$paper_pathologic_stage=="Stage_I"] = 1
X16.c$stage[X16$paper_pathologic_stage=="Stage_II"] = 2
X16.c$stage[X16$paper_pathologic_stage=="Stage_III"] = 3
X16.c$stage[X16$paper_pathologic_stage=="Stage_IV"] = 4
X16.c$ER[X16$estrogen_receptor_status=="Negative"] = 0
X16.c$ER[X16$estrogen_receptor_status=="Positive"] = 1
X16.c$PR[X16$progesterone_receptor_status=="Negative"] = 0
X16.c$PR[X16$progesterone_receptor_status=="Positive"] = 1
X16.c$HER2 = 1
X16.c$HER2[is.na(X16$her2_level_cell_percent_category)] = NA
X16.c$HER2[X16$her2_level_cell_percent_category==""] = NA
X16.c$HER2[X16$her2_level_cell_percent_category=="<10%"] = 0
X16.c$radiation[X16$radiation_therapy=="YES"] = 1
X16.c$radiation[X16$radiation_therapy=="NO"] = 0

X16 = cbind(X16.c, X16)

##### CHECK FOR DUPLICATES #####

Xs = list(X1, X2,X3A,X3B,X4,X5,X6A,X7,X8,X9A,X10,X11,X12,X13,X14,X15)
m=0
for (x in Xs) {
  m = m + dim(x)[1]
  print(dim(x)[1])
}
m #[1] 5518

n= c(
X1$geo_accession,
X2$geo_accession,
X3A$Row.names,
X3B$Row.names,
X5$GSM_ID,
rownames(X6A),
X7$geo_accession,
X8$geo_accession,
X9A$geo_accession,
X10$geo_accession,
X11$geo_accession,
X12$geo_accession,
X13$geo_accession,
X14$geo_accession
)

n_occur <- data.frame(table(n))
n_occur[n_occur$Freq > 1,]
dim(n_occur[n_occur$Freq > 1,])

##### METABRIC #####
MB = get(load("../../../metabric_new_raw_data/data/METABRIC_data.RData"))
dim(MB)

MB.c = empty_df(df_length = dim(MB)[1])

MB.c$age =         MB$Age.at.Diagnosis 
MB.c$chemo =       as.factor(ifelse(MB$Chemotherapy=="YES",1,0))
MB.c$ER =          as.factor(ifelse(MB$ER.Status=="Positive",1,0))
MB.c$PR =          as.factor(ifelse(MB$PR.Status=="Positive",1,0))
MB.c$HER2 =        as.factor(ifelse(MB$HER2.Status=="Positive",1,0))
MB.c$hormonal =    as.factor(ifelse(MB$Hormone.Therapy=="YES",1,0))
MB.c$node_status = as.factor(ifelse(MB$Lymph.nodes.examined.positive==0,0,1))
MB.c$event =       as.factor(ifelse(MB$Overall.Survival.Status=="1:DECEASED",1,0))
MB.c$radiation =   as.factor(ifelse(MB$Radio.Therapy=="YES",1,0))
MB.c$tumor_size =  MB$Tumor.Size/10.0
MB.c$stage =       as.factor(MB$Tumor.Stage)
MB.c$grade   =     as.factor(MB$Neoplasm.Histologic.Grade)
MB.c$survival_time=MB$Overall.Survival..Months.

MB = cbind(MB.c, MB)
sb(MB)
se(MB)
data.MB.clin = MB[,1:51]
data.MB.expr = MB[,52:ncol(MB)]

##### ExprImage ##### 

E =      get(load('../../ExprImage/data/Xraw.Rdata'))
e =      get(load("../../ExprImage/data/X_age75_fu9_tte12.Rdata"))

E$node_status   = as.factor(ifelse((E$pN)==0, 0, 1)) 
E$grade         = E$Grad
E$tumor_size    = E$Tumordurchmesser 
E$survival_time = as.numeric(E$survival_time/31) 
E$event         = E$event_jn 
E$age           = E$Alter 
E$chemo         = as.factor(ifelse(as.character(E$Chemo_jn)=="0", 0, 1)) 
E$hormonal      = E$Antihormon_jn 
E$radiation     = E$Radiatio_jn
E$stage         = E$pT 
E$ER            = E$ER 
E$PR            = E$PR 
E$HER2          = E$Her2neu_jn 

e$node_status   = as.factor(ifelse((e$pN)==0, 0, 1)) 
e$grade         = e$Grad
e$tumor_size    = e$Tumordurchmesser 
e$survival_time = as.numeric(e$survival_time/31) 
e$event         = e$event_jn 
e$age           = e$Alter 
e$chemo         = as.factor(ifelse(as.character(e$Chemo_jn)=="0", 0, 1)) 
e$hormonal      = e$Antihormon_jn 
e$radiation     = e$Radiatio_jn
e$stage         = e$pT 
e$ER            = e$ER 
e$PR            = e$PR 
e$HER2          = e$Her2neu_jn 



data.ExprImage.raw = E
data.ExprImage.fil = e

##### Divide into clinical and expression ######

sb(X1)
c = grep("A1CF", colnames(X1))
data.X1.clin =      X1[, 1:(c-1)]
data.X1.expr =      X1[, c:ncol(X1)]

sb(X2)
c = grep("A1BG", colnames(X2))
data.X2.clin =      X2[, 1:(c-1)]
data.X2.expr =      X2[, c:ncol(X2)]

sb(X3A)
c = grep("LOC400451", colnames(X3A))
data.X3A.clin =      X3A[, 1:(c-1)]
data.X3A.expr =      X3A[, c:ncol(X3A)]

sb(X3B)
c = grep("LOC400451", colnames(X3B))
data.X3B.clin =      X3B[, 1:(c-1)]
data.X3B.expr =      X3B[, c:ncol(X3B)]

sb(X4)
c = grep("A1CF", colnames(X4))
data.X4.clin =      X4[, 1:(c-1)]
data.X4.expr =      X4[, c:ncol(X4)]

names(X5)[1:150]
c = grep("5_8S_rRNA", colnames(X5))
data.X5.clin =      X5[, 1:(c-1)]
data.X5.expr =      X5[, c:ncol(X5)]

names(X6A)[1:150]
c = grep("A1CF", colnames(X6A))
data.X6A.clin =      X6A[, 1:(c-1)]
data.X6A.expr =      X6A[, c:ncol(X6A)]
names(X6B)[1:150]
c = grep("A1BG", colnames(X6B))
data.X6B.clin =      X6B[, 1:(c-1)]
data.X6B.expr =      X6B[, c:ncol(X6B)]

names(X7)[1:150]
c = grep("A1BG", colnames(X7))
data.X7.clin =      X7[, 1:(c-1)]
data.X7.expr =      X7[, c:ncol(X7)]

names(X8)[1:150]
c = grep("A1CF", colnames(X8))
data.X8.clin =      X8[, 1:(c-1)]
data.X8.expr =      X8[, c:ncol(X8)]

names(X9A)[1:150]
c = grep("A1CF", colnames(X9A))
data.X9A.clin =      X9A[, 1:(c-1)]
data.X9A.expr =      X9A[, c:ncol(X9A)]
names(X9B)[1:150]
c = grep("A1BG", colnames(X9B))
data.X9B.clin =      X9B[, 1:(c-1)]
data.X9B.expr =      X9B[, c:ncol(X9B)]

names(X10)[1:150]
c = grep("A1CF", colnames(X10))
data.X10.clin =      X10[, 1:(c-1)]
data.X10.expr =      X10[, c:ncol(X10)]

names(X11)[1:150]
c = grep("MAPK3", colnames(X11))
data.X11.clin =      X11[, 1:(c-1)]
data.X11.expr =      X11[, c:ncol(X11)]

names(X12)[1:150]
c = grep("A1BG", colnames(X12))
data.X12.clin =      X12[, 1:(c-1)]
data.X12.expr =      X12[, c:ncol(X12)]

names(X13)[1:150]
c = grep("RRM1", colnames(X13))
data.X13.clin =      X13[, 1:(c-1)]
data.X13.expr =      X13[, c:ncol(X13)]

names(X14)[1:150]
c = grep("A1CF", colnames(X14))
data.X14.clin =      X14[, 1:(c-1)]
data.X14.expr =      X14[, c:ncol(X14)]

names(X15)[1:150]
c = grep("GREM2", colnames(X15))
data.X15.clin =      X15[, 1:(c-1)]
data.X15.expr =      X15[, c:ncol(X15)]

names(X16)[1:150]
c = grep("TSPAN6", colnames(X16))
data.X16.clin =      X16[, 1:(c-1)]
data.X16.expr =      X16[, c:ncol(X16)]


##### SAVE DATA #####
data.ExprImage.raw$study_id = "ExprImage"
data.ExprImage.fil$study_id = "ExprImage"
data.MB.clin$study_id = "METABRIC"
data.X1.clin$study_id  = "GSE11121"
data.X2.clin$study_id  = "GSE19615"
data.X3A.clin$study_id  = "GSE22226"
data.X3B.clin$study_id  = "GSE22226"
data.X4.clin$study_id  = "GSE26971"
data.X5.clin$study_id  = "GSE96058"
data.X6A.clin$study_id = "GSE1456"
data.X6B.clin$study_id = "GSE1456"
data.X7.clin$study_id  = "GSE42568"
data.X8.clin$study_id  = "GSE45255"
data.X9A.clin$study_id = "GSE4922"
data.X9B.clin$study_id = "GSE4922"
data.X10.clin$study_id = "GSE7390"
data.X11.clin$study_id = "GSE7849"
data.X12.clin$study_id = "GSE9195"
data.X13.clin$study_id = "GSE9893"
data.X14.clin$study_id = "GSE2603"
data.X15.clin$study_id = "NKI"
data.X16.clin$study_id = "TCGA"

save(data.ExprImage.raw,
     data.ExprImage.fil,
     data.MB.clin,
     data.MB.expr,
     data.X1.clin,
     data.X1.expr,
     data.X2.clin,
     data.X2.expr,
     data.X3A.clin,
     data.X3A.expr,
     data.X3B.clin,
     data.X3B.expr,
     data.X4.clin,
     data.X4.expr,
     data.X5.clin,
     data.X5.expr,
     data.X6A.clin,
     data.X6A.expr,
     data.X6B.clin,
     data.X6B.expr,
     data.X7.clin,
     data.X7.expr,
     data.X8.clin,
     data.X8.expr,
     data.X9A.clin,
     data.X9A.expr,
     data.X9B.clin,
     data.X9B.expr,
     data.X10.clin,
     data.X10.expr,
     data.X11.clin,
     data.X11.expr,
     data.X12.clin,
     data.X12.expr,
     data.X13.clin,
     data.X13.expr,
     data.X14.clin,
     data.X14.expr,
     data.X15.clin,
     data.X15.expr,
     data.X16.clin,
     data.X16.expr,
     file = "data.Rdata"
)









