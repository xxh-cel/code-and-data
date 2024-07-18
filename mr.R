library(R.matlab)
setwd("E:\\NSGA-III\\NSGA-III\\NSGA-III")
inipara <- readMat("chromo.mat")
Population=inipara$chromo
inipara <- readMat("p.mat")
p=inipara$dat
inipara <- readMat("M.mat")
M=inipara$f.num
M=M[1]
setwd("F:\\dssa\\code") 
MA=c()
sz=c()
Gen=150#迭代次数
pop=990#种群大小
ptm=proc.time()
setwd("F:\\dssa\\code")
Populationi=matrix(floor(150*runif(66*3)) ,ncol = 3,nrow = 66)
Populationf=matrix(floor(150*runif(66*3)) ,ncol = 3,nrow = 66)
Populationi=Population[,1:3]
Populationf=Population[,4:6]
Population=cbind(Populationi,Populationf)
Population=cbind(Populationi,Populationf)
inputxlsx=paste0("gtest.xlsx")
library(openxlsx)
treatment=read.xlsx(inputxlsx,sheet = 'Treatment',startRow = 4)
simulation=read.xlsx(inputxlsx,sheet = 'Simulation',startRow = 4)
irrigation=read.xlsx(inputxlsx,sheet = 'Irrigation',startRow = 4)
fer=read.xlsx(inputxlsx,sheet = 'Fertilizers',startRow = 4)
initial=read.xlsx(inputxlsx,sheet = 'Initial Conditions',startRow = 3)
addition=read.xlsx(inputxlsx,sheet = 'Addition',startRow = 4)
soilanalysis=read.xlsx(inputxlsx,sheet = 'Soilanalysis',startRow = 4)
envir=read.xlsx(inputxlsx,sheet = 'Environment',startRow = 4)
tillage =read.xlsx(inputxlsx,sheet = 'Tillage',startRow = 4)
harvest =read.xlsx(inputxlsx,sheet = 'Harvest',startRow = 4)
if(M==2){
  k=0
  for (i in 1:(pop/15)) {
    for (j in 1:3) {
      k=k+1
      irrigation$IRVAL[k]=Populationi[i,j]
    }
    for (u in 1:14) {
      for (uu in 1:3) {
        k=k+1
        irrigation$IRVAL[k]=Populationi[i,uu]
      }
    }
  }
  for (i in 1:pop) {
    for (j in 1:3) {
      k=(i-1)*3+j
      irrigation$IDATE[k] =p[[1]][[1]][j]#杩欓噷瀹為檯涓婂綋鎴戜滑鍙墿涓嬩袱涓洰鏍囩殑鏃跺€欙紝灏辩敤ff
      fer$FDATE[k]=p[[1]][[1]][j]
    }
  }
} else{
  k=0
  for (i in 1:(pop/15)) {
    for (j in 1:3) {
      k=k+1
      irrigation$IRVAL[k]=Populationi[i,j]
    }
    for (u in 1:14) {
      for (uu in 1:3) {
        k=k+1
        irrigation$IRVAL[k]=Populationi[i,uu]
      }
    }
  }
  for (i in 1:pop) {
    for (j in 1:3) {
      k=(i-1)*3+j
      irrigation$IDATE[k] =p[[1]][[1]][j]#杩欓噷瀹為檯涓婂綋鎴戜滑鍙墿涓嬩袱涓洰鏍囩殑鏃跺€欙紝灏辩敤ff
      fer$FDATE[k]=p[[1]][[1]][j]
    }
  }
  k=0
  for (i in 1:(pop/15)) {
    for (j in 1:3) {
      k=k+1
      fer$FAMN[k] =Populationf[i,j]
    }
    for (u in 1:14) {
      for (uu in 1:3) {
        k=k+1
        fer$FAMN[k] =Populationf[i,uu]
      }
    }
  }
  
}
sheets = list(" T" =  treatment,"Addition" = addition,"Initial Conditions" =initial ,"Irrigation" =irrigation ,
              "Fertilizers" = fer, "Simulation" = simulation,"Soilanalysis" =soilanalysis,
              "Environment" =envir ,"Tillage" = tillage,"Harvest" =
                harvest)
write.xlsx(sheets,paste0("gtest1.xlsx"))
setwd("F:\\dssa\\code")
source("write dssat file.R")
inputxlsx=paste0("gtest1.xlsx")
inputdata=Preparedata(pop,inputxlsx,IR_FER_all = F)
make_xbuid(pop,inputdata,overwrite = T)
source("make batch file.R")
# makebatchfile("Ri?e",inputdata)
makebatchfile(crop="Maize",inputdata,"DSSBatch.v46")
system("c:/DSSAT46/DSCSM046.exe B F:/dssa/code/DSSBatch.v46")
summ<-readLines("Summary.OUT")
if(M==4){
  nl<-readLines("SoilNiBal.OUT")
  
  nll=grep("!   Leached NO3",nl)
  nle=substr(nl[nll],72,78)
  nle=as.numeric(nle)
  
}
summary(summ)
first=grep("@",summ)
Treatmentname=substr(summ,4,6)

outname=unlist(strsplit(summ[4],"\\s+"))
DailyResult=data.frame(matrix(nrow = 0,ncol =length(outname)+2 ))
colnames(DailyResult)=paste0(outname)
summ[1:3]

Treatment1=summ[first[1]:994]

Treatment1[1]=gsub("#","_",Treatment1[1])

sumtext=read.table(textConnection(Treatment1[2:length(Treatment1)]),header = FALSE,
                   stringsAsFactors = FALSE)

name=read.table(textConnection(as.character(Treatment1[1])),header = FALSE,
                stringsAsFactors = FALSE)

colnames(sumtext)=name[2:83]

sumtext$Treatment=Treatmentname[1]
sumtext$TreatmentNO=1

DailyResult=rbind(DailyResult,sumtext)
head(DailyResult)
write.csv(DailyResult,"results.csv")
SUMMARY=read.csv("results.csv", header=T, na.strings=c("NA"))
FunctionValue=matrix(0,ncol = M,nrow = pop/15)
temp=SUMMARY$HWAM
for(io in 1:M){
  if(io==1){
    for (i in 1:(pop/15)) {
      for (j in 1:15) {
        te=(i-1)*15+j
        FunctionValue[i,io]=FunctionValue[i,io]+temp[te]
      }
      FunctionValue[i,io]=FunctionValue[i,io]/15
      FunctionValue[i,io]=10000-FunctionValue[i,io]
      
    }
  }else  if(io==2){
    for (i in 1:(pop/15)) {
      FunctionValue[i,io]=Populationi[i,1]+Populationi[i,2]+Populationi[i,3]
    }
    
  }else if(io==3){
    
    for (i in 1:(pop/15)) {
      FunctionValue[i,io]=  Populationf[i,1]+Populationf[i,2]+Populationf[i,3]
    }
  }else if(io==4){
    for (i in 1:(pop/15)) {
      for (j in 1:15) {
        te=(i-1)*15+j
        FunctionValue[i,io]=FunctionValue[i,io]+nle[te]
      }          
      FunctionValue[i,io]=FunctionValue[i,io]/15
      
    }
  }
  
}
setwd("E:\\NSGA-III\\NSGA-III\\NSGA-III")
writeMat("FUNV.mat",FunctionValue=FunctionValue)
