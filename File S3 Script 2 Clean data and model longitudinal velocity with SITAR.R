#SITAR
library(sitar)

purify<-function(x,y,id,d,limit=3,count=3){
  d$count<-ave(rep(1,nrow(d)),d$PID,FUN=sum)
  d<-d[d$count >= count,]
  veloutT<-paste0("velout(x=",x,",y=",y,",id=",id,",data=d,limit=",limit,",linearise = TRUE)")
  #velout evaluates in parent, so need to load d into parent
  d<<-d
  outliers<-eval(parse(text=veloutT))
  print(table(outliers$code))
  out46<-sum(outliers$code %in% c(4,6))>0
  
  while(out46) {
    d[,y][outliers$code %in% c(4,6)]<-NA
    d$count<-ave(rep(1,nrow(d)),d$PID,FUN=sum)
    d<-d[d$count >= count,]
    d<<-d
    outliers<-eval(parse(text=veloutT))
    print(table(outliers$code))
    out46<-sum(outliers$code %in% c(4,6))>0
  }
  
  d<<-d
  outliers <- outliers<-eval(parse(text=veloutT))
  print(table(outliers$code))
  return(d)
}


anth<-read.csv("File S1 Tsimane growth data with z-scores.csv",stringsAsFactors = FALSE)

#height
anthH<-anth[!is.na(anth$Height) & abs(anth$HAZ.T)<3,]
anthH<-anthH[order(anthH$PID,anthH$Age),]

#older female height
anthHFO<-anthH[anthH$Sex=="Female" & anthH$Age>5 & anthH$Age<22,]
anthHFO<-purify(x="Age",y="Height",id="PID",d=anthHFO,limit=3)
anthHFO<-anthHFO[!is.na(anthHFO$Height),]
outliers<-velout(Age,Height,PID,anthHFO,limit=3,linearise = TRUE)


plot(Height~I(Age^0.65),data=anthHFO[outliers$code==0,])
mplot(Age,Height,PID,data=anthHFO)

S.FO4<-sitar(x=Age^0.65,y=Height,PID,df=4,data=anthHFO,method="ML",control = nlmeControl(maxIter=50),verbose=TRUE)
plot(S.FO4, y2par=list(col='blue'), apv=TRUE)
# S.FO5 <- update(S.FO4, df=5)
# S.FO6 <- update(S.FO5, df=6)
# S.FO7 <- update(S.FO6, df=7)
# S.FO3 <- update(S.FO4, df=3)
# S.FO2 <- update(S.FO3, df=2)
# 
# varexp(S.FO2,S.FO3,S.FO4,S.FO5,S.FO6,S.FO7)
# BICadj(S.FO2,S.FO3,S.FO4,S.FO5,S.FO6,S.FO7)
#S.F04 appears to have the lowest BIC, though 3 is close
S.FO<-S.FO4 
      
#young female height
anthHFY<-anthH[anthH$Sex=="Female" & anthH$Age<5,]
anthHFY<-purify(x="Age",y="Height",id="PID",d=anthHFY,limit=3)
anthHFY<-anthHFY[!is.na(anthHFY$Height),]
plot(Height~I(Age^0.65),data=anthHFY)
mplot(Age^0.65,Height,PID,data=anthHFY)
S.FY4<-sitar(x=Age^0.75,y=Height,PID,df=4,random="a+c",xoffset=3,data=anthHFY,control = nlmeControl(maxIter=100),verbose=TRUE)
#S.FY3 <- update(S.FY4, df=3)
#S.FY5 <- update(S.FY4, df=5)
#S.FY6 <- update(S.FY5, df=6)
#S.FY7 <- update(S.FY6, df=7)
#S.FY3 <- update(S.FY4, df=4)
#S.FY2 <- update(S.FY3, df=2)
# varexp(S.FY3,S.FY4,S.FY5,S.FY6,S.FY7)
# BICadj(S.FY3,S.FY4,S.FY5,S.FY6,S.FY7)
plot(S.FY4, y2par=list(col='blue'))
S.FY<-S.FY4

#older male height
anthHMO<-anthH[anthH$Sex=="Male" & anthH$Age>5 & anthH$Age<25,]
anthHMO<-purify(x="Age",y="Height",id="PID",d=anthHMO,limit=3)
anthHMO<-anthHMO[!is.na(anthHMO$Height),]
outliers<-velout(Age,Height,PID,anthHMO,limit=3,linearise = TRUE)

plot(Height~I(Age^0.65),data=anthHMO[outliers$code==0,])
mplot(Age,Height,PID,data=anthHMO)
S.MO3<-sitar(x=Age^0.65,y=Height,PID,df=3,data=anthHMO,method="REML",control = nlmeControl(maxIter=100),verbose=TRUE)

#S.MO5 <- update(S.MO4, df=5)
#S.MO6 <- update(S.MO5, df=6)
#S.MO7 <- update(S.MO6, df=7)
#S.MO4 <- update(S.MO3, df=4)
#S.MO2 <- update(S.MO3, df=2) won't converge
# varexp(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7)
# BICadj(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7)
#S.MO3 seems to win
plot(S.MO3, y2par=list(col='blue'),apv=TRUE)
#save(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7,file="S.MO")
S.MO<-S.MO3

#young male height
anthHMY<-anthH[anthH$Sex=="Male" & anthH$Age<7,]
anthHMY<-purify(x="Age",y="Height",id="PID",d=anthHMY,limit=3)
anthHMY<-anthHMY[!is.na(anthHMY$Height),]
outliers<-velout(Age,Height,PID,anthHMY,limit=3,linearise = TRUE)
plot(Height~I(Age^0.6),data=anthHMY)
mplot(Age^0.65,Height,PID,data=anthHMY[1000:1500,])
S.MY4<-sitar(x=Age^0.6,y=Height,PID,random="a+c",xoffset=3,df=4,data=anthHMY,control = nlmeControl(maxIter=100),verbose=TRUE)
plot(S.MY4, y2par=list(col='blue'))
#S.MY3 <- update(S.MY4, df=3)
#S.MY5 <- update(S.MY4, df=5)
#S.MY6 <- update(S.MY5, df=6)
#S.MY7 <- update(S.MY6, df=7)
#S.MY2 <- update(S.MY3, df=2)
#varexp(S.MY3,S.MY4)
#BICadj(S.MY3,S.MY4)
#save(S.MY2,S.MY3,S.MY4,S.MY5,S.MY6,file="S.FM")
#plot(S.MY4, y2par=list(col='blue'))
#4df lowest BIC
S.MY<-S.MY4

save(S.FO,S.FY,S.MO,S.MY,anthHFY,anthHFO,anthHMY,anthHMO,file="sitarModelsH")


#weight
anthW<-anth[!is.na(anth$Weight) & abs(anth$WAZ.T)<3,]
anthW<-anthW[order(anthW$PID,anthW$Age),]

#older female Weight
anthWFO<-anthW[anthW$Sex=="Female" & anthW$Age>5 & anthW$Age<22,]
anthWFO<-purify(x="Age",y="Weight",id="PID",d=anthWFO,limit=3)
anthWFO<-anthWFO[!is.na(anthWFO$Weight),]
outliers<-velout(Age,Weight,PID,anthWFO,limit=3,linearise = TRUE)


plot(Weight~I(Age^0.65),data=anthWFO[outliers$code==0,])
mplot(Age,Weight,PID,data=anthWFO)

SW.FO4<-sitar(x=Age^0.65,y=Weight,PID,df=4,data=anthWFO,method="ML",control = nlmeControl(maxIter=50),verbose=TRUE)
plot(SW.FO4, y2par=list(col='blue'), apv=TRUE)
# S.FO5 <- update(SW.FO4, df=5)
# S.FO6 <- update(S.FO5, df=6)
# S.FO7 <- update(S.FO6, df=7)
# S.FO3 <- update(SW.FO4, df=3)
# S.FO2 <- update(S.FO3, df=2)
# 
# varexp(S.FO4,S.FO5)
# BICadj(S.FO4,S.FO5)
#only 4 and 5 will converge. 4 is best
SW.FO<-SW.FO4 

#young female Weight
anthWFY<-anthW[anthW$Sex=="Female" & anthW$Age<7,]
anthWFY<-purify(x="Age",y="Weight",id="PID",d=anthWFY,limit=3)
anthWFY<-anthWFY[!is.na(anthWFY$Weight),]
plot(Weight~I(Age^0.75),data=anthWFY)
mplot(Age,Weight,PID,data=anthWFY)
SW.FY3<-sitar(x=Age^0.65,y=Weight,PID,random="a+c",xoffset=3,df=3,data=anthWFY,control = nlmeControl(maxIter=50),verbose=TRUE)
plot(SW.FY3, y2par=list(col='blue'))
# SW.FY4 <- update(SW.FY3, df=4)
# SW.FY5 <- update(SW.FY4, df=5)
# SW.FY6 <- update(SW.FY5, df=6)
# SW.FY7 <- update(SW.FY6, df=7)
# SW.FY2 <- update(SW.FY3, df=2)
 #2,3,4, converge
# varexp(S.FY5,S.FY3,S.FY4)
#BICadj(SW.FY5,SW.FY3,SW.FY4)
SW.FY<-SW.FY3

#older male Weight
anthWMO<-anthW[anthW$Sex=="Male" & anthW$Age>5 & anthW$Age<25,]
anthWMO<-purify(x="Age",y="Weight",id="PID",d=anthWMO,limit=3)
anthWMO<-anthWMO[!is.na(anthWMO$Weight),]
outliers<-velout(Age,Weight,PID,anthWMO,limit=3,linearise = TRUE)

plot(Weight~I(Age^0.65),data=anthWMO[outliers$code==0,])
mplot(Age,Weight,PID,data=anthWMO)
SW.MO3<-sitar(x=Age^0.65,y=Weight,PID,df=3,data=anthWMO,control = nlmeControl(maxIter=100),verbose=TRUE)
plot(SW.MO3, y2par=list(col='blue'),apv=TRUE)
# S.MO4 <- update(SW.MO3, df=4)
# S.MO5 <- update(S.MO4, df=5)
# S.MO6 <- update(S.MO5, df=6)
# S.MO7 <- update(S.MO6, df=7)
# 
# #S.MO2 <- update(S.MO3, df=2)
# 
# varexp(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7)
# BICadj(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7)
#S.MO3 seems to win
SW.MO<-SW.MO3

#young male Weight
anthWMY<-anthW[anthW$Sex=="Male" & anthW$Age<7,]
anthWMY<-purify(x="Age",y="Weight",id="PID",d=anthWMY,limit=3)
anthWMY<-anthWMY[!is.na(anthWMY$Weight),]
outliers<-velout(Age,Weight,PID,anthWMY,limit=3,linearise = TRUE)
plot(Weight~I(Age^0.75),data=anthWMY)
mplot(Age^0.75,Weight,PID,data=anthWMY[1000:1500,])
SW.MY4<-sitar(x=Age^0.65,y=Weight,PID,random="a+c",xoffset=2,df=4,data=anthWMY,control = nlmeControl(maxIter=100),verbose=TRUE)
plot(SW.MY4, y2par=list(col='blue'))
# SW.MY3 <- update(SW.MY4, df=3)
# SW.MY2 <- update(SW.MY3, df=2)
# SW.MY5 <- update(SW.MY4, df=5)
# S.MY6 <- update(S.MY5, df=6)
# S.MY7 <- update(S.MY6, df=7)
# 
# 
# varexp(S.MY2,S.MY3)
#BICadj(SW.MY4,SW.MY5)
#6df lowest BIC
SW.MY<-SW.MY4

save(SW.FO,SW.FY,SW.MO,SW.MY,anthWFY,anthWFO,anthWMY,anthWMO,file="sitarModelsW")



#recombine height and wight data frames to create a merged data frame for BMI 
anthC<-merge(anthH[,c("PID","Sex","Age","Height")],anthW[,c("PID","Sex","Age","Weight")],all=TRUE)

#BMI
anthC$BMI<-anthC$Weight/(anthC$Height/100)^2

anthB<-anthC[!is.na(anthC$BMI),]
anthB<-anthB[order(anthB$PID,anthB$Age),]

#older female BMI
anthBFO<-anthB[anthB$Sex=="Female" & anthB$Age>5 & anthB$Age<22,]
anthBFO<-purify(x="Age",y="BMI",id="PID",d=anthBFO,limit=3)
anthBFO<-anthBFO[!is.na(anthBFO$BMI),]
outliers<-velout(Age,BMI,PID,anthBFO,limit=3,linearise = TRUE)


plot(BMI~I(Age^0.65),data=anthBFO)
mplot(Age^0.4,BMI,PID,data=anthBFO)

SB.FO3<-sitar(x=Age^0.65,y=BMI,PID,df=3,data=anthBFO,control = nlmeControl(maxIter=50),verbose=TRUE)
plot(SB.FO3, y2par=list(col='blue'), apv=TRUE)
# SB.FO4 <- update(SB.FO3, df=4)
# S.FO5 <- update(S.FO4, df=5)
# S.FO6 <- update(S.FO5, df=6)
# S.FO7 <- update(S.FO6, df=7)
# S.FO2 <- update(S.FO3, df=2)

# varexp(S.FO4,S.FO5)
# BICadj(S.FO4,S.FO5)
#3 only one that will converge
SB.FO<-SB.FO3

#young female BMI
anthBFY<-anthB[anthB$Sex=="Female" & anthB$Age<7,]
anthBFY<-purify(x="Age",y="BMI",id="PID",d=anthBFY,limit=3)
anthBFY<-anthBFY[!is.na(anthBFY$BMI),]
plot(BMI~I(Age^0.9),data=anthBFY)
mplot(Age^0.9,BMI,PID,data=anthBFY)
SB.FY4<-sitar(x=Age^0.9,y=BMI,PID,random="a+c",xoffset=2,df=4,data=anthBFY,control = nlmeControl(maxIter=50),verbose=TRUE)
plot(SB.FY4, y2par=list(col='blue'))
# SB.FY5 <- update(SB.FY4, df=5)
# SB.FY6 <- update(SB.FY5, df=6)
# SB.FY7 <- update(SB.FY6, df=7)
# SB.FY3 <- update(SB.FY4, df=4)
# SB.FY2 <- update(SB.FY3, df=2)
# varexp(SB.FY3,SB.FY4,SB.FY5)
# BICadj(SB.FY3,SB.FY4,SB.FY5)
#only 4 converges
SB.FY<-SB.FY4

#older male BMI
anthBMO<-anthB[anthB$Sex=="Male" & anthB$Age>5 & anthB$Age<25,]
anthBMO<-purify(x="Age",y="BMI",id="PID",d=anthBMO,limit=3,count=3)
anthBMO<-anthBMO[!is.na(anthBMO$BMI),]
outliers<-velout(Age,BMI,PID,anthBMO,limit=3,linearise = TRUE)

plot(BMI~I(Age^0.65),data=anthBMO[outliers$code==0,])
mplot(Age,BMI,PID,data=anthBMO)
SB.MO4<-sitar(x=Age^0.65,y=BMI,PID,df=4,data=anthBMO,control = nlmeControl(maxIter=100,minScale=0.000001),verbose=TRUE)
#SB.MO6<-sitar(x=Age^0.65,y=BMI,PID,df=6,data=anthBMO,control = nlmeControl(maxIter=100),verbose=TRUE)
plot(SB.MO4, y2par=list(col='blue'),apv=TRUE)
# 
# SB.MO5 <- update(SB.MO4, df=5)
# SB.MO6 <- update(SB.MO5, df=6)
# SB.MO7 <- update(SB.MO6, df=7)
# SB.MO3 <- update(SB.MO4, df=3)
# S.MO2 <- update(S.MO3, df=2) won't converge
# 
# varexp(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7)
# BICadj(SB.MO4,SB.MO6)
#4 only one to converge
#plot(S.MO4, y2par=list(col='blue'),apv=TRUE)
#save(S.MO3,S.MO4,S.MO5,S.MO6,S.MO7,file="S.MO")
SB.MO<-SB.MO4

#young male BMI
anthBMY<-anthB[anthB$Sex=="Male" & anthB$Age<7,]
anthBMY<-purify(x="Age",y="BMI",id="PID",d=anthBMY,limit=3)
anthBMY<-anthBMY[!is.na(anthBMY$BMI),]
outliers<-velout(Age,BMI,PID,anthBMY,limit=3,linearise = TRUE)
plot(BMI~I(Age^0.9),data=anthBMY)
mplot(Age^0.9,BMI,PID,data=anthBMY)
SB.MY4<-sitar(x=Age^0.9,y=BMI,PID,random="a+c",xoffset=0,df=4,data=anthBMY,control = nlmeControl(maxIter=100),verbose=TRUE)
plot(SB.MY4, y2par=list(col='blue'))
# SB.MY5 <- update(SB.MY4, df=5)
# SB.MY6 <- update(SB.MY5, df=6)
# SB.MY7 <- update(SB.MY6, df=7)
# SB.MY3 <- update(SB.MY4, df=3)
# SB.MY2 <- update(SB.MY3, df=2)
# varexp(S.MY4,S.MY5)
# BICadj(SB.MY4,SB.MY5,SB.MY6,SB.MY7)
# only 4 again
# plot(S.MY3, y2par=list(col='blue'))
#3df lowest BIC
SB.MY<-SB.MY4

save(SB.FO,SB.FY,SB.MO,SB.MY,anthBMY,anthBMO,anthBFY,anthBFO,file="sitarModelsB")
