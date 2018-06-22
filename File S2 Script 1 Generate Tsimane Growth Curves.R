library(gamlss)
library(foreach)
library(doParallel)


fulldata<-read.csv("File S1 Tsimane growth data with z-scores.csv",stringsAsFactors = FALSE)


# R code for working with LMS curves



#This function produces centile values for given ages and LMS values.
centilesLMS<-function(xvar,LMS,cent=c(0.02,0.05,0.25,0.50,0.75,0.95,0.98),xname="Age"){
  getcent<-function(la,mu,si){
    sapply(cent,function(x) exp(log(qnorm(x)*la*si+1)/la+log(mu)))
  }
  o<-cbind(xvar,t(apply(LMS,1,function(x) getcent(x[1],x[2],x[3]))))
  o<-data.frame(o)
  names(o)<-c(xname,cent)
  o
}


#This function calculates a z-score from a measure and the appropriate LMS values. 
ZfromLMS<-function(measure,la,mu,si){
  zs <- ((measure / mu)^la - 1)/(la * si)
  zs
}


#function to get z-scores for a given age and value. Sex, age, and value need to be in the same units or with the same coding scheme as is used in the LMS table. Will fail with NAs so remove them first
getZ<-function(Age, Sex, Value, LMSlookuptable){
  useLMS<-apply(data.frame(Age,Sex), 1, function(x) which.min(abs(LMSlookuptable$Age-as.numeric(x[1]))+(x[2]!=LMSlookuptable$Sex)*999))
  LMSmatches<-LMSlookuptable[useLMS,]
  ZfromLMS(Value,LMSmatches$Lambda,LMSmatches$Mu,LMSmatches$Sigma)
}

getZ.WFH<-function(Height, Sex, Value, LMSlookuptable){
  useLMS<-apply(data.frame(Height,Sex), 1, function(x) which.min(abs(LMSlookuptable$Height-as.numeric(x[1]))+(x[2]!=LMSlookuptable$Sex)*999))
  LMSmatches<-LMSlookuptable[useLMS,]
  ZfromLMS(Value,LMSmatches$Lambda,LMSmatches$Mu,LMSmatches$Sigma)
}

#just a wrapper. Uses pnorm to get centiles from the Z-scores in getZ
getCentile<-function(Age, Sex, Value, LMSlookuptable){
  zs<-getZ(Age, Sex, Value, LMSlookuptable)
  cent<-pnorm(zs)
  cent
}

getCentile.WFH<-function(Height, Sex, Value, LMSlookuptable){
  zs<-getZ.WFH(Height, Sex, Value, LMSlookuptable)
  cent<-pnorm(zs)
  cent
}

#Extracts LMS values from a GAMLSS model
getLMS<-function(M1,newdata = data.frame(xvar=seq(0,25,0.1))){
	m<-predict(M1, what = "mu", newdata = newdata)
	s<-predict(M1, what = "sigma", newdata =newdata)
	n<-try(predict(M1, what = "nu", newdata =newdata))
	if (inherits(n, "try-error")) n<-rep(1,nrow(newdata))
	LMS<-data.frame(newdata[,1],Mu=m,Sigma=s,Lambda=n)
	names(LMS)[1]<-names(newdata)[1]
	LMS$Sigma<-exp(LMS$Sigma)
	LMS
}


#Function to run through GAMLSS models to determine appropriate parameter values
getGAMLSSparams<-function(usedata,filen,yvar,xvar="Age",forcel=NA,lrange=seq(0.3,1.0,0.05),cutoff=3){
	SMin<-c(xvar,"Sex",yvar,"PID")
	usedata<-na.omit(usedata[,SMin])
	usedata$yvar<-usedata[,yvar]
	usedata$xvar<-usedata[,xvar]
	usedata<<-usedata
	
	tiff(paste("getGAMLSSparams-",filen,".tif",sep=""),width=234,height=234,units="mm",res=800,pointsize=9,compression="lzw")
	par(mfrow=c(2,2))
	
	#out<-matrix(ncol=2)
	out<-foreach (l=lrange, .packages=c("gamlss"), .combine=rbind) %dopar% {
		MT<-try(gamlss(yvar~pb(I(xvar^l),df=14), sigma.formula=~pb(I(xvar^l),df=3), nu.fo=~1, tau.formula=~1, nu.fix=TRUE, nu.start=1, tau.fix=TRUE, tau.start=2,family=BCPE, data=usedata, control=gamlss.control(n.cyc=50)))
		if (inherits(MT, "try-error")) outp<-c(l,NA) else outp<-c(l,MT$G.deviance)
		outp
	}
	usel<-out[which.min(out[,2]),1]
	plot(out,main=usel)
	if(!is.na(forcel)) usel<-forcel
	#fit low df for removing outliers
	MT<<-gamlss(yvar~pb(I(xvar^usel),df=14), sigma.formula=~pb(I(xvar^usel),df=4), nu.fo=~1, tau.formula=~1, nu.fix=TRUE, nu.start=1, tau.fix=TRUE, tau.start=2, family=BCPE, data=usedata, control=gamlss.control(n.cyc=50))
	centiles(MT,usedata$xvar^usel,cent=c(pnorm(cutoff * -1),pnorm(cutoff))*100)
	Zs<<-centiles.pred(MT,xname="xvar", xvalues=usedata$xvar, yval=usedata$yvar, type="z-scores")
	
	usedata<<-usedata
	
	excl<-paste("Excluded",nrow(usedata[abs(Zs)> cutoff,]),"/",nrow(usedata))
	print(excl)
	#exit
	#remove outliers as determined by either of the two models above
	usedata<-usedata[!(abs(Zs)>cutoff),]
	usedata<<-usedata
  print("2")
  #determine mu and sigma
	dfs<-expand.grid(4:22,4:20)
	out2<-matrix(ncol=7)
	out2<-foreach (i=1:nrow(dfs), .packages=c("gamlss"), .combine=rbind) %dopar% {
		MT<-try(gamlss(yvar~pb(I(xvar^usel),df=dfs[i,1]), sigma.formula=~pb(I(xvar^usel),df=dfs[i,2]), nu.fo=~1, tau.formula=~1, nu.fix=TRUE, nu.start=1, tau.fix=TRUE, tau.start=2, family=BCPE, data=usedata, control=gamlss.control(n.cyc=50)))
		if (inherits(MT, "try-error")) outp<-c(dfs[i,1],dfs[i,2],NA,NA,NA,NA,NA) else outp<-c(dfs[i,1],dfs[i,2],AIC(MT,k=2),GAIC(MT,k=3),GAIC(MT,k=3.5),GAIC(MT,k=4),AIC(MT,k=2,c=TRUE))
		outp
	}
	usedfs<-out2[which.min(out2[,4]),c(1,2)]
	plot(0,0,type="n",xlim=c(1,22),ylim=c(min(out2[,3],na.rm=TRUE),max(out2[,4],na.rm=TRUE)),main=paste(usedfs[1],usedfs[2],sep=","))
	for (i in 4:20) lines(out2[out2[,2]==i,1],out2[out2[,2]==i,4],col=rev(heat.colors(16))[i-1])
print("3")	
	#next model nu
	out3<-matrix(ncol=7)
	#check whether nu should be fixed
	MT2<<-try(gamlss(yvar~pb(I(xvar^usel),df=usedfs[1]), sigma.formula=~pb(I(xvar^usel),df=usedfs[2]), nu.fo=~pb(I(xvar^usel),df=1),nu.start=1,nu.fix=TRUE, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=usedata, control=gamlss.control(n.cyc=50)))
	if (inherits(MT2, "try-error")) out3<-rbind(out3,c(0,NA,NA,NA,NA,NA,1)) else out3<-rbind(out3,c(1,AIC(MT2,k=2),GAIC(MT2,k=3),GAIC(MT2,k=3.5),GAIC(MT2,k=4),AIC(MT2,k=2,c=TRUE),1))
	out3p<-foreach (d=0:14, .packages=c("gamlss"), .combine=rbind) %dopar% {
		MT<-try(gamlss(yvar~pb(I(xvar^usel),df=usedfs[1]), sigma.formula=~pb(I(xvar^usel),df=usedfs[2]), nu.fo=~pb(I(xvar^usel),df=d), tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=usedata, control=gamlss.control(n.cyc=50)))
		if (inherits(MT, "try-error")) outp<-c(d,NA,NA,NA,NA,NA,NA) else outp<-c(d,AIC(MT,k=2),GAIC(MT,k=3),GAIC(MT,k=3.5),GAIC(MT,k=4),AIC(MT,k=2,c=TRUE),0)
		outp
	}
	out3<-rbind(out3,out3p)
	usedfs<-c(usedfs,out3[which.min(out3[,3]),1])
	plot(out3[,1],out3[,3],main=usedfs[3])
	nufix<-out3[which.min(out3[,3]),7]==1
	taufix<-TRUE
	dev.off()
	return(list(usedata=usedata,usel=usel,usedfs=usedfs,nufix=nufix,taufix=taufix,out=out,out2=out2,out3=out3,exclude=excl))	
}

#remove anyone with age==0, since we don't actually have birth weights, so these are rounded ages and seem to mess with the curves.

#get parameters
#adjust this cluster parameter to 1 less than the number of processors in computer
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

H.F3<-getGAMLSSparams(fulldata[fulldata$Sex=="Female" & fulldata$Age>0,],"H.F3","Height",cutoff=3.5)
#new Excluded 58 / 15612
H.M3<-getGAMLSSparams(fulldata[fulldata$Sex=="Male" & fulldata$Age>0,],"H.M3","Height",cutoff=3.5)
#new Excluded 43 / 14561"
W.F3<-getGAMLSSparams(fulldata[fulldata$Sex=="Female" & fulldata$Age>0,],"W.F3","Weight",cutoff=3.5)
#new Excluded 56 / 15640"
W.M3<-getGAMLSSparams(fulldata[fulldata$Sex=="Male" & fulldata$Age>0,],"W.M3","Weight",cutoff=3.5)
#new Excluded 41 / 14578

#WFH also has some big edge effects due to small sample at low end, so need to add some anchor points below lowest observation
add<-data.frame(PID=rep("Fake",100),Age=rep(5,100),Height=rep(38,100),Weight=rnorm(100,2.5,0.5),Sex=c(rep("Female",50),rep("Male",50)))
fulldataAnchor<-merge(fulldata,add,all=TRUE)

WH.F3<-getGAMLSSparams(fulldataAnchor[fulldataAnchor$Sex=="Female" & fulldataAnchor$Age>0 & fulldataAnchor$Age<=20,],"WH.FA","Weight","Height",lrange=seq(0.5,3,0.1),cutoff=3.5)
#Excluded 108 / 15501
#plot(Weight~I(Height^1),data=usedata,xlim=c(30,80)^1,ylim=c(0,15))
WH.M3<-getGAMLSSparams(fulldataAnchor[fulldataAnchor$Sex=="Male" & fulldataAnchor$Age>0 & fulldataAnchor$Age<=20,],"WH.MA","Weight","Height",lrange=seq(0.5,3,0.1),cutoff=3.5)
#Excluded 127 / 14429
fulldata$BMI<-fulldata$Weight/(fulldata$Height/100)^2


B.F3<-getGAMLSSparams(fulldata[fulldata$Sex=="Female" & fulldata$Age>0 & fulldata$Age<30,],"B.F3","BMI",cutoff=3.5,forcel=0.4)
#Excluded 116 / 15500, 205
B.M3<-getGAMLSSparams(fulldata[fulldata$Sex=="Male" & fulldata$Age>0 & fulldata$Age<30,],"B.M3","BMI",cutoff=3.5,forcel=0.4)
#Excluded 132 / 14427, 207


stopImplicitCluster()


params<-rbind(c(length(unique(H.M3$usedata$PID)),nrow(H.M3$usedata),H.M3$usel,H.M3$usedfs, H.M3$nufix),
              c(length(unique(H.F3$usedata$PID)),nrow(H.F3$usedata),H.F3$usel,H.F3$usedfs, H.F3$nufix),
              c(length(unique(W.M3$usedata$PID)),nrow(W.M3$usedata),W.M3$usel,W.M3$usedfs, W.M3$nufix),
              c(length(unique(W.F3$usedata$PID)),nrow(W.F3$usedata),W.F3$usel,W.F3$usedfs, W.F3$nufix),
              c(length(unique(B.M3$usedata$PID)),nrow(B.M3$usedata),B.M3$usel,B.M3$usedfs, B.M3$nufix),
              c(length(unique(B.F3$usedata$PID)),nrow(B.F3$usedata),B.F3$usel,B.F3$usedfs, B.F3$nufix),
              c(length(unique(WH.M3$usedata$PID[WH.M3$usedata$PID!="Fake"])),nrow(WH.M3$usedata[WH.M3$usedata$PID!="Fake",]),WH.M3$usel,WH.M3$usedfs, WH.M3$nufix),
              c(length(unique(WH.F3$usedata$PID[WH.F3$usedata$PID!="Fake"])),nrow(WH.F3$usedata[WH.F3$usedata$PID!="Fake",]),WH.F3$usel,WH.F3$usedfs, WH.F3$nufix))

#running the paramater search takes a long time, so useful to save objects
save(list = c("H.M3", "H.F3","W.F3", "W.M3","B.M3", "B.F3","WH.M3","WH.F3"), file="TsimaneZscore.Rdata")
#can reload and start from here if doing a second time
#load(file="TsimaneZscore.Rdata")

#female height
gH.F3<-gamlss(Height~pb(I((Age)^H.F3$usel),df=H.F3$usedfs[1]), sigma.formula=~pb(I((Age)^H.F3$usel),df=H.F3$usedfs[2]), nu.fo=~pb(I((Age)^H.F3$usel),df=H.F3$usedfs[3]), nu.fix=H.F3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=H.F3$usedata, control=gamlss.control(n.cyc=50))

centiles(gH.F3,(H.F3$usedata$Age)^0.6,legend=F)
centiles(gH.F3,H.F3$usedata$Age,legend=F)
centiles(gH.F3,H.F3$usedata$Age,xlim=c(0,5),legend=F)
LMS.HF2<-getLMS(gH.F3,data.frame(Age=seq(0,25,0.1)))
LMS.HF3<-getLMS(gH.F3,data.frame(Age=seq(0,25,0.01)))
C.HF2<-centilesLMS(LMS.HF2[,1],LMS.HF2[,c(4,2,3)])

#male height
gH.M3<-gamlss(Height~pb(I(Age^H.M3$usel),df=H.M3$usedfs[1]), sigma.formula=~pb(I(Age^H.M3$usel),df=H.M3$usedfs[2]), nu.fo=~pb(I(Age^H.M3$usel),df=H.M3$usedfs[3]), nu.fix=H.M3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=H.M3$usedata, control=gamlss.control(n.cyc=50))
centiles(gH.F3,(H.F3$usedata$Age)^0.6,legend=F)
centiles(gH.M3,H.M3$usedata$Age,legend=F)
centiles(gH.M3,H.M3$usedata$Age,xlim=c(0,5),legend=F)
LMS.HM2<-getLMS(gH.M3,data.frame(Age=seq(0,25,0.1)))
LMS.HM3<-getLMS(gH.M3,data.frame(Age=seq(0,25,0.01)))
C.HM2<-centilesLMS(LMS.HM2[,1],LMS.HM2[,c(4,2,3)])

#female weight
gW.F3<-gamlss(Weight~pb(I(Age^W.F3$usel),df=W.F3$usedfs[1]), sigma.formula=~pb(I(Age^W.F3$usel),df=W.F3$usedfs[2]), nu.fo=~pb(I(Age^W.F3$usel),df=W.F3$usedfs[3]), nu.fix=W.F3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=W.F3$usedata, control=gamlss.control(n.cyc=50))
centiles(gW.F3,W.F3$usedata$Age,legend=F)
centiles(gW.F3,W.F3$usedata$Age,xlim=c(0,5),legend=F,ylim=c(0,40))
LMS.WF2<-getLMS(gW.F3,data.frame(Age=seq(0,25,0.1)))
LMS.WF3<-getLMS(gW.F3,data.frame(Age=seq(0,25,0.01)))
C.WF2<-centilesLMS(LMS.WF2[,1],LMS.WF2[,c(4,2,3)])

#male weight
#the male curves under age 0.1 have an edge effect issue leading to a skewed distribution, likely due to small age errors for some male infants. To minimize this, reducing the lambda degrees of freedom to match the female curve
gW.M3<-gamlss(Weight~pb(I(Age^W.M3$usel),df=W.M3$usedfs[1]), sigma.formula=~pb(I(Age^W.M3$usel),df=W.M3$usedfs[2]), nu.fo=~pb(I(Age^W.M3$usel),df=W.M3$usedfs[3]-2), nu.fix=W.M3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=W.M3$usedata, control=gamlss.control(n.cyc=50))
centiles(gW.M3,W.M3$usedata$Age,legend=F)
centiles(gW.M3,W.M3$usedata$Age,xlim=c(0,5),legend=F,ylim=c(0,40))
LMS.WM2<-getLMS(gW.M3,data.frame(Age=seq(0,25,0.1)))
LMS.WM3<-getLMS(gW.M3,data.frame(Age=seq(0,25,0.01)))
C.WM2<-centilesLMS(LMS.WM2[,1],LMS.WM2[,c(4,2,3)])
head(C.WM2)

#female BMI
#0.55 original
gB.F3<-gamlss(BMI~pb(I(Age^B.F3$usel),df=B.F3$usedfs[1]), sigma.formula=~pb(I(Age^B.F3$usel),df=B.F3$usedfs[2]), nu.fo=~pb(I(Age^B.F3$usel),df=B.F3$usedfs[3]), nu.fix=B.F3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=B.F3$usedata, control=gamlss.control(n.cyc=50))
centiles(gB.F3,B.F3$usedata$Age,legend=F,xlim=c(0,2),ylim=c(5,25))
centiles(gB.F3,B.F3$usedata$Age^0.4,xlim=c(0,5),legend=F)
LMS.BF2<-getLMS(gB.F3,data.frame(Age=seq(0,25,0.1)))
LMS.BF3<-getLMS(gB.F3,data.frame(Age=seq(0,25,0.01)))
C.BF2<-centilesLMS(LMS.BF2[,1],LMS.BF2[,c(4,2,3)])

#male BMI
gB.M3<-gamlss(BMI~pb(I(Age^B.M3$usel),df=B.M3$usedfs[1]), sigma.formula=~pb(I(Age^B.M3$usel),df=B.M3$usedfs[2]), nu.fo=~pb(I(Age^B.M3$usel),df=B.M3$usedfs[3]), nu.fix=B.M3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=B.M3$usedata, control=gamlss.control(n.cyc=50))
centiles(gB.M3,B.M3$usedata$Age,legend=F,xlim=c(0,2),ylim=c(5,25))
centiles(gB.M3,B.M3$usedata$Age^0.4,xlim=c(0,5),legend=F)
LMS.BM2<-getLMS(gB.M3,data.frame(Age=seq(0,25,0.1)))
LMS.BM3<-getLMS(gB.M3,data.frame(Age=seq(0,25,0.01)))
C.BM2<-centilesLMS(LMS.BM2[,1],LMS.BM2[,c(4,2,3)])

#female WFH
gWH.F3<-gamlss(Weight~pb(I(Height^WH.F3$usel),df=WH.F3$usedfs[1]), sigma.formula=~pb(I(Height^WH.F3$usel),df=WH.F3$usedfs[2]), nu.fo=~pb(I(Height^WH.F3$usel),df=WH.F3$usedfs[3]), nu.fix=WH.F3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=WH.F3$usedata, control=gamlss.control(n.cyc=50))
centiles(gWH.F3,WH.F3$usedata$Height,legend=F)
centiles(gWH.F3,WH.F3$usedata$Height^1.5,legend=F)
LMS.WHF2<-getLMS(gWH.F3,data.frame(Height=seq(45,160,0.1)))
LMS.WHF3<-getLMS(gWH.F3,data.frame(Height=seq(45,160,0.01)))
C.WHF2<-centilesLMS(LMS.WHF2[,1],LMS.WHF2[,c(4,2,3)],xname="Height")

#male WFH
gWH.M3<-gamlss(Weight~pb(I(Height^WH.M3$usel),df=WH.M3$usedfs[1]), sigma.formula=~pb(I(Height^WH.M3$usel),df=WH.M3$usedfs[2]), nu.fo=~pb(I(Height^WH.M3$usel),df=WH.M3$usedfs[3]), nu.fix=WH.M3$nufix, nu.start=1, tau.formula=~1, tau.fix=TRUE, tau.start=2, family=BCPE, data=WH.M3$usedata, control=gamlss.control(n.cyc=50))
centiles(gWH.M3,WH.M3$usedata$Height,legend=F)
centiles(gWH.M3,WH.M3$usedata$Height^1.5,legend=F)
LMS.WHM2<-getLMS(gWH.M3,data.frame(Height=seq(45,170,0.1)))
LMS.WHM3<-getLMS(gWH.M3,data.frame(Height=seq(45,170,0.01)))
C.WHM2<-centilesLMS(LMS.WHM2[,1],LMS.WHM2[,c(4,2,3)],xname="Height")

LMS.HF3$Sex<-"Female"
LMS.HM3$Sex<-"Male"
LMS.WF3$Sex<-"Female"
LMS.WM3$Sex<-"Male"
LMS.BF3$Sex<-"Female"
LMS.BM3$Sex<-"Male"
LMS.WHF3$Sex<-"Female"
LMS.WHM3$Sex<-"Male"

HT.LMS<-rbind(LMS.HM3,LMS.HF3)
WT.LMS<-rbind(LMS.WM3,LMS.WF3)
BMI.LMS<-rbind(LMS.BM3,LMS.BF3)
WFH.LMS<-rbind(LMS.WHM3,LMS.WHF3)

write.csv(HT.LMS,"File S4 Tsimane Height By Age LMS.csv")
write.csv(WT.LMS,"File S5 Tsimane Weight By Age LMS.csv")
write.csv(BMI.LMS,"File S6 Tsimane BMI by Age LMS.csv")
write.csv(WFH.LMS,"File S7 Tsimane Weight by Height LMS.csv")

write.csv(C.HF2,"File S8 Tsimane Girls Height Centiles.csv")
write.csv(C.HM2,"File S9 Tsimane Boys Height Centiles.csv")
write.csv(C.WF2,"File S10 Tsimane Girls Weight Centiles.csv")
write.csv(C.WM2,"File S11 Tsimane Boys Weight Centiles.csv")
write.csv(C.BF2,"File S12 Tsimane Girls BMI Centiles.csv")
write.csv(C.BM2,"File S13 Tsimane Boys BMI Centiles.csv")
write.csv(C.WHF2,"File S14 Tsimane Girls WFH Centiles.csv")
write.csv(C.WHM2,"File S15 Tsimane Boys WFH Centiles.csv")

#add z-scores into data file
fulldata$HAZ.T<-getZ(fulldata$Age, fulldata$Sex, fulldata$Height, HT.LMS)
fulldata$WAZ.T<-getZ(fulldata$Age, fulldata$Sex, fulldata$Weight, WT.LMS)
fulldata$BAZ.T<-getZ(fulldata$Age, fulldata$Sex, fulldata$BMI, BMI.LMS)
fulldata$WFHZ.T<-NA
fulldata$WFHZ.T[!is.na(fulldata$Height)]<-getZ.WFH(fulldata$Height[!is.na(fulldata$Height)], fulldata$Sex[!is.na(fulldata$Height)], fulldata$Weight[!is.na(fulldata$Height)], WFH.LMS)
#write.csv(fulldata,"File S1 Tsimane growth data with z-scores.csv",row.names=FALSE)

