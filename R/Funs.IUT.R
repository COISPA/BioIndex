#' @importFrom stats predict AIC nls
#' @importFrom stats pnorm cor.test rnorm dchisq deviance df.residual na.exclude sd time
#' @importFrom mgcv gam
.ChangeYears.func<-function(Year, Indic ,Sd=NULL,alpha=0.05,species,B=199,ytext,CV=0.01,lastn=5,write=FALSE,speciesname,lastligne=FALSE) {
    oldpar <- par(no.readonly = TRUE)
    oldpar$new <- NULL
    on.exit(suppressWarnings(par(oldpar)))

  if (FALSE) {
    Year=years
    Indic=indi
    Sd=rep(NA,length(indi))
    alpha=alpha
    species=system
    B=199
    ytext=indicators[j]
    CV=CVindicators
    lastn=lastn
    write=FALSE
    speciesname=indiname[j]
    lastligne=FALSE

    }

    utils::globalVariables(c("lab_y", "AGE", "LENGTH", "LENGTH_CLASS", "med_bathy", "strata"))


#searches for change in second derivative of index
# Year input data = year vector
# Indic input data = indicator time series estimates as vector
# Sd input data = standard deviation estimates as vector
# species = species name
# if no Sd exists, a CV of 1% is assumed; needed for incorporating uncertainty
# B number of samples for parametric bootstrap (normal distribution)
# ytext : label of indicator name for plot
# lastn : number of recent years (time horizon) for diagnosis of change



####AUXILIARY FUNCTIONS###########
.Ggam.func<-function(G,annee,first=2,df=1){
 #function for fitting Gam to G series
 #first=1 returns df for fit to raw data and p-value for chisquare test of fit
 #first=2 returns fitted values
 if(first==1){
      if(length(annee)>10) fit<-gam(G~s(annee,bs="tp"))
      else fit<-gam(G~s(annee,bs="tp",k=length(annee)-1))
      pchisq<-dchisq(deviance(fit),df.residual(fit))
 }
 else fit<-gam(G~s(annee,k=df+1,fx=TRUE,bs="tp"))
 ifelse(first==1,return(list(df=summary(fit)$edf,pchisq=pchisq)),return(fit$fitted))
}

.Deriv.func<-function(G,annee,df)
{
#uses local difference between fitted function to estimate gradients
#returns first derivative for years with observations
#and  second derivative
 delta <- 1e-7
 data<-data.frame(G=G,annee=annee)
 fit<-gam(G~s(annee,k=df+1,fx=TRUE,bs="tp"),data=data)

 #predict years-delta and years on scale of linear predictor
 data.pred<-data.frame(annee=c(annee-delta,annee))
 Xp1<-predict.gam(fit,data.pred,type="lpmatrix")

 #predict years and years+delta on scale of linear predictor
 data.pred<-data.frame(annee=c(annee,annee+delta))
 Xp2<-predict.gam(fit,data.pred,type="lpmatrix")

 #calculate finite difference Xp(annee+delta) - Xp(annee) of linear predictor
#
 ny<-length(annee)
 X1<-(Xp1[(ny+1):(2*ny),]-Xp1[1:ny,])/delta
 X2<-(Xp2[(ny+1):(2*ny),]-Xp2[1:ny,])/delta

  # GAM model Y= b X = f(years)
 # gradients dY=b dX
 deriv1<-X2%*%coef(fit) ## gradient
 #Vg <- X1%*%vcov(fit)%*%t(X1) ## cov matrix of gradiant

 #2nd derivative is difference between adjacent gradients
 deriv2<-deriv1 - X1%*%coef(fit)

 return(data.frame(deriv1=deriv1,deriv2=deriv2))
}

######
.Bootstrap.func<-function(D)
{
#bootrap using normal distribution
res<-rep(NA,nrow(D))
for(i in 1:nrow(D)) res[i]<-rnorm(1,D[i,1],D[i,2])
res
}

######ChangeSummary.func############
.ChangeSummary.func<-function(D=Smoothed,lastn,Indic,Year,write=write,alpha)
{
if(FALSE){
D=Smoothed
Indic=Indic
Year=Year
lastn=lastn
alpha=alpha
}
#for carrying out individual union tests and final test intersection test
#D= results table of smoothed trends and derivatives with sd
#Indic= raw index estimate
#lastn= time horizon in which minimum should not be, ie 5 most recent years
#lastnmin= time horizon for minimum value
#names(Smoothed)=c("Species","Year","G","Gsd","G1","G1sd","G2","G2sd")

if(is.factor(D$Year)) D$Year<-as.numeric(levels(D$Year)[D$Year])
if(is.factor(D$G)) D$G<-as.numeric(levels(D$G)[D$G])
if(is.factor(D$Gsd)) D$Gsd<-as.numeric(levels(D$Gsd)[D$Gsd])
if(is.factor(D$G1)) D$G1<-as.numeric(levels(D$G1)[D$G1])
if(is.factor(D$G1sd)) D$G1sd<-as.numeric(levels(D$G1sd)[D$G1sd])
if(is.factor(D$G2)) D$G2<-as.numeric(levels(D$G2)[D$G2])
if(is.factor(D$G2sd)) D$G2sd<-as.numeric(levels(D$G2sd)[D$G2sd])

years<-seq(min(D$Year),max(D$Year),1)
res<-matrix(NA,1,ncol=16)

#Intersection-Union test
#There are 3 union tests, ie for each group H0 is true if any one of the tests is true
#In the final test, the global H0 of no change is true if it is true for all 3 groups

#1. test for location of maximum and minimum of series
#H0 no decrease: max in final lastn years: binomial(1,theta), theta => 1-alpha
#H1 decrease: max NOT before final lastnmin years
#H0 no increase: min in final lastn years: binomial(1,theta), theta => 1-alpha
#H1 increase: min NOT before final lastnmin years
T1Dec<-ifelse(is.element(max(D$G),rev(D$G)[1:lastn]),1,0)
T1Inc<-ifelse(is.element(min(D$G),rev(D$G)[1:lastn]),1,0)

#2. test for local slopes in last lastn years
#HO no decrease: One of the annual slopes >0
#H1 decrease: ALL annual slopes <=0
#HO no increase: One of the annual slopes <0
#H1 increase: ALL annual slopes >=0
t2dec=t2inc=rep(0,lastn)
i=1
for(i in 1:lastn){
 testdec<-1-pnorm((rev(as.numeric(D$G1))[i]/rev(as.numeric(D$G1sd))[i]),lower.tail=FALSE)
 t2dec[i]=sum(testdec<=alpha) #=1 if H0 not true
}
T2Dec<-ifelse(sum(t2dec)<lastn,1,0) # =0 if H0 rejected, =1 if H0 accepted

for(i in 1:lastn){
 testinc<-1-pnorm((rev(as.numeric(D$G1))[i]/rev(as.numeric(D$G1sd))[i]),lower.tail=TRUE)
 t2inc[i]=sum(testinc<=alpha) #=1 if H0 not true
}
T2Inc<-ifelse(sum(t2inc)<lastn,1,0) #H0 not true only if ALL tests are significant!

#3. test for second derivatives in lastn years, ignoring the final year estimate,
#H0 no decrease: some second derivatives over lastn years >=0
#H1 decrease: all second derivatives over lastn years <0
#H0 no increase: some second derivatives <=0
#H1 increase: all second derivatives >0

t3dec=t3inc=rep(0,lastn)
for(i in 1:lastn){
 testdec<-pnorm((rev(as.numeric(D$G2))[i]/rev(as.numeric(D$G2sd))[i]),lower.tail=TRUE)
 t3dec[i]=sum(testdec<=alpha) #=1 if H0 not true
}
T3Dec<-ifelse(sum(t3dec,na.rm=TRUE)<lastn,1,0) # =0 if H0 rejected, =1 if H0 accepted

for(i in 1:lastn){
 testinc<-pnorm((rev(as.numeric(D$G2))[i]/rev(as.numeric(D$G2sd))[i]),lower.tail=FALSE)
 t3inc[i]=sum(testinc<=alpha) #=1 if H0 not true
}
T3Inc<-ifelse(sum(t3inc,na.rm=TRUE)<lastn,1,0) #H0 not true only if ALL tests are significant!

#######linear model test#####################
#linear slope over nlastmin years
#print(years)
id<-Year>=(rev(years)[lastn])
SlopeLinRec<- coefficients(lm(Indic[id]~Year[id]))[2]
#p-value of linear slope over lastn years
PLinRec<- coefficients(summary(lm(Indic[id]~Year[id])))[8]
#linear slope over all years
SlopeLinAll<- coefficients(lm(Indic~Year))[2]
#p-value of linear slope over all years
PLinAll<- coefficients(summary(lm(Indic~Year)))[8]

###Mann Kendall test
#print(Indic[id])
#print(Year[id])
PMKinc<-cor.test(scale(Year[id]),Indic[id],method="kendall",alternative="greater")$p.value
PMKdec<-cor.test(scale(Year[id]),Indic[id],method="kendall",alternative="less")$p.value

#create recent diagnose ###################
#code -1 for recent "decreasing" :
#codage 1 for recent "increase" :
#T1Dec etc =1 if H0 true!!!
#Global H0 rejected if ALL H0s are rejected!!!, ie T1+T2+T3=0
#Tx=1 if corresponding H0 is true!
RecentDec<-ifelse((T1Dec+T2Dec)==0,-1,0)
ChangePointDec<-ifelse(T3Dec==1,0,-1)#=0 if change point occured, =-1 if not
RecentInc<-ifelse((T1Inc+T2Inc)==0,1,0)
ChangePointInc<-ifelse(T3Inc==1,0,1)#=0 if change point occured, =1 if not
#print(c(RecentDec,RecentInc))
#print(c("all Tinc",T1Inc, T2Inc,T3Inc,"RecentInc",RecentInc))

#Create diagnose based on linear slope for last nlast (5) years
RecentDiagnosL5<-0
RecentDiagnosL5[SlopeLinRec<0&PLinRec<=alpha]<-(-1)
RecentDiagnosL5[SlopeLinRec>0&PLinRec<=0.05]<-(1)

#create diagnos based on Mann-Kendall test
RecentMannKendallDec<-0
RecentMannKendallDec[PMKdec<=alpha]<-(-1)
RecentMannKendallInc<-0
RecentMannKendallInc[PMKinc<=alpha]<-1

res=c(RecentDec,ChangePointDec,RecentInc,ChangePointInc,RecentDiagnosL5,PLinRec,
RecentMannKendallDec,RecentMannKendallInc,PMKdec,PMKinc)
return(res)
}

#################################
#main
# remove years with missing data
id<-!is.na(Indic);Indic<-Indic[id];Year<-Year[id];Sd<-Sd[id]

#create vector of standar deviations if not provided using given CV
#same if all standard deviations are NA
if(is.null(Sd)|sum(Sd)==0|is.na(sum(Sd)))Sd<-Indic*CV
ind <- scale(1:length(Year), scale = FALSE)

resmax<-res2nd<-matrix(0,1,length(Year))
reslow<-matrix(0,1,4);reshigh<-matrix(0,1,2)
reslast<-matrix(NA,1,6)
dimnames(resmax)<-list(species,as.character(Year))
dimnames(reslow)<-list(species,c("mindata","minGam","lastYear","lastyear"))
dimnames(reshigh)<-list(species,c("maxdata","maxGam"))
dimnames(reslast)<-list(species,c("minGdata","minGGam","lastGdata","lastGGam","maxGdata","maxGGam"))
i<-1

#statistics G
G<-Indic
anneesp<-Year
ida<-(1:length(Year))[is.element(Year,anneesp)]
ida<-ida[!is.na(ida)]
#dflnG<-ceiling(SecondDeriv.func(G,anneesp,first=1))

#degrees of freedom for smooth function, minimum 2
fit<-.Ggam.func(G,anneesp,first=1)
dfG<-ceiling(fit$df)
#pvalue for chisquare test of gam fit (residual deviance)
pchisq<-fit$pchisq
if(dfG==1) dfG<-2
#linear model whole time series
indsp <- scale(1:length(anneesp), scale = FALSE)
Gfit<-lm(G ~ indsp, na.action = na.exclude)
Gpvalue <- summary(Gfit)$coefficients[8]

#smooth model; get predicted values at years with observations
Ggam<-.Ggam.func(G,anneesp,first=2,df=dfG)

Db<-Db.orig<-cbind(Indic,Sd)
Db<-data.frame(Db)
Gboot<-Gderiv1<-Gderiv2<-Ggam<-matrix(NA,nrow=B,ncol=length(Year))
Gboot[1,ida]<-G
Ggam[1,ida]<-.Ggam.func(G,anneesp,first=2,df=dfG)#predicted values
Derivs<-.Deriv.func(G,anneesp,df=dfG)
Gderiv1[1,ida]<-Derivs$deriv1
Gderiv2[1,ida]<-Derivs$deriv2

#bootstrap for standard deviations of derivatives
for(j in 2:B) {
 		  Db$Estimation<-.Bootstrap.func(Db.orig)
		  Gboot[j,ida]<- Db$Estimation
		  Ggam[j,ida]<-.Ggam.func(Gboot[j,ida],anneesp,first=2,df=dfG)
		  Derivs<-.Deriv.func(Gboot[j,ida],anneesp,df=dfG)
		  Gderiv1[j,ida]<-Derivs$deriv1
		  Gderiv2[j,ida]<-Derivs$deriv2
}
G1sd<-apply(Gderiv1,2,sd,na.rm=TRUE)
G2sd<-apply(Gderiv2,2,sd,na.rm=TRUE)
Ggamsd<-apply(Ggam,2,sd,na.rm=TRUE)

#create results matrix:
yr=length(anneesp)
Smoothed<-cbind(rep(species,yr),anneesp,Ggam[1,],Ggamsd,Gderiv1[1,],G1sd,Gderiv2[1,],G2sd)
tab=cbind(rep(species,yr),anneesp,round(pchisq,3),round(Ggam[1,],2),round(Ggamsd,2),round(Gderiv1[1,],2),
round(G1sd,2),Gderiv2[1,],G2sd)
Smoothed=data.frame(Smoothed)
names(Smoothed)=c("Species","Year","G","Gsd","G1","G1sd","G2","G2sd")
tab=data.frame(tab)
names(tab)=c("Species","Year","pchisq","G","Gsd","G1","G1sd","G2","G2sd")
if(write==TRUE)
 write.table(tab,paste(species,ytext,".txt",sep=""),col.names=TRUE,row.names=FALSE,
 sep=";",quote=FALSE)

Diagnos<-.ChangeSummary.func(D=Smoothed,Indic=Indic,Year=Year,lastn=lastn,alpha=alpha)
Diagnos=data.frame(matrix(Diagnos,nrow=1))
names(Diagnos)=c("RecentDec","ChangePointDec","RecentInc","ChangePointInc",
"RecentDiagnoseL5","Plin","RecentMannKendallDec","RecentMannKendallInc",
"PMKdec","PMKinc")
#plot indicator time series with smooths
par(xaxs = "r")
if(Diagnos$RecentDiagnoseL5==0)textdiaglin="stable"
if(Diagnos$RecentDiagnoseL5<0)textdiaglin="decrease"
if(Diagnos$RecentDiagnoseL5>0)textdiaglin="increase"

textregle="unknown"
if(Diagnos$RecentDec==0 &Diagnos$RecentInc==0) textregle="stable"
if(Diagnos$RecentDec<0) textregle="decrease"
#if(Diagnos$ChangePointDec==0) textregle=paste(textregle,"*",sep="")
if(Diagnos$RecentInc>0) textregle="increase"
#if(Diagnos$ChangePointInc==0) textregle=paste(textregle,"*",sep="")

#par(mar=c(2,2,2,0.5))
.error.bar(as.integer(anneesp), G, lower = apply(Gboot[,ida],2,quantile,0.025,na.rm=TRUE),
		upper = apply(Gboot[,ida],2,quantile,0.975,na.rm=TRUE), incr = FALSE,
main = speciesname,xlab ="year", ylab =lab_y,cex.main=1.5,axes=TRUE,bty="o")
axis(2,labels=TRUE)
ifelse(lastligne,axis(1,labels=TRUE),axis(1,labels=TRUE))
box()
#main = paste(species,"\n","Changes in last", lastn,"years", "linear:",textdiaglin,
#"; rule:",textregle),xlab ="Year", ylab = ytext,cex.main=0.7)
#add smooth lines from gam + 2.5 and 97.5 percentiles
lines(anneesp,Ggam[1,ida],lty=1)
lines(anneesp,apply(Ggam[,ida],2,quantile,0.025,na.rm=TRUE),lty=2)
lines(anneesp,apply(Ggam[,ida],2,quantile,0.975,na.rm=TRUE),lty=2)
#savePlot(paste(ytext,species),type="wmf")
#
return(c(species,ytext,pchisq,Diagnos))

}



.error.bar<-function(x, y = NULL, lower, upper, incr = TRUE, bar.ends = FALSE, gap = FALSE, add = FALSE, horizontal = FALSE, ..., xlab = deparse(substitute(x)), xlim, ylim) {
	draw.null.warn <- function(draw, gap)
	{
		if(!any(draw)) {
			warning("Not enough room for a gap.")
			draw <- !draw
			gap <- 0
		}
		invisible(list(draw = draw, gap = gap))
	}
	if(missing(x))
		stop("no data for x or y")
	if(missing(y)) {
		if(missing(xlab))
			xlab <- "Index"
		y <- x
		x <- time(x)
	}
	n <- length(x)
	if(length(y) != n)
		stop("length of y must equal the length of x")
	center <- if(horizontal) x else y
	if(missing(lower))
		stop("you must provide lower")
	if(length(lower) > 1 && length(lower) != n)
		stop("length of lower must be 1 or equal to the length of x")
	#if incr=T lower is assumed >=0
	if(incr) lower <- center - abs(lower) else lower <- rep(lower,
			length = n)
	if(any(lower >= center))
		warning(paste(
			"There are values of 'lower' which are greater or equal to ",
			if(horizontal) "x" else "y"))
	if(missing(upper))
		upper <- 2 * center - lower
	else {
		if(length(upper) > 1 && length(upper) != n)
			stop("length of upper must be 1 or equal to the length of x"
				)
		if(incr)
			upper <- center + upper
		else upper <- rep(upper, length = n)
	}
	if(any(upper <= center))
		warning(paste(
			"There are values of 'upper' which are smaller or equal to ",
			if(horizontal) "x" else "y"))
	if(!add)
		if(horizontal) {
			if(missing(ylim))
				plot(x, y, xlim = if(missing(xlim)) range(
						c(lower, upper), na.rm = TRUE
						) else xlim, xlab = xlab,
					...)
			else plot(x, y, xlim = if(missing(xlim)) range(
						c(lower, upper), na.rm = TRUE
						) else xlim, ylim = ylim,
					xlab = xlab, ...)
		}
		else {
			if(missing(xlim))
				plot(x, y, ylim = if(missing(ylim)) range(
						c(lower, upper), na.rm = TRUE
						) else ylim, xlab = xlab,
					...)
			else plot(x, y, ylim = if(missing(ylim)) range(
						c(lower, upper), na.rm = TRUE
						) else ylim, xlim = xlim,
					xlab = xlab, ...)
		}
	if(horizontal) {
		if(gap)
			gap <- 0.75 * par("cxy")[1]
		draw <- x - lower > gap
		z <- draw.null.warn(draw, gap)
		draw <- z$draw
		gap <- z$gap
		segments(lower[draw], y[draw], x[draw] - gap, y[draw])
		draw <- upper - x > gap
		z <- draw.null.warn(draw, gap)
		draw <- z$draw
		gap <- z$gap
		segments(x[draw] + gap, y[draw], upper[draw], y[draw])
		if(bar.ends) {
			size.bar <- par("cxy")[2]
			segments(lower, y - size.bar, lower, y + size.bar)
			segments(upper, y - size.bar, upper, y + size.bar)
		}
	}
	else {
		if(gap)
			gap <- 0.75 * par("cxy")[2]
		draw <- upper - y > gap
		z <- draw.null.warn(draw, gap)
		draw <- z$draw
		gap <- z$gap
		segments(x[draw], y[draw] + gap, x[draw], upper[draw])
		draw <- y - lower > gap
		z <- draw.null.warn(draw, gap)
		draw <- z$draw
		gap <- z$gap
		segments(x[draw], y[draw] - gap, x[draw], lower[draw])
		if(bar.ends) {
			size.bar <- par("cxy")[1]
			segments(x - size.bar, upper, x + size.bar, upper)
			segments(x - size.bar, lower, x + size.bar, lower)
		}
	}
}



