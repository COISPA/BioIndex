

#' Interception Union Tets
#'
#' @param abundance dataframe of abundance time series as produced by indices_ts function
#' @param biomass dataframe of biomass time series as produced by indices_ts function
#' @param species reference species for the analysis (MEDITS code)
#' @param lastn number of recent years for diagnosis of change
#' @param save boolean. If TRUE results are saved in the output folder
#' @importFrom utils write.table
#' @importFrom grDevices jpeg dev.copy dev.off
#' @export
#' @import mgcv

IUT <- function(abundance, biomass, species, lastn=5, save=TRUE) {

    if (FALSE) {
        GSA=18
        save=TRUE
        depth_range <- c(10,800)

        verbose=TRUE
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex_app\\R_BioIndex_3.3_(in update)"
        ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=TRUE)
        tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=TRUE)
        tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=TRUE)
        country="all"

        species <- "MERLMER"

        m <- merge_TATBTC(ta, tb, tc, species="MERLMER", country="all", wd=wd, verbose=TRUE)
        mTATB <- m[[1]]
        ind <- indices_ts(mTATB, GSA=18, country="all", depth_range=c(10,800), strata_scheme=BioIndex::strata_scheme, stratification=BioIndex::stratification,wd, save=FALSE)

        IUT(abundance=ind[[1]], biomass=ind[[2]], species="MERLMER", lastn=4, save=TRUE)

    }

    write=FALSE
    alpha=0.05
    CVindicators=0.05 	#used in bootstrap for uncertainty
    system<-paste("GSA", GSA, sep="")
    sspp <- species
    #####END PART THAT HAS TO BE ADAPTED BY USER##########


    ##MAIN part for carrying out analysis###########

    #read data table
    df1 <- abundance
    df2 <- biomass
    dat <- data.frame(year = df1$year , ecosystem = paste("GSA",GSA, sep=""), abundance = df1$abundance, biomass = df2$biomass)
    response <- c("Abundance (n/km^2)", "Biomass (kg/km^2)")
    cols_1_2 <- c("year", "ecosystem")

    dat[,3] <- dat[,3]
    dat[,4] <- dat[,4]
    indicators= colnames(dat)[!(colnames(dat) %in% cols_1_2)]  # c("ARITANT")
    indiname=indicators


    #col.names(dat)=c("Year","ecosystem","delta_10","L95_10","10ARISFOLF2D","10PAPE LONF2D","10PAPE LONM2D","10ARISFOLM 2D","10NEPRNORM2D")
    nrows <- length(indicators)/2
    ncols <- 2


    # windows(width=(10),height=(5*nrows))
    par(mfrow=c(nrows,ncols),mar=c(5,5,1.5,0.5),las=1,lab=c(5,3,5))


    #for each indicator run test and store results
    #results matrix
    Res<-matrix(NA,nrow=length(indicators),ncol=14)
    dimnames(Res)<-list(NULL,c("Area","Indicator","Indicators","pchisq","RecentDecrease","ChangePointDec","RecentIncrease","ChangePointInc","RecentLinear","Plinear","DecMK","IncMK","PMKdec","PMKinc"))
    hasdata=0
    j=1
    for(j in 1:length(indicators)){
        lastligne=ifelse(j==length(indicators),TRUE,FALSE)
        indi=dat[,colnames(dat)==indicators[j]]
        id=!is.na(indi)# &dat[,1]<2006
        indi=indi[id]
        years<-dat[id,1]
        Res[j,1]<-system
        Res[j,2:3]<-indicators[j]
        lab_y <- response[j]
        if(length(indi)>=lastn){
            Res[j,2:14]<-unlist(.ChangeYears.func(Year=years,Indic=indi,Sd=rep(NA,length(indi)),species=system,ytext=indicators[j],CV=CVindicators,lastn=lastn,write=write, speciesname=indiname[j],alpha=alpha,lastligne=lastligne))
            hasdata=1
        }
    }


    if(save){
        write.table(Res,paste(wd,"/output/",system,"IUtest",lastn,"_RSS.csv",sep=""),sep=";",quote=FALSE,row.names=FALSE,append=FALSE)

       if(hasdata==1)dev.copy(jpeg,paste(wd,"/output/",sspp,"_",system,"_SmoothedIndicators_RSS.jpg",sep=""),height=5,width=10, units='in', res=300)
        dev.off()
    }
    #summarise results across indicators########

    ni=length(indicators)
    res=matrix(NA,nrow=ni*length(system),ncol=14)
    i=1
    for(i in 1:length(system)){
        res[((i-1)*ni+1):(i*ni),1:14]=as.matrix(Res)
    }
    header=names(data.frame(Res))

    dimnames(res)=dimnames(res)=list(NULL,header)
    res=data.frame(res[,c(1,3,4,5,7)])
    res[,4] <- as.factor(res[,4])
    res[,5] <- as.factor(res[,5])
    j=3
    for(j in 3:5) {res[,j]=as.numeric(as.character(res[,j]))}
    res$Change=res$RecentDecrease+res$RecentIncrease

    #create table with results, one column per indicator, one row per ecosystem
    indi=sort(unique(res[,2]))
    table=data.frame(cbind(res[res[,2]==indi[1],6],res[res[,2]==indi[2],6]))
    names(table)=indi

    table$system=res[res[,2]==indi[1],1]
    table=table[,c(3,1:2)]

    #replace -1 by - 1 by + and 0 by ""
    for (c in 2:3){
        if (table[c] == -1){table[c] = "decrease"}
        if (table[c] == 1){table[c] = "increase"}
        if (table[c] == 0){table[c] = "nochange"}

    }

    if (save){
        write.table(table,paste(wd, "/output/",sspp,"_GSA",GSA,"_IUT_results_",lastn,"years_RSS.csv",sep=""),sep=";",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
return(table)
}
