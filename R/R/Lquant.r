#' Estimation of L50 and L95
#'
#' @param lfd data frame of combined LFD
#' @param wd working directory
#' @param sspp MEDITS code for the selected species
#' @param GSA reference area for the analysis
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE messages are reported in the console
#' @export Lquant



Lquant <- function(lfd, wd=NA, sspp, GSA, save=TRUE, verbose=TRUE) {

    if (FALSE) {
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        # ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        # tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        # tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)

        ta <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TA_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TB_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TC_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        # ta <- ta[ta$YEAR==2017, ]
        # tb <- tb[tb$YEAR==2017, ]
        # tc <- tc[tc$YEAR==2017, ]
        depth_range=c(10,800)
        country="all"
        sex <- "all"
        species <- "MERLMER"
        save=TRUE

        m <- merge_TATBTC(ta, tb, tc, species=species, country=country, wd=wd, verbose=TRUE)
        mTATC <- m[[2]]

        lfd <- LFD(mTATC, sex=sex, GSA=18, country=country, depth_range=c(10,800), strata_scheme, stratification, wd, save=TRUE)
        lfd <- lfd[[1]]

        Lquant(lfd,wd, save=TRUE)
    }

    if (is.na(wd) & save) {
        save =FALSE
        if (verbose){
            message("Missing working directory. Results are not saved in the local folder.")
        }
    }

    Dati=lfd

col_names <- colnames(Dati)
col_years <- col_names[-1]
# col_years <- as.numeric(substr(col_years,2,5))
years_tab <- col_years

table=data.frame(matrix(NA,nrow=6,ncol=(ncol(Dati)-1)))
i=2
for (i in c(2:ncol(Dati))){
    Dati_temp=cbind(Dati[,1],Dati[,i])
    table [1,i] = round(quant(Dati_temp, 0.50)[1],4)
    table [2,i] = round(quant(Dati_temp, 0.50)[2],4) #round(sqrt(sum(Dati_temp[,2]*(Dati_temp[,1]-table[6,i])^2)/sum(Dati_temp[,2])),2)
    table [3,i] = round(quant(Dati_temp, 0.50)[3],4)
    table [4,i] = round(quant(Dati_temp, 0.95)[1],4)
    table [5,i] = round(quant(Dati_temp, 0.95)[2],2)
    table [6,i] = round(quant(Dati_temp, 0.95)[3],4)
    Dati_temp[,] = 0
}
table[,1]=c("50th perc.","SD 50th perc.","CV 50th perc.","95th perc.","SD 95th perc.","CV 95th perc.")
colnames(table)=c("Indices", years_tab)
if (save){
write.table(table,paste(wd, "/output/",sspp,"_GSA",GSA,"_L50_L95_RSS.csv", sep=""),sep=";",row.names=F)
}

#formatting table for plots
par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
t <- data.frame(t(table))
headers <- as.character(t[1,])
t <- t[-1,]
years <- row.names(t)
t <- cbind(years, t)
colnames(t) <- c("years",headers)
row.names(t) <- seq(1,nrow(t),1)
for (i in 1:ncol(t)){
 t[,i] <- as.numeric(t[,i])
}
#plot
rL50 <- range(t$`50th perc.`)
r1 <- (rL50[2]-rL50[1])*20/100
rL50[2] <- rL50[2]+r1
rL95 <- range(t$`95th perc.`)
r2 <- (rL95[2]-rL95[1])*20/100
rL95[2] <- rL95[2]+r2

plot(t[,1],  t[,2], type="b", col="black", pch=16, xlab="year", ylab="L50 (mm)", ylim=rL50) # ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab # ylim=c(0,max_index*1.2)
lines(t[,1], (t[,2]-1.96*t[,3]), type="l",lty=2, col="red" )
lines(t[,1], (t[,2]+1.96*t[,3]), type="l",lty=2, col="red" )
legend("topright", c("time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))

plot(t[,1],  t[,5], type="b", col="black", pch=16, xlab="year", ylab="L95 (mm)", ylim=rL95) # ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab # ylim=c(0,max_index*1.2)
lines(t[,1], (t[,5]-1.96*t[,6]), type="l",lty=2, col="red" )
lines(t[,1], (t[,5]+1.96*t[,6]), type="l",lty=2, col="red" )
legend("topright", c("time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))

if (save){
jpeg(paste(wd,"/output/L50_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
   plot(t[,1],  t[,2], type="b", col="black", pch=16, xlab="year", ylab="L50 (mm)", ylim=rL50) # ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab # ylim=c(0,max_index*1.2)
   lines(t[,1], (t[,2]-1.96*t[,3]), type="l",lty=2, col="red" )
   lines(t[,1], (t[,2]+1.96*t[,3]), type="l",lty=2, col="red" )
   legend("topright", c("L50 time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))
dev.off()


jpeg(paste(wd,"/output/L95_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
   plot(t[,1],  t[,5], type="b", col="black", pch=16, xlab="year", ylab="L95 (mm)", ylim=rL95) # ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab # ylim=c(0,max_index*1.2)
   lines(t[,1], (t[,5]-1.96*t[,6]), type="l",lty=2, col="red" )
   lines(t[,1], (t[,5]+1.96*t[,6]), type="l",lty=2, col="red" )
   legend("topright", c("l95 time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))
dev.off()
}
return(t)

}


