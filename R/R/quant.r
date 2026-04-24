#' Quantile estimation
#'
#' @param weighted LFD data.frame
#' @param qlin reference quantile for the analysis
#' @export
quant <- function(weighted,qlin=0.95){
    weighted <- data.frame(weighted)
    weighted <- weighted[weighted$X2 != 0 , ]
    w <- aggregate(weighted[,2],list(G=weighted[,1]),FUN=sum,na.rm=T)
    x <- as.numeric(as.character(w$G))
    w <- w$x
    if(qlin<0.5) {
        x <- -x
        ql <- 1-qlin
    } else {
        ql <- qlin
    }
    tot <- sum(w)
    frac <- tot*(1-ql)
    w <- w[order(-x)]
    x <- x[order(-x)]
    somme <- w[1]
    i <- 2
    while(somme < frac){
        somme <- somme + w[i]
        i <- i+1
    }
    lq <- round(x[i] + (x[i-1]-x[i])/w[i-1]*(somme-frac),1)
    f1 <- w[i]/tot
    vlq <- ql*(1-ql)/(f1^2)/tot
    SD <- sqrt(vlq)
    cv <- sqrt(vlq)/ifelse(qlin<0.5,-lq,lq)
    c(quantile=ifelse(qlin<0.5,-lq,lq),SD=SD,cv=cv)
}
