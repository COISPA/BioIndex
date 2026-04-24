#' Spearman test for timeseries
#'
#' @param x time series
#' @importFrom stats pt
#' @export
spear<-function (x){

    N = length(x)
    S = rank(x,ties.method ="min") # crank(x)
    R = c(1:length(x))

    r = 1-6*sum((R-S)^2)/N/(N^2-1)
    t=r*sqrt((N-2)/(1-r^2))

    p=2*(1-pt(abs(t),N-2)) #
    results=data.frame(r=r,t=t,p=p)
    return(results)
}
