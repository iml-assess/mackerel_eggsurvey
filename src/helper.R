##' save pngs in one line
##' @param plot plot
##' @param wd working directory in which to plase image
##' @param name file name
##' @param dim dimensions of png in cm
##' @param ... arguments passed to png
##' @importFrom gridExtra grid.arrange
##' @details Saves plot as png in one line
##' @rdname savepng
##' @export
savepng <- function(plot,wd, name,dim,res=300,...){
    png(file=paste0(wd,name,".png"),units="cm",width=dim[1],height=dim[2],res=res,...)
    if("ggplot" %in% class(plot)) grid.arrange(plot) else plot
    dev.off()
}


##' read excel files more easily
##' @param file file
##' @param sheet sheet
##' @import XLConnect
##' @details read excel file
##' @rdname readexcel
##' @export
readexcel <- function(file=file,sheet=sheet){
    #require(XLConnect, pos=4)  #from excel file
    .Workbook <- XLConnect::loadWorkbook(file)
    if(length(sheet)>1){
        l=list()
        for(k in 1:length(sheet)){l[[k]]=XLConnect::readWorksheet(.Workbook, sheet[k])}
    }else{l <- XLConnect::readWorksheet(.Workbook, sheet)}
    remove(.Workbook)
    return(l)
}

##' Maximum numbers
##' @param x vector
##' @param N Numbers of top
##' @details Gives the top N values of a vector
##' @export
maxN <- function(x, N=2){
    x <- x[!is.na(x)]
    len <- length(x)
    if(N>len){
        warning('N greater than length(x).  Setting N=length(x)')
        N <- length(x)
    }
    head(sort(x,decreasing = T),N)
}