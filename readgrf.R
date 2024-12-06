## read grf file

readgrf <- function(srctype) {
    require(stringr)

    if (srctype == "SI") {
        varpick <- c("VE", "AnnInfil", "LeachFlux")
        varlab <- c("ve", "anninfil", "leachflux")
    }
    else {
        varpick <- c("CE", "VE", "SWLoadChem", "CTda", "CTss", "AnnInfil", "Runoff",
                      "LeachFlux")
        varlab <- c("ce", "ve", "swloadchem", "ctda", "ctss", "anninfil",
                    "rflwa","leachflux")
    }

    fp <- file("./hwir/grf/sr.grf", open = "r")
    a <- readLines(fp)
    close(fp)
    removecomma <- function(x) {
        ## remove trailing comma from input string
        ## also remove escaped quotation marks
        x <- substring(x, 1, nchar(x)-1)
        x <- gsub("\"", "", x)
        return(x)
    }
    a <- removecomma(a)
    b <- str_split(a, ",")

    ip <-  which(sapply(b, function(x) is.na(as.numeric(x[1]))))

    ip <- ip[-1]
    varname <- sapply(b[ip], function(x) x[1])

    selvec <- regexpr("NY", varname) != -1
    varpick <- c(varpick, varname[selvec])

    ## set up storage
    nvar <- length(varpick)
    dat <- vector("list", nvar)
    names(dat) <- varpick

    ipp <- rep(NA, times = length(varpick))
    for (i in 1:length(varpick)) {
        isel <- which(varname == varpick[i])
        if (length(isel) == 1) ipp[i] <- ip[isel]
    }

    icount <-0
    for (i in ipp) {
        icount <- icount + 1
        if (!is.na(i)) {
            ## get name and dimensions
            ndim <- as.numeric(b[[i]][2])
            if (ndim == 0) {
                x <- as.numeric(b[[i+1]])
                dat[[icount]] <- x
            }
            if (ndim == 1) {
                ic <- i+1
                x <- as.numeric(b[[ic]])[-1]
                dat[[icount]] <- x
            }
            if (ndim == 2) {
                ic <- i + 2
                d1 <- as.numeric(b[[i+1]])
                d2 <- sapply(b[ic:(ic+d1-1)], function(x) as.numeric(x[1]))
                d2sav <- max(d2)

                x <- matrix(0, nrow = d1, ncol = d2sav)
                for (j in 1:d1)
                    x[j,] <- as.numeric(b[[ic + j - 1]])[-1]
                dat[[icount]] <- x
            }
            if (ndim == 3) {
                ## figure out third dimension of this 3 dimensional array
                d <- as.numeric(b[[i+1]])
                ic <- i+2
                d3sav <- 0
                for (j in 1:d[1]) {
                    d3 <- sapply(b[ic:(ic+d[2]-1)], function(x) as.numeric(x[1]))
                    d3sav <- max(d3, d3sav)
                    ic <- ic + d[2] + 1
                    if (j < d[1]) d[2] <- as.numeric(b[[ic-1]])
                }
                ## define array
                ## 365 days needed if I save the daily data
                ## otherwise, only 2
                x <- rep(0, times = d[1]*2*d3sav)
                dim(x) <- c(d[1], 2, d3sav)
                ## fill the array
                d <- as.numeric(b[[i+1]])
                ic <- i+2
                for (j in 1:d[1]) {
                    d3 <- sapply(b[ic:(ic+d[2]-1)], function(x) as.numeric(x[1]))
                    ifill <- which(d3 > 0)
                    for (ii in ifill) {
                        x[j, ii,1:d3[ii]] <- as.numeric(b[[ic + ii-1]])[-1]
                    }
                    ic <- ic + d[2] + 1
                    if (j < d[1]) d[2] <- as.numeric(b[[ic-1]])
                }
                dat[[icount]] <- x
            }
        }
        else {
            dat[[icount]] <- 0
        }
    }

    return(dat)

}

#readgrf()

