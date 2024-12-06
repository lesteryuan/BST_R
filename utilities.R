## generate scenarios from selections
makescen <- function(siteid, cas, srctype, srctypeid, linertypeid) {

    ncas<- length(cas)
    nsite <- length(siteid)
    nsrctype <- length(srctype)
    grps <- data.frame(siteid = rep(siteid, times = ncas),
                       cas = rep(cas, times = rep(nsite, times = ncas)))

    one <- F
    two <- F
    if ("LAU" %in% srctype) {
        n <- length(srctypeid)
        ngrp <- nrow(grps)
        grps1 <- do.call("rbind", replicate(n, grps, simplify = F))
        grps1$srctype <- "LAU"
        grps1$srctypeid <- rep(srctypeid, times = rep(ngrp, times = n))
        grps1$linertypeid <- 1
        one <- T
    }
    if ("SI" %in% srctype) {
        n <- length(linertypeid)
        ngrp <- nrow(grps)
        grps2 <- do.call("rbind", replicate(n, grps, simplify = F))
        grps2$srctype <- "SI"
        grps2$srctypeid <- 1
        grps2$linertypeid <- rep(linertypeid, times = rep(ngrp, times = n))
        two <- T
    }
    if (one & two) {
        grps.all <- rbind(grps1, grps2)
    }
    else {
        if (one) grps.all <- grps1
        else if (two) grps.all <- grps2
    }
    grps.all$runid <- 1:nrow(grps.all)

    return(grps.all)
}


getval <- function(df, code) {
    if (code %in% df$modelcode) {
        x <- df$value[df$modelcode == code]
        if (! is.na(as.numeric(x))) x <- as.numeric(x)
        return(x)
    }
    else return(NA)
}

## return names vector will all modelcodes
getval.all <- function(dflist) {
    x <- dflist[[1]]$value
    names(x) <- dflist[[1]]$modelcode

    y <- dflist[[2]]$value
    names(y) <- dflist[[2]]$modelcode

    ## check for duplicates with priority to x
    dropsel <- names(y) %in% names(x)
    y <- y[!dropsel]
    x <- c(x,y)

    z <- dflist[[3]]$value
    names(z) <- dflist[[3]]$modelcode
    dropsel <- names(z) %in% names(x)
    z <- z[!dropsel]
    x <- c(x,z)

    z <- dflist[[4]]$value
    names(z) <- dflist[[4]]$modelcode
    dropsel <- names(z) %in% names(x)
    z <- z[!dropsel]
    x <- c(x,z)

    incvec <- !is.na(as.numeric(x))
    y <- as.numeric(x[incvec])
    names(y) <- names(x)[incvec]
    x <- x[!incvec]

    ## drop string values containing 'specific'
    incvec <- regexpr("specific", tolower(x)) != -1
    x <- x[!incvec]

    incvec <- y == -9999
    if (sum(incvec) > 0) y[incvec] <- NA

    return(list(x,y))
}
