removecomma <- function(x) {
    ## remove trailing comma from input string
    ## also remove escaped quotation marks
    x <- substring(x, 1, nchar(x)-1)
    x <- gsub("\"", "", x)
    return(x)
}

## check that valid options are specified
getinputs <- function(chemdb) {
    # choose from 13880 (wet), 94018 (dry), and 94846 (temperate)
    # 1: crop, 2:pasture, 3:reclamation
    require(stringr)
    flink <- file("./input choices.csv", open = "r")
    x <- readLines(flink)
    close(flink)
    b <- str_split(x, ",")
    namesav <- rep("", times = length(b))
    for (i in 1:length(b)) {
        nelem <- length(b[[i]])
        namesav[i] <- b[[i]][1]
        b[[i]] <- b[[i]][-1]
        if (namesav[i] %in% c("cas", "linertypeid"))
            b[[i]] <- sort(as.numeric(b[[i]])) # sort to get rid of NAs
        else {
            b[[i]] <- sort(b[[i]])
            incvec <- b[[i]] != ""
            b[[i]] <- b[[i]][incvec]
        }
    }

    names(b) <- namesav

    if (any(! b$siteid %in% c("wet", "dry", "temperate")))
        stop("Invalid siteid")
    else {
        x <- b$siteid
        b$siteid <- rep(0, times = length(x))
        incvec <- x == "wet"
        b$siteid[incvec] <- 13880
        incvec <- x == "dry"
        b$siteid[incvec] <- 94018
        incvec <- x == "temperate"
        b$siteid[incvec] <- 94846
    }
    if (any(! b$cas %in% chemdb$cas))
        stop("Invalid CAS")
    if (any(! b$srctype %in% c("LAU", "SI")))
        stop("Invalid srctype")
    if (any(! b$srctypeid %in% c("crop", "pasture", "reclamation")))
        stop("Invaid srctypeid")
    else {
        x <- b$srctypeid
        b$srctypeid <- rep(0, times = length(x))
        incvec <- x=="crop"
        b$srctypeid[incvec] <- 1
        incvec <- x=="pasture"
        b$srctypeid[incvec] <- 2
        incvec <- x=="reclamation"
        b$srctypeid[incvec] <- 3
    }
    if (any(! b$linertypeid %in% 0:2))
        stop("Invalid linertypeid")

    return(b)
}

## generate scenarios from selections
makescen <- function(b) {

    ncas<- length(b$cas)
    nsite <- length(b$siteid)
    nsrctype <- length(b$srctype)
    grps <- data.frame(siteid = rep(b$siteid, times = ncas),
                       cas = rep(b$cas, times = rep(nsite, times = ncas)))

    one <- F
    two <- F
    if ("LAU" %in% b$srctype) {
        n <- length(b$srctypeid)
        ngrp <- nrow(grps)
        grps1 <- do.call("rbind", replicate(n, grps, simplify = F))
        grps1$srctype <- "LAU"
        grps1$srctypeid <- rep(b$srctypeid, times = rep(ngrp, times = n))
        grps1$linertypeid <- 1
        one <- T
    }
    if ("SI" %in% b$srctype) {
        n <- length(b$linertypeid)
        ngrp <- nrow(grps)
        grps2 <- do.call("rbind", replicate(n, grps, simplify = F))
        grps2$srctype <- "SI"
        grps2$srctypeid <- 1
        grps2$linertypeid <- rep(b$linertypeid, times = rep(ngrp, times = n))
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
        if (! is.na(suppressWarnings(as.numeric(x)))) x <- as.numeric(x)
        return(x)
    }
    else return(NA)
}

## return names vector with all modelcodes
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

    incvec <- !is.na(suppressWarnings(as.numeric(x)))
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
