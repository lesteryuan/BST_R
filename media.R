media <- function(parm, pnum,
                  dfsource, nyear, siteid, srctypeid, srctype) {
    ## calculate kdsoil
    chemtype <- parm["chemtype"]
    if (is.na(pnum["koc"])) pnum["koc"] <- 0
    if (chemtype == "M" | chemtype == "Hg") {
        kdsoil <- pnum["kd"]
        calc.kdsw.hhrap <- kdsoil
        calc.kdbs.hhrap <- kdsoil
    }
    else {
        if (chemtype == "D" | chemtype == "O") {
            kdsoil <- pnum["koc"]*pnum["foc_soil"]
            calc.kdsw.hhrap <- pnum["koc"]*pnum["foc_sw"]
            calc.kdbs.hhrap <- pnum["koc"]*pnum["foc_bs"]
        }
    }

    imp.area.field <- pnum["area_1"]*pnum["pi_1"]
    imp.area.rws <- pnum["area_4"]*pnum["pi_4"]

    newd <- T
    if (newd) {
        intleachfluxny <- dfsource$LeachFluxNY
        intnyrmet <- dfsource$NyrMet
    ## check all NY values for the maximum number of years
        isel <- which(regexpr("NY", names(dfsource)) != -1)
        numyear <- 0
        for (ii in isel) numyear <- max(numyear, dfsource[[ii]])
    }
    else {
        intleachfluxny <- dfsource$Numeric_Value[dfsource$VarName == "LeachFluxNY"]
        intnyrmet <- dfsource$Numeric_Value[dfsource$VarName == "NyrMet"]
        ## check all NY values for the maximum number of years
        incvec <- regexpr("NY", dfsource$VarName) != -1
        numyear <- max(dfsource$Numeric_Value[incvec])
    }

    if (srctype == "SI") {
        varnames <- c("VE", "AnnInfil", "LeachFlux")
        varlab <- c("ve", "anninfil", "leachflux")
        dim0 <- c(1,2,2)
    }
    else {
        ## make dfsource into matrices
        varnames <- c("CE", "VE", "SWLoadChem", "CTda", "CTss", "AnnInfil", "Runoff",
                      "LeachFlux")
        varlab <- c("ce", "ve", "swloadchem", "ctda", "ctss", "anninfil",
                    "rflwa","leachflux")
        dim0 <- c(1,1,1,2,2,3,3,2)
    }

    src.list <- as.list(rep(NA, times = length(varnames)))
    names(src.list) <- varlab
    if (newd) {
        for (i in 1:length(varnames)) {
            if (! is.null(dfsource[[varnames[i]]])) {
                dim0 <- dim(dfsource[[varnames[i]]])
                if (is.null(dim0)) {
                    x <- dfsource[[varnames[i]]]
                }
                else {
                    if (length(dim0) ==2) {
                        x <- dfsource[[varnames[i]]][1,]
                    }
                    else {
                        if (length(dim0) == 3) {
                            x <- dfsource[[varnames[i]]][1,1:2,]
                            x <- t(x)
                        }
                    }
                }
                ## pad x if it's shorter than numyear
                if (is.null(dim(x))) {
                    if (length(x) < numyear)
                        x <- c(x, rep(0, times = numyear - length(x)))
                    ## otherwise chop it at numyear
                    else x <- x[1:numyear]
                }
                else {
                    if (nrow(x) < numyear) {
                        matadd <- matrix(0, nrow = numyear - nrow(x),
                                         ncol = 2)
                        x <- rbind(x, matadd)
                    }
                    else {
                        x <- x[1:numyear,]
                    }
                }
                src.list[[i]] <- x
            }
        }
    }
    else {
        for (i in 1:length(varnames)) {

            incvec <- dfsource$VarName == varnames[i]

            ## 1 dimension
            if (dim0[i] == 1) {
                x <- rep(0, times = numyear)
                x[dfsource$Index1[incvec]] <- dfsource$Numeric_Value[incvec]
            }
            else {
                if (dim0[i] == 2) {
                    if (all(is.na(dfsource$Index3[incvec]))) {
                        ## this is for leachflux.  not sure
                        ## this format occurs for anything else
                        x <- rep(0, times = numyear)
                        x[dfsource$Index2[incvec]] <- dfsource$Numeric_Value[incvec]
                    }
                    else {
                        nrow <- numyear
                        ncol <- max(dfsource$Index2[incvec])
                        x <- matrix(0, ncol = ncol, nrow = nrow)
                        print(dfsource$Index2[incvec])
                        for (j in 1:nrow(dfsource[incvec,])) {
                            ip <- dfsource$Index3[incvec][j]
                            jp <- dfsource$Index2[incvec][j]
                            if (ip <= numyear) {
                                x[ip, jp] <- dfsource$Numeric_Value[incvec][j]
                            }
                        }
                    }
                }
                else {
                    if (dim0[i] == 3) {
                        ## non-continuous. I *think* the issue is that this
                        ## output might not be written in order so we have
                        ## read the indices, which is what I'm doing anyway
                                        #                    nx <- max(dfsource$Index2[incvec])
                        x <- rep(0, times = numyear)
                        x[dfsource$Index2[incvec]] <- dfsource$Numeric_Value[incvec]
                        if (length(x) > numyear) x <- x[1:numyear]
                    }
                }
            }
            src.list[[i]] <- x
        }
    }


    airpar <- ""
    if (srctype == "SI") {
        src.list <- append(src.list, list(ce = 0))
        airpar <- "si"  # added parameter specification for SI case
    }

    q <- (src.list[["ce"]] + src.list[["ve"]])/86400
    fv <- (src.list[["ve"]]/86400)/q
    fv[q==0] <- 0
    cair <- matrix(0, nrow = numyear, ncol = 5)
    cvapor <- matrix(0, nrow = numyear, ncol= 5)
    dp <- matrix(0, nrow = numyear, ncol = 5)
    dv <- matrix(0, nrow = numyear, ncol = 5)
    for (i in 1:5) {
        ## need to incorporate special case of Hg
        cyv <- pnum[paste("cyv",airpar, "_", i, sep = "")]
        cyp <- pnum[paste("cyp",airpar, "_", i, sep = "")]
        dydp <- pnum[paste("dydp",airpar, "_", i, sep = "")]
        dywp <- pnum[paste("dywp",airpar, "_", i, sep = "")]
        dydv <- pnum[paste("dydv", airpar, "_", i, sep = "")]
        dywv <- pnum[paste("dywv", airpar, "_", i, sep = "")]
        cair[,i] <- q*(fv*cyv+(1-fv)*cyp)*0.001
        cvapor[,i] <- q*fv*cyv
        if (i %in% c(1,2,4)) {
            dp[,i] <- (1000*q)*(1-fv)*(dydp + (pnum["fw"]*dywp))
            dv[,i] <- (1000*q)*fv*(dydv + (pnum["fw"]*dywv))
        }
    }

    cleach <- src.list[["leachflux"]]/src.list[["anninfil"]]
    incvec <- src.list[["anninfil"]] == 0
    cleach[incvec] <- 0

    if (srctype == "SI") cgwater <- cleach/pnum["daf_si"]
    else cgwater <- cleach/pnum["daf"]

    if (srctype != "SI") {

        vvwm_out <- as.list(rep(NA, times = 2))
        names(vvwm_out) <-c("3", "5")
        vvwm_max <- as.list(rep(NA, times = 2))
        names(vvwm_max) <- c("3", "5")
        vvwm_mn <- as.list(rep(NA, times = 2))
        names(vvwm_mn) <- c("3", "5")
        ## calculate inputs and run waterbody model
        for (ii in c(3,5)) {
            if (ii == 3) {
                area.contrib <- pnum["area_2"]
                waterbody.area <- pnum[paste("area", ii, sep = "_")]
            }
            else {
                area.contrib <- pnum["area_1"]
                waterbody.area <- pnum["area_fp"]
            }
            cyv <- pnum[paste("cyv", ii, sep = "_")]
            dydp <- pnum[paste("dydp",ii, sep = "_")]
            dywp <- pnum[paste("dywp",ii, sep = "_")]
            dydv <- pnum[paste("dydv",ii, sep = "_")]
            dywv <- pnum[paste("dywv",ii, sep = "_")]
            dytp <- dydp + dywp
            dytv <- dydv + dywv

            kl <- sqrt(pnum["cd"])*pnum["uw"]*
                sqrt(pnum["rho_air"]/pnum["rho_water"])*pnum["kappa"]^0.33/pnum["lambda_z"]*
                    (pnum["mu_water"]/(pnum["rho_water"]*pnum["dw"]))^(-0.67)*31536000
            kg <- sqrt(pnum["cd"])*pnum["uw"]*pnum["kappa"]^0.33/pnum["lambda_z"]*
                (pnum["mu_air"]/(pnum["rho_air"]*pnum["da"]))^(-0.67)*31536000

            if (is.na(pnum["hlc"])) pnum["hlc"] <- 0
            hprime <- pnum["hlc"]/0.00008205/pnum["twater"]
            tempadjust <- pnum["theta_water"]^(pnum["twater"] - pnum["thlc"])
            kv <- (1/((1/kl) + (1/(kg*hprime))))*tempadjust
            ldep <- q*(fv*dytv + (1-fv)*dytp)*waterbody.area
            if (pnum["hlc"] > 0)
                ldif <- (kv*(q*fv*cyv)*waterbody.area*0.000001)/(pnum["hlc"]/
                                            (0.00008205*pnum["twater"]))
            else ldif <- 0
            ltotal <- ldep + ldif
            ltotal <- ltotal/1000

            ## trim out zeroes in ltotal
            incvec <- ltotal > 0
            ltotal <- ltotal[incvec]

            ## calculate start days for vvwm.
            ## first day of first year is 4/1/2005
            years <- seq(2005, length = length(ltotal), by = 1)
            day0 <- as.numeric(as.POSIXct(paste("4/1/", years, sep = ""),
                                          format = "%m/%d/%Y"))
            origin <- as.numeric(as.POSIXct("1/1/2005", format = "%m/%d/%Y"))
            ind <- round((day0 - origin)/3600/24) + 1

            varloc <- c("koc", "kaer", "temp_ref_aer_all", "kanaer",
                        "temp_ref_anae_all", "kpo", "kh", "hlc", "heat_of_henry")
            if (chemtype == "O") {
                parloc <- pnum[varloc]
                iskoc <- "TRUE"
                parloc["hlc"] <- round(parloc["hlc"]/(0.00008206*298.15),
                                       digits = 8)
            }
            else {
                if (chemtype =="M" | chemtype == "Hg") {
                    parloc <- rep(0, times = length(varloc))
                    names(parloc)<- varloc
                    iskoc <- "FALSE"
                    parloc["koc"] <- pnum["kd"]
                }
            }
            metfilename <- paste(getwd(), 'VVWM', 'MetData', sep = "/")

            d_over_dx <- pnum[paste("d_over_dx_", ii, sep = "")]
            if (ii==3) area_w <- pnum["area_3"]
            else area_w = pnum["area_fp"]
            depth_0 <- pnum[paste("dwc", ii, sep = "_")]
            depth_max <- pnum[paste("dwc", ii, sep = "_")]

            ## hard wired variables
            qt <- 2
            if (ii==3) {
                simtypeflag <- 3
                wbtype <- "Reservoir"
            }
            else {
                if (ii==5) {
                    simtypeflag <- 2
                    wbtype <- "Pond"
                }
            }

            split_path <- function(path) {
                if (dirname(path) %in% c(".", path)) return(basename(path))
                return(c(basename(path), split_path(dirname(path))))
            }

            ## change directory
            dirsav <- getwd()
            setwd(paste(dirsav, "vvwm", sep = "/"))

            ## write outinput file
            zz <- file(paste(getwd(), "outputs/VVWMInputBST.txt",sep = "/"), "w")
            file1 <- split_path(paste(getwd(), "outputs/VVWMInputBST", sep ="/"))
            file1 <- rev(file1)
            cat(file1, sep = "\\", file = zz)
            cat("\n", file = zz)
            cat(parm["chemname"], "\n", file = zz, sep = "")
            cat("1\n", file = zz, sep = "")
            cat(iskoc, "\n", file = zz, sep = "")
            for (i in 1:6) cat(parloc[i], "\n", file = zz,sep = "")
            cat(pnum["sitelatitude"], "\n", file= zz, sep ="")
            cat(parloc["kh"], "\n", file = zz, sep = "")
            for (i in 1:3) cat("\n", file = zz, sep = "")
            cat(pnum["mw"], "\n", file=zz, sep = "")
            for (i in 1:2) cat("0\n", file = zz, sep = "")
            for (i in 1:7) cat("\n", file = zz, sep = "")
            for (i in 8:9) cat(parloc[i], "\n", file = zz, sep ="")

            cat(qt, "\n", file = zz, sep = "")
            ## srctype id defined in main
            cat(c("Crop", "Pasture", "Reclamation")[srctypeid], "\n", file = zz, sep ="")
            ## siteid defined in main
            file1 <- rev(split_path(paste(getwd(), "/metdata/", siteid, ".wea", sep = "")))
            cat(file1, sep = "\\", file = zz)
            cat("\n", file = zz)
            cat(pnum["sitelatitude"], "\n", file = zz, sep = "")
            for (i in 1:2) cat("\n", file = zz, sep = "")
            cat(parm["burialflag"], "\n", file = zz, sep = "")
            for (i in 1:4) cat("\n", file = zz, sep = "")
            cat(d_over_dx, "\n", file = zz, sep = "")
            cat(paste(parm["is_calc_prben"], ", ", pnum["prben"], sep = ""),
                "\n", file = zz)
            varloc <- c("db", "bsp", "bulk_density", "foc_bs", "doc2", "bnmas",
                        "dfac", "tsswc", "chl", "froc1", "doc1", "plmas")

            for (i in 1:length(varloc)) cat(pnum[varloc[i]], "\n", file = zz, sep ="")
            for (i in 1:3) cat("\n", file = zz, sep = "")
            napp <- pnum["nappl"]*pnum["cutoffyr"]
            cat(napp, "\n", file = zz, sep ="")
            cat(ind, sep = ", ", file = zz)
            cat("\n", file = zz, sep ="")
            cat(simtypeflag, "\n", file = zz, sep = "")
            cat(pnum["area_1"], "\n", file = zz, sep = "")
            cat(area_w, "\n", file = zz, sep ="")
            cat(depth_0, "\n", file = zz, sep ="")
            cat(depth_max, "\n", file = zz, sep ="")
            cat(ltotal, sep = ", ", file = zz)
            cat("\n", file = zz)
            cat(pnum["flow_averaging"], "\n", file = zz, sep = "")
            cat(pnum["baseflow"], "\n", file = zz, sep ="")
            cat(pnum["laufrac"], "\n", file = zz, sep = "")
            cat(parm["is_hed_files_made"], "\n", file = zz, sep ="")
            cat(parm["is_add_return_frequency"], ", 0\n", sep = "", file = zz)
            file1 <- rev(split_path(paste(getwd(), "outputs",
                                          paste(wbtype, "daily.csv", sep = "_"),
                                          sep ="/")))
            cat(file1, sep = "\\", file = zz)
            cat("\n", file = zz)
            for (i in 1:2) cat(", ,\n", file = zz, sep = "")
            fname <- paste(getwd(), "outputs",
                           paste(wbtype, "parent.csv", sep = "_"),
                           sep = "/")
            file1 <- rev(split_path(fname))
            cat(file1, sep = "\\", file = zz)
            cat("\n", file = zz)
            for (i in 1:14) cat(", ,\n", file = zz, sep ="")
            close(zz)

            system("./vvwm.exe ./outputs/VVWMInputBST.txt")

            ## load results from vvwm
            ## assuming 150 year simulation
            dfvv <- read.fwf(fname, skip = 19, widths = c(4, rep(11, times = 17)),
                             n = 150)
            names(dfvv)[1] <- "year"

            names(dfvv)[10:18] <- c("cwctotann", "cwcdann", "cwbstotann",
                                    "cwctotmx1","cwctotmax4", "cwctotmx21",
                                    "cwcdmx1", "cwcdmx4", "cwcdmx21")

            numswyear <- pnum["numswyear"]

            vvwm_max[[paste(ii)]] <- apply(dfvv[1:numswyear, 2:ncol(dfvv)], 2, max)
            vvwm_mn[[paste(ii)]] <- apply(dfvv[1:numswyear, 2:ncol(dfvv)], 2, mean)
            vvwm_out[[paste(ii)]] <- dfvv
            setwd(dirsav)
        }

        return(list(ctda = src.list[["ctda"]], cvapor = cvapor,
                    dp = dp, dv = dv, ctss = src.list[["ctss"]],
                    vvwm_out = vvwm_out, vvwm_max = vvwm_max,
                    vvwm_mn = vvwm_mn, numyear = numyear,
                    cgwater = cgwater,
                    cair = cair))
    }
    else {
        return(list(cair = cair, cvapor = cvapor,
                    cgwater = cgwater))
    }

}

