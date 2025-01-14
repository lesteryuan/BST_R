## run eco receptors
eco <- function(cas, mout, foodout, eco_const) {
    tssmax <- apply(mout[["ctss"]], 2, max, na.rm = T)
    foodmax <- sapply(foodout, max, na.rm = T)

    receptors <- read.csv("./data/tbl_ecoreceptors.csv")
    rectype <- read.csv("./data/tbl_ecoreceptortype.csv")

    receptors <- fnames.un(receptors)
    rectype <- fnames.un(rectype)
    receptors <- merge(receptors, rectype, by = "ecorecid")

    append.results <- function(df, ecorecid, intpath, hbid, routeid,
                               cprey, dosepath, risk) {

        dftemp <- data.frame(ecorecid, intpath, hbid, routeid, cprey, dosepath,
                             risk)
        incvec <- dftemp$risk >0
        incvec[is.na(incvec)] <- FALSE
        if (sum(incvec) > 0) {
            df <- rbind(dftemp[incvec,], df)
        }
        return(df)
    }
    ## initialize dataframe for storing resulst
    ecoriskout <- data.frame(ecorecid = numeric(0),
                             path = character(0), hbid = character(0),
                             routeid = numeric(0), cprey = numeric(0),
                             dosepath = numeric(0), risk = numeric(0))

    pathname <- c("Worms","Other soil inverts", "Small mammals",
                  "Herbivorous vertebrates", "Omnivorous vertebrates",
                  "Small birds", "Benthic filter feeders",
                  "T3 fish", "T4 fish", "Aquatic plants",
                  "Exposed fruits", "Exposed vegetables",
                  "Forage", "Grains", "Roots", "Silage",
                  "Soil", "Sediment", "Water",
                  "Small herpetofauna", "Sediment (dissolved)",
                  "Water (dissolved)", "Total Ingestion",
                  "Total food diet")

    ## loop through receptors (all selected for now)
    for (irec in receptors$ecorecid) {
        incvec <- irec == receptors$ecorecid

        ## select correct eco_const
        incvec <- eco_const$ecorecid == irec
        incvec[is.na(incvec)] <- T
        incvec2 <- eco_const$cas == cas
        incvec2[is.na(incvec2)] <- T
        if (sum(incvec & incvec2) > 0) {
            ecop <- eco_const[incvec & incvec2,"value"]
            names(ecop) <- eco_const[incvec & incvec2, "modelcode"]
            if (sum(duplicated(names(ecop))) > 0) stop("more than one code selected")
        }

        crdiet <- ecop["cr_food"]
        if (is.na(crdiet)) crdiet < 0
        crwater <- ecop["cr_water"]

        ## set missing parameters to zero
        plist <- c("bmc_sed", "bmc_soil", "bmc_water", "bmd", "ed_eco",
                   "crfrac_sed", "crfrac_soil")
        for (i in plist) {
            if (is.na(ecop[i])) {
                ecop <- c(0, ecop)
                names(ecop)[1] <- i
            }
        }

        rectype <- receptors$rectype[receptors$ecorecid == irec]

        ## run risk model if there's enough parameter info
        if (rectype == "D" & ecop["bmd"] > 0) {

            orgtype <- c("worms", "soilinvert", "smmammals", "herbvert",
                         "omnvert", "smbirds", "smherp")
            baf <- ecop[paste("baf_", orgtype, sep = "")]
            baf[is.na(baf)] <- 0
            ctype <- baf*tssmax[1]
            names(ctype) <- orgtype
            cprey <- c(ctype[1:6], foodmax[c("c_fish_t3w", "c_fish_t3w",
                                             "c_fish_t4w", "c_fish_t3w",
                                             "p_exfruit", "p_exveg", "p_forage",
                                             "p_grain", "p_root","p_silage")],
                       NA, NA, NA,  ctype[7])
            dfrac <- ecop[paste("df_", 1:20, sep ="")]


            dfrac[is.na(dfrac)] <- 0
            cdiet <- cprey*dfrac
            pathway <- crdiet*cdiet/ecop["bw"]
            rechq <- pathway/ecop["bmd"]

            ## select paths to include
            selvec <- c(rep(T, times = 16), F,F,F,T)
            ecoriskout <- append.results(ecoriskout, irec,
                                         pathname[c(1:16,20)], rep("EC", times = 17),
                                         rep(2, times = 17), cprey[selvec],
                                         pathway[selvec], rechq[selvec])

            dosediet <- sum(pathway, na.rm = T)

            ecoriskout <- append.results(ecoriskout,irec,
                                         pathname[24], "EC", 2, 0, dosediet,
                                         dosediet/ecop["bmd"])
            dosemedia <- NA
            if (ecop["crfrac_sed"] > 0) {
                dosemedia <- crdiet*mout[["vvwm_mn"]][["5"]]["cwbstotann"]*
                    ecop["crfrac_sed"]/ecop["bw"]
                ecoriskout <- append.results(ecoriskout, irec, pathname[18], "EC", 2,
                                            mout[["vvwm_mn"]][["5"]]["cwbstotann"],
                                             dosemedia, dosemedia/ecop["bmd"])

            }
            if (ecop["crfrac_soil"] >0) {
                dosemedia <- crdiet*tssmax[1]*ecop["crfrac_soil"]/ecop["bw"]
                ecoriskout <- append.results(ecoriskout,irec, pathname[17], "EC", 2,
                                             tssmax[1], dosemedia, dosemedia/ecop["bmd"])
            }

            dowater <- ecop["bw"] >0
            dowater[is.na(dowater)] <- F
            if (dowater) {
                dosewater <- crwater*mout[["vvwm_mn"]][["5"]]["cwctotann"]/ecop["bw"]
                ecoriskout <- append.results(ecoriskout, irec, pathname[19], "EC", 2,
                                             mout[["vvwm_mn"]][["5"]]["cwctotann"],
                                             dosewater, dosewater/ecop["bmd"])
            }
            else dosewater <- 0

            dosetotal <- dosediet + dosemedia + dosewater
            if (!is.na(ecop["bmd"])) {
                if (ecop["bmd"] >0) {
                    ecoriskout <- append.results(ecoriskout, irec, pathname[23], "EC", 2,
                                                 0, dosetotal, dosetotal/ecop["bmd"])
                }
            }

        }
        else {
            if (rectype == "C" & ((ecop["bmc_water"] > 0&ecop["ed_eco"] > 0) |
                                      (ecop["bmc_sed"] > 0|ecop["bmc_soil"] > 0))) {
                bmc_soil <- ecop["bmc_soil"] > 0
                bmc_soil[is.na(bmc_soil)] <- F
                if (bmc_soil) {
                    rechq <- tssmax[1]/ecop["bmc_soil"]
                    ecoriskout <- append.results(ecoriskout,irec,
                                                 pathname[17], "EC", 1, 0,tssmax[1],
                                                 rechq)
                }
                bmc_sed <- ecop["bmc_sed"] > 0
                bmc_sed[is.na(bmc_sed)] <- F
                if (bmc_sed) {
                    rechq <- mout[["vvwm_mn"]][["5"]]["cwbstotann"]/
                        ecop["bmc_sed"]
                    ecoriskout <- append.results(ecoriskout, irec,
                                                pathname[18], "EC", 1, 0,
                                                mout[["vvwm_mn"]][["5"]]["cwbstotann"],rechq)
                }
                cmax <- mout[["vvwm_max"]][["5"]]
                if (ecop["bmc_water"] > 0) {
                    if (!is.na(ecop["ed_eco"])) {
                        if (ecop["ed_eco"] <= 1) tempc <- cmax["cwctotmx1"]
                        else {
                            if (ecop["ed_eco"] <= 4) tempc <- cmax["cwctotmax4"]
                            else {
                                if (ecop["ed_eco"] <= 21) tempc <- cmax["cwctotmx21"]
                                else {
                                    tempc <- mout[["vvwm_mn"]][["5"]]["cwctotann"]
                                }
                            }
                        }
                        rechq <- tempc/ecop["bmc_water"]
                        if (rechq >= 1) {
                            ## do short term exceedances
                        }
                        else exceed <- 0
                        ecoriskout <- append.results(ecoriskout, irec, pathname[19], "EC",
                                                     1, tempc, exceed, rechq)
                        if (ecop["ed_eco"] <= 1) tempc <- cmax["cwcdmx1"]
                        else {
                            if (ecop["ed_eco"] <= 4) tempc <- cmax["cwcdmx4"]
                            else {
                                if (ecop["ed_eco"] <= 21) tempc <- cmax["cwcdmx21"]
                                else {
                                    tempc <- mout[["vvwm_mn"]][["5"]]["cwcdann"]
                                }
                            }
                        }
                        rechq <- tempc/ecop["bmc_water"]
                        if (rechq >= 1) {
                            ## do short term exceedances
                        }
                        else exceed <- 0
                        ## strange. Dissolved concentrations are used.
                    }
                }
            }
        }
    }

    ecoriskout <- merge(ecoriskout,
                        unique.data.frame(receptors[, c("ecorecid", "receptor.name")]),
                        by = "ecorecid")


    return(ecoriskout)

}

