## risk module

risk <- function(mout, foodout, pnum, chemtype) {

    getmean<- function(x, istart, istop, ip) {
        if (ip == 0) return(mean(x[istart:istop], na.rm = T)) else
        return(mean(x[istart:istop, ip], na.rm = T))
    }

    sa <- c(1,20) ## starting age, another possibility is 20
    salab <- c("child", "adult")
    cpick <- c(1,4) ## cohort selection for output

    ## get some max values from the media and food models
    tssmax <- max(mout[["ctss"]][,1], na.rm = T)
    tssmax2 <- max(mout[["ctss"]][,2], na.rm = T)

    if (tssmax > tssmax2) {
        idmax <- tssmax == mout[["ctss"]][,1]
    }
    else {
        idmax <- tssmax2 == mout[["ctss"]][,2]
    }
    idmax[is.na(idmax)] <- F

    ctssmaxyear <- which(idmax)
    foodout_max <- sapply(foodout, max, na.rm = T)

    gwatermax <- max(mout[["cgwater"]], na.rm = T)
    swatermax <- max(mout[["vvwm_out"]][["3"]][, "cwcdann"], na.rm = T)
    cair.max <- max(mout[["cair"]][,2], na.rm = T)

    for (kk in 1:length(sa)) { # loop through two starting ages
        intcohort <- as.numeric(cut(sa[kk], c(0,6,12,20,100), right = F))

        ed <- pnum[paste("ed_", intcohort, sep = "")]
        sy <- round(ctssmaxyear - 0.5*ed)
        sy <- max(sy, 1)
        if (sy > 150) sy <- 150-ed

        ey <- sy + ed - 1
        cindoor <- shower(sy, ed, pnum, mout[["cgwater"]], kk, chemtype)

        ## add shower module at some point!
        cair.avg <- getmean(mout[["cair"]], sy, ey, 2)
        cindoor.avg <- getmean(cindoor, sy, ey, 0)

        cohort_ed2 <- function(sy, ed, sy_cohort, ey_cohort) {
            ey <- sy + ed -1
            if ((sy > ey_cohort) | (ey < sy_cohort)) return(0)
            else {
                istart <- ifelse(sy > sy_cohort, sy, sy_cohort)
                iend <- ifelse(ey < ey_cohort, ey, ey_cohort)
                val <- iend - istart + 1
                return(val)
            }
        }

        mlist <- c("cgwater", "ctss")
        vlist <- "cwcdann"
        flist <- c("p_root",
                   "p_exfruit", "p_profruit", "p_exveg",
                   "p_proveg", "a_beef", "a_milk", "c_fish_t3f",
                   "c_fish_t4f")
        mnout <- rep(0, times = length(mlist) + length(vlist) + length(flist))
        names(mnout) <- c(mlist, vlist, flist)
        mnall <- rep(0, times = length(mlist) + length(vlist) + length(flist))
        names(mnall) <- c(mlist, vlist, flist)

        ey <- sy + ed-1
        endyrexp <- sy + ed - 1
        cohort_sy <- sy
        ## intakes stores the intakes for each cohort stage (1-4)
        intakes <- matrix(0, ncol = 11, nrow = 4)
        foodveg <- c("root", "proveg", "exveg", "profruit", "exfruit")
        foodan <- c("beef", "milk")
        dimnames(intakes)[[2]] <- c(foodveg, foodan, "gwater", "swater",
                                    "soil", "fish")
        intakes.max <- matrix(0, ncol = 11, nrow = 4)
        dimnames(intakes.max)[[2]] <- c(foodveg, foodan, "gwater", "swater",
                                        "soil", "fish")

        c_val <- matrix(0, ncol  = 11, nrow = 2)
        dimnames(c_val)<- list(c("avg", "max"),
                               c(foodveg, foodan, "gwater", "swater",
                                 "soil", "fish"))

        ## fname should be sorted in the same order as dimnames for c_avg
        ## and intakes (except for fish)
        fname <- c("p_root", "p_proveg", "p_exveg", "p_profruit", "p_exfruit",
                   "a_beef", "a_milk", "cgwater", "cwcdann", "ctss")

        ## calculate grand means
        ip <- c(0, 1)  # flag to tell getmean whether a single vector
        ## or a matrix is coming for calculating the mean
        for (i in 1:length(mlist)) {
            mnall[mlist[i]] <- getmean(mout[[mlist[i]]], sy, ey, ip[i])
        }
        mnall[3] <- getmean(mout[["vvwm_out"]][["3"]][, vlist], sy, ey, 0)
        for (i in 1:length(flist)) {
            mnall[flist[i]] <- getmean(foodout[[flist[i]]], sy, ey, 0)
        }
        mnall[is.na(mnall)] <- 0

        cohort_edval <- rep(0, times = 4)
        for (cohort in intcohort:4) {
            cohort_sa <- c(1,6,12,20)[cohort]
            cohort_ea <- c(5,11,19,70)[cohort]
            cohort_edval[cohort] <- cohort_ed2(sa[kk], ed, cohort_sa, cohort_ea)
            cohort_ey <- cohort_sy + cohort_edval[cohort] - 1

            if (cohort_edval[cohort] > 0) {
                ip <- c(0, 1)
                for (i in 1:length(mlist)) {
                    mnout[mlist[i]] <- getmean(mout[[mlist[i]]], cohort_sy,
                                               cohort_ey, ip[i])
                }
                mnout[3] <- getmean(mout[["vvwm_out"]][["3"]][, vlist], cohort_sy,
                                    cohort_ey, 0)
                for (i in 1:length(flist)) {
                    mnout[flist[i]] <- getmean(foodout[[flist[i]]], cohort_sy,
                                               cohort_ey, 0)
                }
            }

            mnout[is.na(mnout)] <- 0

            ## veg
            L <- pnum[paste("l_", foodveg, sep = "")]
            CR <- pnum[paste("cr_", foodveg, "_", cohort, sep = "")]
            f <- pnum[paste("f_", foodveg, sep = "")]
            intakes[cohort, 1:5] <- CR*f/1000*mnout[fname[1:5]]*(1-L)
            intakes.max[cohort,1:5] <- CR*f/1000*foodout_max[fname[1:5]]*(1-L)

            ## meat and dairy
            CRi <- pnum[paste("cr_", foodan, "_", cohort, sep = "")]
            Fi <- pnum[paste("f_", foodan, sep = "")]
            L2 <- pnum[paste("l2_", foodan, sep = "")]
            L1 <- pnum[paste("l1_", foodan, sep = "")]
            intakes[cohort, 6:7] <- mnout[fname[6:7]]*CRi/1000*Fi*(1-L1)*(1-L2)
            intakes.max[cohort,6:7] <- foodout_max[fname[6:7]]*CRi/1000*Fi*(1-L1)*(1-L2)

            ## Water
            CR_dw <- pnum[paste("cr_dw_", cohort, sep = "")]
            intakes[cohort, 8:9] <- mnout[fname[8:9]]*CR_dw*pnum["f_dw"]/1000
            intakes.max[cohort, 8] <- gwatermax*CR_dw*pnum["f_dw"]/1000
            intakes.max[cohort, 9] <- swatermax*CR_dw*pnum["f_dw"]/1000

            ## soil
            intakes[cohort, 10] <- mnout[fname[10]]*pnum[paste("crs_", cohort, sep = "")]*pnum["f_soil"]/pnum[paste("bw_", cohort, sep = "")]*0.000001
            intakes.max[cohort, 10] <- tssmax*pnum[paste("crs_", cohort, sep = "")]*pnum["f_soil"]/pnum[paste("bw_", cohort, sep = "")]*0.000001

            ## fish
            fish_avg <- pnum["f_t3"]*mnout["c_fish_t3f"] +
                pnum["f_t4"]*mnout["c_fish_t4f"]
            intakes[cohort,11] <- fish_avg*pnum[paste("cr_fish_", cohort, sep = "")]*pnum["f_fish"]/(1000*pnum[paste("bw_", cohort, sep = "")])

            fish_max <- pnum["f_t3"]*foodout_max["c_fish_t3f"] +
                pnum["f_t4"]*foodout_max["c_fish_t4f"]
            intakes.max[cohort,11] <- fish_max*pnum[paste("cr_fish_", cohort, sep = "")]*pnum["f_fish"]/(1000*pnum[paste("bw_", cohort, sep = "")])

            cohort_sy <- cohort_sy + cohort_edval[cohort]
        }
        ## get mean intake
        intakes.mn <- apply(cohort_edval * intakes, 2, sum,
                            na.rm = T)/sum(cohort_edval)

        ## sort grand means and maxes correctly
        c_val[1,1:10] <- mnall[fname]
        c_val[2,1:7] <- foodout_max[fname[1:7]]
        c_val[2,8:11] <- c(gwatermax, swatermax, tssmax, fish_max)
        fish_avg <- pnum["f_t3"]*mnall["c_fish_t3f"] +
            pnum["f_t4"]*mnall["c_fish_t4f"]
        c_val[1,11] <- fish_avg

        efi <- pnum["ef"]
        ifelse(ed + sa[kk]>70, at <- ed, at <- pnum["at"])
        ladd <- (intakes.mn*ed*efi)/(at*365)

        ## ladd_sum <- sum(ladd)
        ## no surface water included in sum?! strange
        ladd_sum <- sum(ladd[-which(names(ladd) == "swater")])

        ## use pathway specific RFDs if needed
        rfd <- rep(NA, times = 4)
        names(rfd) <- c("fish", "food", "soil", "water")
        rfd <- pnum[paste("rfd_", names(rfd), sep = "")]
        rfd[rfd==0] <- NA
        ## choose pathway RFD here! need to pass as a parameter
        pathwayRFD <- FALSE
        if (! pathwayRFD) {
            rfd <- rep(pnum["rfd"], times = length(rfd))
        }
        names(rfd) <- c("fish", "food", "soil", "water")

        quot.max <- matrix(NA, ncol =11, nrow = 4)
        quot.max[,1:7] <- intakes.max[, 1:7]/rfd["food"]
        quot.max[,8] <- intakes.max[,8]/rfd["water"]
        quot.max[,9] <- intakes.max[,9]/rfd["water"]
        quot.max[,10] <- intakes.max[,10]/rfd["soil"]
        quot.max[,11] <- intakes.max[,11]/rfd["fish"]

        ## exclude surface water when getting total ingestion
        quot.maxsum <- apply(quot.max[, c(rep(T, times = 8), F,T)],
                             1, sum, na.rm = T)

        airquotmax <- cair.max/pnum["rfc"]
        airriskavg <- cair.avg*pnum["iur"]*1000*ed/pnum["at"]

        cindoor.max <- max(cindoor)
        cindoor.quot <- cindoor.max/pnum["rfc"]
        cindoorriskavg <- cindoor.avg*pnum["iur"]*1000*ed/pnum["at"]

        dftemp2 <- data.frame(pathway = dimnames(intakes.max)[[2]],
                             conc =c_val[2,],
                             dose = intakes.max[cpick[kk],],
                             risk = quot.max[cpick[kk],])
        dftemp <- data.frame(pathway = "all",
                             conc = 0,
                             dose = sum(intakes.max[cpick[kk],]),
                             risk = quot.maxsum[cpick[kk]])
        dftemp <- rbind(dftemp, dftemp2)
        dftemp <- rbind(dftemp, data.frame(pathway = c("air","showerair"),
                                           conc = c(cair.max,cindoor.max),
                                           dose = c(0,0),
                                           risk = c(airquotmax,cindoor.quot)))
        dftemp$receptor <- salab[kk]
        dftemp$endpoint <- "Noncancer"
        dftemp <- dftemp[,c("receptor", "pathway", "endpoint",
                           "risk", "dose", "conc")]

        if (kk == 1) dfrisk <- dftemp else dfrisk <- rbind(dfrisk, dftemp)

        if (pnum["csforal"] > 0) {
            dftemp <- data.frame(pathway = names(ladd),
                                 conc = c_val[1,],
                                 dose = as.vector(ladd),
                                 risk = ladd*pnum["csforal"])
            dftemp2 <- data.frame(pathway = "all",
                                 conc = 0,
                                 dose = ladd_sum,
                                 risk = ladd_sum*pnum["csforal"])
            dftemp <- rbind(dftemp2, dftemp)
            dftemp <- rbind(dftemp, data.frame(pathway = c("air","showerair"),
                                               conc = c(cair.avg,cindoor.avg),
                                               dose = c(0,0),
                                               risk = c(airriskavg,cindoorriskavg)))
            dftemp$receptor <- salab[kk]
            dftemp$endpoint <- "Cancer"
            dftemp <- dftemp[, c("receptor", "pathway", "endpoint",
                                 "risk", "dose", "conc")]
            dfrisk <- rbind(dfrisk, dftemp)
        }

    }
    options(scipen =0)
    dfrisk <- dfrisk[order(dfrisk$endpoint, dfrisk$pathway, dfrisk$receptor),]
    rownames(dfrisk) <- paste(1:nrow(dfrisk))

    return(dfrisk)
}
