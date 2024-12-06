## risk module for SI

risk_si<- function(mout, pnum, chemtype) {

    getmean<- function(x, istart, istop, ip) {
        if (ip == 0) return(mean(x[istart:istop])) else
        return(mean(x[istart:istop, ip]))
    }

    sy <- 1
    sa <- c(1,20) ## starting age, another possibility is 20
    salab <- c("child", "adult")
    cpick <- c(1,4) ## cohort selection for output

    ## get some max values from the media and food models
    gwatermax <- max(mout[["cgwater"]])

    for (kk in 1:length(sa)) {
        intcohort <- as.numeric(cut(sa[kk], c(0,6,12,20,100), right = F))

        ed <- pnum[paste("ed_", intcohort, sep = "")]
        ey <- sy + ed - 1

        cindoor <- shower(sy,ed, pnum, mout[["cgwater"]], kk, chemtype)

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

        ey <- sy + ed-1
        endyrexp <- sy + ed - 1
        cohort_sy <- sy

        ## summary stats for gwater
        intakes <- rep(0, times = 4)
        intakes.max <- rep(0, times = 4)

        ## calculate grand means
        mnall <- getmean(mout[["cgwater"]], sy, ey, 0)
        mnall[is.na(mnall)] <- 0
        cair.avg <- getmean(mout[["cair"]], sy, ey, 2)
        cindoor.avg <- getmean(cindoor, sy, ey, 0)

        cohort_edval <- rep(0, times = 4)
        for (cohort in intcohort:4) {
            cohort_sa <- c(1,6,12,20)[cohort]
            cohort_ea <- c(5,11,19,70)[cohort]
            cohort_edval[cohort] <- cohort_ed2(sa[kk], ed, cohort_sa, cohort_ea)
            cohort_ey <- cohort_sy + cohort_edval[cohort] - 1

            if (cohort_edval[cohort] > 0) {
                mnout <- getmean(mout[["cgwater"]], cohort_sy, cohort_ey, 0)
                mnout[is.na(mnout)] <- 0
            }
            ## Water
            CR_dw <- pnum[paste("cr_dw_", cohort, sep = "")]
            intakes[cohort] <- mnout*CR_dw*pnum["f_dw"]/1000
            intakes.max[cohort] <- gwatermax*CR_dw*pnum["f_dw"]/1000

            cohort_sy <- cohort_sy + cohort_edval[cohort]
        }
        ## get mean intake
        intakes.mn <- sum(cohort_edval * intakes)/sum(cohort_edval)

        efi <- pnum["ef"]
        ifelse(ed + sa[kk]>70, at <- ed, at <- pnum["at"])
        ladd <- (intakes.mn*ed*efi)/(at*365)

        ## use pathway specific RFDs if needed. set to false for now
        if (FALSE) {
            rfd <- pnum["rfd_water"]
        }
        else {
            rfd <- pnum["rfd"]
        }
        rfd[rfd==0] <- NA

        quot.max <- intakes.max/rfd

        cindoor.max <- max(cindoor)
        cindoor.quot <- cindoor.max/pnum["rfc"]

        cair.max <- max(mout[["cair"]][sy:endyrexp,2])
        airquotmax <- cair.max/pnum["rfc"]
        airriskavg <- cair.avg*pnum["iur"]*1000*ed/pnum["at"]
        cindoorriskavg <- cindoor.avg*pnum["iur"]*1000*ed/pnum["at"]

        dftemp <- data.frame(pathway = "gwater",
                             conc = gwatermax,
                             dose = intakes.max[intcohort],
                             risk = quot.max[intcohort])
        dftemp2 <- data.frame(pathway = c("air","showerair"),
                              conc = c(cair.max,cindoor.max),
                              dose = c(0,0),
                              risk = c(airquotmax, cindoor.quot))

        dftemp <- rbind(dftemp, dftemp2)
        dftemp$receptor <- salab[kk]
        dftemp$endpoint <- "Noncancer"
        dftemp <- dftemp[,c("receptor", "pathway", "endpoint",
                           "risk", "dose", "conc")]

        if (kk == 1) dfrisk <- dftemp else dfrisk <- rbind(dfrisk, dftemp)

        if (pnum["csforal"] > 0) {
            dftemp <- data.frame(pathway = "gwater",
                                 conc = mnall,
                                 dose = ladd,
                                 risk = ladd*pnum["csforal"])
            dftemp <- rbind(dftemp, data.frame(pathway = c("air","showerair"),
                                               conc = c(cair.avg,cindoor.avg),
                                               dose = c(0,0),
                                               risk = c(airriskavg,
                                                   cindoorriskavg)))

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
