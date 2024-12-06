
## shower module
shower <- function(sy,ed, pnum, cgwater,kk, chemtype) {

    ey <- sy + pnum["ed_4"]-1
    endyrexp <- sy + ed -1

    cin <- cgwater[round(sy):round(endyrexp)]

    if (any(cin > 0) & kk == 2 & chemtype != "M") {

        ## weird. The coefficient is rounded.
        hprime <- round(1/(0.00008205*298))*pnum["hlc"]
        ifelse (hprime > 0,
                kol <- 216/((2.5/pnum["dw"]^0.667) + (1/((pnum["da"]^0.667)*hprime))),
                kol <- 0)
        n <- kol*6/pnum["dropdiam"]*pnum["nozheight"]*100/pnum["dropvel"]
        bsrestime <- pnum["showertime"] + pnum["showerstalltime"] +
            pnum["t_bathroom"]
        showerrestime <- pnum["showertime"] + pnum["showerstalltime"]

        ts <- 0.2
        nc <- length(cin)
        ys <- rep(0, times = nc)
        yb <- rep(0, times = nc)
        ys_sum <- rep(0, times = nc)
        yb_sum <- rep(0, times = nc)
        intshowercount <- rep(0, times = nc)
        intbathcount <- rep(0, times=nc)
        fsat <- rep(0, times = nc)

        for (time in seq(0, bsrestime, by = ts)) {

            if (time < pnum["showertime"]) {
                emax <- (hprime*cin - ys)*pnum["vs"]*1000
                fem <- (1-fsat)*(1-exp(-n))
                et <- cin*pnum["showerrate"]*ts*fem
                incvec <- et > emax
                es <- et
                es[incvec] <-emax[incvec]
            }
            else es <- rep(0, times = nc)

            ys_new<- ys + ((es - pnum["qsb"]*(ys - yb)*ts)/(pnum["vs"]*1000))
            yh <- 0
            yb_new<- yb + (pnum["qsb"]*(ys_new - yb) - pnum["qbh"]*(yb- yh))*ts/
                (pnum["vb"]*1000)
            y_eq <- hprime*cin
            ifelse(y_eq >0, fsat <- ys_new/y_eq, fsat <- 0)

            if (time < showerrestime) {
                ys_avg <- (ys + ys_new)*0.5
                ys_sum <- ys_sum + ys_avg
                intshowercount <- intshowercount + 1
            }
            else {
                if (time >= showerrestime & time < bsrestime) {
                    yb_avg <- (yb + yb_new)*0.5
                    yb_sum <- yb_sum + yb_avg
                    intbathcount <- intbathcount + 1
                }
            }
            ys <- ys_new
            yb <- yb_new
        }
        cair_shower <- ys_sum*1000/intshowercount
        ifelse (intbathcount == 0, cair_bathroom <- 0, cair_bathroom <- yb_sum*1000/intbathcount)
        cindoor <- (cair_shower*showerrestime +
                        cair_bathroom*pnum["t_bathroom"])/1440
    }
    else cindoor <- rep(0, times = length(cin))
    return(cindoor)

}



