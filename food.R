## food module

food <- function(parm, pnum, mout, srctypeid) {

    cfish <- function(fishtype, parm, pnum, cwctot, cwd) {
        if (parm["chemtype"] == "M") {
            cfish <- cwctot*pnum[paste("bcf_",fishtype, sep = "")]
        }
        else {
            cfish <- cwd*pnum[paste("bcf_", fishtype, sep = "")]
        }
        return(cfish)
    }

    calc_acsoil <- function(parm, pnum, product, csoil) {
        acsoil <- pnum[paste("df_soil_", product, sep = "")]*
            pnum["bs"]*csoil
        return(acsoil)
    }
    calc_apveg <- function(parm, pnum, product, vegsource, pveg_wet) {
        pveg_dry <- pveg_wet/((100 - pnum[paste("maf_", vegsource, sep = "")])/100)
        apveg <- pnum[paste("df_", vegsource, "_", product, sep = "")]*
            pveg_dry*pnum[paste("f_", vegsource, sep = "")]
        return(apveg)
    }

    calc_risk_prbg <- function(parm, pnum, csoil) {
        if ((parm["chemtype"] == "M") | (parm["chemtype"] == "Hg")) {
            prbg <- csoil*pnum["brroot"]*(100 - pnum["maf_root"])/100
        }
        else {
            if ((parm["chemtype"] == "O") | (parm["chemtype"] == "D")) {
                kd <- pnum["koc"]*pnum["foc_soil"]
            }
            else kd <- pnum["kd"]
            prbg <- csoil*pnum["rcf"]*pnum["vg_root"]/kd
        }

        return(prbg)
    }

    calc_risk_pveg <- function(vegsource, vegtype,parm, pnum, csoil, cvapor,
                           dp, dv) {
        ## only calculate for 1st column of each of the matrices (L=1)
        if (parm["chemtype"]== "D") {
            if (vegtype == "protected") {
                pdveg <- 0
                pvveg <- 0
            }
            else {
                pvveg = cvapor*pnum["bv"]*
                    pnum[paste("vg_", vegsource, sep = "")]*1000/1200
                ## pd_veg Need to add!!
            }
        }
        else {
            if (vegtype == "protected") {
                pdveg <- 0
                pvveg <- 0
                prveg <- pnum[paste("br", vegsource, sep = "")]*csoil
            }
            else {
                if (vegtype == "exposed") {
                    if (parm["chemtype"] == "O") {
                        pvveg = cvapor*pnum["bv"]*pnum[paste("vg_", vegsource, sep = "")]*1000/1200
                    }
                    else {
                        pvveg <- 0
                    }
                    pdveg <- (dp*pnum[paste("rp_",vegsource, sep = "")])*
                        (1 - exp(-pnum["kppar"]*pnum[paste("tp_", vegsource, sep = "")]))/
                            (pnum[paste("yp_", vegsource, sep = "")]*pnum["kppar"])
                    prveg <- pnum[paste("br", vegsource, sep = "")]*csoil
                }
            }
        }
        pveg <- (pdveg + pvveg + prveg)*
            (100 - pnum[paste("maf_", vegsource, sep = "")])/100
    }


    if ((srctypeid == 1) | (srctypeid == 4)) {
        p_forage_2 <- 0
        p_silage_1 <- 0
        p_grain_1 <- 0
        p_proveg_1 <- calc_risk_pveg("proveg", "protected", parm, pnum,
                                     mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                     mout[["dp"]][,1], mout[["dv"]][,1])
        p_profruit_1 <- calc_risk_pveg("profruit", "protected", parm, pnum,
                                       mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                       mout[["dp"]][,1], mout[["dv"]][,1])
        p_exveg_1 <- calc_risk_pveg("exveg", "exposed", parm, pnum,
                                    mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                    mout[["dp"]][,1], mout[["dv"]][,1])
        p_exfruit_1 <- calc_risk_pveg("exfruit", "exposed", parm, pnum,
                                      mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                      mout[["dp"]][,1], mout[["dv"]][,1])
        p_root_1 <- calc_risk_prbg(parm, pnum, mout[["ctda"]][,1])
        p_forage <- 0
        p_grain <- calc_risk_pveg("grain", "protected", parm, pnum,
                                  mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                  mout[["dp"]][,1], mout[["dv"]][,1])
        p_silage <- calc_risk_pveg("silage", "exposed", parm, pnum,
                                   mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                   mout[["dp"]][,1], mout[["dv"]][,1])
        p_proveg <- calc_risk_pveg("proveg", "protected", parm, pnum,
                                   mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                   mout[["dp"]][,1], mout[["dv"]][,1])
        p_profruit <- calc_risk_pveg("profruit", "protected", parm, pnum,
                                     mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                     mout[["dp"]][,1], mout[["dv"]][,1])
        p_exveg <- calc_risk_pveg("exveg", "exposed", parm, pnum,
                                  mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                  mout[["dp"]][,1], mout[["dv"]][,1])
        p_exfruit <- calc_risk_pveg("exfruit", "exposed", parm, pnum,
                                    mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                    mout[["dp"]][,1], mout[["dv"]][,1])
        p_root <- calc_risk_prbg(parm, pnum, mout[["ctda"]][,1])
        acsoil <- 0
        acforage <- 0
        acgrain <- 0
        acsilage <- 0
        a_beef <- 0
        a_milk <- 0
    }
    else {
        p_forage_2 <- calc_risk_pveg("forage", "exposed",parm, pnum,
                                     mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                     mout[["dp"]][,1], mout[["dv"]][,1])
        p_grain_1 <- 0
        p_silage_1 <- calc_risk_pveg("silage", "exposed", parm, pnum,
                                     mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                     mout[["dp"]][,1], mout[["dv"]][,1])
        p_proveg_1 <- 0
        p_profruit_1 <- 0
        p_exveg_1 <- 0
        p_exfruit_1 <- 0
        p_root_1 <- 0
        p_forage <- calc_risk_pveg("forage", "exposed",parm, pnum,
                                   mout[["ctda"]][,1], mout[["cvapor"]][,1],
                                   mout[["dp"]][,1], mout[["dv"]][,1])
        p_grain <- 0
        p_silage <- 0
        p_proveg <- 0
        p_profruit <- 0
        p_exveg <- 0
        p_exfruit <- 0
        p_root <- 0
        ## beef calculation
        acsoil <- calc_acsoil(parm, pnum, "beef", mout[["ctss"]][,1])
        acforage <- calc_apveg(parm, pnum, "beef", "forage", p_forage_2)
        acgrain <- 0
        acsilage <- calc_apveg(parm, pnum, "beef", "silage", p_silage_1)
        a_beef <- pnum[paste("bcf_", "beef", sep = "")]*
            (acforage + acsilage + acgrain +acsoil)
        ## milk calculation
        acsoil <- calc_acsoil(parm, pnum, "milk", mout[["ctss"]][,1])
        acforage <- calc_apveg(parm, pnum, "milk", "forage", p_forage_2)
        acgrain <- 0
        acsilage <- calc_apveg(parm, pnum, "milk", "silage", p_silage_1)
        a_milk <- pnum[paste("bcf_", "milk", sep = "")]*
            (acforage + acsilage + acgrain + acsoil)
    }

    ## fish concentrations
    c_fish_t3f <- cfish("t3f", parm, pnum,
                        mout[["vvwm_out"]][["5"]][, "cwctotann"],
                        mout[["vvwm_out"]][["5"]][, "cwcdann"])
    c_fish_t4f <- cfish("t4f", parm, pnum,
                        mout[["vvwm_out"]][["5"]][, "cwctotann"],
                        mout[["vvwm_out"]][["5"]][, "cwcdann"])
    c_fish_t3w <- cfish("t3w", parm, pnum,
                        mout[["vvwm_out"]][["5"]][, "cwctotann"],
                        mout[["vvwm_out"]][["5"]][, "cwcdann"])
    c_fish_t4w <- cfish("t4w", parm, pnum,
                        mout[["vvwm_out"]][["5"]][, "cwctotann"],
                        mout[["vvwm_out"]][["5"]][, "cwcdann"])
    c_fish_t3fmax <- cfish("t3f", parm, pnum,
                        mout[["vvwm_max"]][["5"]]["cwctotann"],
                        mout[["vvwm_max"]][["5"]]["cwcdann"])
    c_fish_t4fmax <- cfish("t4f", parm, pnum,
                        mout[["vvwm_max"]][["5"]]["cwctotann"],
                        mout[["vvwm_max"]][["5"]]["cwcdann"])

    foodout <- list(p_forage = p_forage, p_grain = p_grain,
                    p_silage = p_silage, p_proveg = p_proveg,
                    p_profruit = p_profruit, p_exveg = p_exveg,
                    p_exfruit = p_exfruit, p_root = p_root,
                    a_milk = a_milk, a_beef = a_beef,
                    c_fish_t3f = c_fish_t3f,
                    c_fish_t4f = c_fish_t4f,
                    c_fish_t3w = c_fish_t3w,
                    c_fish_t4w = c_fish_t4w)

    ## create teq food concentrations if necessary
    if (parm["chemtype"] == "D") {
        foodout.teq <- as.list(rep(NA, times = length(foodout)))
        names(foodout.teq) <- paste(names(foodout), "teq", sep = "")
        for (i in 1:length(foodout))
            foodout.teq[[i]] <- foodout[[i]]*pnum["tef"]
        foodout <- c(foodout, foodout.teq)
    }

    return(foodout)
}





