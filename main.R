## main driver script

## Things to check: community receptor results, cwbsdannentire: Does this exist?
## need to figure out shortterm exceedances in eco
main <- function() {
    source("utilities.r")
    source("load_data.r")
    source("utilities.r")
    source("runsource.r")
    source("media.r")
    source("food.r")
    source("risk.r")
    source("risk_si.r")
    source("eco.r")
    source("readgrf.R")
    source("shower.r")

    ## load data files
    ## original chemical and eco files
#    dat <- c(load_source_parameters(),
#             loadchem(fchem = "chemical_constant.csv",
#                      feco = "ecological_constant.csv"))
    ## Paul's new files
    dat <- c(load_source_parameters(),
             loadchem(fchem = "chem_const_imp.csv",
                      feco = "eco_const_imp.csv"))

    save(dat, file = "dat.rda")

    ## load and check input parameters
    inputpars <- getinputs(dat$chemical_constant)

    ## generate scenarios as all the permutations of the
    ## options that are selected
    grps <- makescen(inputpars)
    print(grps)

    ## set up storage locations for model output
    riskout <- as.list(rep(NA, times = nrow(grps)))
    ecoout <- as.list(rep(NA, times = nrow(grps)))

    for (ii in 1:nrow(grps)) {

        ## modify srctypeid based on liner type in SI scenario
        if (grps$srctype[ii] == "SI")
            srctypeid.loc <- c(1:3)[grps$linertypeid[ii] + 1]
        else
            srctypeid.loc <- grps$srctypeid[ii]

        parameters <- dat[["parameters"]]
        fate_constant <- dat[["fate_constant"]]
        chemical_constant <- dat[["chemical_constant"]]
        data_dict <- dat[["data_dict"]]
        tbl_model_parm <- dat[["tbl_model_parm"]]
        parm_dim_sum <- dat[["parm_dim_sum"]]
        exposure_constant <- dat[["exposure_constant"]]

        ## 1 parameter in data_dict depends on srctype
        incvec <- data_dict$wmutype == tolower(grps$srctype[ii]) |
            data_dict$wmutype == ""
        data_dict <- data_dict[incvec,]

        ## restrict fate_constant to selected srctype and site
        ## both srctype and site specified
        incvec1<- (fate_constant$siteid == grps$siteid[ii]) &
            (fate_constant$srctypeid  ==  srctypeid.loc)
        notboth <- is.na(incvec1)
        incvec1[is.na(incvec1)] <- F

        incvec2<- fate_constant$srctypeid == srctypeid.loc
        incvec2[is.na(incvec2)] <- F
        incvec2 <- incvec2 & notboth
        incvec3 <- fate_constant$siteid == grps$siteid[ii]
        incvec3[is.na(incvec3)] <- F
        incvec3 <- incvec3 & notboth

        incvec5 <- fate_constant$linertypeid == grps$linertypeid[ii]
        incvec5[is.na(incvec3)] <- F
        incvec4 <- is.na(fate_constant$siteid) & is.na(fate_constant$srctypeid) &
            is.na(fate_constant$linertypeid)

        fate_constant <- fate_constant[incvec1 | incvec2 | incvec3| incvec4 | incvec5,]

        ## restrict chemical_constant by cas number
        incvec <- chemical_constant$cas == grps$cas[ii]
        incvec[is.na(incvec)] <- TRUE
        chemical_constant <- chemical_constant[incvec,]

        nyear <- getval(tbl_model_parm, "numswyear")

        ## combine all parameters into two names vectors: 1 for
        ## strings and 2 for numbers
        parm <- getval.all(list(chemical_constant, fate_constant,
                                exposure_constant, tbl_model_parm))

        pnum <- parm[[2]]
        dfsource <- runsource(parameters, fate_constant, chemical_constant,
                              data_dict, tbl_model_parm, parm_dim_sum,
                              grps$srctype[ii], grps$linertypeid[ii], pnum)
#        save(dfsource, file = "dfsource.rda")
#            load("dfsource.rda")
        mout <- media(parm[[1]], pnum, dfsource, nyear, grps$siteid[ii],
                      srctypeid.loc,
                      grps$srctype[ii])

        if (grps$srctype[ii] == "SI") {
            riskout[[ii]] <-risk_si(mout, pnum, parm[[1]]["chemtype"])
        }
        else {
            foodout <- food(parm[[1]], pnum, mout, srctypeid.loc)
            riskout[[ii]] <- risk(mout,foodout, pnum, parm[[1]]["chemtype"])
            ecoout[[ii]] <- eco(grps$cas[ii], mout, foodout,
                                dat$eco_const)
        }
    }
    return(list(riskout, ecoout, grps))
}

output <- main()
