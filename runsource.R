runsource <- function(parameters, fate_constant, chemical_constant,
                      data_dict, tbl_model_parm, parm_dim_sum,
                      srctype, linertypeid, pnum) {
#    require("odbc")
#    require("DBI")
    options(scipen = 999) # turn off scientific notation

    ## write rows out to ssf file formatted correctly
    writerow <- function(zz, dfloc) {
        if (nrow(dfloc) > 1) {
            dfloc <- dfloc[order(as.numeric(dfloc$value)),]
        }
        if (is.na(dfloc$value[1])) dfloc$value[1] <- " "
            cat('"',dfloc$code[1],'",', dfloc$dimensions[1],',"',
                dfloc$datatype[1],'",0,"',dfloc$units[1],'",\n',
                file = zz, sep = "")
        if (dfloc$dimensions[1]==0) {
            cat(dfloc$value[1], ',\n', file = zz, sep = "")
        }
        else {
            if (dfloc$dimensions[1]==1) {
                cat('1, ',dfloc$value[1], ',\n', file = zz, sep = "")
            }
            else {
                if (dfloc$dimensions[1]==2) {
                    cat('1,\n', file = zz)
                    cat(dfloc$dim2[1],",", file = zz, sep = "")
                    if (nrow(dfloc) == 1) {
                        cat(rep(dfloc$value[1],times = dfloc$dim2[1]), sep = ",", file = zz)
                    }
                    else {
                        cat(c(dfloc$value[2], dfloc$value[1]), sep =",",
                            file = zz)
                    }
                    cat(",\n", file = zz)
                }
                else {
                    if (dfloc$dimensions[1] == 3) {
                        cat("1, 2,\n", file = zz)
                        for (k in 1:2) {
                            cat(dfloc$dim3[1],", ", file = zz, sep = "")
                            cat(rep(dfloc$value[1], times = dfloc$dim3[1]),
                                sep =",", file = zz)
                            cat(",\n", file = zz)
                        }
                    }
                }
            }
        }
    }

    ## area_1 has to be changed to area_si when running SI
    if (srctype == "SI") {
        incvec <- fate_constant$modelcode == "area_1"
        fate_constant$value[incvec] <- pnum["area_si"]
    }

    ## update fate_constant modelcodes to names needed in ssf
    ## Area_1 gets copied over to two parameter spots
    fate_constant <- merge(fate_constant, data_dict, by = "modelcode",
                            all.x = T)
    incvec <- ! is.na(fate_constant[, "source"])
    fate_constant$modelcode[incvec] <- fate_constant[incvec,"source"]

    ## update chemical_constant modelcodes to names needed in ssf
    chemical_constant <- merge(chemical_constant, data_dict,
                                by = "modelcode",
                                all.x = T)
    incvec <- ! is.na(chemical_constant[, "source"])
    chemical_constant$modelcode[incvec] <- chemical_constant[incvec,"source"]

#    dfpnum <- data.frame(modelcode = names(pnum), val = as.vector(pnum))
#    dfpnum <- merge(dfpnum, data_dict, by = "modelcode", all.x = T)
#    incvec <-! is.na(dfpnum[, "source"])
#    dfpnum$modelcode[incvec] <- dfpnum[incvec,"source"]

    ## get input file names for HWIR
    fnames.all <- sort(unique(parameters$input.file.name))
    fnames <- fnames.all[regexpr("ssf", fnames.all) != -1]

    ## loop through fnames and try to find values for them
    for (i in fnames) {

        incvec <- parameters$input.file.name == i
        dftemp <- parameters[incvec,c("code", "dimensions", "units", "datatype")]
        ## merge in dimensions
        dftemp <- merge(dftemp, parm_dim_sum, by.y = "modelcode",
                        by.x = "code", all.x = T)

        dftemp <- merge(dftemp, chemical_constant[, c("modelcode", "value")],
                        by.x = "code", by.y = "modelcode", all.x = T)
        dftemp2 <- fate_constant[, c("modelcode", "value")]
        names(dftemp2)[[2]] <- "value2"
        dftemp <- merge(dftemp, dftemp2,
                        by.x = "code", by.y = "modelcode", all.x = T)
        dftemp3 <- tbl_model_parm[, c("modelcode", "value")]
        names(dftemp3)[[2]] <- "value3"
        dftemp <- merge(dftemp, dftemp3,
                        by.x = "code", by.y = "modelcode",all.x = T)
        ## combine all of the merged in values
        incvec <- is.na(dftemp$value)
        dftemp$value[incvec] <- dftemp$value2[incvec]
        incvec <- is.na(dftemp$value)
        dftemp$value[incvec] <- dftemp$value3[incvec]

        ## get rid of "chemical specific"
        incvec <- regexpr("chemical", tolower(dftemp$value)) != -1
        dftemp$value[incvec] <- ""

#        incvec <- parameters$input.file.name == i
#        dfnew <- parameters[incvec,c("code", "dimensions", "units", "datatype")]
#        dfnew <- merge(dfnew, parm_dim_sum, by.y = "modelcode",
#                       by.x = "code", all.x = T)
#        dfnew <- merge(dfnew, dfpnum, by.y = "modelcode", by.x = "code",
#                       all.x=T)
#        print(dfnew)

        incvec <- dftemp$code == "srctype"
        dftemp$value[incvec] <- srctype
        ## run id doesn't have to change with the data organization
        ## I run now.
        incvec <- dftemp$code == "runid"
        dftemp$value[incvec] <- 1
        incvec <- dftemp$code == "inputdb"
        dftemp$value[incvec] <- paste(getwd(),
                                      "BST_v1-1m-allchems_clean_20240131_NOT_FOR_DISTRIBUTION.mdb", sep = "/")
        incvec <- dftemp$code == "metdir"
        dftemp$value[incvec] <- paste(getwd(), "HWIR/MetData", sep = "/")
        incvec <- dftemp$code == "outputdb"
        dftemp$value[incvec] <- paste(getwd(), "source_results/sourceoutput.mdb", sep = "/")
        incvec <- dftemp$code == "cpdirectory"
        dftemp$value[incvec] <- paste(getwd(), "hwir/cpdata", sep = "/")
        incvec <- dftemp$code == "linertype"
        dftemp$value[incvec] <- linertypeid

        ## some changes to parameters for SI scenario
        if (srctype == "SI") {
            incvec <-dftemp$code == "c_in"
            dftemp$value[incvec] <- pnum["ctpwastedry"]*pnum["tss_in"]
        }

        zz <- file(paste(getwd(), "hwir/ssf", i, sep = "/"), "w")
        cat('1,\n', file = zz, sep = "")
        cat( i, ', "data group",\n', file = zz, sep = '')
        code.u <- unique(dftemp$code)
        cat(length(code.u),",\n", sep = "", file = zz)

        for (j in code.u) {
            incvec <- dftemp$code == j
            writerow(zz,dftemp[incvec,])
        }
        close(zz)
    }

#    require(RODBC)
#    con <- odbcConnectAccess2007(paste(getwd(), "source_results/sourceoutput.mdb", sep = "/"))

    dirsav <- getwd()
    setwd(paste(getwd(), "hwir", sep = "/"))
    if (srctype == "SI")  system2("./si_grf.bat")
    else system2("./la_grf.bat")

#    dfsource <- sqlFetch(con, "OutputData")
#    odbcCloseAll()
    setwd(dirsav)
    ## after loading the results into R, copy the old sourceoutput
    ## file over the new one to avoid accumulating results
#    system("cp ./SOURCE_RESULTS/SourceOutput0.mdb ./SOURCE_RESULTS/SourceOutput.mdb")

    dfsource <- readgrf(srctype)

    return(dfsource)
}

