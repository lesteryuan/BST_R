## load parameters file
## turn flat file into dataframe
make_df <- function(df1, fname.row) {

    df1$modelcode <- tolower(df1$modelcode)

    rownames <- unique(df1[, fname.row])
    nrow <- length(rownames)
    codes <- unique(df1$modelcode)
    ncol <- length(codes)
    a <- matrix(NA, ncol = ncol, nrow = nrow)
    dimnames(a) <- list(rownames, codes)
    for (i in 1:length(rownames)) {
        incvec <- df1[, fname.row]==rownames[i]
        for (j in 1:length(codes)) {
            selvec <- incvec & codes[j] == df1$modelcode
            if (sum(selvec) == 1) {
                a[i,j]<- df1$value[selvec]
            }
            else {
                if (sum(selvec) > 1) {
                    if (length(unique(df1$value[selvec])) > 1) {
                        print("multiple values")
                        print(rownames[i])
                        print(codes[j])
                    }
                    else {
                        a[i,j] <- df1$value[selvec][1]
                    }
                }
            }
        }
    }
    dfa <- data.frame(a)
    for (i in 1:ncol(dfa)) {
        ## set blanks to NA
        incvec <- dfa[,i] == ""
        incvec[is.na(incvec)] <- F
        dfa[incvec,i] <- NA

        x <- suppressWarnings(as.numeric(dfa[,i]))
        if (sum(is.na(x)) == sum(is.na(dfa[,i]))) dfa[,i] <- x
    }
    return(dfa)
}

load_source_parameters <- function() {
    datapath <- paste(getwd(), "data", sep = "/")
    parameters <- read.csv(paste(datapath, "parameters.csv", sep = "/"))
    source("fnames.un.r")
    parameters <- fnames.un(parameters)
    ## make all input.file.names lower case
    parameters$input.file.name <- tolower(parameters$input.file.name)
    parameters$code <- tolower(parameters$code)

    fate_constant <- read.csv(paste(datapath, "fate_constant.csv", sep ="/"))
    fate_constant <- fnames.un(fate_constant)
    fate_constant$modelcode <- tolower(fate_constant$modelcode)

    chemical_constant <- read.csv(paste(datapath, "chemical_constant.csv",
                                        sep = "/"))
    chemical_constant <- fnames.un(chemical_constant)
    chemical_constant$modelcode <- tolower(chemical_constant$modelcode)
    df_chemical_constant <- make_df(chemical_constant, "cas")

    data_dict <- read.csv(paste(datapath, "tbl_data_dictionary_crosstab.csv",
                                sep = "/"), strip.white = T)
    data_dict <- fnames.un(data_dict)
    data_dict$modelcode <- tolower(data_dict$modelcode)
    data_dict[,"source"] <- tolower(data_dict[, "source"])
    data_dict[,"wmutype"] <- tolower(data_dict[, "wmutype"])
    ## area_1 is in data_twice.  drop one that is not SrcArea
#    incvec <- data_dict$modelcode == "area_1" & data_dict[, "source"] != "SrcArea"
#    data_dict <- data_dict[!incvec,]
    data_dict <- unique.data.frame(data_dict[, c("modelcode", "source",
                                                 "index1", "index2", "index3",
                                                 "wmutype")])
    incvec <- ! is.na(data_dict$index2) | ! is.na(data_dict$index3)
    parm_dim <- data_dict[incvec,]

    getmax <- function(x) {
        if (any(! is.na(x))) return(max(x, na.rm = T))
        else return(NA)
    }
    dim2 <- tapply(parm_dim$index2, parm_dim[, "source"], getmax)
    dim3 <- tapply(parm_dim$index3, parm_dim[, "source"], getmax)

    parm_dim_sum <- data.frame(modelcode = names(dim2),
                               dim2= as.vector(dim2),
                               dim3= as.vector(dim3))

    ## only save entries where there is a different name
    data_dict <- unique.data.frame(data_dict[, c("modelcode", "source",
                                                 "wmutype")])

    incvec <- ! data_dict[, "source"]==""
    data_dict <- data_dict[incvec,]
    incvec <- data_dict$modelcode == data_dict[, "source"]
    data_dict <- data_dict[!incvec,]

    ## correct some strange characters
    incvec <- regexpr("da", data_dict$modelcode)!= -1
    data_dict$modelcode[incvec] <- "da"
    incvec <- regexpr("dw", data_dict$modelcode)!= -1
    data_dict$modelcode[incvec] <- "dw"

    tbl_model_parm <- read.csv(paste(datapath, "tbl_model_parameters.csv",
                                     sep ="/"))
    tbl_model_parm <- fnames.un(tbl_model_parm)
    tbl_model_parm$modelcode <- tolower(tbl_model_parm$modelcode)
    #tbl_model_parm$value <- as.numeric(tbl_model_parm$value)
#    incvec <- ! is.na(tbl_model_parm$value)
#    tbl_model_parm <- tbl_model_parm[incvec,]

    exposure_constant <- read.csv(paste(datapath, "exposure_constant.csv",
                                        sep = "/"))
    exposure_constant <- fnames.un(exposure_constant)
    exposure_constant$modelcode <- tolower(exposure_constant$modelcode)
    ## no variations in exposure constant across runs so boiling this
    ## data frame down to means
    mnval <- tapply(exposure_constant$value, exposure_constant$modelcode,
                    mean)
    dftemp.exp <- data.frame(modelcode = names(mnval),
                             value = as.vector(mnval))

    dat <- list(parameters = parameters,
                fate_constant = fate_constant,
                chemical_constant = chemical_constant,
                data_dict = data_dict,
                tbl_model_parm = tbl_model_parm,
                parm_dim_sum = parm_dim_sum,
                exposure_constant = dftemp.exp)
    return(dat)
}


