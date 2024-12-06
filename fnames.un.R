# 3.19.2013
# make field names lower case and remove underscores

fnames.un <- function(df) {
    names0 <- names(df)
    names0 <- tolower(names0)
    names0 <- gsub("_", ".", names0)
    names(df) <- names0
    return(df)
}
