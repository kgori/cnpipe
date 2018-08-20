#' @importFrom "logging" logerror
#' @export
check_file <- function(filename) {
    if (!file.exists(filename)) {
        logerror(paste("File", filename, "not found"))
        stop("ERROR - file error")
    } else {
        return (normalizePath(filename))
    }
}

#' @importFrom "logging" logerror
#' @export
check_dir_and_create <- function(dirname) {
    if (!dir.exists(dirname)) {
        dir.create(dirname)
        if (!dir.exists(dirname)) {
            logerror(paste("Unable to create directory", dirname))
            stop("ERROR - directory error")
        }
    }
    return (normalizePath(dirname))
}
