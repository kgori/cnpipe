magic_number <- function(filename, n = 2) {
    if (!file.exists(filename)) {
        stop("File not found")
    }

    conn <- file(filename, "rb")
    x <- readBin(conn, "raw", n = n, size = 1)
    close(conn)
    x
}

is_gzip <- function(filename) {
    if (!file.exists(filename)) {
        stop("File not found")
    }
    all(magic_number(filename, 2) == charToRaw("\x1f\x8b"))
}

is_bzip2 <- function(filename) {
    if (!file.exists(filename)) {
        stop("File not found")
    }
    (rawToChar(magic_number(filename, 3)) == "BZh")
}
