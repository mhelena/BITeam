# Parses config.sh for the R steps. The path comes from the first argument or $PIPELINE_CONFIG.
read_config <- function() {
  path <- commandArgs(trailingOnly = TRUE)[1]
  if (is.na(path)) path <- Sys.getenv("PIPELINE_CONFIG")
  if (!nzchar(path) || !file.exists(path))
    stop("pass the config file as the first argument or set PIPELINE_CONFIG")
  x <- trimws(sub("#.*$", "", readLines(path)))
  x <- x[nzchar(x)]
  val <- gsub("^[\"']|[\"']$", "", trimws(sub("^[^=]*=", "", x)))
  cfg <- setNames(as.list(path.expand(val)), trimws(sub("=.*$", "", x)))
  function(key) {
    if (is.null(cfg[[key]])) stop("missing ", key, " in ", path)
    cfg[[key]]
  }
}
