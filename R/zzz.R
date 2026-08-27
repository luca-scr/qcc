.onLoad <- function(lib, pkg) {
  defaults <- .qcc_default_options()
  missing <- setdiff(names(defaults), names(options()))
  if(length(missing))
    options(defaults[missing])
  invisible(NULL)
}

# TODO: add .onUnload() to cleanuup .onLoad()

#' Package Startup Message
#'
#' Builds the qcc package startup message with its version and citation
#' reminder. The greeter for interactive sessions was obtained by
#' running `figlet qcc`.
#'
#' @keywords internal
#' @noRd
qccStartupMessage <- function()
{
  msg <- c(paste0(
"  __ _  ___ ___ 
 / _  |/ __/ __|  Quality Control Charts and 
| (_| | (_| (__   Statistical Process Control
 \\__  |\\___\\___|
    |_|           version ", 
packageVersion("qcc")),
"\nType 'citation(\"qcc\")' for citing this R package in publications.")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  msg <- qccStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'qcc' version", packageVersion("qcc"))
  packageStartupMessage(msg)
  invisible()
}
