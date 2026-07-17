#' Package Startup Message
#'
#' Builds the qcc package startup message with its version and citation
#' reminder. The greeter for interactive sessions was obtained by
#' running `figlet qcc`.
#'
#' @keywords internal
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

.onAttach <- function(lib, pkg)
{
  # unlock .qcc.options variable allowing its modification
  unlockBinding(".qcc.options", asNamespace("qcc")) 
  # startup message
  msg <- qccStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'qcc' version", packageVersion("qcc"))
  packageStartupMessage(msg)
  invisible()
}
