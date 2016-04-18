.onAttach <- function(lib, pkg)
{
  unlockBinding(".qcc.options", asNamespace("qcc")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    { # > figlet qcc
      packageStartupMessage(
"  __ _  ___ ___ 
 / _  |/ __/ __|  Quality Control Charts and 
| (_| | (_| (__   Statistical Process Control
 \\__  |\\___\\___|
    |_|           version ", version)
}
else
  { packageStartupMessage("Package 'qcc' version ", version) } 

  packageStartupMessage("Type 'citation(\"qcc\")' for citing this R package in publications.")
  invisible()
}
