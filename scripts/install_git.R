args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

## This script installs packages needed in the workflow that are not available
## in conda channels, but are available in a github repository.
## Here we demonstrate the installation of the package 'tximeta', which is 
## necessary to run the workflow:

suppressPackageStartupMessages(library(devtools))

print(outtxt)

#Change TAR environment variable (for ubuntu 16.04, should be tested in mac)
options(unzip = "internal")

#Package name to be checked/installed
pkg <- c("tximeta") # <--- Add your package name here 

#Check if package is installed
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]

if (length(new.pkg)){ 
  
  #Install package from repository
  t <- install_github("mikelove/tximeta")
  
  if(t) {
  
  #Send output to .txt file  
  sink(outtxt)
  cat(sprintf("Package %s installation successful", pkg))
  sink()
  } else{
    stop("Installation failed")
    
    # If the installation fails, it is most likely because there are dependencies 
    # that have not been added through the environment.yaml file. For this
    # purpose you must verify you have added all the necessary packages to 
    # install your desired package.
  }

} else {
  sink(outtxt)
  print("Package %s already installed")
  sink()
}





