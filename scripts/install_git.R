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

#Repo of package to be checked/installed
pkg <- c("mikelove/tximeta") # <--- Add repo(s) of your package(s) here 

#Check if package is installed
new.pkg <- pkg[!(basename(pkg) %in% installed.packages()[, "Package"])]

if (length(new.pkg)){ 
  
  #Install package from repository
  status <- install_github(new.pkg)
  
  if(status) {
  
  #Send output to .txt file  
  sink(outtxt)
  cat(sprintf("%s installation successful", new.pkg))
  sink()
  } else{
    stop(sprintf("%s installation failed", new.pkg))
    
    # If the installation fails, it is most likely because there are dependencies 
    # that have not been added through the environment.yaml file. For this
    # purpose you must verify you have added all the necessary packages to 
    # install your desired package.
  }

} else {
  sink(outtxt)
  cat(sprintf("%s already installed", new.pkg))
  sink()
}





