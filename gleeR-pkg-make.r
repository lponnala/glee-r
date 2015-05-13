
pkg_name = "gleeR"
pkg_title = "glee for differential protein expression"
version_number = 0.1
license_type = "GPL-2"

# -- create barebones package --
if (file.exists(pkg_name)) { unlink(pkg_name, recursive = TRUE) }
devtools::create(pkg_name, rstudio=FALSE)

# -- assign title, version & license --
L = readLines(paste(pkg_name,"/DESCRIPTION",sep=""))
L[stringr::str_detect(L, "^Title")] = paste("Title: ", pkg_title, sep="")
L[stringr::str_detect(L, "^Version")] = paste("Version: ", version_number, sep="")
L[stringr::str_detect(L, "^License")] = paste("License: ", license_type, sep="")
cat(L, file=paste(pkg_name,"/DESCRIPTION",sep=""), sep="\n")

# -- copy the code --
file.copy(from = paste0(pkg_name,".r"), to = paste0("./",pkg_name,"/R"), overwrite = TRUE)

# -- list all packages that are used in your code --
# these are typically the ones you would specify via library(...) during initial development
pkgs_used = c("dplyr")

# # -- print package dependencies, to select which packages (if any) to mention in Depends --
# # weirdly, if this section is not commented out, the devtools::check step below fails!
# cat("\n", "dependencies for packages used", "\n", sep="")
# pkgs_avail = available.packages()
# pkgs_deps = tools::package_dependencies(pkgs_used, db=pkgs_avail, which="Depends")
# for (k in 1:length(pkgs_deps)) {
	# if (is.null(pkgs_deps[[k]])) { stop("this package is not available: ", names(pkgs_deps)[k]) }
	# if (length(pkgs_deps[[k]]) > 0) { cat(names(pkgs_deps)[k], " : ", paste(pkgs_deps[[k]],collapse=","), "\n", sep="") }
# }

# -- specify Depends --
# use sparingly! the only packages you might mention here are (examples in brackets):
# - the dependencies (zoo) of poorly-written packages (xts, uses as.yearmon instead of zoo::as.yearmon)
#	- packages (xts) that are dependencies of several used packages (PerformanceAnalytics,quantmod,Quandl)
Depends = c()
if (length(Depends) > 0) {
	cat("\n", "specifying Depends: ", paste(Depends,collapse=","), "\n", sep="")
	for (pkg in Depends) { devtools::use_package(pkg, type="Depends", pkg=pkg_name) }
}

# -- specify Imports --
# everything from pkgs_used not in Depends will be specified in Imports
Imports = setdiff(pkgs_used,Depends)
if (length(Imports) > 0) {
	cat("\n", "specifying Imports: ", paste(Imports,collapse=","), "\n", sep="")
	for (pkg in Imports) { devtools::use_package(pkg, type="Imports", pkg=pkg_name) }
}

# -- check --
check_dir = "check_dir"
if (pkg_name %in% as.character(installed.packages()[,1])) { remove.packages(pkg_name) }
if (!file.exists(check_dir)) { dir.create(check_dir) }
devtools::check(pkg_name, document=TRUE, quiet=FALSE, cleanup=FALSE, cran=FALSE, force_suggests=FALSE, check_dir=paste(getwd(),"/",check_dir,sep=""))
# open log file to view (throws error but opens file)
cmd = tryCatch({ edit(file=paste0(check_dir,"/",pkg_name,".Rcheck/00check.log"), editor="notepad++.exe") }, error = function(e) { NA }, warning = function(w) { NA })

# -- build --
# [could ignore this step, but if you build a binary .zip package this way, install it by running "R CMD INSTALL gleeR_0.1.zip" and not via devtools::install_local("gleeR_0.1.zip")]
# devtools::build(pkg_name, binary=TRUE, vignettes=FALSE, manual=FALSE)

# -- install --
devtools::install(pkg_name)
