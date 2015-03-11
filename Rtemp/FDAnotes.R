

require(knitr)
require(markdown)

go <- function(page = "BHSpresentationPlots", open = FALSE, knit = TRUE, purl = TRUE, ...) {

  # create file names
  Rmd <- paste0("B:/Loch_Ard_Felling/Package/rivertemp/R/", page,".Rmd")
  md <-  paste0("B:/Loch_Ard_Felling/", page, ".md")
  html <-  paste0("B:/Loch_Ard_Felling/", page, ".html")
  R <- paste0("B:/Loch_Ard_Felling/Package/rivertemp/R/", page,".R")

  # added this as it was crashing if present
  if (exists(".Random.seed", envir = .GlobalEnv)) rm(.Random.seed, envir = .GlobalEnv)

  if (knit) {
    # knit!
    knit(Rmd, output = md)
    # other options to markdownToHTML??
    markdownToHTML(md, out = html)

    if (open) openpage(page, ...)
  }
  if (purl) {
    knit(Rmd, output = R, tangle = TRUE)
  }

  cat("Loading environment\n")
  load("B:/Loch_Ard_Felling/Package/rivertemp/R/currentImage.Rdata", envir = .GlobalEnv)
}


openpage <- function(page, args = "", ...) {
  system(paste0('"c:/Program Files/Internet Explorer/iexplore.exe" ', args,' B:/Loch_Ard_Felling/', page, '.html'), wait = FALSE)
}


#go(action = "purl")
