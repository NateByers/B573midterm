# read in all of the lines in the helpers file and check to see
# if they have been installed (if not, install), then load
# the libraries

lines <- readLines("helpers.R")
lines <- grep("::", lines, value = TRUE)
libraries <- sapply(lines, function(x) {
  library_ <- strsplit(x, split = "::")[[1]][1]
  if(grepl("\\s{1,}", library_)) {
    library_ <- strsplit(library_, "\\s{1,}")[[1]]
    library_ <- tail(library_, 1)
  }
  library_
  })

libraries <- libraries[libraries != "shiny"]

for(i in unique(libraries)) {
  install.packages(i)
  library(i, character.only = TRUE)
}
