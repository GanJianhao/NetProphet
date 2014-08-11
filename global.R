libs <- c("shiny", "knitr", "R.utils", "markdown", "igraph", "shinyBS")
# check if all packages in libs are available
available <- suppressWarnings(sapply(libs, require, character.only=TRUE))
inst.libs <- libs[available == FALSE]
# if(length(inst.libs) != 0) {
#   install.packages(inst.libs, dependencies = TRUE)
#   suppressWarnings(sapply(inst.libs, require, character.only=TRUE))
# }
# function to render .Rmd files to html on-the-fly
includeRmd <- function(path){
  # shiny:::dependsOnFile(path)
  
  contents <- paste(readLines(path, warn = FALSE), collapse = '\n')
  # do not embed image or add css
  html <- knit2html(text = contents, fragment.only = TRUE, options = "", stylesheet = "www/empty.css")
  Encoding(html) <- 'UTF-8'
  HTML(html)
}
