defaults <- function(x, defaults) {
  c(x, defaults[setdiff(names(defaults), names(x))])
}
