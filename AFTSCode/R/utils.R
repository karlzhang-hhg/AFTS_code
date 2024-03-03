#' Perform a deep-copy on an object.
#'
#' @param object A data series.
#' @export 
deep_copy <- function(object) {
  unserialize(serialize(object, NULL))
}
