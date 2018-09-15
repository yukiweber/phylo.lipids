#' "not in"
#'
#' @name not_in
#' @export
#' @usage %in%
#' @keywords internal

'%!in%' <- function(x,y)!('%in%'(x,y))

