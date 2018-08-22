#' "not in"
#'
#' @name not_in
#' @export
#' @usage %in%
#' @keywords intern

'%!in%' <- function(x,y)!('%in%'(x,y))

