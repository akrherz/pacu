#' Set pacu options
#'
#' @title Setting some options for the package
#' @name pacu_options
#' @description Set settings regarding messages and default behaviors of the package
#' @param suppress.warnings whether to suppress warning messages
#' @param suppress.messages whether to suppress  messages
#' @return as a side effect it modifies the \sQuote{pacu.options} environment.
#' @export
#' @examples
#'\dontrun{
#' names(pacu.options)
#' pacu_options(exe.path = "some-new-path-to-executable")
#' pacu.options$exe.path
#' }

pacu_options <- function(suppress.warnings = FALSE, 
                         suppress.messages = FALSE){

  assign('suppress.warnings', suppress.warnings, pacu.options)
  assign('suppress.messages', suppress.messages, pacu.options)
}

#' Environment which stores PACU options
#'
#' @title Environment which stores PACU options
#' @name pacu.options
#' @description Environment which can store the path to the executable, warning settings and
#'              where examples are located.
#'              Creating an environment avoids the use of global variables or other similar practices
#'              which would have possible undesirable consequences.
#' @return This is an environment, not a function, so nothing is returned.
#' @export
#' @examples
#' \donttest{
#' names(pacu.options)
#' ## to suppress messages
#' pacu_options(suppress.messages = TRUE)
#' }
#'

pacu.options <- new.env(parent = emptyenv())
assign('suppress.warnings', FALSE, pacu.options)
assign('suppress.messages', FALSE, pacu.options)
assign('.run.local.tests', TRUE, pacu.options)
assign('.run.experimental.tests', FALSE, pacu.options)


## Import packages needed for pacu to work correctly
#' @import apsimx gstat httr jsonlite sf stars units XML
#' @importFrom grDevices hcl.colors
## @importFrom utils read.table write.table packageVersion zip unzip
## @importFrom tools file_path_sans_ext file_ext
## @importFrom stats aggregate anova coef cor cov2cor deviance lm optim qt var sd setNames sigma integrate median na.omit weighted.mean
## @importFrom graphics legend, par, points
NULL
