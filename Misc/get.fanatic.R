#' @title Fanatic
#'
#' @description Getting the fanatic achievement
#'
#' @param browser \code{"character"}, the name of the browser
#' @param browser.path \code{"character"}, the path to the browser
#' @param connections \code{"numeric"}, the number of connections to make
#' @param browsing Two \code{"numeric"} for the amount of time to "browse" in seconds (default = \code{c(60, 600)})
#' @param sleeping Two \code{"numeric"} for the amount of time to "sleep" in seconds (default = \code{c(36000, 72000)})
#' 
#' @examples
#' \dontrun{
#' ## Connecting two times for less than 5 seconds every less than 10 seconds
#' get.fanatic("Safari", "/Applications", 2, browsing = c(1, 5), sleeping = c(5, 10))
#' }
#'
#' @seealso
#' 
#' @author Thomas Guillerme
#' @export

get.fanatic <- function(browser, browser.path, connections, browsing = c(60, 600), sleeping = c(36000, 72000)) {

    ## Converting the browsing/sleeping distributions in minutes/seconds
    ## When "browsing" (up to ten minutes)
    browsing.fun <- function(browsing) {sample(browsing[1]:browsing[2], 1)} #*60

    ## When "sleeping" (between 10 and 20 hours)
    sleeping.fun <- function(sleeping) {sample(sleeping[1]:sleeping[2], 1)} #*60*60

    ## Browser's path
    browser_path <- paste0(browser.path,"/",browser, ".app")

    for(connection in 1:connections) {

        cat(paste0(Sys.time(), " - Attempting connection ", connection, ":"))

        ## Open the browser
        system(paste0("open -a ", browser_path, " https://stackoverflow.com/"))

        cat(paste0(" successful.\n"))

        ## "Browsing"
        Sys.sleep(browsing.fun(browsing))

        cat(paste0(Sys.time(), " - Browsing complete :"))
        
        ## Close the browser
        system(paste0("quit ", browser))

        cat(paste0(" successful.\n\n"))

        ## "Sleeping"
        Sys.sleep(sleeping.fun(sleeping))

        ## Connection increment
        connection <- connection + 1
    }
}