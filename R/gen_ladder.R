#' Function to generate a grid of integers
#'
#' @param Ntotal Total number of predictors
#' @param ratio Rate for decrease the grid
#' @param Nmin The minimum number of the grid of integers
#' @returns A vector of integers in decreasing order starting from Ntotal and ending at Nmin.
#'
#' @examples
#' ladder <- gen_ladder(Ntotal = 1000, ratio = 0.9, Nmin = 20)
#'
#' @export
gen_ladder <- function(Ntotal, ratio = 0.9, Nmin = 5){
  x <- x.last <- Ntotal
  while(x.last > Nmin){
    x.next <- round(x.last * ratio)
    if(x.next == x.last){
      x.next <- x.next - 1
    }
    x.last <- max(x.next, Nmin)
    x <- c(x, x.last)
  }
  x
}
