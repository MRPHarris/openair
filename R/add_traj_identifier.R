##' Add a single-digit identifier to each trajectory within a given data frame.
##'
##' This function adds a simple numeric identifier to each trajectory in the
##' given data.frame. It is used within the spatial variance cluster calculations
##' in order to discern between trajectories in a straightforward manner. It operates
##' diffeerently depending on whether the dataset contains trajectories run once per day
##' or multiple times per day, relying on hour and minute initialisation differences in
##' each respective case.
##'
##' @param traj a data.frame containing trajectory data.
##' @param ntraj_1 TRUE or FALSE. Does the dataset contain trajectories run once per day, or more than once?
##'
##' @export
##'
add_traj_identifier <- function(traj, ntraj_1 = TRUE){
  # This function assumes that convert_openair() has been used on the dataset.
  x <- traj # pass var to previous notation
  if(!is.data.frame(x)){
    x <- as.data.frame(x)
  }
  # Compute trajectory number. Starts at 1.
  if(!isTRUE(ntraj_1)){
    x$start.time.minutes <- substr(x[,12],12,19)
    x$start.time.minutes <- 60*24*as.numeric(chron::times(x$start.time.minutes))
    x$trajectory <- cumsum(c(0, as.numeric(diff(x$start.time.minutes)) != 0)) + 1
  } else {
    if(sign(x$hour.inc[1]) == -1){
      x$diffs <- 1
      x$diffs[2:nrow(x)] <- diff(x$hour.inc)
      x$trajectory <- cumsum(c(0, as.numeric(x$diffs[2:nrow(x)]) != 1)) + 1
    } else {
      x$diffs <- 1
      x$diffs[2:nrow(x)] <- diff(x$hour.inc)
      x$trajectory <- cumsum(c(0,as.numeric(x$diffs[2:nrow(x)]) != -1 )) + 1
    }

  }
  #  x_new <<- as.data.frame(x)
  x
  #  rm(x_new, envir = parent.frame())
}
