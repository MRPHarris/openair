## Implementing variance quantification in openair's clustering framework

## SETUP

## Source functions, load all packages the openair package uses.
proj_dir <- "C:/Users/matth/Desktop/University Files/mh phd 2020/software_and_coding/rstudio/MyPackages/openair/"
scripts_dir <- paste0(proj_dir,"R/")
files <- list.files(scripts_dir, full.names = TRUE)
# #file.edit(files)
flist <- vector(mode = "list", length = length(files))
for(i in seq_along(flist)){
  source(files[i])
}

## Load packages
library(pacman)
pacman::p_load(
  grid,dplyr,purrr,tidyr,readr,mgcv,lattice,latticeExtra,lubridate,
cluster,mapproj,hexbin,Rcpp,grDevices,graphics,methods,stats,MASS,utils,KernSmooth,maps,mapdata,quantreg
)


##### Unpacking the cluster function


## A version of the trajectory clustering function with no plotting.
trajCluster_lite <- function(traj, method = "Euclid", n.cluster = 5) {
  # silence R check
  freq <- hour.inc <- default <- NULL

  if (tolower(method) == "euclid") {
    method <- "distEuclid"
  } else {
    method <- "distAngle"
  }

  # remove any missing lat/lon
  traj <- filter(traj, !is.na(lat), !is.na(lon))

  # check to see if all back trajectories are the same length
  traj <- group_by(traj, date) %>%
    mutate(traj_len = length(date))

  if (length(unique(traj$traj_len)) > 1) {
    warning("Trajectory lengths differ, using most common length.")
    ux <- unique(traj$traj_len)
    nmax <- ux[which.max(tabulate(match(traj$traj_len, ux)))]
    traj <- ungroup(traj) %>%
      filter(traj_len == nmax)
  }
  # Calculate clusters
  calcTraj <- function(traj) {

    ## make sure ordered correctly
    traj <- traj[order(traj$date, traj$hour.inc), ]

    ## length of back trajectories
    traj <- group_by(traj, date) %>%
      mutate(len = length(date))

    ## find length of back trajectories
    ## 96-hour back trajectories with origin: length should be 97
    n <- max(abs(traj$hour.inc)) + 1

    traj <- subset(traj, len == n)
    len <- nrow(traj) / n

    ## lat/lon input matrices
    x <- matrix(traj$lon, nrow = n)
    y <- matrix(traj$lat, nrow = n)

    z <- matrix(0, nrow = n, ncol = len)
    res <- matrix(0, nrow = len, ncol = len)

    if (method == "distEuclid") {
      res <- .Call("distEuclid", x, y, res)
    }

    if (method == "distAngle") {
      res <- .Call("distAngle", x, y, res)
    }

    res[is.na(res)] <- 0 ## possible for some to be NA if trajectory does not move between two hours?

    dist.res <- as.dist(res)
    clusters <- pam(dist.res, n.cluster)
    cluster <- rep(clusters$clustering, each = n)
    traj$cluster <- as.character(paste("C", cluster, sep = ""))
    traj
  }
  # return traj
  return(traj)
}
