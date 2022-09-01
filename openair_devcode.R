## Implementing variance quantification in openair's clustering framework

## SETUP

## Source functions, load all packages the openair package uses.
proj_dir <- "C:/Users/matth/Desktop/University Files/mh phd 2020/software_and_coding/rstudio/MyPackages/openair/"
scripts_dir <- paste0(proj_dir,"R/")
files <- list.files(scripts_dir, full.names = TRUE)
# #file.edit(files)
# flist <- vector(mode = "list", length = length(files))
# for(i in seq_along(flist)){
#   source(files[i])
# }

## Load packages
library(pacman)
pacman::p_load(
  grid,dplyr,purrr,tidyr,readr,mgcv,lattice,latticeExtra,lubridate,
cluster,mapproj,hexbin,Rcpp,grDevices,graphics,methods,stats,MASS,utils,KernSmooth,maps,mapdata,quantreg
)

## Get example HYSPLIT trajectory data
prev_proj_dir <- "C:/Users/matth/Desktop/University Files/mh phd 2020/3 thesis/Chapter X Patriot Hills chemistry and transport/"
hysplit_data_dir <- "D:/DATA/HYSPLIT Files/HYSPLIT Runs Apr2022/Collated runs Apr 2022/"
# GDAS1 100m 120hr
PH_GDAS1_100m_120hr <- read.csv(file = paste0(hysplit_data_dir,"PH_100m_120hr_files_collated.csv"))
## Update date formats
PH_GDAS1_100m_120hr_2010 <- PH_GDAS1_100m_120hr %>%
  mutate(date.start = lubridate::ymd_hms(date.start)) %>%
  dplyr::filter(format(date.start, format = "%Y") == "2010") %>%
  dplyr::select(-c(X,diffs,trajectory)) %>%
  rename(date2 = date.inc,
         date = date.start)

## Based on HYSPLIT explanation:
# https://www.ready.noaa.gov/documents/Tutorial/html/traj_cluseqn.html

# Get single cluster
library(openair)
clus1 <- openair::trajCluster(traj = PH_GDAS1_100m_120hr_2010,
                     n.cluster = 2,
                     plot = TRUE, showplot = FALSE)


##### Unpacking the cluster function

trajCluster_dev <- function(traj, method = "Euclid", n.cluster = 5,
                            plot = TRUE, showplot = FALSE, type = "default",
                            cols = "Set1", split.after = FALSE, map.fill = TRUE,
                            map.cols = "grey40", map.alpha = 0.4,
                            projection = "lambert",
                            parameters = c(51, 51), orientation = c(90, 0, 0),
                            by.type = FALSE, origin = TRUE, ...) {

  ## vars for intra-function testing
  traj = PH_GDAS1_100m_120hr_2010
  method = "Euclid"
  n.cluster = 5
  plot = TRUE
  showplot = FALSE
  type = "default"
  cols = "Set1"
  split.after = FALSE
  map.fill = TRUE
  map.cols = "grey40"
  map.alpha = 0.4
  projection = "lambert"
  parameters = c(51, 51)
  orientation = c(90, 0, 0)
  by.type = FALSE
  origin = TRUE



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

  Args <- list(traj = PH_GDAS1_100m_120hr_2010,
               method = "Euclid",
               n.cluster = 5,
               plot = TRUE,
               showplot = FALSE,
               type = "default",
               cols = "Set1",
               split.after = FALSE,
               map.fill = TRUE,
               map.cols = "grey40",
               map.alpha = 0.4,
               projection = "lambert",
               parameters = c(51, 51),
               orientation = c(90, 0, 0),
               by.type = FALSE,
               origin = TRUE)

  ## set graphics
  current.strip <- trellis.par.get("strip.background")
  current.font <- trellis.par.get("fontsize")

  ## reset graphic parameters
  on.exit(trellis.par.set(

    fontsize = current.font
  ))

  ## label controls
  Args$plot.type <- if ("plot.type" %in% names(Args)) {
    Args$plot.type
  } else {
    Args$plot.type <- "l"
  }
  Args$lwd <- if ("lwd" %in% names(Args)) {
    Args$lwd
  } else {
    Args$lwd <- 4
  }

  if ("fontsize" %in% names(Args)) {
    trellis.par.set(fontsize = list(text = Args$fontsize))
  }

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

  ## this bit decides whether to separately calculate trajectories for each level of type

  # if (split.after) {
  #   traj <- group_by(traj, default) %>%
  #     do(calcTraj(.))
  #   traj <- cutData(traj, type)
  # } else {
    traj <- cutData(traj, type)

    traj <- traj %>%
      group_by(across(type)) %>%
      calcTraj()
  # }

  # trajectory origin
  origin_xy <- head(subset(traj, hour.inc == 0), 1) ## origin
  tmp <- mapproject(
    x = origin_xy[["lon"]][1],
    y = origin_xy[["lat"]][1],
    projection = projection,
    parameters = parameters,
    orientation = orientation
  )
  receptor <- c(tmp$x, tmp$y)

  if (plot) {
    ## calculate the mean trajectories by cluster

    vars <- c("lat", "lon", "date", "cluster", "hour.inc", type)
    vars2 <- c("cluster", "hour.inc", type)

    agg <- select(traj, vars) %>%
      group_by(across(vars2)) %>%
      summarise(across(everything(), mean))

    # the data frame we want to return before it is transformed
    resRtn <- agg

    ## proportion of total clusters

    vars <- c(type, "cluster")

    clusters <- traj %>%
      group_by(across(vars)) %>%
      tally() %>%
      mutate(freq = round(100 * n / sum(n), 1))

    ## make each panel add up to 100
    if (by.type) {
      clusters <- clusters %>%
        group_by(across(type)) %>%
        mutate(freq = 100 * freq / sum(freq))

      clusters$freq <- round(clusters$freq, 1)
    }

    ## make sure date is in correct format
    class(agg$date) <- class(traj$date)
    attr(agg$date, "tzone") <- "GMT"

    ## xlim and ylim set by user
    if (!"xlim" %in% names(Args)) {
      Args$xlim <- range(agg$lon)
    }

    if (!"ylim" %in% names(Args)) {
      Args$ylim <- range(agg$lat)
    }

    ## extent of data (or limits set by user) in degrees
    trajLims <- c(Args$xlim, Args$ylim)

    ## need *outline* of boundary for map limits
    Args <- setTrajLims(traj, Args, projection, parameters, orientation)

    ## transform data for map projection
    tmp <- mapproject(
      x = agg[["lon"]],
      y = agg[["lat"]],
      projection = projection,
      parameters = parameters,
      orientation = orientation
    )
    agg[["lon"]] <- tmp$x
    agg[["lat"]] <- tmp$y

    plot.args <- list(
      agg,
      x = "lon", y = "lat", group = "cluster",
      col = cols, type = type, map = TRUE, map.fill = map.fill,
      map.cols = map.cols, map.alpha = map.alpha,
      projection = projection, parameters = parameters,
      orientation = orientation, traj = TRUE, trajLims = trajLims,
      clusters = clusters, receptor = receptor,
      origin = origin
    )

    ## reset for Args
    plot.args <- listUpdate(plot.args, Args)

    ## plot
    if(isTRUE(showplot)){
      plt <- do.call(scatterPlot, plot.args)
    }
    ## create output with plot
    output <- list(plot = plt, data = list(traj = traj, results = resRtn), call = match.call())
    class(output) <- "openair"
    invisible(output)
  } else {
    ## create output without plot
    return(traj)
  }

}
