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

## Calculating intra-cluster spatial variance (SV):
add_traj_identifier <- function(x, ntraj_1 = TRUE){
  # This function assumes that convert_openair() has been used on the dataset.
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


## Extract points for a single cluster.
clus1_pts <- clus1$data$traj %>%
  dplyr::filter(cluster == "C1") %>%
  add_traj_identifier()
clus1_mean <- clus1$data$results %>%
  dplyr::filter(cluster == "C1") %>%
  arrange(hour.inc)

## Now calculate variance
# Each trajectory j gets an SV val.

traj <- PH_GDAS1_100m_120hr_2010 %>%
  add_traj_identifier()
length(unique(test$trajectory))

j = 1

clust_ntraj = max(clus1_pts$trajectory, na.rm = T)

trajCluster_TSV <- function(traj,
                            verbose = T,
                            nclusters = 5,
                            normalise = T){
  TSVlist <- vector('list',length = length(nclusters))
  for(clus in seq_along(TSVlist)){
    if(isTRUE(verbose)){message("Starting calculations for ",nclusters[clus]," clusters.")}
    clus_it <- trajCluster(traj = traj,
                           n.cluster = nclusters[clus],
                           plot = TRUE, showplot = FALSE)
    TSVlist[[clus]] <- cluster_TSV(cluster_res = clus_it,
                          output = 'TSV',
                          verbose = FALSE)
    if(isTRUE(verbose)){message("Calculations complete for ",nclusters[clus]," clusters.")}
  }
  TSVvals <- unlist(TSVlist) %>%
    data.frame() %>%
    mutate(clusters = nclusters) %>%
    'colnames<-'(c('TSV','clusters')) %>%
    dplyr::relocate(clusters) %>%
    mutate(indexc = row_number())
  if(isTRUE(normalise)){
    maxTSV <- max(TSVvals$TSV, na.rm = T)
    TSVvals <- TSVvals %>%
      mutate(TSV = TSV/maxTSV)
  }
  TSVvals
}



PH_GDAS_2010_TSV <- trajCluster_TSV(traj = PH_GDAS1_100m_120hr_2010,
                                    nclusters = seq(2,30,1))

PH_GDAS_2010_TSV_5rpt <- trajCluster_TSV(traj = PH_GDAS1_100m_120hr_2010,
                                    nclusters = rep(5,20))

plot_cluster_TSV <- function(TSVres){
  plotobj <- ggplot() +
    geom_line(data = testTSV, aes(x = clusters, y = TSV)) +
    geom_point(data = testTSV, aes(x = clusters, y = TSV)) +
    theme_cowplot(12)



}

savedir <- 'C:/Users/matth/Desktop/University Files/mh phd 2020/3 thesis/Chapter X Contemporary transport dynamics/trajectory data/TSV calculations/'

saveRDS(PH_GDAS_2010_TSV_5rpt, file = paste0(savedir,"PH_GDAS_2010_100m_TSV_5rpt.rds"))

library(cowplot)

ggplot() +
  geom_line(data = PH_GDAS_2010_TSV_5rpt, aes(x = indexc, y = TSV)) +
  geom_point(data = PH_GDAS_2010_TSV_5rpt, aes(x = indexc, y = TSV)) +
  theme_cowplot(12)





cluster_TSV <- function(cluster_res,
                        output = 'TSV',
                        cluster_sp = NULL,
                        verbose = FALSE){
  ## Check for existence of trajectory identifier column
  if('trajectory' %in% colnames(cluster_res$data$traj)){
    cluster_res$data$traj <- cluster_res$data$traj %>%
      dplyr::select(-c(trajectory))
  }
  ## Two options:
  # TSV across all clusters, or CSV of an individual cluster.
  if(is.null(cluster_sp)){
    # Get n clusters and names of clusters
    clusters <- unique(cluster_res$data$results$cluster)
    nclusters <- length(clusters)
    # Init list for cluster iterations
    CSVlist <- vector('list',length = nclusters)
    for(i in seq_along(CSVlist)){
      if(isTRUE(verbose)){message("Calculating CSV for ",clusters[i])}
      # obtain endpoints for this cluster
      clus_it_pts <- cluster_res$data$traj %>%
        dplyr::filter(cluster == clusters[i]) %>%
        add_traj_identifier()
      # Obtain cluster mean points
      clus_it_mean <- cluster_res$data$results %>%
        dplyr::filter(cluster == clusters[i]) %>%
        arrange(hour.inc)
      # How many trajectories in this cluster?
      clust_ntraj = max(clus_it_pts$trajectory, na.rm = T)
      traj_list <- vector('list',length = clust_ntraj)
      for(j in seq_along(traj_list)){
        # For each trajectory, obtain the haversine distance between each endpoint and the corresponding cluster mean point
        endpt_it <- clus_it_pts %>%
          dplyr::filter(trajectory == j) %>%
          mutate(clusmeanpt_lat = clus_it_mean$lat) %>%
          mutate(clusmeanpt_lon = clus_it_mean$lon) %>%
          mutate(SV = spatialrisk::haversine(lat_from = lat,
                                             lon_from = lon,
                                             lat_to = clusmeanpt_lat,
                                             lon_to = clusmeanpt_lon))
        # obtain Spatial Variance SV for this trajectory
        traj_list[[j]] <- sum(endpt_it$SV)
      }
      # Sum SV for each trajectory to obtain Cluster Spatial Variance CSV
      CSVlist[[i]] <- sum(unlist(traj_list), na.rm = T)
      if(isTRUE(verbose)){message("Cluster ",i," complete. CSV: ",CSVlist[[i]])}
    }
    # Sum CSV for each cluster to obtain Total Spatial Variance TSV for this n clusters
    TSV <- sum(unlist(CSVlist), na.rm = T)
    if(isTRUE(verbose)){message("TSV: ",TSV)}
    if(output == 'TSV'){
      return(TSV)
    } else if(output == 'CSV'){
      names(CSVlist) <- clusters
      return(CSVlist)
    }
  } else {
    ## Individual cluster.
    if(length(cluster_sp) > 1){
      stop("Only a single cluster_sp is supported at this stage.")
    }
    if(isTRUE(verbose)){message("Calculating CSV for ",clusters[i])}
    # obtain endpoints for this cluster
    clus_it_pts <- cluster_res$data$traj %>%
      dplyr::filter(cluster == cluster_sp) %>%
      add_traj_identifier()
    # Obtain cluster mean points
    clus_it_mean <- cluster_res$data$results %>%
      dplyr::filter(cluster == cluster_sp) %>%
      arrange(hour.inc)
    # How many trajectories in this cluster?
    clust_ntraj = max(clus_it_pts$trajectory, na.rm = T)
    traj_list <- vector('list',length = clust_ntraj)
    for(j in seq_along(traj_list)){
      # For each trajectory, obtain the haversine distance between each endpoint and the corresponding cluster mean point
      endpt_it <- clus_it_pts %>%
        dplyr::filter(trajectory == j) %>%
        mutate(clusmeanpt_lat = clus_it_mean$lat) %>%
        mutate(clusmeanpt_lon = clus_it_mean$lon) %>%
        mutate(SV = spatialrisk::haversine(lat_from = lat,
                                           lon_from = lon,
                                           lat_to = clusmeanpt_lat,
                                           lon_to = clusmeanpt_lon))
      # obtain Spatial Variance SV for this trajectory
      traj_list[[j]] <- sum(endpt_it$SV)
    }
    CSV <- sum(unlist(traj_list), na.rm = T)
    return(CSV)
  }
}




Clusters_TSV <- function(traj,
                       verbose = T){
  # This function generates the total spatial variance of a given cluster.


}


## Generate clusters
Cluster_SV <- function(traj,
                       verbose = T,
                       nclust = 'ntraj'){
  if('trajectory' %in% colnames(traj)){
    traj <- traj %>%
      dplyr::select(-c(trajectory)) %>%
      add_traj_identifier()
  } else {
    traj <- traj %>%
      add_traj_identifier()
  }
  # total_clusters =
  if(nclust == 'ntraj'){
    nclusts = max(traj$trajectory, na.rm = T)
    carry_list <- vector('list', length = 2)
  } else {
    nclusts = nclust
  }
  TSV_list <- vector('list', length = nclusts - 1)
  for(cl in seq_along(TSV_list)){
    # following the HYPSLIT guidance.
    # iterations add additional clusters. i = n clusters. At each cluster step, total variance is calculated.
    # But TSV is only calculated for each merge. Thus, we need to know where each trajectory merge occurs.
    if(nclust == 'ntraj'){

      ## For each iteration, determine CSV of new additions.
      if(cl == 1){
        # Generate data as normal.
        clust_it <- trajCluster(traj,
                                n.cluster = ((nclusts-1) - (cl - 1)),
                                plot = TRUE,
                                showplot = FALSE)
        clustlist_it <- unique(clust_it$data$traj$cluster)
        clustdat <- clust_it$data$traj %>%
          add_traj_identifier() %>%
          group_by(cluster) %>%
          distinct(trajectory)
        carry_list[[1]] <- clustdat
        TSV_list[[cl]] <- cluster_TSV(clust_it, verbose = F)
      } else {
        # Generate data as normal.
        clust_it <- trajCluster(traj,
                                n.cluster = ((nclusts-1) - (cl - 1)),
                                plot = TRUE,
                                showplot = FALSE)
        clustlist_it <- unique(clust_it$data$traj$cluster)
        clustdat <- clust_it$data$traj %>%
          add_traj_identifier() %>%
          group_by(cluster) %>%
          distinct(trajectory)
        # carry_list[[1]] <- clustdat
        # Get u

        test <- clustdat %>%
          count() %>%
          arrange(desc(n))

        test2 <- carry_list[[1]] %>%
          count() %>%
          arrange(desc(n)) %>%
          head(-1)

        clusterdiff <- test$cluster[which(test$n != test2$n)]

        compdat <- carry_list[[1]]
        compvals <- compdat %in% clustdat
        TSV_list[[cl]] <- cluster_TSV(clust_it, verbose = F)



      }

      cluster_1st <- 365
      cluster_2nd <- 364



      tnum <- max(clus1_pts$trajectory, na.rm = T)
      clust_it <- trajCluster(traj,
                              n.cluster = tnum - (cl - 1),
                              plot = TRUE,
                              showplot = FALSE)
      clustlist_it <- unique(clust_it$data$traj$cluster)
      clustdat <- clust_it$data$traj %>%
        add_traj_identifier() %>%
        group_by(cluster) %>%
        distinct(trajectory)





    } else {
      clust_it <- trajCluster(traj,
                              n.cluster = nclusts[cl],
                              plot = TRUE,
                              showplot = FALSE)
    }
    TSV_list[[cl]] <- cluster_TSV(clust_it, verbose = F)
    if(isTRUE(verbose)){message("TSV calc for cluster group ",cl," complete")}
  }
  TSVtot <- unlist(TSV_list)
  return(TSVtot)
}

calcclust <- Cluster_SV(PH_GDAS1_100m_120hr_2010)


cluster_res1 <- openair::trajCluster(traj = PH_GDAS1_100m_120hr_2010,
                                    n.cluster = 2,
                                    plot = TRUE, showplot = FALSE)

cluster_res <- clust_it
output = 'CSV'
cluster_sp = NULL
verbose = T


GDAS1_2TSV <- cluster_TSV(cluster_res1)


TSV_curve <- function(traj,
                      n_clusters = seq(1,10,1),
                      verbose = T,
                      plot = FALSE){
  TSVlist <- vector('list', length = length(n_clusters))
  for(c in seq_along(clustlist)){
    if(isTRUE(verbose)){message("Calculating cluster n = ",c," for ",n_clusters[c]," clusters")}
    clust_it <- openair::trajCluster(traj = traj,
                                           n.cluster = n_clusters[c],
                                           plot = TRUE, showplot = FALSE)
    TSVlist[[c]] <- cluster_TSV(clust_it)
  }
  TSVframe <- data.frame(
    c(n_clusters),
    unlist(TSVlist)
  )
  return(TSVframe)
}

TSVsample <- TSV_curve(traj = PH_GDAS1_100m_120hr_2010)

plot(TSVsample)

endpt_it <- clus1_pts %>%
  dplyr::filter(trajectory == j) %>%
  mutate(clusmeanpt_lat = clus1_mean$lat) %>%
  mutate(clusmeanpt_lon = clus1_mean$lon) %>%
  mutate(SV = spatialrisk::haversine(lat_from = lat,
                                     lon_from = lon,
                                     lat_to = clusmeanpt_lat,
                                     lon_to = clusmeanpt_lon))
  # mutate(SVlat = calcSV_lat(endpt_lat = lat,clusmeanpt_lat = clusmeanpt_lat)) %>%
  # mutate(SVlon = calcSV_lon(endpt_lon = lon, clusmeanpt_lon = clusmeanpt_lon)) %>%
  # mutate(SV = SVlat + SVlon)

library(spatialrisk)

haversine_sets <- function(pts1,pts2){

  pts1 <- matrix(c(endpt_it$lon,endpt_it$lat),
                 nrow = nrow(endpt_it), ncol = 2)
  pts2 <- matrix(c(endpt_it$clusmeanpt_lon,endpt_it$clusmeanpt_lat),
                 nrow = nrow(endpt_it), ncol = 2)

  apply(pts1, 1, function(r){distHaversine()})

}

test <- apply(trajectory,1,fn(pt){
  distHaversine(p1)
})

test <- endpt_it %>%
  nest(data = c(lat,lon))


my_points <- matrix(c(77.65465, 91.54323,    # Create longitude/latitude matrix
                      21.35444, 17.65465),
                    nrow = 2)
colnames(my_points) <- c("longitude", "latitude")
rownames(my_points) <- c("pt_1", "pt_2")

ptdists_SV = function(pts_1, pts_2, method = 'haversine'){

pts_1 <- data.frame(endpt_it$lat,
                endpt_it$lon)
pts_2 <- data.frame(endpt_it$clusmeanpt_lat,
                endpt_it$clusmeanpt_lon)




}

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
