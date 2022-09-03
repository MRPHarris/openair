##' Calculate the Total Spatial Variance TSV for a given trajectory cluster dataset.
##'
##' Uses an adaptation of the intra-cluster variance quantification method employed by
##' the HYSPLIT clustering algorithm. See \code{\link{https://www.ready.noaa.gov/documents/Tutorial/html/traj_cluseqn.html}}
##' for the inspiration. It uses the haversine calculation to ensure distances between trajectory and cluster mean endpoints
##' are accurate at higher latitudes.
##'
##' @param cluster_res An output from trajCluster().
##' @param output one of either 'CSV' or 'TSV'. In the former case, a list of CSV values for all clusters is returned.
##' @param cluster_sp NULL or numeric to calculate CSV for a specific cluster.
##' @param verbose TRUE or FALSE to display progress text.
##'
##' @importFrom spatialrisk haversine
##' @importFrom dplyr select
##' @importFrom dplyr filter
##'
##' @export
##'
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
          mutate(SV = haversine(lat_from = lat,
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
        mutate(SV = haversine(lat_from = lat,
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
