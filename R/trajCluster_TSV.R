##' Calculate the TSV values for a set of cluster analyses.
##'
##' Uses \code{\link{trajCluster}} and \code{\link{cluster_TSV}} to generate clusters and
##' their corresponding total spatial variance (TSV), giving an indication as to how well a given
##' clustering attempt reduces relative spatial variance. Generally it is best to calculate the
##' TSV for a series of clusters, e.g. 2 through 20, and then identify the points at which
##' Variance increases suddenly. See e.g. Kassomenos et al., 2010.
##'
##' @param traj A trajectory results data.frame, formatted for openair.
##' @param verbose TRUE or FALSE to display progress text.
##' @param nclusters a numeric vector, containing the number of clusters to run. E.g. seq(2,10,1) will calculate TSV for 2 through 10 clusters.
##' @param normalise TRUE or FALSE to normalise calculated TSV values. Non-normalised values are left as total summed haversine distances.
##'
##' @importFrom dplyr relocate
##' @importFrom dplyr mutate
##'
##' @export
##'
trajCluster_TSV <- function(traj,
                            verbose = TRUE,
                            nclusters = 5,
                            normalise = TRUE){
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
    relocate(clusters) %>%
    mutate(indexc = row_number())
  if(isTRUE(normalise)){
    maxTSV <- max(TSVvals$TSV, na.rm = T)
    TSVvals <- TSVvals %>%
      mutate(TSV = TSV/maxTSV)
  }
  TSVvals
}
