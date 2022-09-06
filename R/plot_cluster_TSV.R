##' Plot calculated TSV and change-in-TSV-% values for a given cluster results grob
##'
##' Takes results from \code{\link{trajCluster_TSV}} and plots the change-in-TSV percentage
##' between different n clusters. Options to set a threshold in % variance change to visualise
##' potentially desirable n cluster specifications
##'
##' @param TSV_res Data.frame, an output from trajCluster_TSV.
##' @param threshold numeric, a percentage change threshold to highlight in the plot.
##' @param max_clusters maximum number of clusters to show on the graph. Will shrink if a smaller-than-specified TSV_res frame is supplied.
##' @param max_change y-axis on the percent change TSV plot.
##'
##' @importFrom dplyr arrange
##' @importFrom dplyr mutate
##' @importFrom cowplot theme_cowplot
##' @importFrom cowplot plot_grid
##'
##' @export
##'
plot_cluster_TSV <- function(TSV_res,
                             threshold = 5,
                             max_clusters = 30,
                             max_change = 10){
  theme_rem_xaxis <- theme(
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
  if(max_clusters > max(TSV_res$clusters, na.rm = T)){
    max_clusters <- max(TSV_res$clusters, na.rm = T)
  }
  if(max_clusters <= 10){
    xlabel_offset = 0.5
  } else if(max_clusters > 10 && max_clusters < 30) {
    xlabel_offset = 0.75
  } else {
    xlabel_offset = 1
  }
  ## formatting
  TSV_res <- TSV_res %>%
    arrange(desc(clusters)) %>%
    mutate(pct_change = as.numeric((TSV/lag(TSV) * 100)-100)) %>%
    mutate(TSV_norm = TSV/max(TSV, na.rm = T))
  TSV_res <- TSV_res %>% # Split these two as having them in the one pipe caused errors
    mutate(flag_threshold = ifelse(lead(TSV_res$pct_change) > threshold,1,0))
  ## Create lines
  xvals_pct_change <- TSV_res$clusters[which(TSV_res$pct_change > threshold)] + 1
  yvals_pct_change <- TSV_res$pct_change[which(TSV_res$pct_change > threshold) - 1]
  lines_pct_change <- data.frame(x = rep(xvals_pct_change,2),y = c(yvals_pct_change, rep(0,length(yvals_pct_change))),
                                 np = rep(seq(1,length(xvals_pct_change),1),2))
  p_TSV <- ggplot() +
    geom_line(data = TSV_res, aes(x = clusters, y = TSV_norm)) +
    geom_point(data = TSV_res, aes(x = clusters, y = TSV_norm)) +
    geom_point(data = TSV_res[which(TSV_res$flag_threshold == 1),], aes(x = clusters, y = TSV_norm), colour = 'red') +
    theme_cowplot(12) +
    scale_x_continuous(expand = c(0,0), limits = c(0,max_clusters + 1)) +
    scale_y_continuous(expand = c(0,0.01)) +
    labs(x = "n clusters", y = "Normalised TSV ") +
    theme_rem_xaxis +
    theme(plot.margin = margin(t = 7, r = 7, b = -2, l = 7))
  p_pctTSV <- ggplot() +
    geom_line(data = TSV_res, aes(x = clusters, y = pct_change)) +
    geom_point(data = TSV_res, aes(x = clusters, y = pct_change)) +
    scale_x_continuous(expand = c(0,0), limits = c(0,max_clusters + 1), breaks = seq(0,max_clusters,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,max_change), breaks = seq(0,max_change,1), position = 'right') +
    ## change points
    geom_line(data = lines_pct_change, aes(x,y, group = np), linetype = 'dashed', colour = 'red') +
    geom_point(data = lines_pct_change[1:(nrow(lines_pct_change)/2),], aes(x,y), colour = 'red') +
    geom_label(data = lines_pct_change[1:(nrow(lines_pct_change)/2),], aes(x + xlabel_offset, y + 1, label = x),
               fill = NA, colour = 'red', label.size = NA) +
    geom_hline(yintercept = as.numeric(threshold), linetype = 'dotted') +
    labs(y = "Change in TSV (%)", x = "Number of clusters") +
    theme_cowplot(12) +
    theme(plot.margin = margin(t = 0, r = 7, b = 7, l = 7),
          axis.text.y = element_text(colour = c("black","NA")),
          axis.text.x = element_text(colour = c("black","NA","NA","NA","NA")))
  plot_grid(p_TSV,p_pctTSV,
                     nrow = 2,ncol = 1,align = 'v')
}
