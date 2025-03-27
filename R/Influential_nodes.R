#' IVI values
#'
#' @param network_ls network list
#' @param group group name
#' @param days days
#'
#' @return data frame with IVI values
#' @export
#'
#' @importFrom influential ivi
IVI_values<- function(network_ls, group, days){

  num_net<-length(network_ls)
  df_ivi<-data.frame(Group=NA, Day=NA)

  index=0
  for(i in 1:num_net){
    index=index+1
    df_ivi[index,1]<-group
    df_ivi[index,2]<-days[i]
    net_ig<-graph_from_adjacency_matrix(network_ls[[i]], mode = "undirected")
    df_ivi[index,3:(dim(network_ls[[1]])[2]+2)]<-ivi(graph = net_ig)
  }
  return(df_ivi)
}
