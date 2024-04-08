##################################
#
# This function is used to add the 
# network score to the TDp
#
# Net_annotation is the dataframe obtained
# from ObtainNetID
#

AddNetworkScore_Schaefer <- function(dat_TDp, net_anna, prefix=1){
  dat_tmp <- dat_TDp
  net_annotation <- net_anna
  network_name <- unique(net_annotation$network) %>% as.character()
  network_pos <- net_annotation$network
  
  dat_tmp <- dat_tmp %>% 
    # Mean
    mutate(td_VIS = rowMeans(across(all_of(prefix + which(network_pos=="VIS"))))) %>%
    mutate(td_SMN = rowMeans(across(all_of(prefix + which(network_pos=="SMN"))))) %>%
    mutate(td_DAN = rowMeans(across(all_of(prefix + which(network_pos=="DAN"))))) %>%
    mutate(td_SAN = rowMeans(across(all_of(prefix + which(network_pos=="SAN"))))) %>%
    mutate(td_LIM = rowMeans(across(all_of(prefix + which(network_pos=="LIM"))))) %>%
    mutate(td_FPN = rowMeans(across(all_of(prefix + which(network_pos=="FPN"))))) %>%
    mutate(td_DMN = rowMeans(across(all_of(prefix + which(network_pos=="DMN")))))
  
  return(dat_tmp)
}