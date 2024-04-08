#######################################
#
# This function is used to obtain the 
# time delay projecton maps from inputs
# while matching the selected participants
#
# This function can load the data of TDp and 
# ETS. Indicates what kinds of data you want to
# load in [type] parameter

ObtainBrainData_individual <- function(file){
  sbj_tmp <- stringr::str_extract(file, pattern = 'sub-[0-9]*')
  td_tmp <- readr::read_csv(file, col_names = F, show_col_types = FALSE)
  dat_td <- data.frame(participant_id = sbj_tmp, td_tmp)
  return(dat_td)
}