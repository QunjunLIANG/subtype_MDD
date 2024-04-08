# color identification --------------------------------------------------------
#
ColorForGroup <- function(subtype){
  if (subtype == "saMDD") {
    return("#ED0000FF")
  }
  if (subtype == "sdMDD") {
    return("#42B540FF")
  }
  if (subtype == "mdMDD") {
    return("#00468BFF")
  }
  return("#0099B4FF")
}

show_col(pal_lancet("lanonc")(4))
