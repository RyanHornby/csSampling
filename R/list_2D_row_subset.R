# list_2D_row_subset: a helper function to convert list of MCMC output from Stan to input for unconstrained transformation
list_2D_row_subset <- function (nmlist, rindex)
{
  temp_list <- list()
  for (k in 1:length(nmlist)) {
    tmpdim <- dim(nmlist[[k]])
    ldim <- length(tmpdim)
    lcommas <- paste(rep(",", ldim - 1), collapse = " ")
    #copy over to new list - drop = FALSE retains ALL dimensions
    eval(parse(text = paste("temp_list$", names(nmlist)[k],
                            " <- ", "(nmlist$", names(nmlist)[k],
                            ")[rindex", lcommas, ",drop = FALSE]", sep = "")))
    #drop only first dimension of array - not others of size 1
    if(ldim > 1){
      eval(parse(text = paste("temp_list$", names(nmlist)[k],
                              " <- ", "array(temp_list$", names(nmlist)[k], ", dim = tmpdim[-1])", sep = "")))
    }
    #if only had 1 dim which is the MCMC draw, make a double (no dim), rather than an array of dim 1 or 0
    if(ldim == 1){
      eval(parse(text = paste("temp_list$", names(nmlist)[k],
                              " <- ", "as.double(temp_list$", names(nmlist)[k], ")", sep = "")))
    }
  }
  return(temp_list)
}

