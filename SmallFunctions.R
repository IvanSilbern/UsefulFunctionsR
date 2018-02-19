#Different small functions

splitName <- function(x, return.position = 1, split = "", fixed = T) {
  #### function calls strsplit function for each of the elements of the character vector x
  #### and returns the n-th element after splitting that is provided by the "return.position" argument
  #### function accepts only a character verctor x and an integer vector of length 1 for "return.position"
  
  
  if(!is.character(x)) stop("vector is not a character")
  if(length(return.position) > 1 | length(return.position) == 0)  stop("length of the \"return.position\" argument is not 1")
  return.position <- as.integer(return.position)
  if(return.position < 1 | is.na(return.position)) stop("\"return.position\" argument cannot be less than 1")
  
  unlist(
    
    lapply(
      
      sapply(x, FUN = strsplit, split = split, fixed = fixed),
      
      "[", return.position
      
    )
    
  )
  
}

loadData <- function(path, name, ...){
  
  df <- read.delim(path, stringsAsFactors = F, ...)
  cat(path, "table dimensions", dim(df))
  df
  
}

deleteFullNa <- function(data){
  
  col.number <- ncol(data)
  sum.na     <- apply(data, 1, function(x) sum(is.na(x)))
  to.delete  <- which(sum.na == col.number) 
  
}

computePercent <- function(data, percent){
  
  apply(data, 2, function(x) quantile(x, probs = percent/100, na.rm = T))
  
}
