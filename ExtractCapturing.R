extractCapturing <- function(pattern, text, n.capture = 2, ...){
  
  unlist(
    
    lapply(
      
      regmatches(text, regexec(pattern, text, perl = T, ...)),
      
      FUN = "[", n.capture
      
    )
  )
  
}
