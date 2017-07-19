disaggregate2 <- function (data) 
{
  n.transitions <- lapply(as.character(data$aggr_Fragment_Annotation), 
                          function(x) strsplit(x, ";"))
  n.transitions2 <- unlist(lapply(n.transitions, function(x) length(unlist(x))))
  n.transitions3 <- lapply(as.character(data$aggr_Peak_Area), 
                           function(x) strsplit(x, ";"))
  n.transitions4 <- unlist(lapply(n.transitions3, function(x) length(unlist(x))))
  
  if (sum(n.transitions2 != n.transitions4) > 0) {
    stop(paste("The number of transitions annotated and measured do not match in the following transitions:\n", 
               paste(unlist(n.transitions[n.transitions2 != n.transitions4]), 
                     collapse = ", ")))
  }
  if (min(n.transitions2) == max(n.transitions2)) {
    message(paste("The library contains", max(n.transitions2), 
                  "transitions per precursor.\n                  
                  \nThe data table was transformed into a table containing one row per transition."))
  }
  if (min(n.transitions2) != max(n.transitions2)) {
    message(paste("The library contains between", min(n.transitions2), 
                  "and", max(n.transitions2), "transitions per precursor.
                  \nThe data table was transformed into a table containing one row per transition."))
  }
  i<-1L
  temp<-as.data.frame(NULL)
  l<-length(n.transitions)
  for(i in seq_along(n.transitions)){
    temp<-rbind(temp, data.frame(FragmentIon = unlist(n.transitions[[i]]),
                                 Intensity = unlist(n.transitions3[[i]]),
                                 stringsAsFactors = F)
    )
  }
  temp2<-data[rep(seq_len(nrow(data)),n.transitions2), ]
  names(temp2)[names(temp2)=="Intensity"]<-"TotalIntensity"
  if(nrow(temp)==nrow(temp2)){
    temp3<-cbind(temp2, temp)
  } else {
    cat("\n Number of rows does not match!")
    break
  }
  
  cols <- colnames(temp3)[colnames(temp3) %in% 
                                  c("transition_group_id", "m_score", "ProteinName", "FullPeptideName", "PeptideSequence", 
                                    "Sequence", "Charge", "PrecursorCharge", "Fragment", 
                                    "FragmentIon", "Intensity", "Condition", "BioReplicate", 
                                    "Run", "RT")]
  
  temp3<-temp3[ ,cols]
  colnames(temp3) <- gsub("FullPeptideName", "PeptideSequence", 
                                    colnames(temp3))
  colnames(temp3) <- gsub("^Charge$", "PrecursorCharge", 
                                    colnames(temp3))
  if ("Sequence" %in% cols) {
    colnames(temp) <- gsub("^Sequence$", "NakedSequence", 
                                      colnames(temp))
  }
  temp3[,"Intensity"]<-as.numeric(temp3$Intensity)
  temp3[,"RT"]<-as.numeric(temp3$RT)
  temp3[,"PrecursorCharge"]<-as.integer(temp3$PrecursorCharge)
  return(temp3)
}