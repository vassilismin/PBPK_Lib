AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
}

# # Example and testing of the faunction
# test_x <- c(0,1,2,3)
# test_y <- c(0,1,2,3)
# 
# # For the given test_x and test_y the AUC=4.5
# AUC(test_x, test_y)