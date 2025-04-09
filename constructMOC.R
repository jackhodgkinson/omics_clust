# constructMOC.R
constructMOC <- function(data              # Input as data frame or list of data frames   
                        )
   {  if (class(data) != "list") {
        data <- as.data.frame(sapply(data, as.factor))
        print(data)
        moc <- do.call(cbind, lapply(data, function(x) model.matrix(~ x - 1)))
        moc <- as.matrix(moc)
        colnames(moc) <- NULL
       }
  
      else {
       moc <- list()
       for (i in 1:length(data)) {
         data[[i]] <- as.data.frame(sapply(data[[i]], as.factor))
         moc[[i]] <- do.call(cbind, lapply(data[[i]], function(x) model.matrix(~ x - 1)))
         moc[[i]] <- as.matrix(moc[[i]])
         colnames(moc[[i]]) <- NULL
         names(moc)[[i]] <- paste0("MOC for Group",i)
        }
       }
   return(moc)
   }