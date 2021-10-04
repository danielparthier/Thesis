###### Parameter post processing

##### extract population parameter for groups
PopulationParameter <- function(Model, ParameterString, PopTransform="none", ParTransform="none", MatrixString, MatrixCol, MatrixRows) {
  PopPar <- posterior::as_draws_matrix(Model$draws(ParameterString))
  switch (ParTransform,
          "none" = {},
          "log" = {PopPar <- log(PopPar)},
          "logit" = {PopPar <- boot::logit(PopPar)
          "neg" = {PopPar <- -PopPar}}
  )
  PopPar <- matrix(data = PopPar, nrow = length(PopPar), ncol = length(MatrixRows)) + posterior::as_draws_matrix(Model$draws(paste0(MatrixString, "[",MatrixRows,",",MatrixCol,"]")))
  switch (PopTransform,
          "none" = PopPar,
          "exp" = exp(PopPar),
          "invLogit" = boot::inv.logit(PopPar)
  )
}

PopulationParameterCond <- function(Model, ParameterString, PopTransform="none", ParTransform="none", MatrixString, MatrixCol, CondMatString, CondMatCol) {
  PopPar <- posterior::as_draws_matrix(Model$draws(ParameterString))
  switch (ParTransform,
          "none" = {},
          "log" = {PopPar <- log(PopPar)},
          "logit" = {PopPar <- boot::logit(PopPar)
          "neg" = {PopPar <- -PopPar}}
  )
  
  MatrixMat <- posterior::as_draws_matrix(Model$draws(MatrixString))
  MatrixMat_name <- colnames(MatrixMat)
  MatrixMat <- MatrixMat[,grepl(pattern = paste0(",",MatrixCol,"]"), x = MatrixMat_name)]
  
  CondMat <- rowMeans(sapply(X = seq_along(CondMatString), function(x) {
    CondMatrixMat <- posterior::as_draws_matrix(Model$draws(CondMatString[x]))
    CondMatrixMat_name <- colnames(CondMatrixMat)
    rowMeans(CondMatrixMat[,grepl(pattern = paste0(",",CondMatCol[x],"]"), x = CondMatrixMat_name)])
  }))
  
  
  PopPar <- matrix(data = PopPar+CondMat, nrow = length(PopPar), ncol = dim(MatrixMat)[2]) + MatrixMat
  switch (PopTransform,
          "none" = PopPar,
          "exp" = exp(PopPar),
          "invLogit" = boot::inv.logit(PopPar)
  )
}

##### compare all groups to each other
DifferenceMat <- function(x, GroupStrings="", Difference="-") {
  dim2 <- dim(x)[2]
  OutputMat <- matrix(data = NaN, nrow = dim(x)[1],ncol = sum(lower.tri(matrix(data = 0, nrow = dim2, ncol = dim2))))
  RowVec  <- matrix(data = 1:dim2, nrow = dim2, ncol = dim2, byrow = F)[lower.tri(matrix(data = 0, nrow = dim2, ncol = dim2))]
  ColVec <- matrix(data = 1:dim2, nrow = dim2, ncol = dim2, byrow = T)[lower.tri(matrix(data = 0, nrow = dim2, ncol = dim2))]
  
  switch(Difference,
         "-" = {
           for(i in   seq_along(RowVec)) {
             OutputMat[,i] <- x[,RowVec[i]]-x[,ColVec[i]]
           }
         },
         "LogOdds" = {
           for(i in   seq_along(RowVec)) {
             OutputMat[,i] <- log((x[,RowVec[i]]/(1-x[,RowVec[i]]))/(x[,ColVec[i]]/(1-x[,ColVec[i]])))
           }
         }
  )
  
  Colnames <- outer(X = paste0(GroupStrings,"-"), Y = GroupStrings, FUN = "paste0")[lower.tri(matrix(data = 0, nrow = dim2, ncol = dim2))]
  colnames(OutputMat) <- Colnames
  OutputMat
}

##### HDI calculation
HDIcalc <- function(x, ci=0.95, TestValue=0) {
  outX <- t(apply(X = x, MARGIN = 2, FUN = function(x) {
    HDIframe <- bayestestR::hdi(x = x, ci=ci)
    c(`2.5%`=HDIframe$CI_low, `50%`=median(x), `97.5%`=HDIframe$CI_high)}))
  OutOfInterval <- apply(X = sign(outX-TestValue), MARGIN = 1, function(x) {max(x)==min(x)})
  Output <- data.table::data.table(cbind(outX,OutOfInterval), keep.rownames = T)
  data.table::setnames(x = Output, "rn", "GroupComparison")
  Output
}

##### IgorPxp Processing
pxp2dt <- function(x, sweepString) {
  SweepNames <- grep(pattern = sweepString, x = names(x), value = T)
  DataTable <- rbindlist(lapply(X = SweepNames, FUN = function(i) {
    Trace <- x[[i]]
    sfA <- unlist(attr(Trace, "WaveHeader")$sfA[1])
    TraceTable <- data.table(amp = as.vector(Trace), sweepName = i)
    TraceTable[,`:=`(index = .I, time=.I*sfA, sweepNr = strtoi(gsub(pattern = "[[:alpha:]]|[[:punct:]]", x = sweepName, replacement = ""))),]
    TraceTable
  }))
  return(DataTable)
}
