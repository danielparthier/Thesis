### Plot Function
TufteBoxPlot <- function(x,
                         LabelX = NULL,
                         LabelY = NULL,
                         RangeY=NULL,
                         BreaksYn = 5,
                         BreaksY=NULL,
                         Grouping=NULL,
                         FacetGroup=F,
                         DataY=NULL,
                         DataGroup=NULL,
                         DataStat="median",
                         Dodge=F,
                         DataJitter=0.04,
                         LineWidth=4,
                         ShortWidth=0.5,
                         DataSpacing=0.25,
                         JitterSeed=NULL,
                         Horizontal=F,
                         Line=NULL,
                         HDI=F,
                         colour="gray90") {
  # browser()
  quantileIn <- c(0.005, 0.025, 0.1, 0.5, 0.9, 0.975, 0.995)
  xDT <- data.table::data.table(x)
  xDT <- data.table::melt.data.table(data = xDT, measure.vars = colnames(xDT))
  if(is.null(RangeY)) {
    RangeY <- c(min(xDT$value)*1.1, max(xDT$value)*1.1)
  }
  if(any(is.na(x))) {
    message("Remove NA")
  }
  if(HDI) {
    xDT <- xDT[,.(bayestestR::hdi(x = value, ci=0.99)$CI_low,bayestestR::hdi(x = value, ci=0.95)$CI_low,bayestestR::hdi(x = value, ci=0.8)$CI_low,median(value),bayestestR::hdi(x = value, ci=0.8)$CI_high, bayestestR::hdi(x = value, ci=0.95)$CI_high, bayestestR::hdi(x = value, ci=0.99)$CI_high),by="variable"]
  }else {
    xDT <- xDT[,lapply(X = quantileIn, function(x) {quantile(x = value, probs = x, na.rm = T)}),by="variable"]
  }
  data.table::setnames(x = xDT, new = c("variable", paste0("Quant_", quantileIn)))
  xDT[,xPos:=.I+0.5,]
  
  if(length(Grouping)>0) {
    GroupingL <- length(Grouping)
  } else {
    GroupingL <- NULL
  }
  
  if(length(Grouping)==xDT[,.N,]) {
    xDT[,Grouping:=Grouping,]
    if(FacetGroup) {
      xDT[,xPos:=1,]
    }
  }
  
  PlotOut <- ggplot(data = xDT, mapping = aes(y=`Quant_0.5`, x=xPos, group=Grouping))
  if(is.numeric(Line)) {
    PlotOut <- PlotOut + geom_hline(yintercept = Line, colour="gray50")
  }
  PlotOut <- PlotOut+
    geom_linerange(mapping = aes(ymin=`Quant_0.975`, ymax=`Quant_0.995`), size=0.25)+
    geom_linerange(mapping = aes(ymin=`Quant_0.025`, ymax=`Quant_0.005`), size=0.25)+
    geom_linerange(mapping = aes(ymin=`Quant_0.1`, ymax=`Quant_0.9`),size=LineWidth, colour=colour)+
    geom_segment(mapping = aes(x = xPos, xend = xPos, y=`Quant_0.5`-grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY), yend=`Quant_0.5`+grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY)), size=LineWidth)+
    coord_cartesian(ylim = RangeY)+
    theme_classic()
  
  if(!is.null(DataY)) {
    if(Dodge) {
      if(DataSpacing<0.5 & DataSpacing>0) {
        DataDodge <- 0.5+DataSpacing
      } else {
        message("Spacing out of range (>0.5 and <0): Reset to default.")
        DataDodge <- 0.75
      }
      
    } else {
      DataDodge <- 0.5
    }
    if(is.null(DataGroup)) {
      DataDT <- data.table::data.table(Data_Y = DataY, xPos = DataDodge+1)
    } else {
      DataDT <- data.table::data.table(Data_Y = DataY, xPos = as.integer(factor(DataGroup))+DataDodge)
      if(length(Grouping)>0) {
        for(i in seq_along(Grouping)) {
          tmpGroup <- levels(factor(Grouping))[i]
          DataDT[xPos==i+Dodge,Grouping:=tmpGroup,] 
        }
      } else {
        DataDT[,`:=`(Grouping=""),]
      }
    }
    if(!is.null(JitterSeed)) {
      set.seed(JitterSeed)
    }
    PlotOut <- PlotOut + geom_jitter(data = DataDT, mapping = aes(y=Data_Y, x=xPos, group=Grouping), height = 0, width = DataJitter, alpha=0.1, shape=19, size=0.5, inherit.aes = F)
    switch(EXPR = DataStat,
           "median" = {StatDT <- DataDT[,.(Data_Y=median(Data_Y)),by=xPos]
           PlotOut <- PlotOut + geom_segment(data = StatDT, mapping = aes(x = xPos, xend = xPos, y=Data_Y-grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY), yend=Data_Y+grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY)), size=LineWidth)},
           "mean" = {StatDT <- DataDT[,.(Data_Y=mean(Data_Y)),by=xPos]
           PlotOut <- PlotOut + geom_segment(data = StatDT, mapping = aes(x = xPos, xend = xPos, y=Data_Y-grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY), yend=Data_Y+grid::convertUnit(x = unit(x = 2, units = "pt"), unitTo = "npc", valueOnly = T)*ShortWidth*diff(RangeY)), size=LineWidth)},
           "none" = {})
  }
  

  if(length(Grouping)>0) {
    if(FacetGroup) {
      Xname <- LabelX
      Xbreaks <- NULL
      Xlabels <- NULL
      Xlimits <- c(0.5,1.5)
      PlotOut <- PlotOut + facet_wrap(~Grouping)     
    } else {
      if(exists(x = "DataDodge")) {
        LabelCenter <- (DataDodge - 0.5)/2
      } else {
        LabelCenter <- 0
      }
      Xname <- LabelX
      Xbreaks <- seq_along(Grouping)+0.5+LabelCenter
      Xlabels <- Grouping
      Xlimits <- c(1, xDT[unique(xPos),.N,]+1)
  #    PlotOut <- PlotOut + scale_x_continuous(name = LabelX, breaks = seq_along(Grouping)+0.5+LabelCenter, labels=Grouping, limits = c(1, xDT[unique(xPos),.N,]+1))
    }
  } else {
    Xname <- LabelX
    Xbreaks <- NULL
    Xlabels <- NULL
    Xlimits <- c(1, xDT[unique(xPos),.N,]+1)
  }
  if(Horizontal) {
    if(is.null(BreaksY)) {
      PlotOut + coord_flip() + scale_x_continuous(name = Xname, breaks = Xbreaks, labels=Xlabels,limits = Xlimits) + scale_y_continuous(name = LabelY,limits = RangeY, n.breaks = BreaksYn) + theme(axis.line = element_blank(), axis.ticks.y = element_blank(),   strip.background = element_blank()) 
    } else {
      PlotOut + coord_flip() + scale_x_continuous(name = Xname, breaks = Xbreaks, labels=Xlabels,limits = Xlimits) + scale_y_continuous(name = LabelY,limits = RangeY,  breaks = BreaksY) + theme(axis.line = element_blank(), axis.ticks.y = element_blank(),   strip.background = element_blank()) 
    }
    
  } else {
    if(is.null(BreaksY)) {
      PlotOut + scale_x_continuous(name = Xname, breaks = Xbreaks, labels=Xlabels,limits = Xlimits) + scale_y_continuous(name = LabelY,limits = RangeY, n.breaks = BreaksYn) + theme(axis.line = element_blank(), axis.ticks.x = element_blank(),   strip.background = element_blank()) 
    } else {
      PlotOut + scale_x_continuous(name = Xname, breaks = Xbreaks, labels=Xlabels,limits = Xlimits) + scale_y_continuous(name = LabelY,limits = RangeY,  breaks = BreaksY) + theme(axis.line = element_blank(), axis.ticks.x = element_blank(),   strip.background = element_blank()) 
    }
  }
}

