###### Spike Analyis
library(data.table)
ReloadFiles <- F
SpikeTableTotal <- readRDS("Data/SpikeTableTotalV2.rds")
# SpikeTableTotal
# SpikeTableTotal[,UniqueID:=as.character(paste(ClusterFac, Session), sep = "_"),] 
# UnitOverview <- SpikeTableTotal[ClusterNr>1&!grepl(pattern = "mua", x = ClusterFac, ignore.case = T),.N,by=.(UniqueID)] 
# SpikeTableTotal[grepl(pattern = "mua",ClusterFac, ignore.case = T),ClusterNr:=1,]
# #SpikeTableTotal[!grepl(pattern = "mua",ClusterFac, ignore.case = T),ClusterNr:=1,]
# SpikeTableTotal[ClusterNr>1,ClusterNr:=strtoi(unlist(strsplit(x = as.character(ClusterFac), split = "_"))[1]),by=ClusterFac][,SpikeIdx:=as.integer(SpikeTime*2e4),]
# SpikeTableTotal[,`:=`(RecDate=strtoi(unlist(strsplit(Session, split = "_"))[2]), AnimalID=unlist(strsplit(Session, split = "_"))[1]),by=Session]
# 
# 
# SpikeTableTotal[UniqueID%in%UnitOverview$UniqueID&hyper_block_frequency<0.5, MeanFreq:=1/mean(diff(SpikeTime), na.rm = T),by=.(UniqueID)]
# SpikeTableTotal[UniqueID%in%UnitOverview$UniqueID, MeanFreq:=mean(MeanFreq, na.rm = T),by=.(UniqueID)]
# SpikeTableTotal[UniqueID%in%UnitOverview$UniqueID,]
# SpikeTableTotal[ClusterNr>1&MeanFreq<1,ClusterNr:=1,by=.(UniqueID)]
# 
# LowOccurenceDropOut <- SpikeTableTotal[ClusterNr>1, .N,by=.(UniqueID)][N<3e3, ][,UniqueID,]  
# SpikeTableTotal[UniqueID %in% LowOccurenceDropOut,ClusterNr:=1,] 


# if(ReloadFiles==T) {
#   FolderList <- list.dirs(path = "Daniel_Data/")
#   UniqueSessions <- SpikeTableTotal[,unique(Session),]
#   OutputFolderList <- grep(pattern = "OutputFolder", FolderList, value = T)
#   for(i in OutputFolderList) {
#     tmpSession <- tail(unlist(strsplit(x = dirname(i), split = "/")), n = 1)
#     if(any(grepl(pattern = tmpSession, x = UniqueSessions))) {
#       tmpStimMat <- readRDS(file = paste0(i,"/StimMatMS.rds"))
#       PulseTimings <- tmpStimMat$PulseMat[tmpStimMat$PulseMat[,16]==1,1]/2e4
#       print("compute")
#       for(pulse in seq_along(PulseTimings)) {
#  #       browser()
#         lowerB <- PulseTimings[pulse]-0.5
#         upperB <- PulseTimings[pulse]+0.5
#         PulseTime <- PulseTimings[pulse]
#         PulseID <- pulse
#         SpikeTableTotal[Session==tmpSession&SpikeTime%between%c(lowerB,upperB),`:=`(SingleStimDiff=SpikeTime-PulseTime, PulseSingleNumber=PulseID),]
#       }
#       print("done")
#     }
#   }
#   saveRDS(SpikeTableTotal, "/alzheimer/Daniel_Data/R/Thesis/Data/SpikeTableTotalV2.rds")
# }

#### isolate ms pv cells ####

MSTestunits <- SpikeTableTotal[ClusterNr>1&frequency==1&UnitLoc=="MS",.N,by=UniqueID][N>1000,UniqueID,] 

SpikeTableTotal[UniqueID%in%MSTestunits,]

SpikeTableMS1Hz[pulse_nr==3305,]

SpikeTableMS1Hz <- SpikeTableTotal[UniqueID %in% MSTestunits & frequency==1&UnitLoc=="MS",] 
SpikeTableMS1Hz <- SpikeTableMS1Hz[!is.na(SingleStimDiff),.(UniqueID, UnitLoc, SingleStimDiff, hyper_block, PulseSingleNumber, AnimalID, RecDate),] 
SpikeTableMS1Hz[,StimRoundPos := round(SingleStimDiff*2, digits = 3)/2,] 
SpikeTableMS1Hz <- SpikeTableMS1Hz[,.(Count=.N),by=.(UniqueID, StimRoundPos)]
SpikeZScoreTable <- rbindlist(l = list(data.table(UniqueID=rep(unique(SpikeTableMS1Hz$UniqueID), 
                                                               each = length(seq(-0.5,0.5,0.0005))),
                                                  StimRoundPos = rep(round(seq(-0.5,0.5,0.0005), digits = 4), times=length(unique(SpikeTableMS1Hz$UniqueID))), 
                                                  Count=0),
                                       SpikeTableMS1Hz))
SpikeZScoreTable <- SpikeZScoreTable[,.(N=sum(Count)),by=.(UniqueID, StimRoundPos)] 
SpikeZScoreTable[,ScaleN:=scale(N), by=UniqueID] 
setorder(SpikeZScoreTable, "StimRoundPos")

OrderSpikeScale <- SpikeZScoreTable[StimRoundPos %between% c(0,0.01), .(ScaleCenter = max(ScaleN), UniqueID),by=UniqueID]
setorder(OrderSpikeScale, ScaleCenter)
OrderSpikeScale[,OrderID:=.I,] 
for(i in OrderSpikeScale[,UniqueID,]){
  SpikeZScoreTable[UniqueID==i, OrderID:=OrderSpikeScale[UniqueID==i,OrderID,],] 
} 

### identify PV cells ####
ScaleCutOff <- OrderSpikeScale[ScaleCenter>8,min(OrderID),] 
PV_Units <- OrderSpikeScale[ScaleCenter>8,UniqueID,]
PV_negUnits <- OrderSpikeScale[ScaleCenter<8,UniqueID,]
SpikeTableTotal[UniqueID%in%PV_Units,PVPos:=TRUE,][UniqueID%in%PV_negUnits,PVPos:=FALSE,] 
### PV-Units per session
SpikeTableTotal[UniqueID %in% PV_Units, unique(UniqueID), by= Session][,.N,by=Session] 


#### make plot for z-score ####
SinglePulseZScorePlot <- ggplot(data = SpikeZScoreTable[StimRoundPos %between% c(-0.035,0.035)]  , aes(x=StimRoundPos*1e3, y=OrderID, fill=ScaleN))+
  geom_raster()+
  geom_vline(xintercept = 0, lty="dashed", size=0.2)+
  scale_y_continuous(expand=c(0,0), name = "Units", seq(0,120,10))+
  scale_x_continuous(expand=c(0,0), name = "Time to Stimulation (ms)", breaks = seq(-40,40,10))+
  geom_hline(yintercept = ScaleCutOff-0.5, lty="dashed", size=0.2, colour="gray")+
  #scale_fill_viridis_c()+
  scale_fill_gradient2(limits = c(-max(SpikeZScoreTable[,max(ScaleN),]),SpikeZScoreTable[,max(ScaleN),]), guide = guide_colorbar(title = "Z-Score"))+
  theme_bw()

SpikeTableTotal[UniqueID%in%MSTestunits,unique(MeanFreq),]

#### Raster plot ####
### good examples 6, 13, 17
RasterPlot <- ggplot(data = SpikeTableTotal[UniqueID%in%PV_Units[8]  &
                                              pulse_nr_block <= 300 & 
                                              hyper_block_frequency==1 & 
                                              SingleStimDiff %between% c(-0.035,0.035),.(pulse_nr_block, SingleStimDiff),],
                     aes(x=SingleStimDiff*1e3, y=pulse_nr_block))+
  geom_rect(data=data.frame(SingleStimDiff=0, pulse_nr_block=0), aes(xmin=0, xmax=2, ymin=0,ymax=Inf), alpha=1, fill="deepskyblue2")+
  geom_point(shape="|", size=0.4)+
  #facet_wrap(facets = ~ UniqueID)+
  scale_y_continuous(expand=c(0,0), name = "Pulse", limits = c(0,300), breaks = seq(0,300,100))+
  scale_x_continuous(expand=c(0,0), name = "Time to Stimulation (ms)", breaks = seq(-40,40,10))+
  theme_classic()+
  theme(strip.background = element_rect(color=NaN, fill=NaN))

PV_Unit_CCF <- rbindlist(lapply(X = SpikeTableTotal[UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS", PaSLoc),unique(UniqueID)], function(x){
  PVCCF <- InVivoR::SpikeCCF(x = SpikeTableTotal[UniqueID%in%x&stimulation_block==0, SpikeTime,])
  data.table(CCF=PVCCF$CCF/max(PVCCF$CCF), Time=PVCCF$xAxis, UniqueID=x)
}))

PV_Unit_CCF[,Time:=Time-0.0005,]
# PV_Unit_CCF_end <- PV_Unit_CCF
# PV_Unit_CCF[,Time:=Time+0.001,]
# PV_Unit_CCF <- rbindlist(l = list(PV_Unit_CCF, PV_Unit_CCF_end))


PVCCFPlotZoom <- ggplot(data = PV_Unit_CCF, mapping = aes(x = Time, y = CCF, group=UniqueID))+
  geom_rect(mapping = aes(xmin = Time, xmax = shift(Time, type = "lead"), 
                          ymin = 0, ymax = CCF))+
  #geom_ribbon(mapping = aes(x = Time, ymax=CCF, ymin=-Inf), fill="grey50")+
  #geom_step(colour="grey10")+
  facet_grid(cols = vars(UniqueID))+
  scale_x_continuous(name= "Time (s)", limits = c(-0.025,0.025))+
  scale_y_continuous(name = "Normalised Rate")+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(hjust = 0.5)) 
PVCCFPlotOut <- ggplot(data = PV_Unit_CCF, mapping = aes(x = Time, y = CCF, group=UniqueID))+
  geom_rect(mapping = aes(xmin = Time, xmax = shift(Time, type = "lead"), 
                          ymin = 0, ymax = CCF))+
  #geom_ribbon(mapping = aes(x = Time, ymax=CCF, ymin=-Inf), fill="grey50")+
  #geom_step(colour="grey10")+
  facet_grid(cols = vars(UniqueID))+
  scale_y_continuous(name = "Normalised Rate")+
  scale_x_continuous(name= "Time (s)", limits = c(-0.5,0.5))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(hjust = 0.5)) 



PVThetaPhase <- ggplot(data = SpikeTableTotal[UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS", PaSLoc),,], aes(x=PaSThetaPhase, group=UniqueID))+
  geom_histogram(aes(y=..count../sum(..count..)))+
  scale_x_continuous(name = expression(theta*" Phase"), breaks = seq(-180,180,90), labels = seq(0,360,90))+
  scale_y_continuous(name = "Density")+
  facet_grid(cols = vars(UniqueID))+
 # facet_wrap(PaSLoc~UniqueID)+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text = element_blank(), plot.title = element_text(hjust = 0.5)) 

ChR2UnitResponse <- RasterPlot + SinglePulseZScorePlot + plot_layout(widths = c(2,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = ChR2UnitResponse, "/alzheimer/Daniel_Data/R/Thesis/Data/ChR2UnitResponse.rds")


UnitDiversity <- PVCCFPlotZoom + PVCCFPlotOut + PVThetaPhase + plot_layout(nrow = 3) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = UnitDiversity, "/alzheimer/Daniel_Data/R/Thesis/Data/UnitDiversity.rds")


PVCCFPlotZoom + PVCCFPlotOut + PVThetaPhase + plot_layout(widths = c(1,2,2)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))

saveRDS(SpikeTableTotal[ClusterNr>1&grepl("MS|PaS|MEC", UnitLoc),.N,by=.(AnimalID, UniqueID, UnitLoc, Session)], "/alzheimer/Daniel_Data/R/Thesis/Data/UnitCountSummary.rds")

######
plot(PVCCF$CCF)
####
####
#### check z score on PaS MEC side
PaSMECTestunits <- SpikeTableTotal[ClusterNr>1&grepl(pattern = "PaS|MEC", UnitLoc)&!is.na(SingleStimDiff),.N,by=UniqueID][N>1000,UniqueID,] 

SpikeTableTotal[UniqueID%in%PaSMECTestunits,]

SpikeTableMS1Hz[pulse_nr==3305,]

SpikeTablePaSMEC1Hz <- SpikeTableTotal[UniqueID %in% PaSMECTestunits & !is.na(SingleStimDiff),] 
SpikeTablePaSMEC1Hz <- SpikeTablePaSMEC1Hz[!is.na(SingleStimDiff),.(UniqueID, UnitLoc, SingleStimDiff, hyper_block, PulseSingleNumber, AnimalID, RecDate),] 
SpikeTablePaSMEC1Hz[,StimRoundPos := round(SingleStimDiff*2, digits = 3)/2,] 
SpikeTablePaSMEC1Hz <- SpikeTablePaSMEC1Hz[,.(Count=.N),by=.(UniqueID, StimRoundPos)]
SpikeZScoreTable <- rbindlist(l = list(data.table(UniqueID=rep(unique(SpikeTablePaSMEC1Hz$UniqueID), 
                                                               each = length(seq(-0.5,0.5,0.0005))),
                                                  StimRoundPos = rep(round(seq(-0.5,0.5,0.0005), digits = 4), times=length(unique(SpikeTablePaSMEC1Hz$UniqueID))), 
                                                  Count=0),
                                       SpikeTablePaSMEC1Hz))
SpikeZScoreTable <- SpikeZScoreTable[,.(N=sum(Count)),by=.(UniqueID, StimRoundPos)] 
SpikeZScoreTable[,ScaleN:=scale(N), by=UniqueID] 
setorder(SpikeZScoreTable, "StimRoundPos")

OrderSpikeScale <- SpikeZScoreTable[StimRoundPos %between% c(0,0.03), .(ScaleCenter = max(ScaleN), UniqueID),by=UniqueID]
setorder(OrderSpikeScale, ScaleCenter)
OrderSpikeScale[,OrderID:=.I,] 
for(i in OrderSpikeScale[,UniqueID,]){
  SpikeZScoreTable[UniqueID==i, OrderID:=OrderSpikeScale[UniqueID==i,OrderID,],] 
} 

### identify PV cells ####
ScaleCutOff <- OrderSpikeScale[ScaleCenter< -8,min(OrderID, na.rm = T),] 
PV_Units <- OrderSpikeScale[ScaleCenter< -8,UniqueID,]
PV_negUnits <- OrderSpikeScale[ScaleCenter> -8,UniqueID,]
#SpikeTableTotal[UniqueID%in%PV_Units,PVPos:=TRUE,][UniqueID%in%PV_negUnits,PVPos:=FALSE,] 
### PV-Units per session
#SpikeTableTotal[UniqueID %in% PV_Units, unique(UniqueID), by= Session][,.N,by=Session] 


#### make plot for z-score ####
SinglePulseZScorePlot <- ggplot(data = SpikeZScoreTable[StimRoundPos %between% c(-0.035,0.035)]  , aes(x=StimRoundPos*1e3, y=OrderID, fill=ScaleN))+
  geom_raster()+
  geom_vline(xintercept = 0, lty="dashed", size=0.2)+
  scale_y_continuous(expand=c(0,0), name = "Units", seq(0,120,10))+
  scale_x_continuous(expand=c(0,0), name = "Time to Stimulation (ms)", breaks = seq(-40,40,10))+
  geom_hline(yintercept = ScaleCutOff-0.5, lty="dashed", size=0.2, colour="gray")+
  #scale_fill_viridis_c()+
  scale_fill_gradient2(limits = c(-max(SpikeZScoreTable[,max(ScaleN),]),SpikeZScoreTable[,max(ScaleN),]), guide = guide_colorbar(title = "Z-Score"))+
  theme_bw()

####
#### test where in wave pi is
Test_hb_transform <- hht::HilbertTransform(sig = ts(data = sin(seq(0,2*pi,0.01)), start = 0, frequency = 1e2))
Test_hb_power <- hht::HilbertEnvelope(Test_hb_transform)
TestData_phase <- (atan2(Re(Test_hb_transform),Im(Test_hb_transform))*180)/pi
plot(TestData_phase)
lines(x = seq(0,2*pi,0.01), y=sin(seq(0,2*pi,0.01))*100)
### 90° peak -90 trough

sin(-180*pi/180+pi)
cos(-180*pi/180+pi)
sin(((-90)+180)/(180)*pi)

plot(seq(-180,180)*pi/180)
sinePhase <- seq(-180,180)*pi/180+pi
plot(sinePhase)
plot(sin(sinePhase)^2)
plot(cos(sinePhase)^2)
plot(sin(sinePhase), cos(sinePhase))

SpikeTableTotal[,`:=`(PaSThetaX=sin(PaSThetaPhase*pi/180+pi), PaSThetaY=cos(PaSThetaPhase*pi/180+pi)),]
PhasePlot <- SpikeTableTotal[(PaSThetaPower/PaSDeltaPower)>3&UniqueID%in%MSTestunits&grepl(pattern = "MEC|PaS", PaSLoc)&stimulation_block==0,.(X=mean(PaSThetaX, na.rm = T), Y=mean(PaSThetaY, na.rm = T) , PV=unique(PVPos)), by=.(UniqueID)]
PhasePlot[PV==TRUE, PVStain:="pos",][PV==F, PVStain:="neg",]
PhaseCirclePlot <- 
  ggplot(data = PhasePlot, mapping = aes(x=X, y=Y, colour=PVStain))+
 # ggforce::geom_circle(aes(x0 = x, y0 = y, r =r, fill = NULL, linetype="dashed"),data=data.frame(x=c(0,0),y=c(0,0),r=c(0.5,1), size=c(0.001,0.001)), inherit.aes = FALSE)+
 # scale_y_continuous(limits = c(-0.8,0.8))+
  #scale_x_continuous(limits = c(-0.8,0.8))+
  annotate("path", x=0+ 1*cos(seq(0,2*pi,length.out=100)), y=0+1*sin(seq(0,2*pi,length.out=100)), size=0.2)+
  annotate("path", x=0+ .5*cos(seq(0,2*pi,length.out=100)), y=0+.5*sin(seq(0,2*pi,length.out=100)), size=0.2)+
  annotate(geom = "segment", x = c(-1,0), xend = c(1,0), y = c(0,-1), yend = c(0,1), size=0.1)+
  # geom_hline(yintercept = 0, size=0.2)+
  # geom_vline(xintercept = 0, size=0.2)+
  geom_point(data = PhasePlot[PV==FALSE,])+
  geom_point(data = PhasePlot[PV==TRUE,])+
#  geom_line(data = data.frame(X=seq(-1,1,length.out = 400), Y=sin(seq(0,4*pi, length.out = 400))/2-2.5), aes(x=X,y=Y), inherit.aes = F)+
    annotate(geom = "text", x = c(3.5), y = c(-3.1)+2.5, label=c("90°"), size=3.5, vjust="outward")+
  annotate(geom = "text", x = c(3.9), y = c(-2.5)+2.5, label=c("180°"), size=3.5, vjust="outward", hjust="inward")+
  annotate(geom = "text", x = c(4.5), y = c(-1.9)+2.5, label=c("270°"), size=3.5, vjust="outward")+
  annotate(geom = "text", x = c(5.1), y = c(-2.5)+2.5, label=c("0°"), size=3.5, vjust="outward", hjust="outward")+
  geom_line(data = data.frame(X=seq(-1,3,length.out = 400)+3, Y=sin(seq(0,4*pi, length.out = 400))/2), aes(x=X,y=Y), inherit.aes = F)+
    
    scale_colour_manual("", values = c("pos" = "red", "neg" = "black"), labels=c(bquote(PV^"-"),bquote(PV^"+")))+
  annotate(geom = "text", x = c(.4,0.75), y = c(-.4,-0.75), label=c("0.5","1"), hjust="inward", vjust="outward")+
  annotate(geom = "text", x = c(1.05), y = c(0), label=c("90°"), hjust="inward", size=5)+
    annotate(geom = "text", x = c(-1.05), y = c(0), label=c("270°"), hjust="outward", size=5)+
    annotate(geom = "text", x = c(0), y = c(-1.1), label=c("180°"), vjust="outward", size=5)+
    annotate(geom = "text", x = c(0), y = c(1.1), label=c("0°"), vjust="outward", size=5)+
    #annotate(geom = "text", x = c(1.05,-1.05,0,0), y = c(0,0,-1.1,1.1), label=c("90°", "270°", "180°", "0°"), hjust="outward", size=5)+
    #  annotate(geom = "text", x = c(1.05,-1.05,0,0), y = c(0,0,-1.1,1.1), label=c("90°", "270°", "180°", "0°"), hjust="outward", size=5)+
#  annotate(geom = "text", x = c(-.25, -0.17,.25,0.6), y = c(-3.1,-2.5,-1.9,-2.5), label=c("90°", "180°","270°", "0°"), size=3.5)+
  #coord_equal(xlim = c(-1.5,1.5), ylim = c(-3,1.1))+
    coord_equal(xlim = c(-1.6,7), ylim = c(-1.2,1.2))+
    #coord_flip()+
  #theme_minimal()
  theme_void()+theme(legend.position = "left", legend.text=element_text(size=10))

PhaseDesignMatrix <- "
A
B
C
D
D
"
UnitDiversity <- PVCCFPlotZoom + PVCCFPlotOut + PVThetaPhase + PhaseCirclePlot + plot_layout(design = PhaseDesignMatrix) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = UnitDiversity, "/alzheimer/Daniel_Data/R/Thesis/Data/UnitDiversity.rds")

PVThetaPhase / PhaseCirclePlot + plot_layout(heights = c(1,2)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))


sin(c(-180,90,0,-90)*pi/180+pi)
cos(c(-180,90,0,-90)*pi/180+pi)

circular::mean.circular(SpikeTableTotal[UniqueID==PV_Units[3]&stimulation_block==0,PaSThetaPhase,])

TestPhasePlotPV <- SpikeTableTotal[(PaSThetaPower/PaSDeltaPower)>3&ClusterNr>1&UnitLoc=="MS"& UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),.(X=mean(PaSThetaX), Y=mean(PaSThetaY)),by=UniqueID]
TestPhasePlotPV[,Rho:=sqrt(Y^2+X^2)]
TestPhasePlot <- SpikeTableTotal[(PaSThetaPower/PaSDeltaPower)>3&ClusterNr>1&UnitLoc=="MS"& !UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),.(X=mean(PaSThetaX), Y=mean(PaSThetaY)),by=UniqueID]
TestPhasePlot[,Rho:=sqrt(Y^2+X^2)]
boxplot(TestPhasePlot$Rho, TestPhasePlotPV$Rho)
#### vm model

VMData <- SpikeTableTotal[!is.na(PVPos)&(PaSThetaPower/PaSDeltaPower)>3&ClusterNr>1&UnitLoc=="MS"&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),][
  ,.(PVPos=unique(PVPos), Phase=circular::mean.circular(PaSThetaPhase/180*pi)),by=UniqueID]
#VMData <- rbindlist(l = list(VMData[PVPos==T,][1:1e3,], VMData[PVPos==F,][1:1e3,]))
stanDataVonMises <- list(
  N=VMData[,.N,],
  y=VMData[,Phase+pi,],
  PVpos=VMData[,PVPos,],
#  UnitN=max(as.integer(factor(VMData[,UniqueID,]))),
#  Unit=as.integer(factor(VMData[,UniqueID,])),
  priorOnly=0
)
rstudioapi::jobRunScript(path = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/R/RunModelVM.R", workingDir = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/", importEnv = T, exportEnv = "R_GlobalEnv")


y_rep <- posterior::as_draws_matrix(fitmodVM$draws("y_rep"))

bayesplot::ppc_dens_overlay_grouped(y = as.vector(stanDataVonMises$y), group = stanDataVonMises$PVpos ,yrep = y_rep[1:200,]%%2*pi)
bayesplot::ppc_intervals(y = as.vector(VMData[,circular::mean.circular(x = PaSThetaPhase),by=UniqueID][,V1,]), yrep = y_rep)
fitmodVM$print()

plot(VMData$Phase[VMData$PVPos==0],col="red", add=T)

library(rstan)
library(circular)

N = 100
y = as.vector(rvonmises(N, -pi, 1.0))
##### write phase lock stan file

plot(TestPhasePlotPV[,.(X,Y),], ylim=c(-0.5,0.5), xlim=c(-0.5,0.5), col="red")
points(TestPhasePlot[,.(X,Y),])
abline(h = 0, v = 0)
ggplot(data = SpikeTableTotal[UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),], aes(x=PaSThetX, group=UniqueID))+
  geom_histogram()+
  facet_wrap(PaSLoc~UniqueID)#+
  geom_vline(xintercept = c(-90,90))

###### theta lock

ggplot(data = SpikeTableTotal[UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),], aes(x=PaSThetaPhase, group=UniqueID))+
  geom_histogram()+
  facet_wrap(PaSLoc~UniqueID)+
  geom_vline(xintercept = c(-90,90))

ggplot(data = SpikeTableTotal[UniqueID%in%PV_Units&stimulation_block==0&grepl(pattern = "PaS|MEC", PaSLoc),], aes(x=MSThetaPhase, group=UniqueID))+
  geom_histogram()+
  facet_wrap(PaSLoc~UniqueID)
SpikeTableTotal[PaSLoc=="PaS"]

#### make ccf graph ####


####
SessionCCFs <- lapply(X = SpikeTableTotal[UniqueID %in% PV_Units,unique(Session),], FUN = function(x){
  InVivoR::SpikeCCF(x = SpikeTableTotal[UnitLoc=="MS"&ClusterNr>1&stimulation_block==0&UniqueID %in% PV_Units & Session==x,SpikeTime,], UnitNr = factor(SpikeTableTotal[UnitLoc=="MS"&ClusterNr>1&stimulation_block==0&UniqueID %in% PV_Units& Session==x,factor(UniqueID),]))
})

SessionCCFs[[2]]

testCCF <- InVivoR::SpikeCCF(x = SpikeTableTotal[stimulation_block==0&UniqueID%in%PV_Units[13] ,SpikeTime,])

plot(testCCF$CCF)

SpikeTableTotal$stimulation_block


CCFmat <- spike_ccf_batch(Time = SpikeData[ClusterNr>1,SpikeTime/2e4,], UnitNr = SpikeData[ClusterNr>1,ClusterNr])
CCFdt <- data.table::data.table(X = rep(CCFmat$xAxis,
                                        times = dim(CCFmat$CcfMatrix)[2]),
                                Units1 = rep(x = CCFmat$Units[,1], each = dim(CCFmat$CcfMatrix)[1]),
                                Units2 = rep(x = CCFmat$Units[,2], each = dim(CCFmat$CcfMatrix)[1]),
                                Count = as.vector(CCFmat$CcfMatrix),
                                CIlow = as.vector(CCFmat$LowerCI),
                                CIup = as.vector(CCFmat$UpperCI))[
                                  ,Group:=paste(Units1,Units2, sep = " - "),]

ggplot(data = CCFdt, aes(x=X, y=Count, group=Group))+
  geom_ribbon(aes(ymax = CIup, ymin = CIlow), alpha=0.1)+
  scale_y_continuous(limits = c(0,NA), expand = c(0,0))+
  geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(facets = Units2~Units1, scales = "free_y", )+
  theme_classic()


CleanUnitsCCF <- SpikeCCFTable(Time = SpikeTableTotal[UniqueID==PV_Units[13]&is.na(onset),SpikeTime,],
                               UnitName = SpikeTableTotal[UniqueID==PV_Units[13]&is.na(onset),UniqueID,])

CCFPlot <- ggplot(data = CleanUnitsCCF[x %between% c(-0.5,0.5), ], aes(x=x, y=NormCount, group=Group))+
  #facet_grid(facets = ~Group, scales = "free_y")+
  scale_y_continuous(limits = c(0,NA), expand = c(0,0), name = "")+
  scale_x_continuous(expand = c(0,0), name = expression(paste(Delta, "Time (s)")))+
  #labs(title = "Cross-Correlograms", subtitle = "Normalised Counts")+
  geom_ribbon(aes(ymax=NormCount, ymin=0), fill="gray50")+
  #geom_line(aes(y=NormCIlow), colour="red", linetype="dashed")+
  #geom_line(aes(y=NormCIup), colour="red", linetype="dashed")+
  #geom_hline(yintercept = unique(CleanUnitsCCF$NormRandom))+
  theme_classic()+
  theme(strip.background = element_rect(color=NaN, fill=NaN))


SinglePulseOverView <- SinglePulseZScorePlot + RasterPlot / (CCFPlot)+ 
  plot_annotation(title = "PV-Unit Identification", subtitle = "ChR2 actviation in MS using 1Hz Pulse", tag_levels = "A") &
  theme(plot.title = element_text(size=20),plot.tag = element_text(size=24)) & theme(plot.tag = element_text(size=24)) 


ggsave(filename = "/alzheimer/Daniel_Data/Analysis/Analysis_output/20200511/SinglePulseOverView_MS.pdf",
       plot = SinglePulseOverView,
       device = "pdf",
       width = 12, 
       height = 7)

PVUnitMS <- readRDS(file = "/alzheimer/Daniel_Data/DSC008154/DSC008154_190611_143146/UnitList")



curve((sin(x-sin(x+10)/2)>0.3)* (sin(x-sin(x)/2)) * sin(x*20)/6, from = 0, to= 50, n = 1e4)
sin(x-sin(x)/2)

xSin <- seq(0,50,0.01)
ySin <- sin(xSin-sin(xSin)/2)+(sin(xSin-sin(xSin+10)/2)>0.4)* sin(xSin) * sin(xSin*20)/6 + sin(xSin/2-1)/3 + sin(xSin/1.3-1.5)/4 + sin(xSin/1.2-1.2)/8*sin(xSin*6)
ySinT <- sin(xSin-sin(xSin)/2)
ySinG <- (sin(xSin-sin(xSin+10)/2)>0.4)* sin(xSin) * sin(xSin*20)/9

plot(ySin, type="l")

WaveData <- data.table(x=rep(xSin, times=3), y=c(ySin, ySinT, ySinG), Group= factor(x = rep(c("LFP", "theta", "gamma"), each=length(xSin)),levels = c("LFP", "theta", "gamma"), labels = c("LFP", "theta*' 4-12Hz'", "gamma*' 30-90Hz'")))
FrequencyBandScheme <- ggplot(data = WaveData, mapping = aes(x=x, y = y, group=Group))+
  geom_line()+
  facet_wrap(~Group, ncol = 1, labeller = label_parsed, strip.position="right")+
  theme_void()+theme(strip.text = element_text(size=15))
saveRDS(FrequencyBandScheme, "/alzheimer/Daniel_Data/R/Thesis/Data/FrequencyBandScheme.rds")
