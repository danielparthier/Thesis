np <- reticulate::import("numpy")
spike_clusters <- np$load("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/spike_clusters.npy")
spike_times <- np$load("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/spike_times.npy")
list.files("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/")
ClusterTable <- fread("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/cluster_groups.csv")
SpikeTable <- data.table(SpikeIdx=as.vector(spike_times), ClusterNr=as.vector(spike_clusters))
SpikeTable[ClusterNr%in%ClusterTable[group=="good",cluster_id],ClusterType:="good",]
SpikeTable[ClusterNr%in%ClusterTable[group=="mua",cluster_id],ClusterType:="mua",]
SpikeTable <- SpikeTable[!is.na(ClusterType),]


testCCF <- InVivoR::SpikeCCF(x = testTable[ClusterFac=="506_4",SpikeTime,], y = testTable[ClusterFac=="506_4",SpikeTime,])
plot(y=testCCF$CCF, testCCF$xAxis)



#######################
np <- reticulate::import("numpy")

Overview_Channel_choice <- data.table(read_excel("/alzheimer/Daniel_Data/Analysis/Overview_Channel_choice.xlsx"))
ChannelLocations <- Overview_Channel_choice[!is.na(session),]
ChannelLocations <- ChannelLocations[,.(...2, MS_Channel, LS_Channel, location_other_MS, PaS_Channel, Other_PaS_Channel, location_other_PaS, ...14),]
ChannelLocations[!is.na(MS_Channel),MSLocation:="MS",][is.na(MS_Channel)&!is.na(LS_Channel),MSLocation:="LS",][is.na(MS_Channel)&is.na(LS_Channel),MSLocation:=...14,]
ChannelLocations[!is.na(PaS_Channel),PaSLocation:="PaS",][!is.na(Other_PaS_Channel),OtherPaSLocation:=location_other_PaS,]
ChannelLocations[grepl(pattern = "MEC", x = PaSLocation), PaSLocation:="MEC"]
ChannelLocations[,.(MS_Channel, LS_Channel, PaS_Channel),]
ChannelLocations[,`:=`(PaS_Channel=PaS_Channel+32, Other_PaS_Channel=Other_PaS_Channel+32),]
#ChannelLocations <- ChannelLocations[-c(1:7),]
#ChannelLocations <- ChannelLocations[-c(1:11),]


FullList <- list.files(path = Overview_Channel_choice$...2, pattern = "spike_times.npy|clu", recursive = T, full.names = T)
Overview_Channel_choice <- Overview_Channel_choice[!is.na(...2),]
SortedData <- sapply(Overview_Channel_choice$...2, function(x) {
  any(grepl(pattern = x, FullList))
})
SortedData <- Overview_Channel_choice$...2[SortedData]

###### Spike####
ExtractSpikes <- function(RefLoc, LFPRefChannels, RecordingFolder, DataPath){
  # browser()
  # for(RecordingFolder in SortedData) {
  
  OutputFolder <- paste0(RecordingFolder, "OutputFolder/")
  SessionName <- tail(unlist(strsplit(RecordingFolder, split = "/")), n = 1)
  if(!dir.exists(OutputFolder)) {
    dir.create(path = OutputFolder)
  }
  
  NDManager <- ifelse(test = length(list.files(path = RecordingFolder,
                                               pattern = paste0(tail(unlist(strsplit(RecordingFolder, split = "/")), n=1), ".xml"),
  ))==1,
  yes = T, no = F)
  
  ProbeLocationFile <- list.files(path = RecordingFolder, pattern = ".xlsx")
  
  #### load location files and adjust channel number so that they are continous
  ChannelMax <- 0
  ShankMax <- 0
  for(i in ProbeOrder){
    if(ChannelMax==0){
      LocationTable <- as.data.table(read_excel(paste0(RecordingFolder, ProbeLocationFile[grep(x = ProbeLocationFile, pattern = i)])))
      ChannelMax <- max(LocationTable$ChNr)
      ShankMax <- max(LocationTable$Shank)
    } else {
      tmp <- as.data.table(read_excel(paste0(RecordingFolder, ProbeLocationFile[grep(x = ProbeLocationFile, pattern = i)])))
      tmp[,`:=`(ChNr=ChNr+ChannelMax, Shank=Shank+ShankMax),]
      LocationTable <- rbindlist(l = list(LocationTable, tmp))
      ChannelMax <- max(LocationTable$ChNr)
      ShankMax <- max(LocationTable$Shank)
    }  
  }
  LocationTable[,loc:=as.factor(loc)] 
  setorder(x = LocationTable, ChNr)
  message("loaded location files")
  
  
  if(NDManager){
    AmpFilePath <- paste0(RecordingFolder, tail(unlist(strsplit(x = RecordingFolder, split = "/")), n = 1), ".dat")
    SpikeDirectory <- RecordingFolder
    RecInfo <- MetaData(FileName = paste0(RecordingFolder, tail(unlist(strsplit(x = RecordingFolder, split = "/")), n = 1), ".xml"))
    SamplingRate <- RecInfo$SamplingRate
    ChannelNr <- RecInfo$ChannelNr
    
    
    #### load spike files NDmanager/Klusers ####
    CluFiles <- list.files(path = SpikeDirectory, pattern = ".clu", full.names = T)
    ResFiles <- list.files(path = SpikeDirectory, pattern = ".res", full.names = T)
    FetFiles <- list.files(path = SpikeDirectory, pattern = ".fet", full.names = T)
    
    for(i in 1:length(CluFiles)) {
      if(length(count.fields(CluFiles[i]))>1) {
        if(i == 1) {
          SpikeTable <- fread(input = CluFiles[1], skip = 1, col.names = "ClusterNr")
          SpikeTable[,SpikeIdx := fread(input = ResFiles[1], skip = 0)][,SpikeGroup:=1,]
        } else {
          tmpTable <- fread(input = CluFiles[i], skip = 1, col.names = "ClusterNr")
          tmpTable[,SpikeIdx := fread(input = ResFiles[i], skip = 0)][,SpikeGroup:=i,]
          SpikeTable <- data.table::rbindlist(l = list(SpikeTable, tmpTable))
          rm(tmpTable)
        }
      }
    }
  } else{
    SamplingRate <- 2e4
    ChannelNr <- 64
    np <- reticulate::import("numpy")
    for(i in seq_along(ProbeOrder)){
      ClusterList <- list.files(path = paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/"),
                                pattern = "cluster_group",
                                full.names = T)
      if(length(ClusterList)!=1) {
        ClusterList <- tcltk::tk_choose.dir(default = paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/"), caption = "Missing Clusterfile")
        ClusterTable <- fread(ClusterList)
      } else {
        ClusterTable <- fread(ClusterList)
      }
      ClusterTable$cluster_id <- ClusterTable$cluster_id+1
      
      if(i==1){
        SpikeTable<- data.table(SpikeIdx=as.vector(np$load(paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/spike_times.npy"))),
                                ClusterNr=as.vector(np$load(paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/spike_clusters.npy")))+1,
                                SpikeGroup=i)
        
        SpikeTable[ClusterNr%in%ClusterTable[group=="noise",cluster_id],ClusterNr:=0,]
        SpikeTable[ClusterNr%in%ClusterTable[group=="mua",cluster_id],ClusterNr:=1,]
        
      } else {
        SpikeTableTmp <- data.table(SpikeIdx=as.vector(np$load(paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/spike_times.npy"))),
                                    ClusterNr=as.vector(np$load(paste0(RecordingFolder,paste0("amplifier_", ProbeOrder[i]),"/spike_clusters.npy")))+1,
                                    SpikeGroup=i)
        
        SpikeTableTmp[ClusterNr%in%ClusterTable[group=="noise",cluster_id],ClusterNr:=0,]
        SpikeTableTmp[ClusterNr%in%ClusterTable[group=="mua",cluster_id],ClusterNr:=1,]
        SpikeTable <- rbindlist(l = list(SpikeTable, SpikeTableTmp))
      }  
    }
    SpikeTable[,ClusterNr:=as.integer(ClusterNr),] 
  } 

  SpikeTable[,SpikeTime := SpikeIdx/SamplingRate,][,ClusterFac:=as.factor(paste(ClusterNr, SpikeGroup, sep = "_")),]
  SpikeTable[ClusterNr>1, .N, by=ClusterFac]
  #### remove noise and clean labels ####
  SpikeTable[ClusterNr==1, ClusterFac:=as.factor(paste0("MUA_",SpikeGroup))]
  UnitFilter <- SpikeTable[ClusterNr>0,1/mean(diff(SpikeTime)),by=ClusterFac][V1>1,ClusterFac,] 
  SpikeTable <- SpikeTable[(ClusterNr==1|ClusterFac %in% UnitFilter), ]
  
  for(UnitSet in SpikeTable[ClusterNr!=1, unique(ClusterFac)]){
    SpikeTable[ClusterFac==UnitSet,SpikeNr:=.I,]
  } 
  
  ColumsnToCarry <- names(SpikeTable)
  StimMatDir <- list.files(path = OutputFolder, pattern = "StimMat", recursive = T, full.names = T)
  if(length(StimMatDir)>1) {
    ### load stim mat files  
    StimMatMS <- readRDS(file = grep(pattern = "StimMatMS.rds",x = StimMatDir, value = T))
    StimMatPaS <- readRDS(file = grep(pattern = "StimMatPaS.rds",x = StimMatDir, value = T))
    
    ### Add MS stim
    SpikeTable <- data.table(SpikeTable ,SpikeStimProperties(SpikeIdx = SpikeTable$SpikeIdx, StimMat = StimMatMS$PulseMat, BlockMat = StimMatMS$BlockMat, SamplingRate = 2e4, Isolated = F))
    
    for(PulseTiming in StimMatMS$PulseMat[StimMatMS$PulseMat[,16]==1,1]) {
      PulseSeconds <- PulseTiming/2e4
      SpikeTable[SpikeTime %between%c(PulseSeconds-0.5, PulseSeconds+0.5),SingleStimDiff:=SpikeTime-PulseSeconds,]
    }
    SpikeTableTmp1 <- SpikeTable[hyper_block_frequency!=0,]
    SpikeTableTmp1[,StimLocation:="MS",]
    
    ### Add PaS stim
    SpikeTableTmp2 <- SpikeTable[hyper_block_frequency==0,..ColumsnToCarry]
    SpikeTableTmp2 <- data.table(SpikeTableTmp2, SpikeStimProperties(SpikeIdx = SpikeTableTmp2$SpikeIdx, StimMat = StimMatPaS$PulseMat, BlockMat = StimMatPaS$BlockMat, SamplingRate = 2e4, Isolated = F))
    SpikeTableTmp2[,SpikeTime:=SpikeIdx/2e4,]
    SpikeTableTmp2[,SingleStimDiff:=NA,]
    for(PulseTiming in StimMatPaS$PulseMat[StimMatPaS$PulseMat[,16]==1,1]) {
      PulseSeconds <- PulseTiming/2e4
      SpikeTableTmp2[SpikeTime %between%c(PulseSeconds-0.5, PulseSeconds+0.5),SingleStimDiff:=SpikeTime-PulseSeconds,]
    }
    SpikeTableTmp2[,StimLocation:="PaS",]
    SpikeTable <- rbindlist(l = list(SpikeTableTmp1, SpikeTableTmp2))
    SpikeTable[,time_stamp:=NULL,]
  } else if(length(StimMatDir)==1) {
    StimMat <- readRDS(file = StimMatDir)
    ### Add MS stim
    SpikeTable <- data.table(SpikeTable ,SpikeStimProperties(SpikeIdx = SpikeTable$SpikeIdx, StimMat = StimMatMS$PulseMat, BlockMat = StimMatMS$BlockMat, SamplingRate = 2e4, Isolated = F))
    
    for(PulseTiming in StimMat$PulseMat[StimMat$PulseMat[,16]==1,1]) {
      PulseSeconds <- PulseTiming/2e4
      SpikeTable[SpikeTime %between%c(PulseSeconds-0.5, PulseSeconds+0.5),SingleStimDiff:=SpikeTime-PulseSeconds,]
    }
    locString <- ifelse(grepl(pattern = "MS", x = StimMatDir), "MS", "PaS")
    SpikeTable[,StimLocation:=locString,]
    SpikeTable[,time_stamp:=NULL,]
  }
  ##### load Amp file
  
  if(NDManager){
    AmpFile <- AmpFileRead(FILENAME = AmpFilePath, ChannelNumber = ChannelNr)
  } else {
    AmpFile <- rbind(AmpFileRead(FILENAME = paste(RecordingFolder,"amplifier_MS/amplifier_MS.dat", sep = "/"), ChannelNumber = 32),
                     AmpFileRead(FILENAME = paste(RecordingFolder,"amplifier_PaS/amplifier_PaS.dat", sep = "/"), ChannelNumber = 32))
  }  
  
  
  message("Extract Reference Channels")
  LFPRefSignal <- AmpFile[LFPRefChannels,] 
  LFPRefSignalTheta <- AmpFile[LFPRefChannels,] 
  LFPRefSignalDelta <- AmpFile[LFPRefChannels,] 
  
  LFPDownSample <- matrix(data = 0, nrow = dim(LFPRefSignal)[1], ncol = floor(dim(LFPRefSignal)[2]/20)) 
  filter_width <- 0.5e3/(SamplingRate/2)
  fir_filter501 <- signal::fir1(501, filter_width, type = "low")
  
  filter_width_Theta <- c(4/500,12/500)
  fir_filter_Theta <- signal::fir1(4001, filter_width_Theta, "DC-0")
  filter_width_Delta <- c(0.5/500,4/500)
  fir_filter_Delta <- signal::fir1(4001, filter_width_Delta, "DC-0")
  
  
  for(RefNr in 1:dim(LFPRefSignal)[1]){
    LFPDownSample[RefNr,] <- InVivoR::decimate(SIGNAL = LFPRefSignal[RefNr,], FIR_FILTER = fir_filter501, M = 20, CORES = Cores)
    message("Reference Channel ", RefNr)
  }  
  
  
  message("importing amplifier file complete")
  UnitLocation <- UnitChannel(SpikeIdx = SpikeTable[ClusterNr!=1 & (SpikeIdx<dim(AmpFile)[2]-21), SpikeIdx],
                              Units = SpikeTable[ClusterNr!=1 & (SpikeIdx<dim(AmpFile)[2]-21), ClusterFac],
                              AmpMatrix = AmpFile,
                              WINDOW = 20)
  saveRDS(UnitLocation, paste0(OutputFolder, "UnitLocation.rds"))
  
  if(!is.null(LocationTable)) {
    UnitLocation$Location <- LocationTable[UnitLocation$ChannelNr,loc]
  }
  
  #saveRDS(object = UnitLocation, file = paste0(RecordingFolder, "UnitList"))
  
  #### assign location according to excel location and channel ####
  for(i in seq_along(UnitLocation$UnitNr)){
    if(is.null(LocationTable)) {
      SpikeTable[as.integer(ClusterFac)==UnitLocation$UnitNr[i],
                 `:=`(UnitChan=UnitLocation$ChannelNr[i],
                      UnitAmp=UnitLocation$Amplitude[i])] 
    } else {
      SpikeTable[as.integer(ClusterFac)==UnitLocation$UnitNr[i],
                 `:=`(UnitChan=UnitLocation$ChannelNr[i],
                      UnitAmp=UnitLocation$Amplitude[i],
                      UnitLoc=LocationTable[ChNr==UnitLocation$ChannelNr[i],loc])]
    }
    
  } 
  # }
  
  SpikeTableSessionIndex <- SpikeTable[,floor(SpikeIdx/20),]
  
  for(i in seq_len(dim(LFPDownSample)[1])) {
    ThetaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Theta, FiltFilt = T, CORES = 10)
    DeltaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Delta, FiltFilt = T, CORES = 10)
    
    ThetaData_hb_transform <- hht::HilbertTransform(sig = ts(data = ThetaData, start = 0, frequency = 1e3))
    Theta_hb_power <- hht::HilbertEnvelope(ThetaData_hb_transform)
    ThetaData_phase <- (atan2(Re(ThetaData_hb_transform),Im(ThetaData_hb_transform))*180)/pi
    
    DeltaData_hb_transform <- hht::HilbertTransform(sig = ts(data = DeltaData, start = 0, frequency = 1e3))
    Delta_hb_power <- hht::HilbertEnvelope(DeltaData_hb_transform)
    DeltaData_phase <- (atan2(Re(DeltaData_hb_transform),Im(DeltaData_hb_transform))*180)/pi
    
    if(i==1) {
      SpikeTable[,`:=`(MSThetaPower=Theta_hb_power[SpikeTableSessionIndex], MSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], MSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], MSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], MSLoc=RefLoc[i]),]
    } else {
      SpikeTable[,`:=`(PaSThetaPower=Theta_hb_power[SpikeTableSessionIndex], PaSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], PaSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], PaSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], PaSLoc=RefLoc[i]),]
    }
  }
  rm(AmpFile)
  SpikeTable[,Session:=SessionName,]
  saveRDS(object = SpikeTable, file = paste0(OutputFolder, "SpikeTableThetaDelta.rds"))
  saveRDS(object = LFPDownSample, file = paste0(OutputFolder, "LFPDownSample.rds"))
  rm(LFPDownSample, SpikeTable)
  gc()
}

library(data.table)
library(readxl)
library(ggplot2)
library(InVivoR)

ProbeOrder <- c("MS", "PaS")

##### Theta phase to Table
Overview_Channel_choice <- data.table(read_excel("/alzheimer/Daniel_Data/Analysis/Overview_Channel_choice.xlsx"))
ChannelLocations <- Overview_Channel_choice[!is.na(session),]
ChannelLocations <- ChannelLocations[,.(...2, MS_Channel, LS_Channel, location_other_MS, PaS_Channel, Other_PaS_Channel, location_other_PaS, ...14),]
ChannelLocations[!is.na(MS_Channel),MSLocation:="MS",][is.na(MS_Channel)&!is.na(LS_Channel),MSLocation:="LS",][is.na(MS_Channel)&is.na(LS_Channel),MSLocation:=...14,]
ChannelLocations[!is.na(PaS_Channel),PaSLocation:="PaS",][!is.na(Other_PaS_Channel),OtherPaSLocation:=location_other_PaS,]
ChannelLocations[grepl(pattern = "MEC", x = PaSLocation), PaSLocation:="MEC"]
ChannelLocations[,.(MS_Channel, LS_Channel, PaS_Channel),]
ChannelLocations[,`:=`(PaS_Channel=PaS_Channel+32, Other_PaS_Channel=Other_PaS_Channel+32),]
#ChannelLocations <- ChannelLocations[-c(1:7),]
#ChannelLocations <- ChannelLocations[-c(1:11),]

#ChannelLocationsFilter <- ChannelLocations[grepl(pattern = paste(unique(SpikeTableTotal[,SessionID,]), collapse = "|"), x = ...2)]
FullList <- list.files(path = Overview_Channel_choice$...2, pattern = "spike_times.npy|clu", recursive = T, full.names = T)
Cores <- 20
Overview_Channel_choice <- ChannelLocations[!is.na(...2),]
SortedData <- sapply(ChannelLocations$...2, function(x) {
  any(grepl(pattern = x, FullList))
})
#SortedData <- Overview_Channel_choice$...2[SortedData]

ChannelLocationsFilter <- ChannelLocations[SortedData,]

DataPaths <- lapply(X = ChannelLocationsFilter$...2, function(x) {
  amplifierFile <- (paste0(gsub(pattern = "/",x = gsub(pattern = dirname(x), replacement = "", x = x), replacement = ""), ".dat"))
  amplifierFile <- list.files(path = x, pattern = amplifierFile, full.names = T)
  if(length(amplifierFile)==0) {
    amplifierFile <- list.files(path = x, pattern = "amplifier.*.dat", recursive = T, full.names = T)
    if(length(amplifierFile)>2) {
      amplifierFile <- grep(pattern = "amplifier.dat", x = amplifierFile, value = T)
    }
  }
  
  digitalinPath <- list.files(path = x, pattern = "digitalin.dat", full.names = T)
  analoginPath <- list.files(path = x, pattern = "analogin.dat", full.names = T)
  
  list(AmplifierFile=amplifierFile,DigitalIn=digitalinPath, AnalogIn=analoginPath)
})


CompletedPaths <- dirname(dirname(list.files(ChannelLocationsFilter$...2, pattern = "SpikeTableThetaDelta.rds", recursive = T, full.names = T)))
ChannelLocationsFilterMissing <- ChannelLocationsFilter[!ChannelLocationsFilter$...2%in%paste0(CompletedPaths,"/"),]


for(FolderNr in 18:19){#seq_along(ChannelLocationsFilter$...2)) {
  RefLoc <- as.vector(na.omit(unlist(ChannelLocationsFilter[FolderNr,9:11])))
  LFPRefChannels <- as.vector(na.omit(unlist(ChannelLocationsFilter[FolderNr,2:6])))
  print(data.frame(RefLoc, LFPRefChannels))
  RecordingFolder <- ChannelLocationsFilter$...2[FolderNr]
  ExtractSpikes(RefLoc = RefLoc[1:2], LFPRefChannels = LFPRefChannels[1:2], RecordingFolder = RecordingFolder, DataPath = DataPaths[[FolderNr]])
gc()
  }



ChannelLocationsFilter[21,]


testTable <- readRDS("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/OutputFolder/SpikeTableThetaDelta.rds")
testTable[,MeanFiring:=mean(1/diff(SpikeTime)),by=ClusterFac]
testTable <- testTable[MeanFiring>1,]
testTable[,ThetaDeltaRatio:=PaSThetaPower/PaSDeltaPower,]
for(i in testTable[,unique(ClusterFac),]) {
  tmpCCF <- InVivoR::SpikeCCF(testTable[ClusterFac==i,SpikeTime,])
  testTable[ClusterFac==i, Clean20:=max(tmpCCF$CCF[tmpCCF$xAxis<0.002&tmpCCF$xAxis>-0.002])<tmpCCF$RandomBinCount*0.1,]
}
TestUnits <- testTable[Clean20==T,unique(ClusterFac),]
CellsCCF <- InVivoR::SpikeCCF(testTable[ClusterFac%in%TestUnits,SpikeTime,], UnitNr = testTable[ClusterFac%in%TestUnits,ClusterFac,])
#tmpCCF <- InVivoR::SpikeCCF(testTable[ClusterFac=="466_2",SpikeTime,])
#tmpCCF$CCF[tmpCCF$xAxis==0]<tmpCCF$RandomBinCount*0.2

UnitSelection <- TestUnits
CCFTable <- rbindlist(lapply(X = seq_along(CellsCCF$Units[,1]), function(x) {
  data.table(Cluster1=CellsCCF$Units[x,1], Cluster2=CellsCCF$Units[x,2], Time=CellsCCF$xAxis, CCF=CellsCCF$CcfMatrix[,x], RandomBin=CellsCCF$RandomBinCount[x], Lower=CellsCCF$LowerCI[,x], Upper=CellsCCF$UpperCI[,x])
}))
for(i in seq_along(UnitSelection)) {
  CCFTable[Cluster1==as.integer(UnitSelection)[i],Cluster1Name:=UnitSelection[i],]
  CCFTable[Cluster2==as.integer(UnitSelection)[i],Cluster2Name:=UnitSelection[i],]
}
CCFTable[,Comparison:=paste(Cluster1Name,Cluster2Name, sep = "-"),]

ggplot(data = CCFTable[Time%between%c(-0.01,0.01)], aes(y=CCF, x=Time-unique(diff(CellsCCF$xAxis))/2, group=Comparison))+
  geom_step()+
  geom_vline(xintercept = 0)+
  facet_wrap(~Comparison, scales = "free_y")


ggplot(data = testTable[Clean20==T&ThetaDeltaRatio>2,,], mapping = aes(x = PaSThetaPhase, group=ClusterFac))+
  geom_histogram()+
  facet_wrap(~ClusterFac, scales = "free_y")


ggplot(data = CCFTable[Time %between% c(-0.5,0.5), ], aes(x=Time, y=CCF, group=Comparison))+
  #facet_grid(facets = ~Group, scales = "free_y")+
  scale_y_continuous(limits = c(0,NA), expand = c(0,0), name = "")+
  scale_x_continuous(expand = c(0,0), name = expression(paste(Delta, "Time (s)")))+
  #labs(title = "Cross-Correlograms", subtitle = "Normalised Counts")+
  geom_ribbon(aes(ymax=CCF, ymin=0), fill="gray50")+
  #geom_line(aes(y=NormCIlow), colour="red", linetype="dashed")+
  #geom_line(aes(y=NormCIup), colour="red", linetype="dashed")+
  #geom_hline(yintercept = unique(CleanUnitsCCF$NormRandom))+
  theme_classic()+
  theme(strip.background = element_rect(color=NaN, fill=NaN))

hist(testTable$ThetaDeltaRatio, breaks = seq(0,1e4, 1), xlim = c(0,10))

seq_along(CellsCCF$Units[,1])[CellsCCF$Units[,1]==CellsCCF$Units[,2]]
plot(x = CellsCCF$xAxis,y = CellsCCF$CcfMatrix[,205], type="l", xlim = c(-0.01,0.01))


setkey(testTable, "SpikeTime")

tmpCCF <- InVivoR::SpikeCCF(testTable[ClusterFac%in%c("394_2", "587_2"),SpikeTime,])
plot(y = tmpCCF$CCF, x = tmpCCF$xAxis, type="l", xlim=c(-0.1,0.1))


FileList <- sapply(ChannelLocationsFilter$...2, function(x){list.files(path = x, pattern = "SpikeTableThetaDelta.rds",recursive = T)})
ChannelLocationsFilter$...2[1:5]

for(i in "/alzheimer/Daniel_Data/DSC008154/DSC008154_190612_133329/"){
  filter_width_Delta <- c(0.5/500,4/500)
  fir_filter_Delta <- signal::fir1(4001, filter_width_Delta, "DC-0")
  if(0!=length(list.files(path = i, pattern = "SpikeTableThetaDelta.rds", recursive = T))) {
    message(i)
    SpikeTableThetaDelta <- readRDS(paste0(i, "OutputFolder/SpikeTableThetaDelta.rds")) 
    SpikeTableSessionIndex <- SpikeTableThetaDelta[,floor(SpikeIdx/20),]
    LFPDownSample <- readRDS(paste0(i, "OutputFolder/LFPDownSample.rds")) 
    message("loading complete")
    DeltaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[1,], FIR_FILTER = fir_filter_Delta, FiltFilt = T, CORES = 10)
    DeltaData_hb_transform <- hht::HilbertTransform(sig = ts(data = DeltaData, start = 0, frequency = 1e3))
    Delta_hb_power <- hht::HilbertEnvelope(DeltaData_hb_transform)
    
    SpikeTableThetaDelta[,`:=`(MSDeltaPower=Delta_hb_power[SpikeTableSessionIndex]),]
    message("1 complete")
    
    DeltaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[2,], FIR_FILTER = fir_filter_Delta, FiltFilt = T, CORES = 10)
    DeltaData_hb_transform <- hht::HilbertTransform(sig = ts(data = DeltaData, start = 0, frequency = 1e3))
    Delta_hb_power <- hht::HilbertEnvelope(DeltaData_hb_transform)
    
    SpikeTableThetaDelta[,`:=`(PaSDeltaPower=Delta_hb_power[SpikeTableSessionIndex]),]
    message("2 complete")
    saveRDS(SpikeTableThetaDelta, paste0(i, "OutputFolder/SpikeTableThetaDelta.rds")) 
    
  }
}

for(i in ChannelLocationsFilter$...2){
  if(0!=length(list.files(path = i, pattern = "SpikeTableThetaDelta.rds"))) {
  }
}



ThetaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Theta, FiltFilt = T, CORES = 10)
DeltaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Delta, FiltFilt = T, CORES = 10)

ThetaData_hb_transform <- hht::HilbertTransform(sig = ts(data = ThetaData, start = 0, frequency = 1e3))
Theta_hb_power <- hht::HilbertEnvelope(ThetaData_hb_transform)
ThetaData_phase <- (atan2(Re(ThetaData_hb_transform),Im(ThetaData_hb_transform))*180)/pi

DeltaData_hb_transform <- hht::HilbertTransform(sig = ts(data = DeltaData, start = 0, frequency = 1e3))
Delta_hb_power <- hht::HilbertEnvelope(DeltaData_hb_transform)
DeltaData_phase <- (atan2(Re(DeltaData_hb_transform),Im(DeltaData_hb_transform))*180)/pi

if(i==1) {
  SpikeTable[,`:=`(MSThetaPower=Theta_hb_power[SpikeTableSessionIndex], MSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], MSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], MSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], MSLoc=RefLoc[i]),]
} else {
  SpikeTable[,`:=`(PaSThetaPower=Theta_hb_power[SpikeTableSessionIndex], PaSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], PaSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], PaSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], PaSLoc=RefLoc[i]),]
}









#
