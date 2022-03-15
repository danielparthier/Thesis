library(data.table)
library(readxl)
library(ggplot2)
library(ggh4x)
library(InVivoR)
library(patchwork)

source(file = "/alzheimer/Daniel_Data/R/In-Vitro-ConnectivityAnalysis/R/ParPostProcessing.R")
source(file = "/alzheimer/Daniel_Data/R/In-Vitro-ConnectivityAnalysis/R/PlotFunctions.R")

ERPExtraction <- function(RefLoc, LFPRefChannels, RecordingFolder, DataPath) {
  SamplingRate <- 2e4
  OutputFolder <- paste0(RecordingFolder, "OutputFolder/")
  if(!dir.exists(OutputFolder)) {
    dir.create(path = OutputFolder)
  }
  
  Cores <- 10
  
  message("load Amplifier File")
  if(length(DataPath$AmplifierFile)==2) {
    AmpFile <- rbind(InVivoR::AmpFileRead(FILENAME = grep(pattern = "amplifier_MS", x = DataPath$AmplifierFile, value = T), ChannelNumber = 32),
                     InVivoR::AmpFileRead(FILENAME = grep(pattern = "amplifier_PaS", x = DataPath$AmplifierFile, value = T), ChannelNumber = 32))
  } else {
    AmpFile <- InVivoR::AmpFileRead(FILENAME = DataPath$AmplifierFile, ChannelNumber = 64)
  }
  message("File Loaded")
  message("Calculate StimMat")
  ### MS
  message("MS StimMat")
  StimTraceMS <- InVivoR::StimFileRead(FILENAME = DataPath$DigitalIn, digital = T)$Output
  StimMatMS <- InVivoR::StimulusSequence(raw = StimTraceMS, sampling_frequency = SamplingRate, threshold =0.5*max(StimTraceMS), max_time_gap = 10, digital = T, CORES = Cores)
  saveRDS(object = StimMatMS, file = paste0(OutputFolder, "StimMatMS.rds"), compress = F)
  ### PaS
  message("PaS StimMat")
  StimTracePaS <- InVivoR::StimFileRead(FILENAME = DataPath$AnalogIn, digital = F)$Output

  StimTracePaS <- StimTracePaS > max(StimTracePaS)/2
  if(RecordingFolder=="/alzheimer/Daniel_Data/DSC008155/DSC008155_190613_105423/") {
    StimTracePaS[1:2e6] <- 0
  }
  StimMatPaS <- InVivoR::StimulusSequence(raw = StimTracePaS, sampling_frequency = SamplingRate, threshold =0.5*max(StimTracePaS), max_time_gap = 10, digital = T, CORES = Cores)
  saveRDS(object = StimMatPaS, file = paste0(OutputFolder, "StimMatPaS.rds"), compress = F)
  
  message("Extract Reference Channels")
  LFPRefSignal <- AmpFile[LFPRefChannels,] 
  LFPDownSample <- matrix(data = 0, nrow = dim(LFPRefSignal)[1], ncol = floor(dim(LFPRefSignal)[2]/20)) 
  filter_width <- 0.5e3/(SamplingRate/2)
  fir_filter501 <- signal::fir1(501, filter_width, type = "low")
  for(RefNr in 1:dim(LFPRefSignal)[1]){
    LFPDownSample[RefNr,] <- InVivoR::decimate(SIGNAL = LFPRefSignal[RefNr,], FIR_FILTER = fir_filter501, M = 20, CORES = Cores)
    message("Reference Channel ", RefNr)
  }  
  
  message("ERP extraction")
  for(i in seq_len(dim(LFPDownSample)[1])) {
    #### MS Stim
    ERPListMS <- InVivoR::ERPList(Trace = LFPDownSample[i,],
                                  BlockMat = StimMatMS$BlockMat,
                                  SamplingFreqStim = SamplingRate,
                                  SamplingFreqTrace = 1e3,
                                  FixStartLength = 5,
                                  WindowLength = 15)
    saveRDS(object = ERPListMS, file = paste0(OutputFolder,"StimMSrec", RefLoc[i],"_StimERP.rds"), compress = F)
    
    ### Process MS ERP
    getProtocols <- grep(pattern = "^2.0_pulse|^4.0_pulse|^8.0_pulse|^16.0_pulse|^32.0_pulse",x = ERPListMS$ProtocolNames)
    
    ERPWTListMSOut <- lapply(X = getProtocols, FUN = function(x) {
      StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = ERPListMS$ProtocolNames[x], split = "_"))[1]))
      if(grepl(pattern = "burst", x = ERPListMS$ProtocolNames[x])) {
        BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = ERPListMS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
      } else {
        BurstFreq <- 0
      }
      
      ERPWT <- InVivoR::WTbatch(ERPMat = ERPListMS$ERP[[x]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F, PhaseAnalysis = F)
      # PhaseRho <- InVivoR::PhaseListAnalysis(x = atan2(Im(ERPWT$Raw), Re(ERPWT$Raw)))
      # PowerWT <- abs(ERPWT$Raw)^2
      # PowerWTFreqCorr <- PowerWT*array(data = rep(seq(0.25,80,0.25), each=dim(PowerWT)[1]), dim = dim(PowerWT))
      # WTPSD <- array(data = 0, dim = dim(PowerWTFreqCorr))
      # for(sliceNr in 1:dim(PowerWT)[3]) {
      #   WTPSD[,,sliceNr] <- PowerWTFreqCorr[,,sliceNr]/sum(PowerWTFreqCorr[,,sliceNr])*length(PowerWTFreqCorr[,,sliceNr])
      # }
      
      # list(PhaseRho=PhaseRho, PowerWT=PowerWT, PowerWTFreqCorr=PowerWTFreqCorr, WTPSD=WTPSD,Frequencies=seq(0.25,80,0.25), StimulationFrequency=StimulationFrequency, BurstFreq=BurstFreq, N=dim(ERPWT$Raw)[3], ProtName=ERPListMS$ProtocolNames[x])
      ERPWT
      })
    StimFrequency <- sapply(X = getProtocols, FUN = function(x){as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = ERPListMS$ProtocolNames[x], split = "_"))[1]))})
    BurstFrequency <- sapply(X = getProtocols, FUN = function(x){
      if(grepl(pattern = "burst", x = ERPListMS$ProtocolNames[x])) {
        BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = ERPListMS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
      } else {
        BurstFreq <- 0
      }
      BurstFreq
    })
    saveRDS(object = list(WT=ERPWTListMSOut, Protocols=getProtocols, StimFrequency=StimFrequency, BurstFrequency=BurstFrequency, StimLoc = "MS", RecLoc=RefLoc[i], StimFrequency=seq(0.25,80,0.25)), file = paste0(OutputFolder,"StimMSrec", RefLoc[i],"_WTs.rds"), compress = F)
    
    #### PaS Stim
    ERPListPaS <- InVivoR::ERPList(Trace = LFPDownSample[i,],
                                   BlockMat = StimMatPaS$BlockMat,
                                   SamplingFreqStim = SamplingRate,
                                   SamplingFreqTrace = 1e3,
                                   FixStartLength = 5,
                                   WindowLength = 15)
    saveRDS(object = ERPListPaS, file = paste0(OutputFolder,"StimPaSrec", RefLoc[i],"_StimERP.rds"), compress = F)
    
    ### Process PaS ERP
    ### Process MS ERP
    getProtocols <- grep(pattern = "^2.0_pulse|^4.0_pulse|^8.0_pulse|^16.0_pulse|^32.0_pulse",x = ERPListPaS$ProtocolNames)
    
    ERPWTListPaSOut <- lapply(X = getProtocols, FUN = function(x) {
      StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = ERPListPaS$ProtocolNames[x], split = "_"))[1]))
      if(grepl(pattern = "burst", x = ERPListPaS$ProtocolNames[x])) {
        BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = ERPListPaS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
      } else {
        BurstFreq <- 0
      }
      
      ERPWT <- InVivoR::WTbatch(ERPMat = ERPListPaS$ERP[[x]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F, PhaseAnalysis = F)
      ERPWT
    })
    StimFrequency <- sapply(X = getProtocols, FUN = function(x){as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = ERPListPaS$ProtocolNames[x], split = "_"))[1]))})
    BurstFrequency <- sapply(X = getProtocols, FUN = function(x){
      if(grepl(pattern = "burst", x = ERPListPaS$ProtocolNames[x])) {
        BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = ERPListPaS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
      } else {
        BurstFreq <- 0
      }
      BurstFreq
    })
    saveRDS(object = list(WT=ERPWTListPaSOut, Protocols=getProtocols, StimFrequency=StimFrequency, BurstFrequency=BurstFrequency, StimLoc = "PaS", RecLoc=RefLoc[i], StimFrequency=seq(0.25,80,0.25)), file = paste0(OutputFolder,"StimPaSrec", RefLoc[i],"_WTs.rds"), compress = F)
    
    rm(list = c("ERPWTListPaSOut", "ERPWTListMSOut"))
  }
  rm(list = c("AmpFile", "StimTraceMS", "StimTracePaS"))
  gc()
}


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

DataPaths <- lapply(X = ChannelLocations$...2, function(x) {
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

for(FolderNr in seq_along(ChannelLocations$...2)) {
  RefLoc <- as.vector(na.omit(unlist(ChannelLocations[FolderNr,9:11])))
  LFPRefChannels <- as.vector(na.omit(unlist(ChannelLocations[FolderNr,2:6])))
  print(data.frame(RefLoc, LFPRefChannels))
  RecordingFolder <- ChannelLocations$...2[FolderNr]
  ERPExtraction(RefLoc = RefLoc, LFPRefChannels = LFPRefChannels, RecordingFolder = RecordingFolder, DataPath = DataPaths[[FolderNr]])
}

##########
##########
FileList <- grep(pattern = "ERP", x = list.files(path = paste0(ChannelLocations$...2, "OutputFolder"), full.names = T), value = T)
FileList <- FileList[81:108]
SessionNames <- dirname(dirname(FileList))
Frequency <- seq(0.25,80,0.25)

for(FileCount in seq_along(FileList)) {

  if(file.exists(gsub(pattern = "StimERP.rds",replacement = "PowerFrequency.rds", x = FileList[FileCount]))) {
    message(gsub(pattern = "StimERP.rds",replacement = "PowerFrequency.rds", x = FileList[FileCount]), ": File exists")
    next
  }
  DirName <- dirname(dirname(FileList[FileCount]))
  AnimalID <- unlist(strsplit(x = FileList[FileCount], split = "/"))[4]
  Date <- strtoi(tail(unlist(strsplit(x = dirname(FileList[FileCount]), split = "_")), n = 2)[1])
  #tmpWT <- readRDS(file = FileList[FileCount])
  tmpERP <- readRDS(file = FileList[FileCount])
  getProtocols <- grep(pattern = "^2.0_pulse|^4.0_pulse|^8.0_pulse|^16.0_pulse|^32.0_pulse",x = tmpERP$ProtocolNames)
  ####
  tmpWT <- lapply(X = getProtocols, FUN = function(x) {
    StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = tmpERP$ProtocolNames[x], split = "_"))[1]))
    if(grepl(pattern = "burst", x = tmpERP$ProtocolNames[x])) {
      BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = tmpERP$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
    } else {
      BurstFreq <- 0
    }
    
    ERPWT <- InVivoR::WTbatch(ERPMat = tmpERP$ERP[[x]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F, PhaseAnalysis = F)
    list(WT=ERPWT, StimFrequency=StimulationFrequency, BurstFrequency=BurstFreq)
  })
  getProtocols <- tmpERP$ProtocolNames[getProtocols]
  tmpWTTable <- rbindlist(lapply(X = seq_len(length(tmpWT)), FUN = function(x) {
    Power <- abs(tmpWT[[x]]$WT$Raw)
    PowerFreqCorrected <- (Power*array(data = rep(seq(0.25,80,0.25), each=dim(Power)[1]), dim = dim(Power)))^2
    DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
    WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd
    Power <- Power^2
    tmpTable <- rbindlist(list(data.table(Power=as.vector(colMeans(Power[1:5000,,])),
               PowerFreqCorrected=as.vector(colMeans(PowerFreqCorrected[1:5000,,])),
               PSD=as.vector(colMeans(WTPSD[1:5000,,])),
               Stimulation="Pre",
               Frequency=rep(Frequency, times=dim(Power)[3]),
               Trial=rep(seq_len(dim(Power)[3]), each=length(Frequency))),
    data.table(Power=as.vector(colMeans(Power[5001:10000,,])),
               PowerFreqCorrected=as.vector(colMeans(PowerFreqCorrected[5001:10000,,])),
               PSD=as.vector(colMeans(WTPSD[5001:10000,,])),
               Stimulation="Stim",
               Frequency=rep(Frequency, times=dim(Power)[3]),
               Trial=rep(seq_len(dim(Power)[3]), each=length(Frequency))),
    data.table(Power=as.vector(colMeans(Power[10001:15000,,])),
               PowerFreqCorrected=as.vector(colMeans(PowerFreqCorrected[10001:15000,,])),
               PSD=as.vector(colMeans(WTPSD[10001:15000,,])),
               Stimulation="Post",
               Frequency=rep(Frequency, times=dim(Power)[3]),
               Trial=rep(seq_len(dim(Power)[3]), each=length(Frequency)))))
    tmpTable[,`:=`(StimulationFrequency = tmpWT[[x]]$StimFrequency, BurstFrequency=tmpWT[[x]]$BurstFrequency, Protocol=getProtocols[x], AnimalID=AnimalID, Date=Date, DirName=DirName),]
  }))
  Locations <- gsub(pattern = "Stim|_StimERP.rds", replacement = "", x = unlist(strsplit(x = gsub(pattern = paste0(dirname(FileList[FileCount]), "/"), replacement = "", x = FileList[FileCount]), split = "rec")))
  
  PhaseAnalysis <- lapply(X = seq_len(length(tmpWT)), FUN = function(x) {
    Phase <- InVivoR::PhaseListAnalysis(x = atan2(Im(tmpWT[[x]]$WT$Raw), Re(tmpWT[[x]]$WT$Raw)), CORES = 10)
    list(PhaseAnalysis=Phase, StimulationFrequency = tmpWT[[x]]$StimFrequency, BurstFrequency=tmpWT[[x]]$BurstFrequency, Protocol=getProtocols[x], AnimalID=AnimalID, Date=Date, DirName=DirName, Count=dim(tmpWT[[x]]$WT$Raw)[3])
  })
  
  if(Locations[1]=="MS") {
    tmpWTTable[,StimLoc:=ChannelLocations[...2==paste0(SessionNames[FileCount], "/"),MSLocation,],]
  } else if(!is.na(ChannelLocations[...2==paste0(SessionNames[FileCount], "/"),PaS_Channel,])) {
    tmpWTTable[,StimLoc:="PaS",]
  } else {
    tmpWTTable[,StimLoc:=ChannelLocations[...2==paste0(SessionNames[FileCount], "/"),OtherPaSLocation,],]
  }
  tmpWTTable[,RecLoc:=Locations[2],]
  message("saving file: ", FileCount)
  saveRDS(object = tmpWT, file = gsub(pattern = "StimERP.rds",replacement = "PowerFrequency.rds", x = FileList[FileCount]), compress = F)
  saveRDS(object = PhaseAnalysis, file = gsub(pattern = "StimERP.rds",replacement = "PhaseAnalysis.rds", x = FileList[FileCount]), compress = F)
  message("saved", gsub(pattern = "StimERP.rds",replacement = "PowerFrequency.rds", x = FileList[FileCount]))
  gc()
}




FileListPhaseAnalysis <- grep(pattern = "PhaseAnalysis", x = list.files(path = paste0(ChannelLocations$...2, "OutputFolder"), full.names = T), value = T)
PhaseTable <- rbindlist(lapply(X = FileListPhaseAnalysis, function(x) {
  PhaseList <- readRDS(x)
  Frequency <- seq(0.25,80,0.25)
  rbindlist(lapply(X = PhaseList, function(Phase) {
    tmpPhaseTable <- rbindlist(list(data.table(Phase=as.vector(colMeans(Phase$PhaseAnalysis$Rho[1:5000,])),
                                                    Stimulation="Pre",
                                                    Frequency=Frequency),
                                         data.table(Phase=as.vector(colMeans(Phase$PhaseAnalysis$Rho[5001:10000,])),
                                                    Stimulation="Stim",
                                                    Frequency=Frequency),
                                         data.table(Phase=as.vector(colMeans(Phase$PhaseAnalysis$Rho[10001:15000,])),
                                                    Stimulation="Post",
                                                    Frequency=Frequency)))
    tmpCount <- Phase$Count
    tmpStimulationFrequency <- Phase$StimulationFrequency
    tmpBurstFrequency <- Phase$BurstFrequency
    tmpProtocol <- Phase$Protocol
    tmpAnimalID <- Phase$AnimalID
    tmpDirName <- Phase$DirName
    tmpDate <- Phase$Date
    tmpPhaseTable[,`:=`(Count=tmpCount, StimulationFrequency=tmpStimulationFrequency, BurstFrequency=tmpBurstFrequency, Protocol=tmpProtocol, AnimalID=tmpAnimalID, DirName=tmpDirName, Date=tmpDate),]
    Locations <- gsub(pattern = "Stim|_PhaseAnalysis.rds", replacement = "", x = unlist(strsplit(x = gsub(pattern = paste0(dirname(x), "/"), replacement = "", x = x), split = "rec")))
    SessionName <- dirname(dirname(x))
    if(Locations[1]=="MS") {
      tmpPhaseTable[,StimLoc:=ChannelLocations[...2==paste0(SessionName, "/"),MSLocation,],]
    } else if(!is.na(ChannelLocations[...2==paste0(SessionName, "/"),PaS_Channel,])) {
      tmpPhaseTable[,StimLoc:="PaS",]
    } else {
      tmpPhaseTable[,StimLoc:=ChannelLocations[...2==paste0(SessionName, "/"),OtherPaSLocation,],]
    }
    tmpPhaseTable[,RecLoc:=Locations[2],]
    }))
  
}))

#### Stimulation Lines ####
justLoad <- T
if(!justLoad) {
  StimulationLine <- data.table(StimulationFrequency=c(2,4,8,16,32), Frequency=c(2,4,8,16,32))
StimulationLine[,StimulationFrequencyString:=paste0(StimulationFrequency, "Hz"),][
  ,StimulationFrequencyString:=factor(StimulationFrequencyString, levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
saveRDS(object = StimulationLine, file = "/alzheimer/Daniel_Data/R/Thesis/Data/StimulationLine.rds")
} else {
StimulationLine <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/StimulationLine.rds")
}
#### Phase Table ####
if(!justLoad) {
PhaseTable[,RandomRho:=1/sqrt(Count),][,RhoLogOdds:=log((Phase/(1-Phase))/(RandomRho/(1-RandomRho))),]
PhaseTable[,Stimulation:=factor(Stimulation, levels = c("Pre", "Stim", "Post")),]
PhaseTable[RecLoc=="through", RecLoc:="PaS",][StimLoc=="through", StimLoc:="PaS",]
PhaseTable[grepl(x=RecLoc, pattern ="MEC"), RecLoc:="MEC",][grepl(x=StimLoc, pattern ="MEC"), StimLoc:="MEC",][,StimulationFrequencyString:=paste0(StimulationFrequency, "Hz"),][,StimulationFrequencyString:=factor(StimulationFrequencyString, levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
PhaseTable[,StimulationLocation:="Posterior",][StimLoc %in% c("MS", "LS"),StimulationLocation:="Anterior",][StimLoc %in% c("SI", "AON"), StimulationLocation:="Off",][,RecLoc:=factor(x = RecLoc, levels = c("MS", "LS", "PaS", "MEC", "PrS", "Sub", "DG-mo", "alv", "AON", "SI")),]
saveRDS(object = PhaseTable, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PhaseTable.rds")
} else {
  PhaseTable <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/PhaseTable.rds")
}

#### Average Phase Table ####
if(!justLoad) {
PhaseTableAvg <- PhaseTable[,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][
  ,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][
    ,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), RhoLogOddssd=sd(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T), Phasesd=sd(Phase, na.rm = T), Count=.N),by=.(Frequency, StimulationFrequency, Stimulation, StimulationLocation, RecLoc)][
      ,StimulationFrequencyString:=paste0(StimulationFrequency, "Hz"),][
      ,StimulationFrequencyString:=factor(StimulationFrequencyString, levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
saveRDS(object = PhaseTableAvg, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PhaseTableAvg.rds")
} else {
  PhaseTableAvg <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/PhaseTableAvg.rds")
}

#### Plot Rho Plots ####
### RhoLogOdd
if(!justLoad) {
LFPRhoLogOddsPlotAnterior <- ggplot(data = PhaseTableAvg[StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = Frequency, y = RhoLogOdds, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  ggtitle("MS Stimulation")+
#  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
  facet_grid(cols = vars(RecLoc), rows = vars(StimulationFrequencyString))+
  scale_y_continuous(expand = c(0,0), name = expression(rho*" Log Odds Ratio"), limits = c(-1,2.5))+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(0:6), limits = c(1,80))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = LFPRhoLogOddsPlotAnterior, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoLogOddsPlotAnterior.rds")
} else {
  LFPRhoLogOddsPlotAnterior <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoLogOddsPlotAnterior.rds")
}

if(!justLoad) {
LFPRhoLogOddsPlotPosterior <- ggplot(data = PhaseTableAvg[StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = Frequency, y = RhoLogOdds, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  ggtitle("Fibre Stimulation")+
  #  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
  facet_grid(cols = vars(RecLoc), rows = vars(StimulationFrequencyString))+
  scale_y_continuous(expand = c(0,0), name = expression(rho*" Log Odds Ratio"), limits = c(-1,2.5))+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(0:6), limits = c(1,80))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = LFPRhoLogOddsPlotPosterior, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoLogOddsPlotPosterior.rds")
} else {
  LFPRhoLogOddsPlotPosterior <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoLogOddsPlotPosterior.rds")
}

### Rho
if(!justLoad) {
LFPRhoRawPlotAnterior <- ggplot(data = PhaseTableAvg[StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = Frequency, y = Phase, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  ggtitle("MS Stimulation")+
  #  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
  facet_grid(cols = vars(RecLoc), rows = vars(StimulationFrequencyString))+
  scale_y_continuous(limits = c(0,1), name = expression(rho))+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(0:6), limits = c(1,80))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = LFPRhoRawPlotAnterior, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoRawPlotAnterior.rds")
} else {
  LFPRhoRawPlotAnterior <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoRawPlotAnterior.rds")
}

if(!justLoad) {
LFPRhoRawPlotPosterior <- ggplot(data = PhaseTableAvg[StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = Frequency, y = Phase, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  ggtitle("Fibre Stimulation")+
  #  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
  facet_grid(cols = vars(RecLoc), rows = vars(StimulationFrequencyString))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0), name = expression(rho))+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(0:6), limits = c(1,80))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0)) 
saveRDS(object = LFPRhoRawPlotPosterior, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoRawPlotPosterior.rds")
} else {
  LFPRhoRawPlotPosterior <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPRhoRawPlotPosterior.rds")
}


if(!justLoad) {
LFPFrequencyRhoCrunchAnteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = Phase, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
#  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = Phase, group=interaction(Stimulation, Date, AnimalID)), alpha=0.2)+
  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(0,1), name = expression(rho))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoCrunchAnteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoCrunchAnteriorPlot.rds")
} else {
  LFPFrequencyRhoCrunchAnteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoCrunchAnteriorPlot.rds")
}

if(!justLoad) {
LFPFrequencyRhoSessionAnteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = Phase, group=Stimulation, colour=Stimulation))+
  #geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = Phase, group=interaction(Stimulation, Date, AnimalID)))+
  #geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(0,1), name = expression(rho))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoSessionAnteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoSessionAnteriorPlot.rds")
} else {
  LFPFrequencyRhoSessionAnteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoSessionAnteriorPlot.rds")
}

if(!justLoad) {
LFPFrequencyRhoCrunchPosteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = Phase, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
 # geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = Phase, group=interaction(Stimulation, Date, AnimalID)), alpha=0.2)+
  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(0,1), name = expression(rho))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoCrunchPosteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoCrunchPosteriorPlot.rds")
} else {
  LFPFrequencyRhoCrunchPosteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoCrunchPosteriorPlot.rds")
}

if(!justLoad) {
LFPFrequencyRhoSessionPosteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = Phase, group=Stimulation, colour=Stimulation))+
#  geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = Phase, group=interaction(Stimulation, Date, AnimalID)))+
#  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(0,1), name = expression(rho))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoSessionPosteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoSessionPosteriorPlot.rds")
} else {
  LFPFrequencyRhoSessionPosteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoSessionPosteriorPlot.rds")
}

#### Frequency Log-Odds Plot
if(!justLoad) {
LFPFrequencyRhoLogOddsCrunchAnteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
#  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=interaction(Stimulation, Date, AnimalID)), alpha=0.2)+
  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(-1,5), name = "Log-Odds-Ratio")+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoLogOddsCrunchAnteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsCrunchAnteriorPlot.rds")
} else {
  LFPFrequencyRhoLogOddsCrunchAnteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsCrunchAnteriorPlot.rds")
}


if(!justLoad) {
LFPFrequencyRhoLogOddsSessionAnteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=Stimulation, colour=Stimulation))+
#  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=interaction(Stimulation, Date, AnimalID)))+
#  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(-1,5), name = "Log-Odds-Ratio")+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoLogOddsSessionAnteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsSessionAnteriorPlot.rds")
} else {
  LFPFrequencyRhoLogOddsSessionAnteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsSessionAnteriorPlot.rds")
}

if(!justLoad) {
LFPFrequencyRhoLogOddsCrunchPosteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
#  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=interaction(Stimulation, Date, AnimalID)))+
  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(-1,5), expand = c(0,0), name = "Log-Odds-Ratio")+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoLogOddsCrunchPosteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsCrunchPosteriorPlot.rds")
} else {
  LFPFrequencyRhoLogOddsCrunchPosteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsCrunchPosteriorPlot.rds")
}

if(!justLoad) {
LFPFrequencyRhoLogOddsSessionPosteriorPlot <- ggplot(data = PhaseTableAvg[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=Stimulation, colour=Stimulation))+
#  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line(data = PhaseTable[StimulationFrequency==Frequency&StimulationLocation=="Posterior"&RecLoc%in%c("MEC", "PaS", "PrS"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T), Phase=mean(Phase, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation, StimulationLocation, StimLoc, RecLoc)], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=interaction(Stimulation, Date,AnimalID)))+
#  geom_line()+
  facet_wrap(~RecLoc)+
  scale_y_continuous(limits = c(-1,5), expand = c(0,0), name = "Log-Odds-Ratio")+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(1:5), limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 
saveRDS(object = LFPFrequencyRhoLogOddsSessionPosteriorPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsSessionPosteriorPlot.rds")
} else {
  LFPFrequencyRhoLogOddsSessionPosteriorPlot <- readRDS(file = "/alzheimer/Daniel_Data/R/Thesis/Data/LFPFrequencyRhoLogOddsSessionPosteriorPlot.rds")
}
#############

##### Rho Model - Asymptotic regression model ####
#### chance level of rho is 1/sqrt(N)
LogOddsRho <-   PhaseTable[Count>1&StimulationFrequency==Frequency&RecLoc%in%c("MEC", "PaS")&StimulationLocation%in%c("Anterior", "Posterior"),]
LogOddsRho[,Trial:=paste0(AnimalID, Date, StimulationFrequency, RecLoc ,BurstFrequency, StimulationLocation),][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")][
  ,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation")),]

stanDataPhaseRhoLogOdds <- list(
  N = LogOddsRho[,.N,],
  RhoLogOdds = LogOddsRho[,RhoLogOdds,],
  TrialN = max(as.integer(factor(LogOddsRho[,Trial,]))),
  Trial = as.integer(factor(LogOddsRho[,Trial,])),
  StimulationInput = model.matrix(object = RhoLogOdds ~ 0 + Stimulation, data = LogOddsRho),
  StimulationFrequency = LogOddsRho[,StimulationFrequency,],
  Stimulation = as.integer(factor(LogOddsRho[,Stimulation,])),
  StimulationN = max(as.integer(factor(LogOddsRho[,Stimulation,]))),
  StimLoc = as.integer(factor(LogOddsRho[,StimulationLocation,])),
  StimLocN = max(as.integer(factor(LogOddsRho[,StimulationLocation,]))),
  Location = as.integer(factor(LogOddsRho[,RecLoc,])),
  LocationN = max(as.integer(factor(LogOddsRho[,RecLoc,]))),
  DateN = max(as.integer(factor(LogOddsRho[,Date,]))),
  Date = as.integer(factor(LogOddsRho[,Date,])),
  AnimalN = max(as.integer(factor(LogOddsRho[,AnimalID,]))),
  Animal = as.integer(factor(LogOddsRho[,AnimalID,])),
  BurstFrequencyN = max(as.integer(factor(LogOddsRho[,BurstFrequency,]))),
  BurstFrequency = as.integer(factor(LogOddsRho[,BurstFrequency,])),
  priorOnly = 1
) 

if(!justLoad) {
  rstudioapi::jobRunScript(path = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/R/RunModelPrior.R", workingDir = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/", importEnv = T, exportEnv = "R_GlobalEnv")
  rstudioapi::jobRunScript(path = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/R/RunModel.R", workingDir = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/", importEnv = T, exportEnv = "R_GlobalEnv")
} else {
  fitmodLogOdds <- readRDS("/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/output/ModelLogOdds.rds")
  fitmodLogOddsPrior <- readRDS("/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/output/ModelLogOddsPrior.rds")
  
}
fitmodLogOdds$print()
#fitmodLogOdds <- fitmodLogOddsPrior
#fitmodLogOdds <- readRDS(file = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/output/ModelLogOdds.rds")


np_lin <- bayesplot::nuts_params(object = fitmodLogOddsPrior)
fitmodLogOddsPrior$print()
bayesplot::mcmc_trace(x =fitmodLogOddsPrior$draws(c("sigma", "betaSlope", "betaMax", "beta0")), np = np_lin)
y_rep <- posterior::as_draws_matrix(x = fitmodLogOddsPrior$draws("y_rep"))
bayesplot::ppc_intervals(y = LogOddsRho$RhoLogOdds, yrep = y_rep)

np_lin <- bayesplot::nuts_params(object = fitmodLogOdds)
fitmodLogOdds$print()
bayesplot::mcmc_trace(x =fitmodLogOdds$draws(c("sigma", "betaSlope", "betaMax", "beta0")), np = np_lin)
y_rep <- posterior::as_draws_matrix(x = fitmodLogOdds$draws("y_rep"))
bayesplot::ppc_intervals(y = LogOddsRho$RhoLogOdds, yrep = y_rep)
bayesplot::ppc_dens_overlay_grouped(y = LogOddsRho$RhoLogOdds, group = LogOddsRho[,paste0(RecLoc,Stimulation),], yrep = y_rep[1:200,])


bayesplot::mcmc_pairs(x =fitmodLogOdds$draws(c("sigma", "betaSlope", "betaMax", "beta0")), np = np_lin)

loo_result <- fitmodLogOdds$loo(cores = 1, r_eff = T, variables = "y_lik", save_psis = TRUE)
print(loo_result)
plot(loo_result)

LogOddsRho[loo_result$diagnostics$pareto_k>0.6,,]

bayesplot::ppc_loo_pit_overlay(
  y =  LogOddsRho$RhoLogOdds,
  yrep = y_rep,
  lw = weights(loo_result$psis_object),
  samples = 500,
  alpha = 0.7
)
      

#fitmodLogOdds <- fitmodLogOddsPrior
betaMaxPreA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,1]")),ncol = 2, nrow = 1e4)
betaMaxStimA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,1]")),ncol = 2, nrow = 1e4)
betaMaxPostA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1))+ matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,1]")),ncol = 2, nrow = 1e4)
betaMaxPreP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,1]")),ncol = 2, nrow = 1e4)
betaMaxStimP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,1]")),ncol = 2, nrow = 1e4)
betaMaxPostP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 1, CondMatCol = c(1,1,1,1))+ matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,1]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,1]")),ncol = 2, nrow = 1e4)

beta0PreA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,2]")),ncol = 2, nrow = 1e4) 
beta0StimA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,2]")),ncol = 2, nrow = 1e4) 
beta0PostA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,2]")),ncol = 2, nrow = 1e4) 
beta0PreP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,2]")),ncol = 2, nrow = 1e4) 
beta0StimP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,2]")),ncol = 2, nrow = 1e4) 
beta0PostP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "beta0", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 2, CondMatCol = c(2,2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,2]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,2]")),ncol = 2, nrow = 1e4) 

betaSlopePreA <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,3]")),ncol = 2, nrow = 1e4))
betaSlopeStimA <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,3]")),ncol = 2, nrow = 1e4))
betaSlopePostA <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,3]")),ncol = 2, nrow = 1e4))
betaSlopePreP <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,3]")),ncol = 2, nrow = 1e4))
betaSlopeStimP <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,3]")),ncol = 2, nrow = 1e4))
betaSlopePostP <- exp(PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log",MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,3]")),ncol = 2, nrow = 1e4))
# 
# betaSlopePreA <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,3]")),ncol = 2, nrow = 1e4))
# betaSlopeStimA <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,3]")),ncol = 2, nrow = 1e4))
# betaSlopePostA <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,3]")),ncol = 2, nrow = 1e4))
# betaSlopePreP <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[1,3]")),ncol = 2, nrow = 1e4))
# betaSlopeStimP <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[2,3]")),ncol = 2, nrow = 1e4))
# betaSlopePostP <- (PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", MatrixString = "l", CondMatString = c("a", "d", "b", "t"), MatrixCol = 3, CondMatCol = c(3,3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("s[3,3]")),ncol = 2, nrow = 1e4))
# 



BetaMaxRhoLogOddsOverview <- cbind(betaMaxPreA, betaMaxStimA, betaMaxPostA,betaMaxPreP, betaMaxStimP, betaMaxPostP)
colnames(BetaMaxRhoLogOddsOverview) <- paste(rep(c("PaS", "MEC"), times = 6),rep(c("MS", "fibre"), each = 6), rep(rep(c("Pre", "Stim", "Post"), each=2), times=2), sep = "-")
saveRDS(BetaMaxRhoLogOddsOverview, "/alzheimer/Daniel_Data/R/Thesis/Data/BetaMaxRhoLogOddsOverview.rds")
saveRDS(HDIcalc(x = BetaMaxRhoLogOddsOverview), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaMaxRhoLogOddsSummary.rds")
saveRDS(HDIcalc(DifferenceMat(BetaMaxRhoLogOddsOverview, GroupStrings = colnames(BetaMaxRhoLogOddsOverview))), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaMaxRhoLogOddsOverviewDiff.rds")

BetaSlopeRhoLogOddsOverview <- cbind(betaSlopePreA, betaSlopeStimA,betaSlopePostA, betaSlopePreP, betaSlopeStimP, betaSlopePostP)
colnames(BetaSlopeRhoLogOddsOverview) <- paste(rep(c("PaS", "MEC"), times = 6),rep(c("MS", "fibre"), each = 6), rep(rep(c("Pre", "Stim", "Post"), each=2), times=2), sep = "-")
saveRDS(BetaSlopeRhoLogOddsOverview, "/alzheimer/Daniel_Data/R/Thesis/Data/BetaSlopeRhoLogOddsOverview.rds")
saveRDS(HDIcalc(x = BetaSlopeRhoLogOddsOverview), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaSlopeRhoLogOddsSummary.rds")
saveRDS(HDIcalc(DifferenceMat(BetaSlopeRhoLogOddsOverview, GroupStrings = colnames(BetaSlopeRhoLogOddsOverview))), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaSlopeRhoLogOddsOverviewDiff.rds")




curve(expr = 1-(1-0)*exp(-0.5*x), from = 0, to = 40)
      
##### Rho Log Odds
RhoLogStim <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxStimA[,1])-(betaMaxStimA[,1]-beta0StimA[,1]) * exp(-betaSlopeStimA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxStimA[,2])-(betaMaxStimA[,2]-beta0StimA[,2]) * exp(-betaSlopeStimA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxStimP[,1])-(betaMaxStimP[,1]-beta0StimP[,1]) * exp(-betaSlopeStimP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxStimP[,2])-(betaMaxStimP[,2]-beta0StimP[,2]) * exp(-betaSlopeStimP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Stim", StimulationFrequency=x),]
}))

RhoLogPre <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPreA[,1])-(betaMaxPreA[,1]-beta0PreA[,1]) * exp(-betaSlopePreA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPreA[,2])-(betaMaxPreA[,2]-beta0PreA[,2]) * exp(-betaSlopePreA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPreP[,1])-(betaMaxPreP[,1]-beta0PreP[,1]) * exp(-betaSlopePreP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPreP[,2])-(betaMaxPreP[,2]-beta0PreP[,2]) * exp(-betaSlopePreP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Pre", StimulationFrequency=x),]
  
}))
RhoLogPost <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostA[,1])-(betaMaxPostA[,1]-beta0PostA[,1]) * exp(-betaSlopePostA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostA[,2])-(betaMaxPostA[,2]-beta0PostA[,2]) * exp(-betaSlopePostA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostP[,1])-(betaMaxPostP[,1]-beta0PostP[,1]) * exp(-betaSlopePostP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostP[,2])-(betaMaxPostP[,2]-beta0PostP[,2]) * exp(-betaSlopePostP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Post", StimulationFrequency=x),]
}))



##### no difference between Pre and Post
# RhoLogPrePost <- rbindlist(lapply(c(1:32), function(x) {
#   HDIcalc((betaMaxPre[,1])-(betaMaxPre[,1]-beta0[,1]) * exp(-betaSlope[,1]*x)-(betaMaxPost[,1])-(betaMaxPost[,1]-beta0[,1]) * exp(-betaSlope[,1]*x))
# }))

RhoLogTable <- rbindlist(l = list(RhoLogPre,RhoLogStim,RhoLogPost))

RhoLogStim <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxStimA[,1])-(betaMaxStimA[,1]-beta0StimA[,1]) * exp(-betaSlopeStimA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxStimA[,2])-(betaMaxStimA[,2]-beta0StimA[,2]) * exp(-betaSlopeStimA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxStimP[,1])-(betaMaxStimP[,1]-beta0StimP[,1]) * exp(-betaSlopeStimP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxStimP[,2])-(betaMaxStimP[,2]-beta0StimP[,2]) * exp(-betaSlopeStimP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Stim", StimulationFrequency=x),]
}))

RhoLogPre <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxPreA[,1])-(betaMaxPreA[,1]-beta0PreA[,1]) * exp(-betaSlopePreA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxPreA[,2])-(betaMaxPreA[,2]-beta0PreA[,2]) * exp(-betaSlopePreA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxPreP[,1])-(betaMaxPreP[,1]-beta0PreP[,1]) * exp(-betaSlopePreP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (pi/sqrt(3)*0.1),x = (betaMaxPreP[,2])-(betaMaxPreP[,2]-beta0PreP[,2]) * exp(-betaSlopePreP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Pre", StimulationFrequency=x),]
  
}))
RhoLogPost <- rbindlist(lapply(c(1:32), function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostA[,1])-(betaMaxPostA[,1]-beta0PostA[,1]) * exp(-betaSlopePostA[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostA[,2])-(betaMaxPostA[,2]-beta0PostA[,2]) * exp(-betaSlopePostA[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostP[,1])-(betaMaxPostP[,1]-beta0PostP[,1]) * exp(-betaSlopePostP[,1]*x)))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(TestValue = (-pi/sqrt(3)*0.1),x = (betaMaxPostP[,2])-(betaMaxPostP[,2]-beta0PostP[,2]) * exp(-betaSlopePostP[,2]*x)))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Post", StimulationFrequency=x),]
}))

RhoLogTable$OutOfInterval <- RhoLogTable$OutOfInterval*rbindlist(l = list(RhoLogPre,RhoLogStim,RhoLogPost))$OutOfInterval

RhoLogTable[,`:=`(Stimulation=factor(Stimulation, levels = c("Pre", "Stim", "Post")), RecLoc=factor(RecLoc, levels=c("PaS", "MEC"))),]

RhoLogTable[StimulationFrequency%in%c(2,4,8,16,32),]
RhoLogTable[StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation"),][
  ,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation")),]
saveRDS(object = RhoLogTable, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogTable.rds")


MSStimulationRhoLogOdds <- ggplot(data = RhoLogTable[StimulationLocation=="Anterior"], aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=LogOddsRho[StimulationLocation=="Anterior"], mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
 # geom_line(data=LogOddsRho[StimulationLocation=="Anterior"], mapping = aes(x=StimulationFrequency, y=RhoLogOdds, group=DirName))+
 # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  #scale_y_continuous(name = expression(rho*" Log Odds Ratio"), limits = c(-1.5,5), breaks = seq(-1,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  ggtitle("MS Stimulation")+
 # geom_hline(yintercept = 1)+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = FibreStimulationRhoLogOdds, file = "/alzheimer/Daniel_Data/R/Thesis/Data/MSStimulationRhoLogOdds.rds")


FibreStimulationRhoLogOdds <- ggplot(data = RhoLogTable[StimulationLocation=="Posterior"], aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=LogOddsRho[StimulationLocation=="Posterior"], mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  scale_y_continuous(name = expression(rho*" Log Odds Ratio"), limits = c(-1.5,5), breaks = seq(-1,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  ggtitle("Fibre Stimulation")+
  # geom_hline(yintercept = 1)+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = FibreStimulationRhoLogOdds, file = "/alzheimer/Daniel_Data/R/Thesis/Data/FibreStimulationRhoLogOdds.rds")

RhoLogOddsFreqPlotOverview <- MSStimulationRhoLogOdds + FibreStimulationRhoLogOdds + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = RhoLogOddsFreqPlotOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsFreqPlotOverview.rds")

RhoLogOddsFreqSpectrumPlotOverview <- LFPRhoLogOddsPlotAnterior + LFPRhoLogOddsPlotPosterior + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = RhoLogOddsFreqSpectrumPlotOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsFreqSpectrumPlotOverview.rds")


RhoLogOddsModelFrequencyPlot <- ggplot(data = RhoLogTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=LogOddsRho, mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc + Stimulation, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = expression(rho*" Log Odds Ratio"), limits = c(-1.5,5), breaks = seq(-1,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = RhoLogOddsModelFrequencyPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsModelFrequencyPlot.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsModelFrequencyPlot.pdf", plot = RhoLogOddsModelFrequencyPlot, device = "pdf")

RhoLogOddsModelFrequencyFlipPlot <- ggplot(data = RhoLogTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=LogOddsRho, mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc + Stimulation, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = expression(rho*" Log Odds Ratio"), limits = c(-1.5,5), breaks = seq(-1,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = RhoLogOddsModelFrequencyFlipPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsModelFrequencyFlipPlot.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsModelFrequencyFlipPlot.pdf", plot = RhoLogOddsModelFrequencyFlipPlot, device = "pdf")

RhoLogOddsN <- merge.data.table(x = LogOddsRho[,.N,by=.(StimulationLocationString, RecLoc, AnimalID)][,.(AnimalN=.N),by=.(StimulationLocationString, RecLoc)], LogOddsRho[,.N,by=.(StimulationLocationString, RecLoc, DirName)][,.(SessionN=.N),by=.(StimulationLocationString, RecLoc)])
saveRDS(object = RhoLogOddsN, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsN.rds")


PhaseTableAvg[StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")][,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation")),]
RhoLogOddsFreqSpectrumPlotFacetOverview <- ggplot(data = PhaseTableAvg[StimulationLocation!="Off"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = Frequency, y = RhoLogOdds, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin = RhoLogOdds-RhoLogOddssd/sqrt(Count), ymax = RhoLogOdds+RhoLogOddssd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  #  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
  facet_nested(StimulationFrequencyString ~ RecLoc + StimulationLocationString, nest_line = element_line(colour = "black"))+
  scale_y_continuous(expand = c(0,0), name = expression(rho*" Log Odds Ratio"), limits = c(-1,2.5))+
  scale_x_continuous(name = "Frequency (Hz)", trans = "log2", breaks = 2^(0:6), limits = c(1,80))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"), guide = guide_legend(keywidth = 0.5, title = ""))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"), guide=F)+
  theme_classic() + theme(panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = RhoLogOddsFreqSpectrumPlotFacetOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsFreqSpectrumPlotFacetOverview.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoLogOddsFreqSpectrumPlotFacetOverview.pdf", plot = RhoLogOddsFreqSpectrumPlotFacetOverview, device = "pdf")



# ggplot(data = RhoLogTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
#   geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
#   geom_point(data=LogOddsRho, mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
#   # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
#   facet_nested(RecLoc+StimulationLocationString ~ Stimulation)+
#   scale_y_continuous(name = expression(rho*" Log Odds Ratio"), limits = c(-1.5,5), breaks = seq(-1,5,1))+
#   scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(2, 4, 8, 16, 32), labels = c("2", "4", "8", "16", "32"), trans = "log2")+
#   scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#   scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#   geom_line()+
#   ggtitle("Fibre Stimulation")+
#   # geom_hline(yintercept = 1)+
#   theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 

ORComparison <- exp(cbind(betaMaxStimA-beta0StimA, betaMaxStimP-beta0StimP))
colnames(x = ORComparison) <- c("MS-PaS", "MS-MEC", "Fibre-PaS", "Fibre-MEC")
ORComparisonPlot <- TufteBoxPlot(x = ORComparison, RangeY = c(0,15), Grouping = colnames(ORComparison), BreaksY = seq(0,15,2.5), Line = 1, LabelY = expression(rho * " Odds Ratio"))
RhoModelOverview <- RhoLogOddsModelFrequencyPlot + ORComparisonPlot + plot_layout(widths = c(3,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = RhoModelOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoModelOverview.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoModelOverview.pdf", plot = RhoModelOverview, device = "pdf", width = 8, height = 4, scale = 2)

ORComparisonTable <- HDIcalc(ORComparison, TestValue = 1)
saveRDS(object = ORComparisonTable, file = "/alzheimer/Daniel_Data/R/Thesis/Data/ORComparisonTable.rds")



##### Immuno In-Vivo ####
MSImage <- tiff::readTIFF(source = "/alzheimer/Daniel_Data/R/Thesis/Data/DSC007849_PV-Cre_MS.tif", info = T, all = T, indexed = T)
MSImageDT <- data.table(reshape2::melt(MSImage))
setnames(x = MSImageDT, old = c("Var1", "Var2"), new = c("y", "x"))
MSImageDT[,`:=`(x=x/attributes(MSImage[[1]])$x.resolution, y=y/attributes(MSImage[[1]])$y.resolution),][,L1:=factor(L1, levels = c(2,1)),]
ScaleBarMSImage <- data.table(x=2500, xend=3500, y=4700, yend=4700, L1=factor(x = c(2), levels = c(2,1)))
StainingTextMSImage <- data.table(x=c(100,100,3000), y=c(100,100,4600), L1=factor(x = c(1,2,2), levels = c(2,1)), text=c("YFP","DiI", "1mm"))
MSCoronalProbe <- ggplot(data = MSImageDT, mapping = aes(x = x, y = y, fill=value))+
  geom_raster()+
  scale_fill_gradient(low = "black", high = "white")+
  scale_y_reverse()+
  coord_fixed(expand = 0)+
  facet_grid(rows = vars(L1))+
  geom_segment(inherit.aes = F, data = ScaleBarMSImage, mapping = aes(x=x, xend=xend, y=y, yend=yend), colour="white", size=1)+
  geom_text(inherit.aes = F, data = StainingTextMSImage[1:2,], aes(x = x, y = y, label=text), colour="white", size=4, vjust = "inward", hjust = "inward")+
  geom_text(inherit.aes = F, data = StainingTextMSImage[3,], aes(x = x, y = y, label=text), colour="white", size=2, vjust="inward")+
  theme_void()+
  theme(legend.position="none", panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(0, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text = element_blank()) 
saveRDS(object = MSCoronalProbe, file = "/alzheimer/Daniel_Data/R/Thesis/Data/MSCoronalProbe.rds")

PaSImage <- tiff::readTIFF(source = "/alzheimer/Daniel_Data/R/Thesis/Data/DSC007849_PV-Cre_PaS.tif", info = T, all = T, indexed = T)
PaSImageDT <- data.table(reshape2::melt(PaSImage))
setnames(x = PaSImageDT, old = c("Var1", "Var2"), new = c("y", "x"))
PaSImageDT[,`:=`(x=x/attributes(PaSImage[[1]])$x.resolution, y=y/attributes(PaSImage[[1]])$y.resolution),][,L1:=factor(L1, levels = c(3,2,1)),]

ScaleBarPaSImage <- data.table(x=2500, xend=3500, y=4700, yend=4700, L1=factor(x = c(3), levels = c(3,2,1)))
StainingTextPaSImage <- data.table(x=c(100,100,100,3000), y=c(100,100,100,4600), L1=factor(x = c(1,2,3,3), levels = c(3,2,1)), text=c("YFP","DiI", "WFS1", "1mm"))
PaSSagittalProbe <- ggplot(data = PaSImageDT, mapping = aes(x = x, y = y, fill=value))+
  geom_raster()+
  scale_fill_gradient(low = "black", high = "white")+
  scale_y_reverse()+
  coord_fixed(expand = 0)+
  facet_grid(rows = vars(L1))+
  geom_segment(inherit.aes = F, data = ScaleBarPaSImage, mapping = aes(x=x, xend=xend, y=y, yend=yend), colour="white", size=1)+
  geom_text(inherit.aes = F, data = StainingTextPaSImage[1:3,], aes(x = x, y = y, label=text), colour="white", size=4, vjust = "inward", hjust = "inward")+
  geom_text(inherit.aes = F, data = StainingTextPaSImage[4,], aes(x = x, y = y, label=text), colour="white", size=2, vjust="inward")+
  theme_void()+
  theme(legend.position="none", panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(0, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text = element_blank()) 
saveRDS(object = PaSSagittalProbe, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PaSSagittalProbe.rds")




####### Wavelet examples and Raster plots ####

StimMS_PaS <- readRDS(file = "/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190409_132939/OutputFolder/StimMSrecPaS_StimERP.rds")
grep(pattern = "32.0", x = StimMS_PaS$ProtocolNames)

StimMS_PaS <- readRDS(file = "/alzheimer/Daniel_Data/DSC008159/DSC008159_190619_150144/OutputFolder/StimMSrecPaS_StimERP.rds")
StimMS_PaS$ProtocolNames
ExampleWT <- InVivoR::WTbatch(ERPMat = StimMS_PaS$ERP[[2]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F)
###
Power <- abs(ExampleWT$Raw)
PowerFreqCorrected <- (Power*array(data = rep(seq(0.25,80,0.25), each=dim(Power)[1]), dim = dim(Power)))^2
DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd

#plot(x = StimMS_PaS$ERP[[22]][1,], type="l")
#image(InVivoR::PowerMat(Power), useRaster=T)

TraceTable <- data.table(Amplitude=StimMS_PaS$ERP[[2]][23,], Time=seq_along(StimMS_PaS$ERP[[6]][23,])/1e3)

LightStimSegmentTraceDT <- data.table(TimeMin=5, TimeMax=10, AmplitudeMin=500, AmplitudeMax = 500)
ScaleBarTraceDT <- data.table(TimeMin=c(13,0.5), TimeMax=c(15,0.5), AmplitudeMin=c(-600,-550), AmplitudeMax = c(-600,-350), Label=c("2s", "200µV"))

TraceExamplePlot <- ggplot(data = TraceTable, aes(x = Time, y = Amplitude))+
  geom_line(size=0.2)+
  scale_y_continuous(limits = c(-800, 600))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  geom_segment(inherit.aes = F, data = LightStimSegmentTraceDT, mapping = aes(x=TimeMin, xend=TimeMax, y=AmplitudeMin, yend=AmplitudeMax), colour="deepskyblue2", size=2)+
  geom_segment(inherit.aes = F, data = ScaleBarTraceDT, mapping = aes(x=TimeMin, xend=TimeMax, y=AmplitudeMin, yend=AmplitudeMax), colour="black", size=1)+
  geom_text(inherit.aes = F, data = ScaleBarTraceDT[1,], mapping = aes(x=TimeMin+1, y=AmplitudeMin, label=Label), vjust=1.4, size=3)+
  geom_text(inherit.aes = F, data = ScaleBarTraceDT[2,], mapping = aes(x=TimeMin, y=AmplitudeMin+100, label=Label), hjust=-0.2, size=3)+
  theme_void()
saveRDS(object = TraceExamplePlot, file = "Data/TraceExamplePlot.rds")


Single8HzTable <- data.table(reshape2::melt(Power[,,23]))
setnames(Single8HzTable, old = c("Var1", "Var2", "value"), new = c("Time", "Frequency", "Power"))
Single8HzTable[,`:=`(WTPSD=as.vector(WTPSD[,,23]), Frequency=Frequency/4, Time=Time/1e3),]
Single8HzTable[,Time:=round(x = Time, digits = 2),]
Single8HzTable <- Single8HzTable[,.(Power=mean(Power), WTPSD=mean(WTPSD)),by=.(Time, Frequency)]

LightStimSingle8HzSegmentDT <- data.table(TimeMin=5, TimeMax=10, FrequencyMin=63, FrequencyMax = 63)

Single8HzPowerPlot <- ggplot(data = Single8HzTable[Frequency<60,], aes(x = Time, y = Frequency, fill=Power^2/1e6))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  geom_segment(inherit.aes = F, data = LightStimSingle8HzSegmentDT, mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(option = "B", limits=c(0,6), breaks=seq(0,6,2), begin = 0, guide = guide_colourbar(title = "Power\n(A.U.)", barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = Single8HzPowerPlot, file = "Data/Single8HzPowerPlot.rds")


Single8HzWTPSDPlot <- ggplot(data = Single8HzTable[Frequency<60,], aes(x = Time, y = Frequency, fill=WTPSD))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  geom_segment(inherit.aes = F, data = LightStimSingle8HzSegmentDT, mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(option = "B", limits=c(0,25),begin = 0, guide = guide_colourbar(title = "Power\nDensity", barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = Single8HzWTPSDPlot, file = "Data/Single8HzWTPSDPlot.rds")


# DesignMatrix <- "
# AABBB
# AABBB
# AABBB
# CCCCC
# CCCCC
# DDDDD
# DDDDD
# EEEEE
# EEEEE
# "
DesignMatrix <- "
ACC
ADD
BDD
BEE
BEE
"
# DesignMatrix <- "
# ABCC
# ABDD
# #BEE
# "

WTSingleExampleOverviewPlot <- MSCoronalProbe + PaSSagittalProbe + TraceExamplePlot + Single8HzPowerPlot + Single8HzWTPSDPlot + plot_layout(design = DesignMatrix, widths = c(1,1,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = WTSingleExampleOverviewPlot, file="Data/WTSingleExampleOverviewPlot.rds")
ggsave(filename = "Data/WTSingleExampleOverviewPlot.pdf", plot = WTSingleExampleOverviewPlot, device = "pdf")


StimMS_PaS <- readRDS(file = "/alzheimer/Daniel_Data/DSC008159/DSC008159_190619_150144/OutputFolder/StimMSrecPaS_StimERP.rds")
StimMS_PaS$ProtocolNames

getProtocolsValue <- grep(pattern = "16|2|4|8|32",x = StimMS_PaS$ProtocolNames, value = T)
getProtocols <- grep(pattern = "16|2|4|8|32",x = StimMS_PaS$ProtocolNames)
FrequencySelect <- seq(0.25,80,0.25)
WTExamplePaS <- rbindlist(lapply(X = getProtocols, function(x) {
  StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = StimMS_PaS$ProtocolNames[x], split = "_"))[1]))
  if(grepl(pattern = "burst", x = StimMS_PaS$ProtocolNames[x])) {
    BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = StimMS_PaS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
  } else {
    BurstFreq <- 0
  }
  ExampleWT <- InVivoR::WTbatch(ERPMat = StimMS_PaS$ERP[[x]], frequencies = FrequencySelect, SamplingRate = 1e3, CORES = 10, compression = F)
  Phase <- InVivoR::PhaseListAnalysis(x = atan2(Im(ExampleWT$Raw), Re(ExampleWT$Raw)), CORES = 10)
  Power <- abs(ExampleWT$Raw)
  PowerFreqCorrected <- (Power*array(data = rep(FrequencySelect, each=dim(Power)[1]), dim = dim(Power)))^2
  DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
  WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd
  Power <- Power^2
  Power <- InVivoR::PowerMat(x = Power)
  tmpTable <- data.table(reshape2::melt(Power))
  setnames(tmpTable, old = c("Var1", "Var2", "value"), new = c("Time", "Frequency", "Power"))
  WTPSD <- InVivoR::PowerMat(x = WTPSD)
  tmpTable[,`:=`(WTPSD=as.vector(WTPSD), Rho=as.vector(Phase$Rho), Frequency=Frequency/4, Time=Time/1e3, BurstFreq=BurstFreq, StimulationFrequency=StimulationFrequency),]
  tmpTable
  }))


WTExamplePaS[,StimulationFrequency:=paste0(StimulationFrequency, "Hz"),][,StimulationFrequency:=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
WTExamplePaS[,Time:=round(Time, digits = 2),]
WTExamplePaS <- WTExamplePaS[,.(WTPSD=mean(WTPSD), Rho=mean(Rho), Power=mean(Power)),by=.(Time, Frequency, StimulationFrequency)]
WTExamplePaS[,`:=`(RecLoc="PaS", StimulationLocationString="MS Stimulation"),]

StimMS_MEC <- readRDS(file = "/alzheimer/Daniel_Data//DSC008159/DSC008159_190619_150144/OutputFolder/StimMSrecMEC-deep_StimERP.rds")

WTExampleMEC <- rbindlist(lapply(X = getProtocols, function(x) {
  StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = StimMS_MEC$ProtocolNames[x], split = "_"))[1]))
  if(grepl(pattern = "burst", x = StimMS_MEC$ProtocolNames[x])) {
    BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = StimMS_MEC$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
  } else {
    BurstFreq <- 0
  }
  ExampleWT <- InVivoR::WTbatch(ERPMat = StimMS_MEC$ERP[[x]], frequencies = FrequencySelect, SamplingRate = 1e3, CORES = 10, compression = F)
  Phase <- InVivoR::PhaseListAnalysis(x = atan2(Im(ExampleWT$Raw), Re(ExampleWT$Raw)), CORES = 10)
  Power <- abs(ExampleWT$Raw)
  PowerFreqCorrected <- (Power*array(data = rep(FrequencySelect, each=dim(Power)[1]), dim = dim(Power)))^2
  DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
  WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd
  Power <- Power^2
  Power <- InVivoR::PowerMat(x = Power)
  tmpTable <- data.table(reshape2::melt(Power))
  setnames(tmpTable, old = c("Var1", "Var2", "value"), new = c("Time", "Frequency", "Power"))
  WTPSD <- InVivoR::PowerMat(x = WTPSD)
  tmpTable[,`:=`(WTPSD=as.vector(WTPSD), Rho=as.vector(Phase$Rho), Frequency=Frequency/4, Time=Time/1e3, BurstFreq=BurstFreq, StimulationFrequency=StimulationFrequency),]
  tmpTable
}))

WTExampleMEC[,StimulationFrequency:=paste0(StimulationFrequency, "Hz"),][,StimulationFrequency:=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
WTExampleMEC[,Time:=round(Time, digits = 2),]
WTExampleMEC <- WTExampleMEC[,.(WTPSD=mean(WTPSD), Rho=mean(Rho), Power=mean(Power)),by=.(Time, Frequency, StimulationFrequency)]
WTExampleMEC[,`:=`(RecLoc="MEC", StimulationLocationString="MS Stimulation"),]


StimPaS_PaS <- readRDS(file = "/alzheimer/Daniel_Data//DSC008159/DSC008159_190619_150144/OutputFolder/StimPaSrecPaS_StimERP.rds")
getProtocolsValue <- grep(pattern = "16|2|4|8|32",x = StimPaS_PaS$ProtocolNames, value = T)
getProtocols <- grep(pattern = "16|2|4|8|32",x = StimPaS_PaS$ProtocolNames)

WTExamplePaSFibre <- rbindlist(lapply(X = getProtocols, function(x) {
  StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = StimPaS_PaS$ProtocolNames[x], split = "_"))[1]))
  if(grepl(pattern = "burst", x = StimPaS_PaS$ProtocolNames[x])) {
    BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = StimPaS_PaS$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
  } else {
    BurstFreq <- 0
  }
  ExampleWT <- InVivoR::WTbatch(ERPMat = StimPaS_PaS$ERP[[x]], frequencies = FrequencySelect, SamplingRate = 1e3, CORES = 10, compression = F)
  Phase <- InVivoR::PhaseListAnalysis(x = atan2(Im(ExampleWT$Raw), Re(ExampleWT$Raw)), CORES = 10)
  Power <- abs(ExampleWT$Raw)
  PowerFreqCorrected <- (Power*array(data = rep(FrequencySelect, each=dim(Power)[1]), dim = dim(Power)))^2
  DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
  WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd
  Power <- Power^2
  Power <- InVivoR::PowerMat(x = Power)
  tmpTable <- data.table(reshape2::melt(Power))
  setnames(tmpTable, old = c("Var1", "Var2", "value"), new = c("Time", "Frequency", "Power"))
  WTPSD <- InVivoR::PowerMat(x = WTPSD)
  tmpTable[,`:=`(WTPSD=as.vector(WTPSD), Rho=as.vector(Phase$Rho), Frequency=Frequency/4, Time=Time/1e3, BurstFreq=BurstFreq, StimulationFrequency=StimulationFrequency),]
  tmpTable
}))

WTExamplePaSFibre[,StimulationFrequency:=paste0(StimulationFrequency, "Hz"),][,StimulationFrequency:=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
WTExamplePaSFibre[,Time:=round(Time, digits = 2),]
WTExamplePaSFibre <- WTExamplePaSFibre[,.(WTPSD=mean(WTPSD), Rho=mean(Rho), Power=mean(Power)),by=.(Time, Frequency, StimulationFrequency)]
WTExamplePaSFibre[,`:=`(RecLoc="PaS", StimulationLocationString="Fibre Stimulation"),]

StimPaS_MEC <- readRDS(file = "/alzheimer/Daniel_Data//DSC008159/DSC008159_190619_150144/OutputFolder/StimPaSrecMEC-deep_StimERP.rds")


WTExampleMECFibre <- rbindlist(lapply(X = getProtocols, function(x) {
  StimulationFrequency <- as.integer(gsub(pattern = "[.]0", replacement = "", x = unlist(strsplit(x = StimPaS_MEC$ProtocolNames[x], split = "_"))[1]))
  if(grepl(pattern = "burst", x = StimPaS_MEC$ProtocolNames[x])) {
    BurstFreq <- as.integer(unlist(strsplit(x = unlist(strsplit(x = StimPaS_MEC$ProtocolNames[x], split = "burst_"))[2], split = "[.]"))[1])
  } else {
    BurstFreq <- 0
  }
  ExampleWT <- InVivoR::WTbatch(ERPMat = StimPaS_MEC$ERP[[x]], frequencies = FrequencySelect, SamplingRate = 1e3, CORES = 10, compression = F)
  Phase <- InVivoR::PhaseListAnalysis(x = atan2(Im(ExampleWT$Raw), Re(ExampleWT$Raw)), CORES = 10)
  Power <- abs(ExampleWT$Raw)
  PowerFreqCorrected <- (Power*array(data = rep(FrequencySelect, each=dim(Power)[1]), dim = dim(Power)))^2
  DimensionProd <- prod(dim(PowerFreqCorrected)[1:2])
  WTPSD <- PowerFreqCorrected/rep(x = colSums(matrix(data = PowerFreqCorrected, nrow = DimensionProd, ncol = dim(PowerFreqCorrected)[3])), each=DimensionProd)*DimensionProd
  Power <- Power^2
  Power <- InVivoR::PowerMat(x = Power)
  tmpTable <- data.table(reshape2::melt(Power))
  setnames(tmpTable, old = c("Var1", "Var2", "value"), new = c("Time", "Frequency", "Power"))
  WTPSD <- InVivoR::PowerMat(x = WTPSD)
  tmpTable[,`:=`(WTPSD=as.vector(WTPSD), Rho=as.vector(Phase$Rho), Frequency=Frequency/4, Time=Time/1e3, BurstFreq=BurstFreq, StimulationFrequency=StimulationFrequency),]
  tmpTable
}))

WTExampleMECFibre[,StimulationFrequency:=paste0(StimulationFrequency, "Hz"),][,StimulationFrequency:=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]
WTExampleMECFibre[,Time:=round(Time, digits = 2),]
WTExampleMECFibre <- WTExampleMECFibre[,.(WTPSD=mean(WTPSD), Rho=mean(Rho), Power=mean(Power)),by=.(Time, Frequency, StimulationFrequency)]
WTExampleMECFibre[,`:=`(RecLoc="MEC", StimulationLocationString="Fibre Stimulation"),]

WTExampleLarge <- rbindlist(l=list(WTExamplePaS, WTExampleMEC, WTExampleMECFibre, WTExamplePaSFibre))
rm(list = c("WTExamplePaS", "WTExampleMEC", "WTExampleMECFibre", "WTExamplePaSFibre"))
WTExampleLarge[,`:=`(StimulationFrequency=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")), RecLoc=factor(x = RecLoc, c("PaS", "MEC")), StimulationLocationString=factor(StimulationLocationString, levels = c("MS Stimulation","Fibre Stimulation"))),]

LightStimSegmentDT <- data.table(TimeMin=5, TimeMax=10, FrequencyMin=63, FrequencyMax = 63, StimulationFrequency = rep(c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz"), each=4), RecLoc=rep(c("PaS", "MEC"), times=10), StimulationLocationString = rep(c("MS Stimulation", "MS Stimulation","Fibre Stimulation", "Fibre Stimulation"), times=5))
LightStimSegmentDT[,`:=`(StimulationFrequency=factor(x = StimulationFrequency, c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")), RecLoc=factor(x = RecLoc, c("PaS", "MEC")), StimulationLocationString=factor(StimulationLocationString, levels = c("MS Stimulation","Fibre Stimulation"))),]
WTPSDExampleRasterPaS <- ggplot(data = WTExampleLarge[Frequency<60&RecLoc=="PaS",], aes(x = Time, y = Frequency, fill=WTPSD))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10), labels = c("0", "", "20", "","40", "","60"))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  facet_grid(cols = vars(StimulationFrequency), rows = vars(StimulationLocationString))+
  geom_segment(inherit.aes = F, data = LightStimSegmentDT[RecLoc=="PaS",], mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(limits=c(0,8),option = "B", guide = guide_colourbar(title = "Power\nDensity", barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 

RhoExampleRasterPaS <- ggplot(data = WTExampleLarge[Frequency<60&RecLoc=="PaS",], aes(x = Time, y = Frequency, fill=Rho))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10), labels = c("0", "", "20", "","40", "","60"))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  facet_grid(cols = vars(StimulationFrequency), rows = vars(StimulationLocationString))+
  geom_segment(inherit.aes = F, data = LightStimSegmentDT[RecLoc=="PaS",], mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(limits = c(0,1), option = "D", guide = guide_colourbar(title = expression(rho), barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 

WaveletPaSOverview <- WTPSDExampleRasterPaS / RhoExampleRasterPaS+ plot_layout(widths = c(3,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))

saveRDS(object = WaveletPaSOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/WaveletPaSOverview.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/WaveletPaSOverview.pdf", plot = WaveletPaSOverview, device = "pdf")

WTPSDExampleRaster <- ggplot(data = WTExampleLarge[Frequency<60,], aes(x = Time, y = Frequency, fill=WTPSD))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10), labels = c("0", "", "20", "","40", "","60"))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  facet_nested(StimulationFrequency ~ RecLoc + StimulationLocationString, nest_line = element_line(colour="black"))+
  geom_segment(inherit.aes = F, data = LightStimSegmentDT, mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(limits=c(0,8),option = "B", guide = guide_colourbar(title = "Power\nDensity", barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = WTPSDExampleRaster, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDWaveletOverview.rds")

RhoExampleRaster <- ggplot(data = WTExampleLarge[Frequency<60,], aes(x = Time, y = Frequency, fill=Rho))+
  geom_raster()+
  scale_y_continuous(name = "Frequency (Hz)", limits = c(0, 65), expand = c(0,0), breaks = seq(0,60,10), labels = c("0", "", "20", "","40", "","60"))+
  scale_x_continuous(name = "Time (s)", expand = c(0,0))+
  facet_nested(StimulationFrequency ~ RecLoc + StimulationLocationString, nest_line = element_line(colour="black"))+
  geom_segment(inherit.aes = F, data = LightStimSegmentDT, mapping = aes(x=TimeMin, xend=TimeMax, y=FrequencyMin, yend=FrequencyMax), colour="deepskyblue2", size=2)+
  scale_fill_viridis_c(limits = c(0,1), option = "D", guide = guide_colourbar(title = expression(rho), barwidth = grid::unit(x = 1, units = "lines")))+
  theme_classic() + theme(panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = RhoExampleRaster, file = "/alzheimer/Daniel_Data/R/Thesis/Data/RhoWaveletOverview.rds")


WaveletOverview <- WTPSDExampleRaster / RhoExampleRaster+ plot_layout(widths = c(3,1)) + plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(size=16))
saveRDS(object = WaveletOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/WaveletOverview.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/WaveletOverview.pdf", plot = WaveletOverview, device = "pdf")
###########

BetaStimDiff <- DifferenceMat(x = cbind(betaMaxStimP, betaMaxStimA), GroupStrings = c("PaS", "MEC"))
TufteBoxPlot(x = BetaStimDiff, Grouping = colnames(BetaStimDiff), Line = 0, HDI = T, ShortWidth = 0.2)

test16 <- cbind((betaMaxStimA[,1])-(betaMaxStimA[,1]-beta0StimA[,1]) * exp(-betaSlopeStimA[,1]*4), (betaMaxPreA[,1])-(betaMaxPreA[,1]-beta0PreA[,1]) * exp(-betaSlopePreA[,1]*4))

TufteBoxPlot(x = test16, Grouping = c("PaS", "MEC"), DataGroup = LogOddsRho[RecLoc%in%c("PaS","MEC")&StimulationLocation=="Anterior"&StimulationFrequency==4&Stimulation=="Stim",RecLoc,],DataY = LogOddsRho[RecLoc%in%c("PaS","MEC")&StimulationLocation=="Anterior"&StimulationFrequency==4&Stimulation=="Stim",RhoLogOdds,], Dodge = T, RangeY = c(-1,5), HDI = T, DataStat = "mean")
LogOddsRho[RecLoc%in%c("PaS","MEC")&StimulationLocation=="Anterior",,]



# ggplot(data = PhaseTable[Count>1&StimulationFrequency==Frequency&StimulationLocation=="Anterior"&RecLoc%in%c("MEC", "PaS"),], mapping = aes(x = StimulationFrequency, y = RhoLogOdds, group=interaction(AnimalID, Date, BurstFrequency, Stimulation)))+
#   #geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
#   #geom_ribbon(mapping = aes(ymin = Phase-Phasesd/sqrt(Count), ymax = Phase+Phasesd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
#   geom_line()+
#   #  geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
#   facet_wrap(RecLoc~Stimulation)+
#   #facet_grid(cols = vars(RecLoc), rows = vars(StimulationFrequency))+
#  # scale_y_continuous(limits = c(0,1), expand = c(0,0))+
#  # scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
# #  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#   theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 
# 

# ggplot(data = PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y = RhoLogOdds, group=interaction(AnimalID, Date, StimulationFrequency, Stimulation), colour=Stimulation))+
#   geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
#   geom_line(alpha=0.1)+
#   geom_line(inherit.aes = F, PhaseTable[StimulationLocation=="Posterior"&grepl(x=RecLoc, pattern ="MEC"),.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(RhoLogOdds=mean(RhoLogOdds, na.rm = T)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  RhoLogOdds, colour=Stimulation), size=1)+
#   facet_wrap(~Stimulation)+
#   scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#   theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 

##### Power analysis ####
FileListPowerFreq <- grep(pattern = "PowerFrequency", x = list.files(path = paste0(ChannelLocations$...2, "OutputFolder"), full.names = T), value = T)

PowerFrequency <- rbindlist(lapply(X = FileListPowerFreq, function(x) {
  readRDS(file = x)
}))
PowerFrequency[,Stimulation:=factor(Stimulation, levels = c("Pre", "Stim", "Post")),]
PowerFrequency[RecLoc=="through", RecLoc:="PaS",][StimLoc=="through", StimLoc:="PaS",]
PowerFrequency[grepl(x=RecLoc, pattern ="MEC"), RecLoc:="MEC",][grepl(x=StimLoc, pattern ="MEC"), StimLoc:="MEC",]
PowerFrequency[,StimulationLocation:="Posterior",][StimLoc %in% c("MS", "LS"),StimulationLocation:="Anterior",][StimLoc %in% c("SI", "AON"), StimulationLocation:="Off",][,RecLoc:=factor(x = RecLoc, levels = c("MS", "LS", "PaS", "MEC", "PrS", "Sub", "DG-mo", "alv", "AON", "SI")),]

PowerFrequency[,PrePSD:=0,][Stimulation=="Pre",PrePSD:=PSD,][,PrePSD:=max(PrePSD, na.rm = T),.(AnimalID, Date, DirName, Protocol, Frequency, Trial, StimLoc, RecLoc)][,PSDRatio:=PSD/PrePSD,]
#PowerFrequency[Frequency==StimulationFrequency,][1:10,.(Stimulation, PSD, PrePSD, PSDRatio, Trial, StimLoc, RecLoc)]
PowerFrequency[,StimulationFrequencyString:=paste0(StimulationFrequency, "Hz"),][
  ,StimulationFrequencyString:=factor(StimulationFrequencyString, levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]

#PowerFrequencyAvg <- PowerFrequency[,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T),PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T)),by=.(AnimalID, Date, StimulationFrequency, Stimulation, Frequency, StimLoc, RecLoc)][,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected), na.rm = T),by=.(AnimalID, StimulationFrequency, Stimulation, Frequency, StimLoc, RecLoc)][,.(PSD=mean(PSD, na.rm = T),PSDsd=sd(PSD, na.rm = T), Power=mean(Power, na.rm = T), Powersd=sd(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PowerFreqCorrectedsd=sd(PowerFreqCorrected, na.rm = T), Count=.N),by=.(StimulationFrequency, Stimulation, Frequency, StimLoc, RecLoc)]
PowerFrequencyAvg <- PowerFrequency[,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T),PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, Date, StimulationFrequency, Stimulation, Frequency, StimulationLocation, RecLoc)][,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, StimulationFrequency, Stimulation, Frequency, StimulationLocation, RecLoc)][,.(PSD=mean(PSD, na.rm = T),PSDsd=sd(PSD, na.rm = T), Power=mean(Power, na.rm = T), Powersd=sd(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PowerFreqCorrectedsd=sd(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T), Count=.N),by=.(StimulationFrequency, Stimulation, Frequency, StimulationLocation, RecLoc)]
PowerFrequencyAvg[,StimulationFrequencyString:=paste0(StimulationFrequency, "Hz"),][
  ,StimulationFrequencyString:=factor(StimulationFrequencyString, levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),]


#PowerFrequency[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected), PSDRatio=mean(PSDRatio)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)]

##### PSD Model 
PSDSession <- PowerFrequency[StimulationLocation!="Off"&StimulationFrequency==Frequency&RecLoc%in%c("PaS","MEC"),
               .(Count=.N, Trial=paste0(DirName, Protocol), PSD=mean(PSD), PSDRatio=mean(PSDRatio), StimulationFrequencyString=factor(paste0(StimulationFrequency, "Hz"), levels = c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz"))),
               by=.(Protocol, AnimalID, DirName, RecLoc, StimLoc, Stimulation, StimulationFrequency, Date, AnimalID, BurstFrequency, StimulationLocation)][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")][,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation")),]
PSDSession[,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation"))]
saveRDS(PSDSession, "/alzheimer/Daniel_Data/R/Thesis/Data/PSDSession.rds")
PSDN <- merge.data.table(x = PSDSession[,.N,by=.(StimulationLocationString, RecLoc, AnimalID)][,.(AnimalN=.N),by=.(StimulationLocationString, RecLoc)], PSDSession[,.N,by=.(StimulationLocationString, RecLoc, DirName)][,.(SessionN=.N),by=.(StimulationLocationString, RecLoc)])
saveRDS(object = PSDN, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDN.rds")

PSDSession[Stimulation=="Pre", ]


stanDataPSD <- list(
  N = PSDSession[,.N,],
  PSD = PSDSession[,PSD,],
  TrialN = max(as.integer(factor(PSDSession[,Trial,]))),
  Trial = as.integer(factor(PSDSession[,Trial,])),
  StimulationInput = model.matrix(object = PSD ~ 0 + Stimulation, data = PSDSession),
  StimulationStim = as.integer(PSDSession[,Stimulation=="Stim",]),
  StimulationStimN = 2,
  StimulationFreq = as.integer(factor(PSDSession[,StimulationFrequency,])),
  StimulationFreqN = max(as.integer(factor(PSDSession[,StimulationFrequency,]))),
  StimulationFrequency = PSDSession[,StimulationFrequency,],
  Stimulation = as.integer(factor(PSDSession[,Stimulation,])),
  StimulationN = max(as.integer(factor(PSDSession[,Stimulation,]))),
  StimLoc = as.integer(factor(PSDSession[,StimulationLocation,])),
  StimLocN = max(as.integer(factor(PSDSession[,StimulationLocation,]))),
  Location = as.integer(factor(PSDSession[,RecLoc,])),
  LocationN = max(as.integer(factor(PSDSession[,RecLoc,]))),
  DateN = max(as.integer(factor(PSDSession[,Date,]))),
  Date = as.integer(factor(PSDSession[,Date,])),
  AnimalN = max(as.integer(factor(PSDSession[,AnimalID,]))),
  Animal = as.integer(factor(PSDSession[,AnimalID,])),
  BurstFrequencyN = max(as.integer(factor(PSDSession[,BurstFrequency,]))),
  BurstFrequency = as.integer(factor(PSDSession[,BurstFrequency,])),
  priorOnly = 1
) 

if(!justLoad) {
  rstudioapi::jobRunScript(path = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/R/RunModelPSDPrior.R", workingDir = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/", importEnv = T, exportEnv = "R_GlobalEnv")
  rstudioapi::jobRunScript(path = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/R/RunModelPSD.R", workingDir = "/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/", importEnv = T, exportEnv = "R_GlobalEnv")
} else {
  fitmodPSD <- readRDS("/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/output/ModelPSD.rds")
  fitmodPSDPrior <- readRDS("/alzheimer/Daniel_Data/R/VRSpikeAnalysis2021/VRSpikeAnalysis2021/output/ModelPSDPrior.rds")
}


fitmodPSDPrior$print()
y_rep <- posterior::as_draws_matrix(fitmodPSDPrior$draws("y_rep"))
np_PSDPrior <- bayesplot::nuts_params(fitmodPSDPrior)
bayesplot::mcmc_trace(x = fitmodPSDPrior$draws(c("beta0", "betaSlope", "betaMax", "sigma")), np = np_PSDPrior)
bayesplot::ppc_intervals(y = stanDataPSD$PSD, yrep = y_rep)
bayesplot::mcmc_trace(x = fitmodPSDPrior$draws(c("beta0", "betaSlope", "betaMax", "sigma")), np = np_PSDPrior)
bayesplot::mcmc_pairs(x =fitmodPSDPrior$draws(c("sigma", "betaSlope", "betaMax", "beta0")), np = np_PSDPrior)


fitmodPSD$print()
summary(fitmodPSD$sampler_diagnostics())
np_PSD <- bayesplot::nuts_params(fitmodPSD)
y_rep <- posterior::as_draws_matrix(fitmodPSD$draws("y_rep"))
bayesplot::mcmc_trace(x =fitmodPSD$draws(c("sigma", "betaSlope", "betaMax")), np = np_PSD)

bayesplot::ppc_intervals(y = stanDataPSD$PSD, yrep = y_rep)
#bayesplot::ppc_intervals(y = stanDataPSD$PSD[stanDataPSD$StimulationFreq==2&stanDataPSD$StimulationInput[,1]==1], yrep = y_rep[,stanDataPSD$StimulationFreq==2&stanDataPSD$StimulationInput[,1]==1])

bayesplot::mcmc_pairs(x =fitmodPSD$draws(c("sigma", "betaSlope", "betaMax")), np = np_PSD)

bayesplot::ppc_dens_overlay_grouped(y =  stanDataPSD$PSD, group = PSDSession[,paste0(RecLoc,Stimulation,StimulationLocationString),], yrep = y_rep[1:200,])

loo_result <-  fitmodPSD$loo(cores = 1, r_eff = T, variables = "y_lik", save_psis = TRUE)
print(loo_result)
plot(loo_result)

#LogOddsRho[loo_result$diagnostics$pareto_k>0.6,,]

bayesplot::ppc_loo_pit_overlay(
  y =  stanDataPSD$PSD,
  yrep = y_rep,
  lw = weights(loo_result$psis_object),
  samples = 500,
  alpha = 0.7
)

bayesplot::ppc_intervals(y =  stanDataPSD$PSD[loo_result$diagnostics$n_eff<3000], yrep = y_rep[,loo_result$diagnostics$n_eff<3000])
PSDSession[loo_result$diagnostics$n_eff<3000,]



#fitmodlogNormal<- fitmodPSD
#fitmodPSD <- fitmodPSDGamma
#  Anterior                    mu = (exp(betaMaxVec) -(exp(betaMaxVec)) .* exp(-StimulationFrequency .* exp(betaSlopeVec))
betaMaxPSDA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 1, CondMatCol = c(1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,1]")),ncol = 2, nrow = 1e4)
betaSlopePSDA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaSlope", ParTransform = "log", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 5, CondMatCol = c(5,5,5)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,5]")),ncol = 2, nrow = 1e4)
betaStimPSDPreA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[1]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 2, CondMatCol = c(2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,2]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDPreA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[1]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 6, CondMatCol = c(6,6,6)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,6]")),ncol = 2, nrow = 1e4)
betaStimPSDStimA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[2]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 3, CondMatCol = c(3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,3]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDStimA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[2]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 7, CondMatCol = c(7,7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,7]")),ncol = 2, nrow = 1e4)
betaStimPSDPostA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[3]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 4, CondMatCol = c(4,4, 4)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,4]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDPostA <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[3]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 8, CondMatCol = c(8,8,8)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[1,8]")),ncol = 2, nrow = 1e4)


## Posterior
betaMaxPSDP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 1, CondMatCol = c(1,1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,1]")),ncol = 2, nrow = 1e4)
betaSlopePSDP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaSlope", ParTransform = "log", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 5, CondMatCol = c(5,5,5)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,5]")),ncol = 2, nrow = 1e4)
betaStimPSDPreP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[1]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 2, CondMatCol = c(2,2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,2]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDPreP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[1]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 6, CondMatCol = c(6,6,6)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,6]")),ncol = 2, nrow = 1e4)
betaStimPSDStimP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[2]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 3, CondMatCol = c(3,3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,3]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDStimP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[2]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 7, CondMatCol = c(7,7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,7]")),ncol = 2, nrow = 1e4)
betaStimPSDPostP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStim[3]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 4, CondMatCol = c(4,4,4)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,4]")),ncol = 2, nrow = 1e4)
betaStimSlopePSDPostP <- PopulationParameterCond(Model = fitmodPSD, ParameterString = "betaStimSlope[3]", MatrixString = "l", CondMatString = c("a", "d", "t"), MatrixCol = 7, CondMatCol = c(8,8,8)) + matrix(data = posterior::as_draws_matrix(x = fitmodPSD$draws("sl[2,8]")),ncol = 2, nrow = 1e4)


BetaMaxOverview <- cbind(exp(betaMaxPSDA+betaStimPSDPreA), exp(betaMaxPSDA+betaStimPSDStimA), exp(betaMaxPSDA+betaStimPSDPostA), exp(betaMaxPSDP+betaStimPSDPreP), exp(betaMaxPSDP+betaStimPSDStimP), exp(betaMaxPSDP+betaStimPSDPostP))
colnames(BetaMaxOverview) <- paste(rep(c("PaS", "MEC"), times = 6),rep(c("MS", "fibre"), each = 6), rep(rep(c("Pre", "Stim", "Post"), each=2), times=2), sep = "-")
saveRDS(BetaMaxOverview, "/alzheimer/Daniel_Data/R/Thesis/Data/BetaMaxOverview.rds")
saveRDS(HDIcalc(x = BetaMaxOverview), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaMaxSummary.rds")
saveRDS(HDIcalc(DifferenceMat(BetaMaxOverview, GroupStrings = colnames(BetaMaxOverview))), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaMaxDiff.rds")

BetaSlopeOverview <- cbind(exp(betaSlopePSDA + betaStimSlopePSDPreA), exp(betaSlopePSDA + betaStimSlopePSDStimA), exp(betaSlopePSDA + betaStimSlopePSDPostA), exp(betaSlopePSDP + betaStimSlopePSDPreP), exp(betaSlopePSDP + betaStimSlopePSDStimP), exp(betaSlopePSDP + betaStimSlopePSDPostP))
colnames(BetaSlopeOverview) <- paste(rep(c("PaS", "MEC"), times = 6),rep(c("MS", "fibre"), each = 6), rep(rep(c("Pre", "Stim", "Post"), each=2), times=2), sep = "-")
saveRDS(BetaSlopeOverview, "/alzheimer/Daniel_Data/R/Thesis/Data/BetaSlopeOverview.rds")
saveRDS(HDIcalc(x = BetaSlopeOverview), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaSlopeSummary.rds")
saveRDS(HDIcalc(DifferenceMat(BetaSlopeOverview, GroupStrings = colnames(BetaSlopeOverview))), "/alzheimer/Daniel_Data/R/Thesis/Data/HDIBetaSlopeDiff.rds")


summary(exp(betaMaxPSDA + betaStimPSDPreA) - exp(betaMaxPSDA + betaStimPSDPreA)*exp(-16*exp(betaSlopePSDA + betaStimSlopePSDPreA)))
summary(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-16*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1])))
summary(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-16*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1])))

##### PSD gamma model
PSDStim <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDStimA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDStimP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDStimP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Stim", StimulationFrequency=x),]
}))

PSDPre <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDPreA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDPreA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDPreP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDPreP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Pre", StimulationFrequency=x),]
}))

PSDPost <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxPSDA[,1] + betaStimPSDPostA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDPostA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDPostA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDA[,2] + betaStimPSDPostA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDPostA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDPostA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,1] + betaStimPSDPostP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDPostP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDPostP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxPSDP[,2] + betaStimPSDPostP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDPostP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDPostP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Post", StimulationFrequency=x),]
}))

PSDRatio <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(1/((exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDPreA[,1])))/(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDPreA[,2])))/(exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDStimA[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDPreP[,1])))/(exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDStimP[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDPreP[,2])))/(exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDStimP[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(StimulationFrequency=x),]
  }))

PSDDiff <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(1/((exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDPreA[,1])))-(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDPreA[,2])))-(exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2]) - exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2])*exp(-x*exp(betaSlopePSDA[,2] + betaStimSlopePSDStimA[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDPreP[,1])))-(exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1]) - exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1])*exp(-x*exp(betaSlopePSDP[,1] + betaStimSlopePSDStimP[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDPreP[,2])))-(exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2]) - exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2])*exp(-x*exp(betaSlopePSDP[,2] + betaStimSlopePSDStimP[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(StimulationFrequency=x),]
}))


#testtmp <- rbindlist(lapply(X = 1:32, function(x){HDIcalc(1/(exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDPreA[,1])))/(exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1]) - exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])*exp(-x*exp(betaSlopePSDA[,1] + betaStimSlopePSDStimA[,1]))), TestValue = 1)}))
PSDDiff[,`:=`(RecLoc=factor(RecLoc, levels=c("PaS", "MEC"))),][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")]
PSDDiff[,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation"))]


PSDRatio[,`:=`(RecLoc=factor(RecLoc, levels=c("PaS", "MEC"))),][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")]
PSDRatio[,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation"))]

saveRDS(object = PSDRatio, "/alzheimer/Daniel_Data/R/Thesis/Data/PSDRatio.rds")

PSDTable <- rbindlist(l = list(PSDPre,PSDStim,PSDPost))
PSDTable[,`:=`(Stimulation=factor(Stimulation, levels = c("Pre", "Stim", "Post")), RecLoc=factor(RecLoc, levels=c("PaS", "MEC"))),][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")]
PSDTable[,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation"))]


AsymPrePaSA <- exp(betaMaxPSDA[,1] + betaStimPSDPreA[,1])
AsymPreMECA <- exp(betaMaxPSDA[,2] + betaStimPSDPreA[,2])
AsymStimPaSA <- exp(betaMaxPSDA[,1] + betaStimPSDStimA[,1])
AsymStimMECA <- exp(betaMaxPSDA[,2] + betaStimPSDStimA[,2])
AsymPostPaSA <- exp(betaMaxPSDA[,1] + betaStimPSDPostA[,1])
AsymPostMECA <- exp(betaMaxPSDA[,2] + betaStimPSDPostA[,2])

AsymPrePaSP <- exp(betaMaxPSDP[,1] + betaStimPSDPreP[,1])
AsymPreMECP <- exp(betaMaxPSDP[,2] + betaStimPSDPreP[,2])
AsymStimPaSP <- exp(betaMaxPSDP[,1] + betaStimPSDStimP[,1])
AsymStimMECP <- exp(betaMaxPSDP[,2] + betaStimPSDStimP[,2])
AsymPostPaSP <- exp(betaMaxPSDP[,1] + betaStimPSDPostP[,1])
AsymPostMECP <- exp(betaMaxPSDP[,2] + betaStimPSDPostP[,2])

CombinedAsymA <- cbind(AsymPrePaSA, AsymPreMECA, AsymStimPaSA, AsymStimMECA, AsymPostPaSA, AsymPostMECA)
colnames(CombinedAsymA) <- c("PaS Pre", "MEC Pre", "PaS Stim", "MEC Stim", "PaS Post", "MEC Post")
CombinedAsymP <- cbind(AsymPrePaSP, AsymPreMECP, AsymStimPaSP, AsymStimMECP, AsymPostPaSP, AsymPostMECP)
colnames(CombinedAsymA) <- c("PaS Pre", "MEC Pre", "PaS Stim", "MEC Stim", "PaS Post", "MEC Post")

CombinedAsym <- cbind(CombinedAsymA[,c(3,4)], CombinedAsymP[,c(3,4)])
colnames(CombinedAsym) <- c("MS-PaS", "MS-MEC", "Fibre-PaS", "Firbe-MEC")

TufteBoxPlot(CombinedAsymA, Grouping = colnames(CombinedAsymA), colour = rep(c("grey30", "deepskyblue2","grey50"), each=2), HDI = T)
PaSAsym <- TufteBoxPlot(LabelY = "Maximum PSD",CombinedAsymA[,c(1,3,5)], RangeY = c(0,3.5), BreaksY = seq(0,3,0.5), Grouping = c("Pre", "Stim", "Post"), colour = c("grey30", "deepskyblue2","grey50"), HDI = T) + ggtitle("PaS") + theme(plot.title = element_text(hjust = 0.5)) + annotate("segment", x = 1.25, xend = 3.75, y = 3.5, yend =  3.5)
MECAsym <- TufteBoxPlot(CombinedAsymA[,c(2,4,6)], RangeY = c(0,3.5), BreaksY = seq(0,3,0.5),Grouping = c("Pre", "Stim", "Post"), colour = c("grey30", "deepskyblue2","grey50"), HDI = T, BreaksYn = NULL) + ggtitle("MEC") + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+ annotate("segment", x = 1.25, xend = 3.75, y = 3.5 ,yend =  3.5)
MaximumPSD <- PaSAsym + MECAsym

AsymA <- TufteBoxPlot(LabelY = "Maximum PSD",CombinedAsymA[,c(3,4)], RangeY = c(0,3.5), BreaksY = seq(0,3,0.5), Grouping = c("MS-PaS", "MS-MEC"), HDI = T) + ggtitle("PaS") + theme(plot.title = element_text(hjust = 0.5)) + annotate("segment", x = 1.25, xend = 3.75, y = 3.5, yend =  3.5)
AsymP <- TufteBoxPlot(CombinedAsymA[,c(2,4,6)], RangeY = c(0,3.5), BreaksY = seq(0,3,0.5),Grouping = c("Pre", "Stim", "Post"), colour = c("grey30", "deepskyblue2","grey50"), HDI = T, BreaksYn = NULL) + ggtitle("MEC") + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+ annotate("segment", x = 1.25, xend = 3.75, y = 3.5 ,yend =  3.5)

AsymAP <- TufteBoxPlot(LabelY = "Maximum PSD",CombinedAsym, RangeY = c(0,10), BreaksY = seq(0,3,0.5), Grouping = c("MS-PaS", "MS-MEC", "Fibre-PaS", "Fibre-MEC"), HDI = T) 


FrequencyStimulationPSDModelOverlap <- ggplot(data = PSDTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSD), size=0.5)+
  # geom_line(data=LogOddsRho[StimulationLocation=="Anterior"], mapping = aes(x=StimulationFrequency, y=PSD, group=DirName))+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested(RecLoc ~ StimulationLocationString)+
  scale_y_continuous(name = expression("Power Density"),breaks = seq(0,5,1), limits = c(0,5))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(2, 4, 8, 16, 32), labels = c("2", "4", "8", "16", "32"), trans = "log2", limits = c(2,32))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  # geom_hline(yintercept = 1)+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = FrequencyStimulationPSDModelOverlap, file = "/alzheimer/Daniel_Data/R/Thesis/Data/FrequencyStimulationPSDModelOverlap.rds")


PSDModelFrequencyPlot <- ggplot(data = PSDTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSD), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc + Stimulation, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density", limits = c(0,5), breaks = seq(0,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PSDModelFrequencyPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDModelFrequencyPlot.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDModelFrequencyPlot.pdf", plot = PSDModelFrequencyPlot, device = "pdf")

PSDModelFrequencyPlot + AsymAP + plot_layout(widths = c(2,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))


PSDRatioModelFrequencyPlot <- ggplot(data = PSDRatio, aes(x = StimulationFrequency, y = `50%`))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`), alpha=0.1, colour=NaN)+
 # geom_point(inherit.aes = F, data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSDRatio), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density Ratio (Stim/Pre)", limits = c(0,4), breaks = seq(0,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
#  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PSDRatioModelFrequencyPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDRatioModelFrequencyPlot.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDRatioModelFrequencyPlot.pdf", plot = PSDRatioModelFrequencyPlot, device = "pdf")

PSDDiffModelFrequencyPlot <- ggplot(data = PSDDiff[StimulationFrequency>2&StimulationLocation=="Anterior",], aes(x = StimulationFrequency, y = `50%`))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`), alpha=0.1, colour=NaN)+
  # geom_point(inherit.aes = F, data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSDRatio), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density Ratio (Stim/Pre)")+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  #  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  #  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PSDRatioModelFrequencyPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDRatioModelFrequencyPlot.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDRatioModelFrequencyPlot.pdf", plot = PSDRatioModelFrequencyPlot, device = "pdf")


PSDModelFrequencyPlot + PSDRatioModelFrequencyPlot + plot_layout(widths = c(3,1)) + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))





###### Spike Analyis
SpikeTableTotal <- readRDS("/alzheimer/Daniel_Data/R/Thesis/Data/SpikeTableTotalV1")


ggplot(data = PSDRatio, aes(x = StimulationFrequency, y = `50%`))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`), alpha=0.1, colour=NaN)+
  geom_point(inherit.aes = F, data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSDRatio), size=0.5)+
  geom_smooth(inherit.aes = F, data=PSDSession, mapping = aes(x=StimulationFrequency, y=PSDRatio))+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density Ratio (Stim/Pre)", limits = c(0,4), breaks = seq(0,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  #  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  #  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 


#### test rho logodds


PSDSession[Stimulation=="Pre",PSDRatioPre:=PSD,]
PSDSession[,PSDRatio:=PSD,]
PSDSession[,PSDRatioPre:=mean(PSDRatioPre, na.rm = T),by=.(Trial, BurstFrequency, Count, RecLoc, StimLoc, DirName, StimulationFrequency, Date, BurstFrequency)]
PSDSession[,PSDRatio:=PSDRatio/PSDRatioPre,]

betaMaxRhoLogOddsA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 1, CondMatCol = c(1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,1]")),ncol = 2, nrow = 1e4)
betaSlopeRhoLogOddsA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 5, CondMatCol = c(5,5)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,5]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsPreA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[1]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 2, CondMatCol = c(2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,2]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsPreA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[1]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 6, CondMatCol = c(6,6)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,6]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsStimA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[2]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 3, CondMatCol = c(3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,3]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsStimA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[2]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 7, CondMatCol = c(7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,7]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsPostA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[3]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 4, CondMatCol = c(4,4)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,4]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsPostA <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[3]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 7, CondMatCol = c(7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[1,8]")),ncol = 2, nrow = 1e4)
## Posterior
betaMaxRhoLogOddsP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaMax", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 1, CondMatCol = c(1,1)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,1]")),ncol = 2, nrow = 1e4)
betaSlopeRhoLogOddsP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaSlope", ParTransform = "log", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 5, CondMatCol = c(5,5)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,5]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsPreP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[1]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 2, CondMatCol = c(2,2)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,2]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsPreP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[1]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 6, CondMatCol = c(6,6)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,6]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsStimP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[2]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 3, CondMatCol = c(3,3)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,3]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsStimP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[2]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 7, CondMatCol = c(7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,7]")),ncol = 2, nrow = 1e4)
betaStimRhoLogOddsPostP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStim[3]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 4, CondMatCol = c(4,4)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,4]")),ncol = 2, nrow = 1e4)
betaStimSlopeRhoLogOddsPostP <- PopulationParameterCond(Model = fitmodLogOdds, ParameterString = "betaStimSlope[3]", MatrixString = "l", CondMatString = c("a", "d"), MatrixCol = 7, CondMatCol = c(7,7)) + matrix(data = posterior::as_draws_matrix(x = fitmodLogOdds$draws("sl[2,8]")),ncol = 2, nrow = 1e4)
##### RhoLogOdds gamma model
RhoLogOddsStim <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsStimA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsStimA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsStimP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsStimP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Stim", StimulationFrequency=x),]
}))

RhoLogOddsPre <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsPreA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsPreA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsPreP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsPreP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Pre", StimulationFrequency=x),]
}))

RhoLogOddsPost <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPostA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPostA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsPostA[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPostA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPostA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsPostA[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPostP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPostP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsPostP[,1]))))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPostP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPostP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsPostP[,2]))))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(Stimulation="Post", StimulationFrequency=x),]
}))

RhoLogOddsRatio <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsPreA[,1])))/(exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsStimA[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsPreA[,2])))/(exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsStimA[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsPreP[,1])))/(exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsStimP[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsPreP[,2])))/(exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsStimP[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(StimulationFrequency=x),]
}))

RhoLogOddsDiff <- rbindlist(lapply(1:32, function(x) {
  tmp <- rbindlist(l = list(data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsPreA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsPreA[,1])))-(exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1]) - exp(betaMaxRhoLogOddsA[,1] + betaStimRhoLogOddsStimA[,1])*exp(-x*exp(betaSlopeRhoLogOddsA[,1] + betaStimSlopeRhoLogOddsStimA[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsPreA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsPreA[,2])))-(exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2]) - exp(betaMaxRhoLogOddsA[,2] + betaStimRhoLogOddsStimA[,2])*exp(-x*exp(betaSlopeRhoLogOddsA[,2] + betaStimSlopeRhoLogOddsStimA[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Anterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsPreP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsPreP[,1])))-(exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1]) - exp(betaMaxRhoLogOddsP[,1] + betaStimRhoLogOddsStimP[,1])*exp(-x*exp(betaSlopeRhoLogOddsP[,1] + betaStimSlopeRhoLogOddsStimP[,1])))), TestValue = 1))[,`:=`(RecLoc="PaS", StimulationLocation="Posterior"),],
                            data.table(HDIcalc(1/((exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsPreP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsPreP[,2])))-(exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2]) - exp(betaMaxRhoLogOddsP[,2] + betaStimRhoLogOddsStimP[,2])*exp(-x*exp(betaSlopeRhoLogOddsP[,2] + betaStimSlopeRhoLogOddsStimP[,2])))), TestValue = 1))[,`:=`(RecLoc="MEC", StimulationLocation="Posterior"),]))
  tmp[,`:=`(StimulationFrequency=x),]
}))



RhoLogTable <- rbindlist(l = list(RhoLogOddsPre,RhoLogOddsStim,RhoLogOddsPost))
RhoLogTable[,`:=`(Stimulation=factor(Stimulation, levels = c("Pre", "Stim", "Post")), RecLoc=factor(RecLoc, levels=c("PaS", "MEC"))),][StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")]
RhoLogTable[,StimulationLocationString:=factor(StimulationLocationString, levels=c("MS Stimulation", "Fibre Stimulation")),]

RhoLogOddsModelFrequencyPlot <- ggplot(data = RhoLogTable, aes(x = StimulationFrequency, y = `50%`, group=Stimulation, colour=Stimulation))+
  geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_point(data=LogOddsRho, mapping = aes(x=StimulationFrequency, y=RhoLogOdds), size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_nested( StimulationLocationString~ RecLoc + Stimulation, nest_line = element_line(colour = "black"))+
 # scale_y_continuous(name = "Power Density", limits = c(0,5), breaks = seq(0,5,1))+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(1, 2, 4, 8, 16, 32), labels = c("1", "2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  geom_line()+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), panel.spacing.y=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 270), plot.title = element_text(hjust = 0.5)) 

####
##### Power density average over frequency #####
PowerFrequencyAvg[StimulationLocation%in%c("Anterior", "Posterior"),StimulationLocationString:=ifelse(StimulationLocation=="Posterior", "Fibre Stimulation", "MS Stimulation")][
  ,StimulationLocationString:=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation")),]

PowerDensitySpectralAnteriorPlot <- ggplot(data = PowerFrequencyAvg[grepl(x=StimulationLocation, pattern ="Anterior")&RecLoc%in%c("PaS", "MEC"),], mapping = aes(x = Frequency, y = PSD, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin=PSD-PSDsd/sqrt(Count), ymax=PSD+PSDsd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  # geom_line(inherit.aes = F, PowerFrequencyAvg[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),], mapping = aes(x = Frequency, y =  PSD, colour=Stimulation), size=1)+
  facet_grid(StimulationFrequencyString ~ RecLoc)+
  scale_x_continuous(trans = "log2", breaks = 2^(0:6), limits = c(1,80), name = "Frequency (Hz)")+
  scale_y_continuous(name = "Power Density", limits = c(0,3.5))+
  ggtitle("Fibre Stimulation")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PowerDensitySpectralPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PowerDensitySpectralPlot.rds")

PowerDensitySpectralPosteriorPlot <- ggplot(data = PowerFrequencyAvg[grepl(x=StimulationLocation, pattern ="Posterior")&RecLoc%in%c("PaS", "MEC"),], mapping = aes(x = Frequency, y = PSD, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin=PSD-PSDsd/sqrt(Count), ymax=PSD+PSDsd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  # geom_line(inherit.aes = F, PowerFrequencyAvg[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),], mapping = aes(x = Frequency, y =  PSD, colour=Stimulation), size=1)+
  facet_grid(StimulationFrequencyString ~ RecLoc)+
  scale_x_continuous(trans = "log2", breaks = 2^(0:6), limits = c(1,80), name = "Frequency (Hz)")+
  scale_y_continuous(name = "Power Density", limits = c(0,3.5))+
  ggtitle("MS Stimulation")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PowerDensitySpectralPlot, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PowerDensitySpectralPlot.rds")


PowerDensitySpectralPlotOverview <- PowerDensitySpectralAnteriorPlot + PowerDensitySpectralPosteriorPlot + plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size=16))
saveRDS(object = PowerDensitySpectralPlotOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PowerDensitySpectralPlotOverview.rds")

PowerDensitySpectralPlotNestedOverview <- ggplot(data = PowerFrequencyAvg[StimulationLocation%in%c("Anterior", "Posterior")&RecLoc%in%c("PaS", "MEC"),], mapping = aes(x = Frequency, y = PSD, group=interaction(StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_ribbon(mapping = aes(ymin=PSD-PSDsd/sqrt(Count), ymax=PSD+PSDsd/sqrt(Count), group=Stimulation, fill=Stimulation), colour=NaN, alpha=0.2)+
  geom_line()+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"), guide=FALSE)+
  # geom_line(inherit.aes = F, PowerFrequencyAvg[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),], mapping = aes(x = Frequency, y =  PSD, colour=Stimulation), size=1)+
  facet_nested(StimulationFrequencyString ~ RecLoc+StimulationLocationString, nest_line = element_line(colour="black"))+
  scale_x_continuous(trans = "log2", breaks = 2^(0:6), limits = c(1,80), name = "Frequency (Hz)")+
  scale_y_continuous(name = "Power Density", limits = c(0,3.5))+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"), guide = guide_legend(keywidth = 0.5, title = ""))+
  theme_classic() + theme(panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PowerDensitySpectralPlotNestedOverview, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PowerDensitySpectralPlotNestedOverview.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PowerDensitySpectralPlotNestedOverview.pdf", plot = PowerDensitySpectralPlotNestedOverview, device = "pdf",width = 7, height = 9)


##### Power density line plots for each frequency ####
LightDF <- data.table(RecLoc = rep(c("PaS", "MEC"), times=10), StimulationLocationString=rep(c("MS Stimulation", "Fibre Stimulation"), each=10), PSD=rep_len(x = PSDSession[,mean(PSD),], length.out = 20), StimulationFrequencyString=rep(rep(c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz"), each=2), times=2),Stimulation="Stim", AnimalID="", Date="", Protocol="")
LightDF[,`:=`(RecLoc=factor(RecLoc, levels = c("PaS", "MEC")),
              Stimulation=factor(Stimulation, levels = c("Pre", "Stim", "Post"), ordered = T),
              StimulationFrequencyString=factor(StimulationFrequencyString, levels=c("2Hz", "4Hz", "8Hz", "16Hz", "32Hz")),
              StimulationLocationString=factor(StimulationLocationString, levels = c("MS Stimulation", "Fibre Stimulation"))),]

PSDStimLinelog10PlotFacet <- ggplot(data = PSDSession, mapping = aes(x = Stimulation, y = PSD, group=interaction(AnimalID, Date, StimulationLocationString, Protocol)))+
  #geom_vline(xintercept = "Stim", colour="deepskyblue2", size=3)+
  geom_tile(data = LightDF, aes(width=0.5, height=Inf, x = Stimulation, y = PSD), fill="deepskyblue2")+
  geom_line(alpha=0.5, size=0.5)+
  # geom_line(inherit.aes = F, PowerFrequency[StimulationLocation=="Anterior"&RecLoc%in%c("MEC","PaS")&StimulationFrequency==Frequency,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequencyString, Stimulation, RecLoc)][,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequencyString, Stimulation, RecLoc)], mapping = aes(x = Stimulation, y = PSDRatio, group=AnimalID), size=0.5)+
  facet_nested(RecLoc ~ StimulationLocationString+StimulationFrequencyString, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density", expand = c(0,0), trans="log10")+
  scale_x_discrete(name="",drop = FALSE, guide = guide_axis(n.dodge=2)) +
  #  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PSDStimLinelog10PlotFacet, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDStimLinelog10PlotFacet.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDStimLinelog10PlotFacet.pdf", plot = PSDStimLinelog10PlotFacet, device = "pdf",width = 7, height = 9)

PSDStimLinePlotFacet <- ggplot(data = PSDSession, mapping = aes(x = Stimulation, y = PSD, group=interaction(AnimalID, Date, StimulationLocationString, Protocol)))+
  #geom_vline(xintercept = "Stim", colour="deepskyblue2", size=3)+
  geom_tile(data = LightDF, aes(width=0.5, height=Inf, x = Stimulation, y = PSD), fill="deepskyblue2")+
  geom_line(alpha=0.5, size=0.5)+
  # geom_line(inherit.aes = F, PowerFrequency[StimulationLocation=="Anterior"&RecLoc%in%c("MEC","PaS")&StimulationFrequency==Frequency,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, Date, Frequency, StimulationFrequencyString, Stimulation, RecLoc)][,.(PSD=mean(PSD, na.rm = T), Power=mean(Power, na.rm = T), PowerFreqCorrected=mean(PowerFreqCorrected, na.rm = T), PSDRatio=mean(PSDRatio, na.rm = T)),by=.(AnimalID, Frequency, StimulationFrequencyString, Stimulation, RecLoc)], mapping = aes(x = Stimulation, y = PSDRatio, group=AnimalID), size=0.5)+
  facet_nested(RecLoc ~ StimulationLocationString+StimulationFrequencyString, nest_line = element_line(colour = "black"))+
  scale_y_continuous(name = "Power Density", expand = c(0,0))+
  scale_x_discrete(name="",drop = FALSE, guide = guide_axis(n.dodge=2)) +
  #  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(1.8, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 
saveRDS(object = PSDStimLinePlotFacet, file = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDStimLinePlotFacet.rds")
ggsave(filename = "/alzheimer/Daniel_Data/R/Thesis/Data/PSDStimLinePlotFacet.pdf", plot = PSDStimLinePlotFacet, device = "pdf",width = 7, height = 9)

###########


### Single lines
ggplot(data = PowerFrequency[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected), PSDRatio=mean(PSDRatio)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y = PSD, group=interaction(AnimalID, Date, StimulationFrequency, Stimulation), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_line(alpha=0.1)+
  geom_line(inherit.aes = F, PowerFrequency[grepl(x=StimLoc, pattern ="MS")&grepl(x=RecLoc, pattern ="PaS"),.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected), PSDRatio=mean(PSDRatio)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected), PSDRatio=mean(PSDRatio)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)][,.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected),PSDRatio=mean(PSDRatio)),by=.(Frequency, StimulationFrequency, Stimulation)], mapping = aes(x = Frequency, y =  PSD, colour=Stimulation), size=1)+
  facet_wrap(~StimulationFrequency)+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank())


ggplot(data = PowerFrequencyAvg[StimulationLocation=="Anterior"&RecLoc%in%c("MEC","PaS")&StimulationFrequency==Frequency], aes(x = StimulationFrequency, y = PSD, group=Stimulation, colour=Stimulation))+
 # geom_ribbon(mapping = aes(x = StimulationFrequency, ymin=`2.5%`, ymax=`97.5%`, fill=Stimulation), alpha=0.1, colour=NaN)+
  geom_line(size=0.5)+
  # facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
  facet_grid(rows = vars(RecLoc), cols = vars(Stimulation))+
 # scale_y_continuous(name = expression(rho*" Log Odds Ratio"), trans = "log")+
  scale_x_continuous(name = "Stimulation Frequency (Hz)", breaks = c(2, 4, 8, 16, 32), labels = c("2", "4", "8", "16", "32"), trans = "log2")+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  scale_fill_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
#  geom_line()+
  ggtitle("MS Stimulation")+
  # geom_hline(yintercept = 1)+
  theme_classic() + theme(legend.position="none", panel.spacing.x=unit(1, "lines"), axis.line = element_blank(), strip.background = element_blank(), strip.text.y = element_text(angle = 0), plot.title = element_text(hjust = 0.5)) 



###
PowerFrequency[StimLoc=="MS"&grepl(x=RecLoc, pattern ="PaS")&StimulationFrequency==Frequency,.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected)),by=.(AnimalID, Date, Frequency, StimulationFrequency, Stimulation)][,.(PSD=mean(PSD), Power=mean(Power), PowerFreqCorrected=mean(PowerFreqCorrected)),by=.(AnimalID, Frequency, StimulationFrequency, Stimulation)]

WTList <- lapply(seq_len(length(StimPaSrecPaS_WTs)), function(x) {sqrt(StimPaSrecPaS_WTs[[x]]$PowerWT)})
WTList
Frequency <- seq(0.25,80,0.25)
StimulationFrequency <- sapply(seq_len(length(StimPaSrecPaS_WTs)),FUN =  function(x){StimPaSrecPaS_WTs[[x]]$StimulationFrequency})
BurstFrequency <- sapply(seq_len(length(StimPaSrecPaS_WTs)),FUN =  function(x){StimPaSrecPaS_WTs[[x]]$BurstFreq})
ProtName <- "test"
PowerTable <- rbindlist(lapply(X = WTList, FUN = function(x) {
  rownames(x = x) <- 1:dim(x)[1]/1e3
  colnames(x = x) <- Frequency
  outTable <- data.table(reshape2::melt(x))
  outTable[,`:=`(StimulationFrequency=StimulationFrequency, BurstFreq=BurstFrequency,Count=dim(x)[3], ProtName=ProtName),]
  setnames(x = outTable, old= c("Var1", "Var2", "Var3","value"), new = c("Time", "Frequency", "Trial","Power"))
  outTable
}))

PowerTable[,NormalisedPower:=(Power*Frequency)^2,]
PowerTable[]
object.size(PowerTable)/1e6/1e3
PowerFreqCorrected <- lapply(X = WTList, FUN = function(x){
  # PhaseRho <- InVivoR::PhaseListAnalysis(x = atan2(Im(ERPWT$Raw), Re(ERPWT$Raw)))
   #PowerWT <- abs(x)^2
   PowerWTFreqCorr <- (x*array(data = rep(seq(0.25,80,0.25), each=dim(x)[1]), dim = dim(x)))^2
  # WTPSD <- array(data = 0, dim = dim(PowerWTFreqCorr))
  # for(sliceNr in 1:dim(PowerWT)[3]) {
  #   WTPSD[,,sliceNr] <- PowerWTFreqCorr[,,sliceNr]/sum(PowerWTFreqCorr[,,sliceNr])*length(PowerWTFreqCorr[,,sliceNr])
  # }
})
WTPSD <- lapply(X = PowerFreqCorrected, function(x) {
  WTPSD <- array(data = 0, dim = dim(x))
  for(sliceNr in 1:dim(x)[3]) {
   WTPSD[,,sliceNr] <- x[,,sliceNr]/sum(x[,,sliceNr])*length(x[,,sliceNr])
  }
  WTPSD
})



ggplot(data = testBind, mapping = aes(x = Frequency, y = PSD, group=interaction(Trial, Stimulation, BurstFrequency), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_line(alpha=0.1)+
  geom_line(inherit.aes = F, testBind[,.(PSD=mean(PSD)),by=.(Stimulation, StimulationFrequency, Frequency)], mapping = aes(x = Frequency, y = PSD, colour=Stimulation), size=1)+
  facet_wrap(~StimulationFrequency)+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 


ggplot(data = StimMSrecPaS_PowerFrequency, mapping = aes(x = Frequency, y = PSD, group=interaction(Trial, Stimulation, BurstFrequency), colour=Stimulation))+
  geom_vline(data = StimulationLine, aes(xintercept = Frequency), linetype="dashed")+
  geom_line(alpha=0.1)+
  geom_line(inherit.aes = F, StimMSrecPaS_PowerFrequency[,.(PSD=mean(PSD)),by=.(Stimulation, StimulationFrequency, Frequency)], mapping = aes(x = Frequency, y = PSD, colour=Stimulation), size=1)+
  facet_wrap(~StimulationFrequency)+
  scale_colour_manual(values = c("Pre"="black", "Stim"="deepskyblue2", "Post"="grey50"))+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 


ggplot(data = testBind[StimulationFrequency==Frequency,], mapping = aes(x = Stimulation, y = PSD, group=interaction(Trial, Protocol)))+
  geom_line(alpha=0.1)+
  geom_line(inherit.aes = F, testBind[,.(PSD=mean(PSD)),by=.(Stimulation, StimulationFrequency, Frequency)][StimulationFrequency==Frequency,], mapping = aes(x = Stimulation, y=PSD, group=StimulationFrequency), size=1)+
  facet_wrap(~StimulationFrequency)+
  theme_classic() + theme(legend.position="none", panel.spacing.y=unit(0.1, "lines")) + theme(axis.line = element_blank(), strip.background = element_blank()) 


image(x = InVivoR::PowerMat(WTPSD[[4]]), useRaster = T) 
PowerAvg <- InVivoR::PowerMat(WTPSD[[4]])
rownames(x = PowerAvg) <- 1:dim(PowerAvg)[1]/1e3
colnames(x = PowerAvg) <- Frequency
outTable <- data.table(reshape2::melt(PowerAvg))
outTable[,`:=`(StimulationFrequency=36, BurstFreq=0),]
setnames(x = outTable, old= c("Var1", "Var2","value"), new = c("Time", "Frequency","PSD"))

library(ggplot2)
ggplot(data = outTable, mapping = aes(x = Time, y = Frequency, fill=PSD))+
  geom_raster()+
  scale_fill_viridis_c(begin = 0)



str(WTPSD)
image(InVivoR::PowerMat(PowerFreqCorrected[[4]]), useRaster = T)
StimMSrecPaS_StimERP$ProtocolNames
getProtocols <- grep(pattern = "^2.0_pulse|^4.0_pulse|^8.0_pulse|^16.0_pulse|^32.0_pulse",x = StimMSrecPaS_StimERP$ProtocolNames)

WTSD <- (InVivoR::PowerMat(x = (sqrt(StimPaSrecPaS_WTs[[4]]$PowerWT)*array(data = rep(seq(0.25,80,0.25), each=dim(StimPaSrecPaS_WTs[[4]]$PowerWT)[1]), dim = dim(StimPaSrecPaS_WTs[[4]]$PowerWT)))^2))
StimMSrecPaS_WTs[[4]]$StimulationFrequency
StimMSrecPaS_WTs[[4]]$N
image(WTSD, useRaster = T)

StimPaS <- readRDS(file = "/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/StimPaSrecPaS_StimERP")
StimMS <- readRDS(file = "/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/StimMSrecPaS_StimERP")

StimMSrecPaS_StimERP <- readRDS(file = "/alzheimer/Daniel_Data/DSC008159/DSC008159_190619_150144/OutputFolder/StimMSrecPaS_StimERP.rds")
StimPaSrecPaS_WTs <- readRDS(file = "/alzheimer/Daniel_Data/DSC008159/DSC008159_190619_150144/OutputFolder/StimMSrecPaS_WTs.rds")
StimMSrecPaS_WTs <- readRDS(file = "/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190409_132939/OutputFolder/StimMSrecPaS_WTs.rds")

image(StimMSrecPaS_WTs[[1]]$PowerWT[,,1], useRaster = T)

WTPSDCollab <- InVivoR::PowerMat(x = StimPaSrecPaS_WTs[[4]]$WTPSD)
image(WTPSDCollab, useRaster = T)
plot(y=colMeans(WTPSDCollab[1:5000,]), x = StimPaSrecPaS_WTs[[1]]$Frequencies, type="l", ylim = c(0,10))
lines(y=colMeans(WTPSDCollab[5001:10000,]), x = StimPaSrecPaS_WTs[[1]]$Frequencies, type="l", col="deepskyblue2")
lines(y=colMeans(WTPSDCollab[10001:15000,]), x = StimPaSrecPaS_WTs[[1]]$Frequencies, type="l")
abline(v = StimPaSrecPaS_WTs[[3]]$StimulationFrequency)

WTPSDCollab <- InVivoR::PowerMat(x = StimMSrecPaS_WTs[[4]]$PowerWTFreqCorr)
image(WTPSDCollab, useRaster = T)
plot(y=colMeans(WTPSDCollab[1:5000,]), x = StimMSrecPaS_WTs[[1]]$Frequencies, type="l", ylim = c(0,7e5))
lines(y=colMeans(WTPSDCollab[5001:10000,]), x = StimMSrecPaS_WTs[[1]]$Frequencies, type="l", col="deepskyblue2")
lines(y=colMeans(WTPSDCollab[10001:15000,]), x = StimMSrecPaS_WTs[[1]]$Frequencies, type="l")
abline(v = StimMSrecPaS_WTs[[4]]$StimulationFrequency)


ProtocolNames <- unlist(lapply(StimMSrecPaS_WTs, function(x) {
  x$ProtName
}))

PowerTable <- rbindlist(lapply(X = StimMSrecPaS_WTs, FUN = function(x) {
  rownames(x = x$PowerWT) <- 1:dim(x$PowerWT)[1]/1e3
  colnames(x = x$PowerWT) <- x$Frequencies
  outTable <- data.table(reshape2::melt(x$PowerWT))
  outTable[,`:=`(PowerWTFreqCorr=x$PowerWTFreqCorr, WTPSD=x$WTPSD, StimulationFrequency=x$StimulationFrequency, BurstFreq=x$BurstFreq,Count=x$N, ProtName=x$ProtName),]
  setnames(x = outTable, old= c("Var1", "Var2", "Var3","value"), new = c("Time", "Frequency", "Trial","Power"))
  outTable
}))
PowerTable[Time%between%c(0,5),Stimulation:="Pre",][Time%between%c(10,15),Stimulation:="Post",][Time%between%c(5,10),Stimulation:="Stim",]

PowerFrequency <- PowerTable[,.(WTPSD=mean(WTPSD),Power=mean(Power),PowerWTFreqCorr=mean(PowerWTFreqCorr)),by=.(ProtName, Stimulation, Frequency, StimulationFrequency, BurstFreq, Trial)]

WTSPSDavgTable <- PowerTable[,.(WTPSD=mean(WTPSD),Power=mean(Power),PowerWTFreqCorr=mean(PowerWTFreqCorr)),by=.(ProtName, Time, Frequency)]

ggplot(data = PowerFrequency[StimulationFrequency==32,], mapping = aes(x = Frequency, y=WTPSD))+
  geom_line()+
  geom_vline(xintercept = 16)+
 # geom_vline(xintercept = 4)+
  facet_wrap(facets = ~Stimulation)

ggplot(data = PowerFrequency[StimulationFrequency==Frequency&Stimulation%in%c("Pre","Stim"),], mapping = aes(x = StimulationFrequency, y=PowerRatio, colour=Stimulation))+
  geom_point()+
  facet_wrap(~Stimulation)

PowerFrequency[Stimulation=="Pre",PrePower:=Power,]
PowerFrequency[,PrePower:=max(PrePower, na.rm = T),by=.(ProtName, Frequency)]
PowerFrequency[,PowerRatio:=Power/PrePower,]

StimMSrecPaS_WTs[[4]]$StimulationFrequency

plot(StimMSrecMS_StimERP$ERP[[1]][1,], type="l")


testWT <- InVivoR::WTbatch(ERPMat = StimPaS$ERP[[2]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F, PhaseAnalysis = F)
WTz <- InVivoR::PowerMat(x = abs(testWT$Raw), ZScore = F)
dim(WTz)
PowerWTFreqCorr <- testWT1$Raw^2*array(data = rep(seq(0.25,80,0.25), each=dim(testWT1$Raw)[1]), dim = dim(testWT1$Raw))
WTPowerZ <- InVivoR::PowerMat(x = abs(PowerWTFreqCorr), ZScore = F)
image(x = WTPowerZ, useRaster = T)

WTz <- t(t(WTz)*testWT$Frequencies)
u <- matrix(data = colMeans(WTz[1:5000,]), nrow = dim(WTz)[2], ncol = dim(WTz)[1])
sd <- matrix(data = apply(X = WTz[1:5000,], MARGIN = 2, FUN = sd), nrow = dim(WTz)[2], ncol = dim(WTz)[1])
dim(sd)
ZscoreWT <- t((t(WTz)-u)/sd)

StimMS$ProtocolNames
testWT1 <- InVivoR::WTbatch(ERPMat = StimMS$ERP[[2]], frequencies = seq(0.25,80,0.25), SamplingRate = 1e3, CORES = 10, compression = F, PhaseAnalysis = F)
WTz1 <- InVivoR::PowerMat(x = abs(testWT1$Raw)^2, ZScore = F)
WTz1 <- t(t(WTz1)*testWT1$Frequencies)
u1 <- matrix(data = colMeans(abs(testWT1$Raw[1:5000,,1])^2), nrow = dim(WTz1)[2], ncol = dim(WTz1)[1])
sd1 <- matrix(data = apply(X = abs(testWT1$Raw[1:5000,,1])^2, MARGIN = 2, FUN = sd), nrow = dim(WTz1)[2], ncol = dim(WTz1)[1])
dim(sd)
ZscoreWT1 <- t((t(abs(testWT1$Raw[,,1])^2)-u1)/sd1)

plot(x =testWT$Frequencies, y = colMeans(ZscoreWT1[5001:10000,]), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = colMeans(ZscoreWT1[1:5000,]), type="l")
lines(x =testWT$Frequencies, y = colMeans(ZscoreWT1[10000:15000,]), type="l")

lines(x =testWT$Frequencies, y = colMeans(ZscoreWT[5001:10000,]), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = colMeans(ZscoreWT[1:5000,]), type="l")
lines(x =testWT$Frequencies, y = colMeans(ZscoreWT[10000:15000,]), type="l")
image(x = ZscoreWT1, useRaster = T)
image(x = ZscoreWT1, useRaster = T)

WTz1 <- t(t(WTz1)*testWT1$Frequencies)
u1 <- matrix(data = colMeans(WTPowerZ[1:5000,]^2), nrow = dim(WTPowerZ)[2], ncol = dim(WTPowerZ)[1])
sd1 <- matrix(data = apply(X = WTPowerZ[1:5000,]^2, MARGIN = 2, FUN = sd), nrow = dim(WTPowerZ)[2], ncol = dim(WTPowerZ)[1])

dim(sd)

ZscoreWT1 <- t((t(WTPowerZ^2)-u1)/sd1)
image(ZscoreWT1, useRaster = T)
image(WTPowerZ/sum(WTPowerZ)*length(WTPowerZ), useRaster = T)
plot(x =testWT$Frequencies, y = colMeans(WTPowerZ[5001:10000,])/sum(colMeans(WTPowerZ[5001:10000,]))*length(colMeans(WTPowerZ[5001:10000,])), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = colMeans(WTPowerZ[1:5000,])/sum(colMeans(WTPowerZ[1:5000,]))*length(colMeans(WTPowerZ[5001:10000,])), type="l")
lines(x =testWT$Frequencies, y = colMeans(WTPowerZ[10000:15000,])/sum((WTPowerZ[10000:15000,]))*length(WTPowerZ[5001:10000,]), type="l")

WTPSD <- WTPowerZ/sum(WTPowerZ)*length(WTPowerZ)
plot(x =testWT$Frequencies, y = colMeans(WTPSD[5001:10000,]), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = colMeans(WTPSD[1:5000,]), type="l")
lines(x =testWT$Frequencies, y = colMeans(WTPSD[10000:15000,]), type="l")


plot(x =testWT$Frequencies, y = 10*log(colMeans(WTPowerZ[5001:10000,])), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = 10*log(colMeans(WTPowerZ[1:5000,])), type="l")
lines(x =testWT$Frequencies, y = 10*log(colMeans(WTPowerZ[10000:15000,])), type="l")


StimMS$ERP[[2]][1,]
testBspec <- bspec::welchPSD(x = ts(data = StimMS$ERP[[2]][1,1:5000], frequency = 1e3), seglength = 2)
testBspec1 <- bspec::welchPSD(x = ts(data = StimMS$ERP[[2]][1,5001:10000], frequency = 1e3), seglength = 2)
testBspec2 <- bspec::welchPSD(x = ts(data = StimMS$ERP[[2]][1,10001:15000], frequency = 1e3), seglength = 2)

plot(x = testBspec1$frequency, y = testBspec1$power, type="l", xlim = c(0,30), col="deepskyblue2")
lines(x = testBspec$frequency, y = testBspec$power, type="l", xlim = c(0,30))
lines(x = testBspec$frequency, y = testBspec2$power, type="l", xlim = c(0,30))


PowerWT <- abs(testWT1$Raw)^2
PowerWTFreqCorr <- PowerWT*array(data = rep(seq(0.25,80,0.25), each=dim(PowerWT)[1]), dim = dim(PowerWT))
WTPSD <- array(data = 0, dim = dim(PowerWTFreqCorr))
for(sliceNr in 1:dim(PowerWT)[3]) {
  WTPSD[,,sliceNr] <- PowerWTFreqCorr[,,sliceNr]/sum(PowerWTFreqCorr[,,sliceNr])*length(PowerWTFreqCorr[,,sliceNr])
}
image(WTPSD[,,2], useRaster = T)

PSDavg <- InVivoR::PowerMat(x = WTPSD)
plot(x =testWT$Frequencies, y = 10*log(colMeans(PSDavg[5001:10000,])), type="l", col="deepskyblue2")
lines(x =testWT$Frequencies, y = 10*log(colMeans(PSDavg[1:5000,])), type="l")
lines(x =testWT$Frequencies, y = 10*log(colMeans(PSDavg[10000:15000,])), type="l")
plot(10*log(1/(1:80)))


StimMSrecPaS_WTs <- readRDS(file = "/alzheimer/Daniel_Data/DSC009444/DSC009444_191125_152005/OutputFolder/StimMSrecPaS_WTs")
StimMSrecPaS_ERP <- readRDS(file = "/alzheimer/Daniel_Data/DSC009444/DSC009444_191125_152005/OutputFolder/WT")
StimMSrecPaS_ERP <- readRDS(file = "/alzheimer/Daniel_Data/DSC009444/DSC009444_191125_152005/OutputFolder/StimPaSrecMS_WTs.rds")

str(StimMSrecPaS_ERP)

image(StimMSrecPaS_ERP$ERP[[1]], useRaster=T)
StimMSrecPaS_ERP$ERP[[1]][1,]
system.time(saveRDS(object = StimMSrecPaS_ERP, file = "slow.rds"))
system.time(saveRDS(object = StimMSrecPaS_ERP, file = "fast.rds", compress = F))



##### !!!!! Spike analysis ####
#SpikeTableTotal <- readRDS("/alzheimer/Daniel_Data/Analysis/Analysis_output/20200511/SpikeTableTotalWithFiringRate.rds")
SpikeTableTotal <- readRDS("Data/SpikeTableTotalV2.rds")
LowOccurenceDropOut <- SpikeTableTotal[ClusterNr>1, .N,by=.(UniqueID)][N<3e3, ][,UniqueID,]  
SpikeTableTotal[UniqueID %in% LowOccurenceDropOut,ClusterNr:=1,] 

for(Unit in SpikeTableTotal[ClusterNr>1,unique(UniqueID),]) {
  testCCF <- InVivoR::SpikeCCF(x = SpikeTableTotal[UniqueID==Unit, SpikeTime,], y = SpikeTableTotal[UniqueID==Unit, SpikeTime,])
  SpikeTableTotal[UniqueID==Unit, Clean20:=testCCF$CCF[testCCF$xAxis==0]<testCCF$RandomBinCount*0.2,]
}

#### isolate ms pv cells ####
MSTestunits <- SpikeTableTotal[Clean20==T&frequency==1&UnitLoc=="MS",.N,by=UniqueID][N>1000,UniqueID,] 

SpikeTableMS1Hz <- SpikeTableTotal[Clean20==T&UniqueID %in% MSTestunits & frequency==1&UnitLoc=="MS",] 
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
SpikeTableTotal[UniqueID %in% PV_Units, unique(UniqueID), by= SessionID][,.N,by=SessionID] 


#### make plot for z-score ####
SinglePulseZScorePlot <- ggplot(data = SpikeZScoreTable[StimRoundPos %between% c(-0.035,0.035)]  , aes(x=StimRoundPos*1e3, y=OrderID, fill=ScaleN))+
  geom_raster()+
  geom_vline(xintercept = 0, lty="dashed", size=0.2)+
  scale_y_continuous(expand=c(0,0), name = "Units", seq(0,80,10))+
  scale_x_continuous(expand=c(0,0), name = "Time to Stimulation (ms)", breaks = seq(-40,40,10))+
  geom_hline(yintercept = ScaleCutOff-0.5, lty="dashed", size=0.2, colour="gray")+
  #scale_fill_viridis_c()+
  scale_fill_gradient2(limits = c(-max(SpikeZScoreTable[,max(ScaleN),]),SpikeZScoreTable[,max(ScaleN),]), guide = guide_colorbar(title = "Z-Score"))+
  theme_bw()





#### Raster plot ####

# RasterDataSet <- SpikeTableTotal[UniqueID%in%PV_Units[13]  &
#                                    PulseSingleNumber <= 100 & 
#                                    hyper_block_frequency==1 & 
#                                    SingleStimDiff %between% c(-0.035,0.035), ,] 
# RasterDataSet <- RasterDataSet[,.(SingleStimDiff=SingleStimDiff*1e3, PulseSingleNumber),]
# RasterDataSet[,`:=`(y=PulseSingleNumber-0.5, yend=PulseSingleNumber+0.5),]

PV_Units
RasterDataSet <- SpikeTableTotal[UniqueID=="516_1 DSC009444_191127_163613"  &
                                   PulseSingleNumber <= 300 & 
                                   hyper_block_frequency==1 & 
                                   SingleStimDiff %between% c(-0.035,0.035), ,] 
RasterDataSet <- RasterDataSet[,.(SingleStimDiff=SingleStimDiff*1e3, PulseSingleNumber),]
RasterDataSet[,`:=`(y=PulseSingleNumber-0.5, yend=PulseSingleNumber+0.5),]

RasterPlot <- ggplot(data = RasterDataSet, aes(x=SingleStimDiff, y=y))+
  geom_rect(data=data.frame(SingleStimDiff=0, y=0), aes(xmin=0, xmax=2, ymin=0,ymax=Inf), alpha=1, fill="deepskyblue2")+
  geom_segment(aes(y=y, yend=yend, x=SingleStimDiff, xend=SingleStimDiff), colour="black")+
  geom_point(shape="|", size=0.4)+
  #facet_wrap(facets = ~ UniqueID)+
  scale_y_continuous(expand=c(0,0), name = "Pulse", limits = c(0,301), breaks = seq(0,100,10))+
  scale_x_continuous(expand=c(0,0), name = "Time to Stimulation (ms)", breaks = seq(-40,40,10))+
  theme_classic()+
  theme(strip.background = element_rect(color=NaN, fill=NaN))

xhist <- cowplot::axis_canvas(RasterPlot, axis = "x") + 
  stat_bin(data = RasterDataSet,
           aes(x = SingleStimDiff),
           breaks = seq(-35,35,1), fill = 'grey50', colour="white")
           
  geom_histogram(data = RasterDataSet,
                 aes(x = SingleStimDiff),
                 binwidth = 1,
                 fill = 'grey50', colour="white")

RasterWithHistogram <- RasterPlot %>%
  cowplot::insert_xaxis_grob(xhist, grid::unit(1, "in"), position = "top") %>%
  cowplot::ggdraw()

SinglePulseZScorePlot + RasterPlot

##### 

setkey(SpikeTableTotal, SpikeTime)
testCCF <- InVivoR::SpikeCCF(x = SpikeTableTotal[SessionID=="DSC009443_191127_135954"&PVPos==T,SpikeTime,], UnitNr = factor(SpikeTableTotal[SessionID=="DSC009443_191127_135954"&PVPos==T,UniqueID,]))
levels(factor(SpikeTableTotal[SessionID=="DSC009443_191127_135954"&PVPos==T,UniqueID,]))
plot(testCCF$CcfMatrix[,4])




testCCF <- InVivoR::SpikeCCF(x = SpikeTableTotal[UniqueID=="539_1 PV_DSC-007849_190408_142258", SpikeTime,], y = SpikeTableTotal[UniqueID=="539_1 PV_DSC-007849_190408_142258", SpikeTime,])

testData <- cumsum(rgamma(n = 1e4, shape = 10, rate = 90))
hist(testData)
testCCF <- InVivoR::SpikeCCF(x = cumsum(rgamma(n = 1e4, shape = 10, rate = 90)), y =cumsum(rgamma(n = 1e4, shape = 10, rate = 90)))


plot(y=testCCF$CCF, testCCF$xAxis, xlim = c(-0.02,0.02))
SpikeTableTotal[UniqueID=="506_1 DSC009443_191127_135954", .N,]

##### Theta phase to Table
SpikeTableTotal <- readRDS("/alzheimer/Daniel_Data/Analysis/Analysis_output/20200511/SpikeTableTotalWithFiringRate.rds")
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

ChannelLocationsFilter <- ChannelLocations[grepl(pattern = paste(unique(SpikeTableTotal[,SessionID,]), collapse = "|"), x = ...2)]

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


for(FolderNr in seq_along(ChannelLocationsFilter$...2)) {
  RefLoc <- as.vector(na.omit(unlist(ChannelLocationsFilter[FolderNr,9:11])))
  LFPRefChannels <- as.vector(na.omit(unlist(ChannelLocationsFilter[FolderNr,2:6])))
  print(data.frame(RefLoc, LFPRefChannels))
  RecordingFolder <- ChannelLocationsFilter$...2[FolderNr]
  ThetaExtraction(RefLoc = RefLoc[1:2], LFPRefChannels = LFPRefChannels[1:2], RecordingFolder = RecordingFolder, DataPath = DataPaths[[FolderNr]])
}

ThetaExtraction <- function(RefLoc, LFPRefChannels, RecordingFolder, DataPath) {
  SamplingRate <- 2e4
  OutputFolder <- paste0(RecordingFolder, "OutputFolder/")
  SessionName <- tail(unlist(strsplit(RecordingFolder, split = "/")), n = 1)
  if(!dir.exists(OutputFolder)) {
    dir.create(path = OutputFolder)
  }
  
  Cores <- 10
  
  message("load Amplifier File")
  if(length(DataPath$AmplifierFile)==2) {
    AmpFile <- rbind(InVivoR::AmpFileRead(FILENAME = grep(pattern = "amplifier_MS", x = DataPath$AmplifierFile, value = T), ChannelNumber = 32),
                     InVivoR::AmpFileRead(FILENAME = grep(pattern = "amplifier_PaS", x = DataPath$AmplifierFile, value = T), ChannelNumber = 32))
  } else {
    AmpFile <- InVivoR::AmpFileRead(FILENAME = DataPath$AmplifierFile, ChannelNumber = 64)
  }
  message("File Loaded")

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
  
  SpikeTableSessionThetaDelta <- SpikeTableTotal[grepl(pattern = SessionName, x = SessionID, ignore.case = T),]
  SpikeTableSessionIndex <- SpikeTableSessionThetaDelta[,floor(SpikeIdx/20),]
  
  for(i in seq_len(dim(LFPDownSample)[1])) {
    ThetaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Theta, FiltFilt = T, CORES = 10)
    DeltaData <- InVivoR::FirFiltering(SIGNAL = LFPDownSample[i,], FIR_FILTER = fir_filter_Delta, FiltFilt = T, CORES = 10)
    
    ThetaData_hb_transform <- hht::HilbertTransform(sig = ts(data = ThetaData, start = 0, frequency = 1e3))
    Theta_hb_power <- hht::HilbertEnvelope(ThetaData_hb_transform)
    ThetaData_phase <- (atan2(Re(ThetaData_hb_transform),Im(ThetaData_hb_transform))*180)/pi
    
    DeltaData_hb_transform <- hht::HilbertTransform(sig = ts(data = DeltaData, start = 0, frequency = 1e3))
    Delta_hb_power <- hht::HilbertEnvelope(ThetaData_hb_transform)
    DeltaData_phase <- (atan2(Re(DeltaData_hb_transform),Im(DeltaData_hb_transform))*180)/pi
    
    if(i==1) {
      SpikeTableSessionThetaDelta[,`:=`(MSThetaPower=Theta_hb_power[SpikeTableSessionIndex], MSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], MSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], MSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], MSLoc=RefLoc[i]),]
    } else {
      SpikeTableSessionThetaDelta[,`:=`(PaSThetaPower=Theta_hb_power[SpikeTableSessionIndex], PaSDeltaPower=Delta_hb_power[SpikeTableSessionIndex], PaSThetaPhase=ThetaData_phase[SpikeTableSessionIndex], PaSDeltaPhase=DeltaData_phase[SpikeTableSessionIndex], PaSLoc=RefLoc[i]),]
    }
  }
  saveRDS(object = LFPDownSample, file = paste0(OutputFolder, "LFPDownSample.rds"))
  saveRDS(object = SpikeTableSessionThetaDelta, file = paste0(OutputFolder, "SpikeTableSessionThetaDelta.rds"))
  rm(list = c("AmpFile", "SpikeTableSessionThetaDelta", "SpikeTableSessionIndex"))
  gc()
}

np <- reticulate::import("numpy")
spike_clusters <- np$load("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/spike_clusters.npy")
spike_times <- np$load("/alzheimer/Daniel_Data/DSC007849/PV_DSC-007849_190408_142258/amplifier_MS/spike_times.npy")

testData <- as.vector(spike_times)[spike_clusters==606]/2e4
