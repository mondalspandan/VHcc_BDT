#!/usr/bin/env python
#Last Updated: 14/04/2020
from ROOT import TMVA, TFile, TCut, TChain
import argparse,os
import time
from random import random

# Setup argparse
parser = argparse.ArgumentParser(description = "TMVA Classification", usage = "./ClassificationKeras.py -s [] -l []")
parser.add_argument("-s", "--sample", type=str, default = "2LL", choices=["0L", "1L", "2LL", "2LH"], help = "Sample to pick up from 0L, 1L, 2LL, 2LH")
parser.add_argument("-l", "--location", default = "Desy", choices = ["rwth", "hpc", "Tier2", "Desy"], help = "Location where the program is executed")
opt = parser.parse_args()

samp = opt.sample

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
output = TFile.Open("TMVASolution_%s.root"%samp, "RECREATE")
factory = TMVA.Factory("TMVAClassification", output, "!V:!Silent:!Color:DrawProgressBar:Transformations=I:AnalysisType=Classification")
dataloader = TMVA.DataLoader("dataset_%s"%samp)

RunningTimeStart = time.time()

#Define the paths:
if opt.location == "hpc": path = "/rwthfs/rz/cluster/hpcwork/sc532682/DataFiles2017/VHccNtuple_Dec2018_DeepJet_newSF-pTincl-noLowStat_skimmed_"
if opt.location == "rwth": path = "/net/scratch_cms3a/rocamora/VHccNtuple_Dec2018_DeepJet_newSF-pTincl-noLowStat_skimmed_"
if opt.location == "Tier2": path = "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/rocamora/DataFiles2017/VHccNtuple_Dec2018_DeepJet_newSF-pTincl-noLowStat_skimmed_"
if opt.location == "Desy": path = "/nfs/dust/cms/user/spmondal/vhccRun2/BDTTraning/200413/"


#Adding the Cuts
if samp == "2LL":
    myCut = TCut("(isZmm || isZee) && JetPt_1>20 && JetPt_2>20 && DeepJet_CvsL_1>0.225 && DeepJet_CvsB_1>0.40 && DeepJet_CvsL_2>0.0 && DeepJet_CvsB_2>0.0 && controlSample==0 && H_mass>50 && H_mass<200 && V_mass >75 && V_mass<105 && V_pt>=50 && V_pt<150 && HVdPhi>2.5")

if samp == "2LH":
    myCut = TCut("(isZmm || isZee) && JetPt_1>20 && JetPt_2>20 && DeepJet_CvsL_1>0.225 && DeepJet_CvsB_1>0.40 && DeepJet_CvsL_2>0.0 && DeepJet_CvsB_2>0.0 && controlSample==0 && H_mass>50 && H_mass<200 && V_mass >75 && V_mass<105 && V_pt>150 && HVdPhi>2.5")

if samp == "1L":
    myCut = TCut("(isWmunu || isWenu) && JetPt_1>25 && JetPt_2>25 && (MET_Pt/sqrt(htJet30))>4 && H_pt>100 && DeepJet_CvsL_1>0.225 && DeepJet_CvsB_1>0.40 && DeepJet_CvsL_2>0.00 && DeepJet_CvsB_2>0.0 && controlSample==0 && H_mass>50 && H_mass<200 && V_pt>100 && HVdPhi>2.5 && nAddJetsFSRsub302p5_puid<1.5")

if samp == "0L":
    myCut = TCut("isZnn && nJetsCloseToMET==0 && dPhi_MET_TkMET < 0.5 && max(JetPt_1,JetPt_2) > 60 && min(JetPt_1,JetPt_2) > 35 && MET_Pt>170 && H_pt>120 && controlSample==0 && H_mass>50 && H_mass<200 && DeepJet_CvsL_1>0.40 && DeepJet_CvsB_1>0.20 && DeepJet_CvsL_2>0.00 && DeepJet_CvsB_2>0.00 && HVdPhi>2.0")


# Load the Data
signalTrainChain = TChain("tree")
backgroundTrainChain = TChain("tree")
# signalTestChain = TChain("tree")
# backgroundTestChain = TChain("tree")

#####################################       Signal
if samp == "2LH" or samp == "2LL":
    sigs = ["_ZH125ToCC_ZLL", "ggZH125ToCC_ZLL"]

if samp == "1L":
    sigs = ["WminusH125ToCC_WLNu","WplusH125ToCC_WLNu"]

if samp == "0L":
    sigs = ["_ZH125ToCC_ZNuNu","ggZH125ToCC_ZNuNu"]


#####################################       Background

bkgs = ["DYToLL","ST_","TT_","WJets","WW","WZ","ZZ","WminusH125_powheg","WplusH125_powheg"]

if samp == "0L": 
    bkgs.extend(["ZJetsToNuNu","_ZH125_ZNuNu","ggZH125_ZNuNu"])

if samp == "2LH" or samp == "2LL": 
    bkgs.extend(["ggZH125_ZLL","_ZH125_ZLL"])

sigfllist = []
bkgfllist = []
for ind,ifl in enumerate(os.listdir(path)):
    for sig in sigs:
        if sig in ifl and ifl not in sigfllist:
            signalTrainChain.Add("%s/%s"%(path,ifl))
            sigfllist.append(ifl)
            # signalTestChain.Add("%s/%s"%(path,ifl))
    for bkg in bkgs:
        if bkg in ifl and ifl not in bkgfllist:
            backgroundTrainChain.Add("%s/%s"%(path,ifl))
            # backgroundTestChain.Add("%s/%s"%(path,ifl))
            bkgfllist.append(ifl)

# print sigfllist
# print bkgfllist

# print signalTrainChain.GetEntries()
# print backgroundTrainChain.GetEntries()

# dataloader.AddSignalTree(signalTrainChain, 1.0, TMVA.Types.kTraining)
# dataloader.AddSignalTree(signalTrainChain, 1.0, TMVA.Types.kTesting)
# dataloader.AddBackgroundTree(backgroundTrainChain, 1.0, TMVA.Types.kTraining)
# dataloader.AddBackgroundTree(backgroundTrainChain, 1.0, TMVA.Types.kTesting)

dataloader.AddTree(signalTrainChain, "Signal", 1.0, myCut+TCut("event%2 == 0"),"train")
dataloader.AddTree(signalTrainChain, "Signal", 1.0, myCut+TCut("event%2 == 1"),"test")
dataloader.AddTree(backgroundTrainChain, "Background", 1.0, myCut+TCut("event%2 == 0"),"train")
dataloader.AddTree(backgroundTrainChain, "Background", 1.0, myCut+TCut("event%2 == 1"),"test")


dataloader.SetSignalWeightExpression("weight")
dataloader.SetBackgroundWeightExpression("weight")

# Addint the variables
if samp == "2LL" or samp == "2LH":
    varNames = ["H_mass", "H_pt", "DeepJet_CvsL_1", "DeepJet_CvsL_2", "DeepJet_CvsB_1", "DeepJet_CvsB_2", "V_pt", "nAddJetsFSRsub302p5_puid", "SA5",
        "JetPt_1", "JetPt_2", "jjVPtRatio", "HJ1_HJ2_dEta", "HJ1_HJ2_dR", "HVdPhi", "Lep2_HJ1_dPhi", "Lep2_HJ2_dPhi", 
        "Lep1_Lep2_dPhi", "Lep1_Lep2_dEta"]

if samp == "1L":
    varNames = ["H_mass", "H_pt", "jjVPtRatio", "V_pt", "DeepJet_CvsL_1", "DeepJet_CvsL_2", "DeepJet_CvsB_1", "DeepJet_CvsB_2", "nAddJetsFSRsub302p5_puid", "HVdPhi", 
        "SA5", "lepMetDPhi", "V_mt", "MET_Pt", "JetPt_1", "JetPt_2", "HJ1_HJ2_dPhi", "HJ1_HJ2_dR", "HJ1_HJ2_dEta", 
        "Lep_HJ1_dPhi"]

if samp == "0L":
    varNames = ["H_mass", "H_pt", "HVdPhi", "V_pt", "HJ1_HJ2_dEta", "DeepJet_CvsL_1", "DeepJet_CvsL_2", "DeepJet_CvsB_1", "DeepJet_CvsB_2",
        "SA5", "HJ1_HJ2_dPhi", "JetPt_1", "JetPt_2", "nAddJetsFSRsub302p5_puid", "MET_Pt", "jjVPtRatio"]

for name in varNames:
    if name == "SA5":
        dataloader.AddVariable("%s" % name, "I")
    else:
        dataloader.AddVariable("%s" % name, "F")


dataloader.PrepareTrainingAndTestTree(myCut, "SplitMode=Random:NormMode=NumEvents:V")


if samp == "2LH" or samp == "2LL":
    factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDTG', "!H:!V:NTrees=750:MinNodeSize=7.5%:VarTransform=None:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4")
if samp == "1L":
    factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=7.5%:VarTransform=None:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3")
if samp == "0L":
    factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=500:MinNodeSize=5%:VarTransform=None:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3")

# Run training, test and evaluation
TrainingTimeStart = time.time()
factory.TrainAllMethods()
TrainingTimeStop = time.time()

factory.TestAllMethods()
factory.EvaluateAllMethods()

RunningTimeStop = time.time()

TrainingTime = time.strftime("%H:%M:%S", time.gmtime(TrainingTimeStop-TrainingTimeStart))
RunningTime = time.strftime("%H:%M:%S", time.gmtime(RunningTimeStop-RunningTimeStart))

print "\n\nTime to run all the code = %s s" % RunningTime
print "Time to train = %s s\n" % TrainingTime


dataSetInfo = dataloader.GetDataSetInfo()
dataSet = dataSetInfo.GetDataSet()

NSgnTrainEvents = dataSet.GetNEvtSigTrain()
NSgnTestEvents = dataSet.GetNEvtSigTest()
NBkgTrainEvents = dataSet.GetNEvtBkgdTrain()
NBkgTestEvents = dataSet.GetNEvtBkgdTest()

print "Dataset information:"
print "\tBackground -- training events            : %i" % int(NBkgTrainEvents)
print "\tBackground -- testing events             : %i" % int(NBkgTestEvents)
print "\tSignal     -- training events            : %i" % int(NSgnTrainEvents)
print "\tSignal     -- testing events             : %i\n" % int(NSgnTestEvents)
print "\tTotal Background Events                  : %i" % (int(NBkgTrainEvents)+int(NBkgTestEvents))
print "\tTotal Signal Events                      : %i" % (int(NSgnTrainEvents)+int(NSgnTestEvents))