#!/usr/bin/env python

import ROOT, os
from os import listdir, makedirs, path, system
import sys
from subprocess import check_output

import argparse
parser = argparse.ArgumentParser(description='Run Skim and add trees for a sample', usage="./SkimTrees.py sample_name")
parser.add_argument('-j', '--ncpu', type=int, default=1, help="Number of CPUs to use with RDataFrame")
parser.add_argument('-o', '--overlap', type=str, default=None, choices=['0L','1L','2L', '2L_Hbb','1L_Hbb','0L_Hbb'], 
                    help="Skim ntuples for Overlap check. Selection in SR is applyed corresponding to 0L, 2L or 2L channels.\
                    Suffix _Hbb allows to apply selection of the VHbb analysis.")
parser.add_argument('--lxplus', action='store_true', default=False, help="Run on LXPLUS")
parser.add_argument("sample")

opt = parser.parse_args()


rootver = ROOT.gROOT.GetVersion()
if '/' in rootver: rootver = rootver.split('/')[0]
if float(rootver) >= 6.14:
    TH1DModel = ROOT.RDF.TH1DModel
    ROOTDataFrame = ROOT.RDataFrame
    print "Detected ROOT %s. Successfully loaded RDataFrame."%rootver
    usingRDF = True
    ROOT.ROOT.EnableImplicitMT(opt.ncpu)

else:
    TH1DModel = ROOT.Experimental.TDF.TH1DModel
    ROOTDataFrame = ROOT.Experimental.TDataFrame
    print "***** WARNING *****: Detected ROOT %s. Loaded TDataFrame instead of RDataFrame."%rootver


#in_path = "/pnfs/desy.de/cms/tier2/store/user/andrey/VHccAnalysisNtuples/2019_07_VJets_NLO/"
if opt.lxplus:
    in_path = "/eos/cms/store/group/phys_higgs/hbb/ntuples/VHcc/VHccNtuple_Jan2020_CvsL_225_CvsB_4/"
    out_path = "/eos/cms/store/group/phys_higgs/hbb/ntuples/VHcc/VHccNtuple_Jan2020_CvsL_225_CvsB_4/skimmed/"
else:
    #in_path  = "/nfs/dust/cms/user/andreypz/VHccAnalysisNtuples/2020_03_10_VHcc_VJetsNLO/"
    #in_path = "/nfs/dust/cms/user/lmastrol/VHccAnalysisNtuples/VHccNtuple_Test2017_DeepJet"
    #out_path = "/nfs/dust/cms/user/andreypz/VHccAnalysisNtuples/2020_03_23_VHcc_VJets_Luca_skimmed/"
    # in_path = "/pnfs/desy.de/cms/tier2/store/user/andrey/VHccAnalysisNtuples/"
    # in_path = "/nfs/dust/cms/user/lmastrol/VHbbAnalysisNtuples/VHccNtuple_March16_forApproval_fixZnnHF/haddjobs/"
    in_path = "/nfs/dust/cms/user/lmastrol/VHccAnalysisNtuples/VHccNtuple_Test2017_DeepJet_LO/haddjobs"
    out_path = "/nfs/dust/cms/user/spmondal/VHccRun2/BDTTraning/200413/"
    

if opt.overlap!=None:
    out_path += '_'+opt.overlap


if opt.sample!='all':
    onlyfiles = []
    if 'haddjobs' in in_path:
        onlyfiles = [f for f in listdir(in_path) if path.isfile(path.join(in_path, f))]
    elif opt.lxplus:
        for d in check_output(['eos', 'root://eoscms.cern.ch/', 'ls', in_path]).splitlines():
            if d=="haddjobs": continue
            if d!=opt.sample: continue
            print "Scanning path", d
            new_path = path.join(in_path, d)
            this_sample_files = [path.join(d,f) for f in check_output(['eos', 'root://eoscms.cern.ch/', 'ls', new_path]).splitlines() if '.root' in f]
            onlyfiles.extend(this_sample_files)
    else:
        for d in listdir(in_path):
            if d=="haddjobs": continue
            if d!=opt.sample: continue
            print "Scanning path", d
            new_path = path.join(in_path, d)
            this_sample_files = [path.join(d,f) for f in listdir(new_path) if path.isfile(path.join(new_path, f))]
            onlyfiles.extend(this_sample_files)

    # print onlyfiles
    print "Total number of files:", len(onlyfiles)
    if len(onlyfiles)==0:
        print "** Error: something is wrong. The number of files for that sample is zero..."
        sys.exit(0)
        

def createDir(myDir):
    if not path.exists(myDir):
        makedirs(myDir)

def doSample(sample):
    print 'Doing sample:', sample
    fChain = ROOT.TChain("Events")

    for ff in onlyfiles:
        #print ff
        if sample in ff:
            if opt.lxplus:
                add_file = 'root://eoscms.cern.ch/' + in_path + ff
                fChain.Add(add_file)
            else:
                add_file = in_path+'/'+ff
                print add_file
                if path.exists(add_file):
                    print 'Loading file', add_file
                    fChain.Add(add_file);

    #fChain.SetBranchStatus("*",0)

    d = ROOTDataFrame(fChain)

    # Pre-selection for 2-Lep SR 
    sel_2L = "controlSample==0 && (isZmm || isZee) && JetPt_1>20 && JetPt_2>20 && H_mass>50 && H_mass<200 && V_mass >75 && V_mass<105 && CvsL_1>0.4 && CvsB_1>0.2 && CvsL_2>0.0 && CvsB_2>0.0 && HVdPhi>2.5 && V_pt>=50"

    # Pre-selection for 1-Lep SR 
    sel_1L = "controlSample==0 && (isWenu || isWmunu) && JetPt_1>25 && JetPt_2>25 && (MET_Pt/sqrt(htJet30))>4 && H_pt>100 && CvsL_1>0.40 && CvsB_1>0.20 && CvsL_2>0.00 && CvsB_2>0.0 && H_mass>50 && H_mass<200 && HVdPhi>2.5 && nAddJetsFSRsub302p5_puid<1.5 && V_pt>=100"

    # Pre-selection for 0-Lep SR 
    sel_0L = "controlSample==0 && isZnn && max(JetPt_1, JetPt_2)>60 &&  min(JetPt_1, JetPt_2)>35 && H_pt>120 && H_mass>50 && H_mass<200 && MET_Pt>170 && CvsL_1>0.40 && CvsB_1>0.20 && CvsL_2>0.00 && CvsB_2>0.0 && HVdPhi>2.0 && nJetsCloseToMET==0 && dPhi_MET_TkMET < 0.5"

    
    # SR preselections in VHHBB analysis
    # Pre-selection for 2-Lep SR 
    sel_2L_Hbb = "controlSample==0 && (isZmm || isZee) && hJets_btagged_0>0.1522 && hJets_btagged_1 > 0.1522 && twoResolvedJets && V_pt>=75  && H_mass_fit_fallback > 90 && H_mass_fit_fallback < 150"
    sel_1L_Hbb = "controlSample==0 && (isWenu || isWmunu) && hJets_btagged_0>0.8001 && hJets_btagged_1 > 0.1522 && twoResolvedJets && V_pt>=150 && H_mass>90 && H_mass<150"
    sel_0L_Hbb = "controlSample==0 && isZnn && hJets_btagged_0>0.8001 && hJets_btagged_1>0.1522 && twoResolvedJets && H_mass>60 && H_mass<160 && V_pt>250"

    # No pre-selection:
    sel = "2*2==4"
    
    if opt.overlap=='2L':
        sel = sel_2L
    if opt.overlap=='1L':
        sel = sel_1L
    if opt.overlap=='0L':
        sel = sel_0L
    if opt.overlap=='2L_Hbb':
        sel = sel_2L_Hbb
    if opt.overlap=='1L_Hbb':
        sel = sel_1L_Hbb
    if opt.overlap=='0L_Hbb':
        sel = sel_0L_Hbb

    filt_0 = d.Filter(sel)
    # Use this to redefine variable names (V_pt and H_mass):
    # dOut = dOut.Define("NewName","OldName")

    dOut = filt_0
    # This is list of variables for Varial and BDT
    varsToSave = ['event', 'run', 'luminosityBlock', 'sampleIndex',
                  'weight', 'cTagWeight', 'recoWReWeight', 'recoZReWeight', 'weight_ptEWK', 'intWeight', 'weight_PU', 'WJetNLOWeight',
                  'controlSample', 'isZee', 'isZmm', 'isZnn', 'isWmunu', 'isWenu',
                  'Pass_nominal', 'cutFlow', 'twoResolvedJets',
                  'V_pt', 'V_mass', 'H_pt', 'H_mass',
                  'DeepJet_CvsL_1', 'DeepJet_CvsL_2', 'DeepJet_CvsB_1', 'DeepJet_CvsB_2',
                  'CvsL_1', 'CvsL_2', 'CvsB_1', 'CvsB_2',
                  'nJet', 'JetPt_1', 'JetPt_2',
                  'HJ1_pt', 'HJ2_pt', 'HJ1_HJ2_dEta', 'HJ1_HJ2_dPhi', 'HJ1_HJ2_dR', 'htJet30',
                  'nMuon', 'nElectron', 'lepInd1', 'lepInd2', 'Electron_pt', 'Electron_eta', 'Muon_pt', 'Muon_eta',
                  'HVdPhi', 'MET_Pt', 'lepMetDPhi', 'V_mt', 'dPhi_MET_TkMET',
                  'Lep_HJ1_dPhi', 'Lep1_HJ1_dPhi', 'Lep1_HJ2_dPhi',  'Lep2_HJ1_dPhi', 'Lep2_HJ2_dPhi', 'Lep1_Lep2_dPhi', 'Lep1_Lep2_dEta', 
                  'SA5', 'nJetsCloseToMET', 'nAddJets302p5_puid', 'nAddJetsFSRsub302p5_puid', 'jjVPtRatio',
                  'CMS_vzcc_BDTG_Wln_13TeV', 'CMS_vzcc_BDTG_Zll_13TeV', 'CMS_vzcc_BDTG_Znn_13TeV',
                  'CMS_vhcc_BDTG_Wln_13TeV', 'CMS_vhcc_BDTG_Zll_LowPT_13TeV', 'CMS_vhcc_BDTG_Zll_HighPT_13TeV', 'CMS_vhcc_BDTG_Znn_13TeV','BDTJetOrderScore'
                  #'CMS_vhcc_BDTG_NLO_Wln_13TeV', 'CMS_vhcc_BDTG_NLO_Zll_LowPT_13TeV', 'CMS_vhcc_BDTG_NLO_Zll_HighPT_13TeV', 'CMS_vhcc_BDTG_NLO_Znn_13TeV'
                  ]
    if 'Run' not  in sample: # only add from MC, not the data
        varsToSave.extend(['LHE_Vpt', 'LHE_HT', 'genWeight'])
        
    if opt.overlap!=None:
        # This is a short list of variables for Event overlap check
        varsToSave = ['event', 'run', 'luminosityBlock', 'V_pt', 'H_mass']

    print "Variables to be saved in the output tree:\n", varsToSave

    outVars = ROOT.std.vector('string')()
    for v in varsToSave:
        outVars.push_back(v)

    if opt.lxplus:
        dOut.Snapshot("tree", "sum_"+sample+".root", outVars)
        #dOut.Snapshot("tree", "sum_"+sample+".root", outVars)
        #system("xrdcp -f sum_%s.root %s" % (sample, out_path) )
        #system("rm sum_%s.root" % sample)
    else:
        dOut.Snapshot("tree", out_path+"/sum_"+sample+".root", outVars)

if __name__ == "__main__":
    print "This is the __main__ part"


    import time
    start = time.time()

    createDir(out_path)

    if opt.sample == 'all':
        print 'Run over all samples'
        print 'Will submit on condor'
        
        createDir('outCondor')

        with open("skimmer_runscript.sh",'w') as runf:
            
            cond_run = '''#!/bin/bash

export ORIG_DIR=$PWD
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_0/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $ORIG_DIR
ls
echo "** Running parameters: "
echo $1 $2 $3
echo "Will start skimming sample $1"
./SkimTrees.py $1 -j $2 $3
'''
            if opt.lxplus:
                cond_run += """
xrdcp -f sum_$1.root root://eoscms.cern.ch/%s
rm sum_$1.root
""" % out_path
            runf.write(cond_run)

        with open("skimmer_config.sub",'w') as subf:
            
            if opt.lxplus:
                cond_conf = '''universe = vanilla
Executable     =  skimmer_runscript.sh
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
Notification     = never
transfer_input_files = SkimTrees.py
requirements = OpSysAndVer == "CentOS7"
request_cpus = %NCPU%
transfer_output_files=""
Output = outCondor/job_$(Cluster)_$(Process).out
Error  = outCondor/job_$(Cluster)_$(Process).err
Log    = outCondor/job_$(Cluster)_$(Process).log
+MaxRuntime = 21600
Queue Arguments from (
%ARGUMENTS%)
'''            
            else:
                cond_conf = '''universe = vanilla
Executable     =  skimmer_runscript.sh
Should_Transfer_Files     = YES
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
Notification     = never
transfer_input_files = SkimTrees.py
WhenToTransferOutput=On_Exit
requirements = OpSysAndVer == "CentOS7"
request_cpus = %NCPU%
Output = outCondor/job_$(Cluster)_$(Process).out
Error  = outCondor/job_$(Cluster)_$(Process).err
Log    = outCondor/job_$(Cluster)_$(Process).log
##+RequestRuntime = 100000
Queue Arguments from (
%ARGUMENTS%)
'''            
            args_for_cond = ''
            for samp in listdir(in_path):
                if 'hadd' in samp: continue
                if 'skim' in samp: continue
                if len([i for i in os.listdir(out_path) if samp in i])>0: continue
                if opt.lxplus: args_for_cond += '%s %i --lxplus\n'% (samp, opt.ncpu)
                else: args_for_cond += '%s %i\n'% (samp, opt.ncpu)
                # print samp
            # print args_for_cond
            cond_conf = cond_conf.replace("%ARGUMENTS%", args_for_cond).replace("%NCPU%", str(opt.ncpu))
            print cond_conf
            subf.write(cond_conf)

            print 'Submition file', subf.name, 'is created'
            print 'Run it with condor_submit'
            # print 'Doing', samp
            # doSample(samp)

    else:
        doSample(opt.sample)


    end = time.time()

    print 'The output will be stored at:\n\t ', out_path
    print "Run time:", (end - start), 'sec'
