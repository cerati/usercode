void doData()
{
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libCMS2NtupleMacrosCORE.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");
  gSystem->SetIncludePath(Form("%s -I./Tools", gSystem->GetIncludePath()));

  gROOT->LoadMacro("myBabyMaker.C++");
  TChain *ch_data = new TChain( "Events" );
  ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-01-02/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-01-02_merged/V04-01-02/merged_ntuple_16*.root") ;
  ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-01-03/DoubleElectron_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/merged_ntuple_16*.root") ;

  myBabyMaker * baby = new myBabyMaker();
  baby->ScanChain(ch_data, "baby_JunkRun2011A.root", true);

}
