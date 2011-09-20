void doData()
{
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libCMS2NtupleMacrosCORE.so");

  //gSystem->Load("./Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");
  //gSystem->SetIncludePath(Form("%s -I./Tools", gSystem->GetIncludePath()));

  gROOT->LoadMacro("myBabyMaker.C++");
  TChain *ch_data = new TChain( "Events" );
  ch_data->Add("/nfs-6/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v6_AOD/V04-02-30/DoubleElectronTriggerSkim/skimmed_ntuple_173*.root");
  myBabyMaker * baby = new myBabyMaker();
  baby->ScanChain(ch_data, "baby_data_v6.root", true);

}
