void doData()
{
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libCMS2NtupleMacrosCORE.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");
  gSystem->SetIncludePath(Form("%s -I./Tools", gSystem->GetIncludePath()));

  gROOT->LoadMacro("myBabyMaker.C++");
  TChain *ch_data = new TChain( "Events" );

  //ch_data->Add("/nfs-3/userdata/cms2/EGMonitor_Run2010B-Nov4ReReco_v2_RECO/V03-06-17/merged_ntuple_78.root");
  ch_data->Add("/nfs-3/userdata/cms2/EGMonitor_Run2010B-Nov4ReReco_v2_RECO/V03-06-17/merged_ntuple*.root");
  //ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-Nov4ReReco_v1_RECO/V03-06-16/merged_ntuple*.root");

  myBabyMaker * baby = new myBabyMaker();
  //baby->ScanChain(ch_data, "baby_EG2010A.root", true);
  baby->ScanChain(ch_data, "baby_EGMon2010B.root", true);

}
