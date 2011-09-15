void doMC()
{
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libCMS2NtupleMacrosCORE.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");
  gSystem->SetIncludePath(Form("%s -I./Tools", gSystem->GetIncludePath()));

  gROOT->LoadMacro("myBabyMaker.C++");

  myBabyMaker * baby = new myBabyMaker();

  // HWW 130 with PU
  TChain *ch_hww130 = new TChain( "Events" );
  cout << "babifying HWW 130" << endl;
  ch_hww130->Add("/nfs-3/userdata/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged_ntuple*.root");
  baby->ScanChain(ch_hww130, 	"baby_GluGluToHToWWTo2L2Nu_M-130_7TeV.root", false);

}
