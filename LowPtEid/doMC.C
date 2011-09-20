void doMC()
{
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  gSystem->Load("./Tools/EgammaAnalysisTools/lib/libCMS2NtupleMacrosCORE.so");

  //gSystem->Load("./Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");
  //gSystem->SetIncludePath(Form("%s -I./Tools", gSystem->GetIncludePath()));

  gROOT->LoadMacro("myBabyMaker.C++");

  myBabyMaker * baby = new myBabyMaker();

//   TChain *ch_dyee = new TChain( "Events" );
//   cout << "babifying DYEE" << endl;
//   ch_dyee->Add("/nfs-7/userdata/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*.root");
//   baby->ScanChain(ch_dyee, 	"baby_dyee.root", false);

  TChain *ch_ww2l = new TChain( "Events" );
  cout << "babifying WW2L" << endl;
  ch_ww2l->Add("/nfs-7/userdata/cms2/WWTo2L2Nu_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged_ntuple*.root");
  baby->ScanChain(ch_ww2l, 	"baby_ww2l.root", false);

}
