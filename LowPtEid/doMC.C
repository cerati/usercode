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
  ch_hww130->Add("/nfs-3/userdata/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/merged_ntuple*.root");
  baby->ScanChain(ch_hww130, 	"baby_GluGluToHToWWTo2L2Nu_M-130_7TeV_LOPU.root", false);

//   // HWW 130 no PU
//   TChain *ch_hww130nopu = new TChain( "Events" );
//   cout << "babifying HWW 130 NO PU" << endl;
//   ch_hww130nopu->Add("/nfs-3/userdata/cms2/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-START38_V12-v1/V03-06-18/merged_ntuple*.root");
//   baby->ScanChain(ch_hww130nopu, "baby_GluGluToHToWWTo2L2Nu_M-130_7TeV_NOPU.root", false);

//   // QCD5080 no PU
//   TChain *ch_qcd5080nopu = new TChain( "Events" );
//   cout << "babifying QCD5080 NO PU" << endl;
//   ch_qcd5080nopu->Add("/nfs-3/userdata/cms2/QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root");
//   baby->ScanChain(ch_qcd5080nopu, "baby_QCD_Pt_50to80_NOPU.root", false);


//   // QCD 50-80
//   TChain *ch_qcd5080 = new TChain( "Events" );
//   ch_qcd5080->Add("/nfs-3/userdata/cms2/QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root");
//   cout << "babifying QCD 50-80" << endl;
//   baby->ScanChain(ch_qcd5080, 	"baby_QCD_Pt_50to80.root", false);

//   // QCD 30-50
//   TChain *ch_qcd3050 = new TChain( "Events" );
//   ch_qcd3050->Add("/nfs-3/userdata/cms2/QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10-START38_V12-v1/V03-06-14/merged_ntuple*.root");
//   cout << "babifying QCD 30-50" << endl;
//   baby->ScanChain(ch_qcd3050, 	"baby_QCD_Pt_30to50.root", false);
}
