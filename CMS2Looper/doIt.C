{

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");
  gSystem->Load("libGenVector.so");
  gSystem->Load("libCMS2NtupleMacrosCORE.so");
  gSystem->Load("liblooper.so");
 
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  chain->Add("/tas/cms2/VVJetsTo4L_TuneD6T_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-12/wwfilter/merged_ntuple.root") ;
  l->ScanChain(chain,"test",0);



}
