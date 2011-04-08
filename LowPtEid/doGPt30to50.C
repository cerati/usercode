{

  gSystem->Load("libCMS2NtupleMacrosCORE.so");

  gROOT->ProcessLine(".L ScanChain.C++");
 
  gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");

  TChain *chain = new TChain("Events");
  chain->Add("ntuple_GPt30to50.root") ;
  ScanChain(chain);



}
