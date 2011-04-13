{

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");
  gSystem->Load("libCMS2NtupleMacrosCORE.so");
  gSystem->Load("liblooper.so");
 
  gSystem->Load("./Tools/MiniFWLite/libMiniFWLite.so");

  looper *l = new looper();
  TChain *chain = new TChain("Events");
  chain->Add("/tas/cms2/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/ntuple_1.root") ;
  l->ScanChain(chain,"test",1);



}
