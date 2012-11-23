{

  gStyle->SetOptStat(0);

  TCanvas c1;

  TFile *_file_G = TFile::Open("hwwof_0j.input_8TeV.root");
  TFile *_file_g = TFile::Open("hwwof_0j.input.root");

  TFile *_file_new = TFile::Open("hwwof_0j.input_mixtest.root","RECREATE");
  _file_new->cd();

  TString discriminator = "MVAMETResBounding";

  TIter next_G(_file_G->GetListOfKeys());
  TKey *key_G;
  while ((key_G = (TKey*)next_G())) {
    TClass *cl = gROOT->GetClass(key_G->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_G = (TH1F*)key_G->ReadObj();
    if (TString(h_G->GetName()).Contains(discriminator.Data())==0) continue;
    TH1F* h_G_clone = new TH1F();
    h_G->Copy(*h_G_clone);
    h_G_clone.Write();
  }

  TIter next_g(_file_g->GetListOfKeys());
  TKey *key_g;
  while ((key_g = (TKey*)next_g())) {
    TClass *cl = gROOT->GetClass(key_g->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_g = (TH1F*)key_g->ReadObj();
    if (TString(h_g->GetName()).Contains(discriminator.Data())) continue;
    TH1F* h_g_clone = new TH1F();
    h_g->Copy(*h_g_clone);
    h_g_clone.Write();
  }

  _file_new->Close();

}
