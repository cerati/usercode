// 
// RUN BY root -b -q plotNuis.C
// 
void plotNuis() {
  int mH = 125; 
  TString inj = "125";
  TString dir_result = "cards_inj_statsyst78_new/";
  gSystem->Exec(Form("mkdir -p %s/plots",dir_result.Data()));
  
  TString ana = "hww";
  int ntoys = 1002; 
  for ( int njet = 0; njet <= 0; njet ++ ) {
    plotNuisSingle(inj,njet,mH,dir_result,ana,ntoys);
  }
}

void createFillTH1F(TDirectory* fout, const char* name,const char* title,int nbins,float minx,float maxx,const char* xtitle,const char* ytitle,float value) {

  TH1F* h = 0;
  h = (TH1F*) fout->Get(name);
  if (h==0) {
    fout->cd();
    //cout << "create new hist" << endl;
    h = new TH1F(name,title,nbins,minx,maxx);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->SetDirectory(fout);
  }
  h->Fill(value);

}

void plotNuisSingle(TString inj,int jet, int mH, TString dir, TString ana, int ntoys) {
  
  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  bool debug = false;

  TFile* fout = TFile::Open("dummy.root","RECREATE");

  for(int i=0; i<ntoys; i++) {
    TString fitresults= Form("%s/logsNorm/%i/mlfit_injm%s_m%i_%sof_%ij_id%i.root", dir.Data(), mH, inj.Data(), mH, ana.Data(), jet, i); 
    if ( debug ) std::cout << "Opening " << fitresults << "\n";
    TFile *File = TFile::Open(fitresults, "READ");
    if ( File == 0x0  ) { continue; }
    RooFitResult *fit_s = (RooFitResult*) File->Get("fit_s");
    if( fit_s == 0x0 )  { File->Close(); continue; }
    if(fit_s->status() != 0) { delete fit_s; File->Close(); continue; } // fit status == 0 : requires fit quality
    
    RooArgList parlist = fit_s->floatParsFinal();
    // 
    // Loop over RooArgList and store the fit results 
    // in the above plots
    // 
    for(int j=0; j<parlist.getSize(); j++) {
      //create and make the plots
      TString name_tstr = parlist[j].GetName(); 
      float pull = ((RooRealVar&)parlist[j]).getVal();
      if (name_tstr=="r")  pull=pull-1.;
      float uncert = ((RooRealVar&)parlist[j]).getError();
      createFillTH1F(fout,TString(name_tstr+"_pull").Data(),TString(name_tstr+"_pull").Data(),500,-5,5,"(nuis_{fit} - nuis_{in})/#sigma_{nuis}","toys/bin",pull/uncert);
      createFillTH1F(fout,TString(name_tstr+"_uncert").Data(),TString(name_tstr+"_uncert").Data(),200,0,2,"#sigma_{nuis}","toys/bin",uncert);
    }
    if ( debug ) std::cout << "now closing" << std::endl;
    delete fit_s;
    File->Close();
  }
  
  fout->Write();

  if ( debug ) std::cout << "now printing" << std::endl;
  //print the plots
  fout->cd();
  TCanvas* c1 = new TCanvas();
  TIter next_G(fout->GetListOfKeys());
  TKey *key_G;
  while ((key_G = (TKey*)next_G())) {
    TClass *cl = gROOT->GetClass(key_G->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_G = (TH1F*)key_G->ReadObj();
    if ( debug ) std::cout << h_G->GetTitle() << std::endl;
    h_G->Draw();
    h_G->Fit("gaus");
    c1->SaveAs(Form("%s/plots/%s.png",dir.Data(),h_G->GetTitle()));
  }

  fout->Close();
  gSystem->Exec("rm dummy.root");
  
}
 
/*

*/
