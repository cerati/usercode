// 
// RUN BY root -b -q plotNorm.C\(\"cards_inj_stat\"\)
// 
void plotNorm(TString dir_result) {
  int mH = 125; 
  TString inj = "125";
  gSystem->Exec(Form("mkdir -p %s/plots",dir_result.Data()));
  
  TString ana = "hww";
  int ntoys = 10000; 
  for ( int njet = 0; njet <= 1; njet ++ ) {
    plotNormSingle(inj,njet,mH,dir_result,ana,ntoys);
  }
}

void plotNormSingle(TString inj,int jet, int mH, TString dir, TString ana, int ntoys) {
  
  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  TGaxis::SetMaxDigits(1);
  bool debug = false;

  // prefit yields
  float yield_ZH(0.), yield_WH(0.), yield_qqH(0.), yield_ggH(0.);
  float yield_qqWW(0.), yield_ggWW(0.), yield_VV(0.), yield_Top(0.0), yield_Zjets(0.0), yield_WjetsE(0.), yield_Wgamma(0.0), yield_Wg3l(0.0), yield_Ztt(0.), yield_WjetsM(0.); 
  
  TString dir_cards="cards_def";
  TString datacardName = Form("%s/%i/%sof_%ij_shape_8TeV.txt",dir_cards.Data(),mH,ana.Data(),jet);
  std::cout << datacardName << "\n";
  ifstream incard (datacardName);
  string line;
  
  if (incard.is_open()) {
    while ( incard.good() ) {
      getline (incard,line);
      size_t found=line.find("rate");
      if (found!=string::npos) {
	std::cout << line << std::endl;
	stringstream stream(line);
	string rate;
	stream >> rate ;
	stream >> yield_ZH; 
	stream >> yield_WH;
	stream >> yield_qqH;
	stream >> yield_ggH;
	stream >> yield_qqWW;
	stream >> yield_ggWW;
	stream >> yield_VV;    
	stream >> yield_Top;
	stream >> yield_Zjets;
	stream >> yield_WjetsE;
	stream >> yield_Wgamma;
	stream >> yield_Wg3l;
	stream >> yield_Ztt;
	stream >> yield_WjetsM;
      }
    }
    incard.close();
  }
  
  // after fit yields

  // for PullMu
  TH1F* h_PullMu_s = new TH1F("PullMu_sig_fit","PullMu",50,-4,4);
  TH1F* h_ggH_s = new TH1F("ggH_sig_fit","ggH",50,-2,2);
  TH1F* h_qqWW_s = new TH1F("qqWW_sig_fit","qqWW",50,-0.5,0.5);
  TH1F* h_ggWW_s = new TH1F("ggWW_sig_fit","ggWW",50,-0.5,0.5);
  TH1F* h_Top_s = new TH1F("Top_sig_fit","Top",50,-1,1);
  TH1F* h_WjetsE_s = new TH1F("WjetsE_sig_fit","WjetsE",50,-1,1);
  TH1F* h_WjetsM_s = new TH1F("WjetsM_sig_fit","WjetsM",50,-1,1);
  TH1F* h_Wgamma_s = new TH1F("Wgamma_sig_fit","Wgamma",50,-1,1);
  TH1F* h_Wg3l_s = new TH1F("Wg3l_sig_fit","Wg3l",50,-1,1);
  TH1F* h_VV_s = new TH1F("VV_sig_fit","VV",50,-0.2,0.2);
  TH1F* h_Ztt_s = new TH1F("Ztt_sig_fit","Ztt",50,-0.2,0.2);
  
  Float_t yield_ggH_sfit, yield_qqWW_sfit, yield_ggWW_sfit, yield_Top_sfit, yield_WjetsE_sfit, yield_WjetsM_sfit, yield_Wgamma_sfit, yield_Wg3l_sfit, yield_VV_sfit, yield_Ztt_sfit;
  
  for(int i=0; i<ntoys; i++) {
    TString fitresults= Form("%s/logsNorm/%i/mlfit_injm%s_m%i_%sof_%ij_id%i.root", dir.Data(), mH, inj.Data(), mH, ana.Data(), jet, i); 
    if ( debug ) 
      std::cout << "Opening " << fitresults << "\n";
    TFile *File = TFile::Open(fitresults, "READ");
    if ( File == 0x0  ) { continue; }
    RooFitResult *fit_s = (RooFitResult*) File->Get("fit_s");
    if( fit_s == 0x0 )  { File->Close(); continue; }
    if(fit_s->status() != 0) { delete fit_s; File->Close(); continue; } // fit status == 0 : requires fit quality
    
    //
    // get signal fit strength
    // 
    RooRealVar *r = (RooRealVar*) fit_s->floatParsFinal().find("r");
    h_PullMu_s->Fill( ( r->getVal() - 1. ) / r->getError() ); 
    if ( debug ) 
      std::cout <<  r->getVal() << " +/- " <<  r->getError() << " ==> " <<  ( r->getVal() - 1. ) / r->getError() << "\n";
    
    // 
    // getting the signal fit yield
    // 
    RooArgSet *norm_fit_s = (RooArgSet*) File->Get("norm_fit_s");
    if ( norm_fit_s ) {
      
      // ggH 
      RooRealVar* yield_sfit_ggH; 
      yield_sfit_ggH = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_ggH", jet)); 
      if ( yield_sfit_ggH != 0x0 ) {
	yield_ggH_sfit = yield_sfit_ggH->getVal();  
	if ( debug ) 
	  std::cout << "yield_ggH_sfit = " << yield_ggH_sfit << "\n";
	h_ggH_s->Fill( (yield_ggH_sfit-yield_ggH)/yield_ggH );
      }
      
      // qqWW 
      RooRealVar* yield_sfit_qqWW; 
      yield_sfit_qqWW = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_qqWW", jet)); 
      if ( yield_sfit_qqWW != 0x0 ) {
	yield_qqWW_sfit = yield_sfit_qqWW->getVal();  
	if ( debug ) 
	  std::cout << "yield_qqWW_sfit = " << yield_qqWW_sfit << "\n";
	h_qqWW_s->Fill( (yield_qqWW_sfit-yield_qqWW)/yield_qqWW );
      }
      
      // ggWW 
      RooRealVar* yield_sfit_ggWW; 
      yield_sfit_ggWW = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_ggWW", jet)); 
      if ( yield_sfit_ggWW != 0x0 ) {
	yield_ggWW_sfit = yield_sfit_ggWW->getVal();  
	if ( debug ) 
	  std::cout << "yield_ggWW_sfit = " << yield_ggWW_sfit << "\n";
	h_ggWW_s->Fill( (yield_ggWW_sfit-yield_ggWW)/yield_ggWW );
      }

      // Top 
      RooRealVar* yield_sfit_Top; 
      yield_sfit_Top = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_Top", jet)); 
      if ( yield_sfit_Top != 0x0 ) {
	yield_Top_sfit = yield_sfit_Top->getVal();  
	if ( debug ) 
	  std::cout << "yield_Top_sfit = " << yield_Top_sfit << "\n";
	h_Top_s->Fill( (yield_Top_sfit-yield_Top)/yield_Top );
      }

      // WjetsE 
      RooRealVar* yield_sfit_WjetsE; 
      yield_sfit_WjetsE = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_WjetsE", jet)); 
      if ( yield_sfit_WjetsE != 0x0 ) {
	yield_WjetsE_sfit = yield_sfit_WjetsE->getVal();  
	if ( debug ) 
	  std::cout << "yield_WjetsE_sfit = " << yield_WjetsE_sfit << "\n";
	h_WjetsE_s->Fill( (yield_WjetsE_sfit-yield_WjetsE)/yield_WjetsE );
      }
      
      // WjetsM 
      RooRealVar* yield_sfit_WjetsM; 
      yield_sfit_WjetsM = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_WjetsM", jet)); 
      if ( yield_sfit_WjetsM != 0x0 ) {
	yield_WjetsM_sfit = yield_sfit_WjetsM->getVal();  
	if ( debug ) 
	  std::cout << "yield_WjetsM_sfit = " << yield_WjetsM_sfit << "\n";
	h_WjetsM_s->Fill( (yield_WjetsM_sfit-yield_WjetsM)/yield_WjetsM );
      }

      // Wg3l 
      RooRealVar* yield_sfit_Wg3l; 
      yield_sfit_Wg3l = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_Wg3l", jet)); 
      if ( yield_sfit_Wg3l != 0x0 ) {
	yield_Wg3l_sfit = yield_sfit_Wg3l->getVal();  
	if ( debug ) 
	  std::cout << "yield_Wg3l_sfit = " << yield_Wg3l_sfit << "\n";
	h_Wg3l_s->Fill( (yield_Wg3l_sfit-yield_Wg3l)/yield_Wg3l );
      }

      // Wgamma 
      RooRealVar* yield_sfit_Wgamma; 
      yield_sfit_Wgamma = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_Wgamma", jet)); 
      if ( yield_sfit_Wgamma != 0x0 ) {
	yield_Wgamma_sfit = yield_sfit_Wgamma->getVal();  
	if ( debug ) 
	  std::cout << "yield_Wgamma_sfit = " << yield_Wgamma_sfit << "\n";
	h_Wgamma_s->Fill( (yield_Wgamma_sfit-yield_Wgamma)/yield_Wgamma );
      }

      // VV 
      RooRealVar* yield_sfit_VV; 
      yield_sfit_VV = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_VV", jet)); 
      if ( yield_sfit_VV != 0x0 ) {
	yield_VV_sfit = yield_sfit_VV->getVal();  
	if ( debug ) 
	  std::cout << "yield_VV_sfit = " << yield_VV_sfit << "\n";
	h_VV_s->Fill( (yield_VV_sfit-yield_VV)/yield_VV );
      }
      
      // Ztt 
      RooRealVar* yield_sfit_Ztt; 
      yield_sfit_Ztt = (RooRealVar*) norm_fit_s->find(Form("n_exp_final_binj%iof_proc_Ztt", jet)); 
      if ( yield_sfit_Ztt != 0x0 ) {
	yield_Ztt_sfit = yield_sfit_Ztt->getVal();  
	if ( debug ) 
	  std::cout << "yield_Ztt_sfit = " << yield_Ztt_sfit << "\n";
	h_Ztt_s->Fill( (yield_Ztt_sfit-yield_Ztt)/yield_Ztt );
      }
      
    }
    delete fit_s;
    File->Close();
  }
  
  TCanvas* c1 = new TCanvas();
  // plotting PullMu
  h_PullMu_s->SetTitle("");
  h_PullMu_s->GetXaxis()->SetTitle("(fitted #mu - 1)/#sigma_{fitted #mu}");
  h_PullMu_s->GetYaxis()->SetTitle("toys/bin");
  h_PullMu_s->Draw();
  h_PullMu_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_PullMu_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_PullMu_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));


  // plot ggH yield 
  h_ggH_s->SetTitle("");
  h_ggH_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_ggH_s->GetYaxis()->SetTitle("toys/bin");
  h_ggH_s->Draw();
  h_ggH_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_ggH_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_ggH_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot qqWW yield 
  h_qqWW_s->SetTitle("");
  h_qqWW_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_qqWW_s->GetYaxis()->SetTitle("toys/bin");
  h_qqWW_s->Draw();
  h_qqWW_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_qqWW_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_qqWW_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  

  // plot ggWW yield 
  h_ggWW_s->SetTitle("");
  h_ggWW_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_ggWW_s->GetYaxis()->SetTitle("toys/bin");
  h_ggWW_s->Draw();
  h_ggWW_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_ggWW_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_ggWW_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  
  // plot Top yield 
  h_Top_s->SetTitle("");
  h_Top_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_Top_s->GetYaxis()->SetTitle("toys/bin");
  h_Top_s->Draw();
  h_Top_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Top_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Top_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot WjetsE yield 
  h_WjetsE_s->SetTitle("");
  h_WjetsE_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_WjetsE_s->GetYaxis()->SetTitle("toys/bin");
  h_WjetsE_s->Draw();
  h_WjetsE_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_WjetsE_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_WjetsE_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));


  // plot WjetsM yield 
  h_WjetsM_s->SetTitle("");
  h_WjetsM_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_WjetsM_s->GetYaxis()->SetTitle("toys/bin");
  h_WjetsM_s->Draw();
  h_WjetsM_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_WjetsM_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_WjetsM_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot Wgamma yield 
  h_Wgamma_s->SetTitle("");
  h_Wgamma_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_Wgamma_s->GetYaxis()->SetTitle("toys/bin");
  h_Wgamma_s->Draw();
  h_Wgamma_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Wgamma_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Wgamma_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot Wg3l yield 
  h_Wg3l_s->SetTitle("");
  h_Wg3l_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_Wg3l_s->GetYaxis()->SetTitle("toys/bin");
  h_Wg3l_s->Draw();
  h_Wg3l_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Wg3l_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Wg3l_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot Ztt yield 
  h_Ztt_s->SetTitle("");
  h_Ztt_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_Ztt_s->GetYaxis()->SetTitle("toys/bin");
  h_Ztt_s->Draw();
  h_Ztt_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Ztt_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_Ztt_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // plot VV yield 
  h_VV_s->SetTitle("");
  h_VV_s->GetXaxis()->SetTitle("(N_{fit,s}-N_{in})/N_{in}");
  h_VV_s->GetYaxis()->SetTitle("toys/bin");
  h_VV_s->Draw();
  h_VV_s->Fit("gaus");
  gPad->Update();
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_VV_%s.png",dir.Data(),inj.Data(),jet,mH,ana.Data()));
  c1->SaveAs(Form("%s/plots/norm_inj%s_%ij_%i_sfit_VV_%s.eps",dir.Data(),inj.Data(),jet,mH,ana.Data()));

  // delete 
  delete h_PullMu_s;
  delete h_ggH_s;
  delete h_qqWW_s;
  delete h_ggWW_s;
  delete h_Top_s;
  delete h_WjetsE_s;
  delete h_WjetsM_s;
  delete h_Wgamma_s;
  delete h_Wg3l_s;
  delete h_VV_s;
  delete h_Ztt_s;
  delete c1; 
  
}
 
/*

*/
